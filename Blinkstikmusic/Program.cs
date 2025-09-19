using System;
using System.Numerics;
using System.Threading;
using NAudio.CoreAudioApi;
using NAudio.Wave;
using BlinkStickDotNet;

class Program
{
    // Configuration
    const int SampleRate = 48000;       // WASAPI default on many systems
    const int BufferMs = 20;            // analysis/update interval
    const int FftSize = 1024;           // power of two
    const int LedCount = 8;             // BlinkStick Square has 8 LEDs

    static volatile bool _running = true;

    [STAThread]
    static void Main()
    {
        Console.CancelKeyPress += (_, e) => { e.Cancel = true; _running = false; };

        Console.WriteLine("BlinkStick music visualizer starting... Press Ctrl+C to stop.");

        // Find BlinkStick device
        var blinkstick = BlinkStick.FindFirst();
        if (blinkstick == null)
        {
            Console.WriteLine("No BlinkStick found. Plug in your BlinkStick Square and try again.");
            return;
        }

        // Try to open device
        try { blinkstick.OpenDevice(); } catch { }

        // Ensure all off on start
        SetAllLeds(blinkstick, (0, 0, 0));

        using var device = new MMDeviceEnumerator().GetDefaultAudioEndpoint(DataFlow.Render, Role.Multimedia);
        using var capture = new WasapiLoopbackCapture(device);

        // Convert captured audio to 32-bit float mono for analysis
        var desiredFormat = WaveFormat.CreateIeeeFloatWaveFormat(SampleRate, 2); // keep stereo for resampler
        IWaveProvider provider = new WasapiCaptureProvider(capture);
        var resampler = new MediaFoundationResampler(provider, desiredFormat) { ResamplerQuality = 60 };
        ISampleProvider sampleProvider = resampler.ToSampleProvider();
        if (sampleProvider.WaveFormat.Channels > 1)
        {
            var toMono = new NAudio.Wave.SampleProviders.StereoToMonoSampleProvider(sampleProvider)
            {
                LeftVolume = 0.5f,
                RightVolume = 0.5f
            };
            sampleProvider = toMono;
        }

        int frameSize = (SampleRate * BufferMs) / 1000; // samples per frame
        float[] frameBuffer = new float[frameSize];
        float[] window = HannWindow(frameSize);

        // Frequency bands for 8 LEDs (approx): 63, 125, 250, 500, 1k, 2k, 4k, 8k Hz
        float[] bandEdges = new float[] { 31.5f, 90f, 180f, 350f, 700f, 1400f, 2800f, 5600f, 11200f };

        var readThread = new Thread(() =>
        {
            var reader = sampleProvider;
            var sampleReadBuffer = new float[frameSize];
            while (_running)
            {
                int read = reader.Read(sampleReadBuffer, 0, frameSize);
                if (read <= 0)
                {
                    Thread.Sleep(5);
                    continue;
                }

                // Copy into analysis buffer
                Array.Copy(sampleReadBuffer, 0, frameBuffer, 0, read);
                for (int i = read; i < frameSize; i++) frameBuffer[i] = 0f;

                // Apply window
                for (int i = 0; i < frameSize; i++) frameBuffer[i] *= window[i];

                // Magnitude per band using FFT
                var spectrum = MagnitudeSpectrum(frameBuffer, FftSize);
                var bands = MapBands(spectrum, SampleRate, bandEdges);

                // Normalize and colorize
                var colors = BandsToColors(bands);
                try
                {
                    for (int i = 0; i < LedCount; i++)
                    {
                        var (r, g, b) = colors[Math.Min(i, colors.Length - 1)];
                        blinkstick.SetColor(0, (byte)i, r, g, b); // channel 0, indexed LED
                    }
                }
                catch (Exception ex)
                {
                    Console.WriteLine($"BlinkStick write error: {ex.Message}");
                }
            }
        }) { IsBackground = true, Name = "AudioReader" };

        capture.StartRecording();
        readThread.Start();

        while (_running) Thread.Sleep(50);

        // Cleanup
        try
        {
            SetAllLeds(blinkstick, (0, 0, 0));
            capture.StopRecording();
        }
        catch { }

        Console.WriteLine("Stopped.");
    }

    static void SetAllLeds(BlinkStick device, (byte r, byte g, byte b) color)
    {
        for (int i = 0; i < LedCount; i++)
        {
            device.SetColor(0, (byte)i, color.r, color.g, color.b);
        }
    }

    static float[] HannWindow(int n)
    {
        var w = new float[n];
        for (int i = 0; i < n; i++) w[i] = 0.5f - 0.5f * (float)Math.Cos(2 * Math.PI * i / (n - 1));
        return w;
    }

    static float[] MagnitudeSpectrum(float[] samples, int fftSize)
    {
        // zero-pad to fftSize          
        System.Numerics.Complex[] buf = new System.Numerics.Complex[fftSize];
        int copy = Math.Min(samples.Length, fftSize);
        for (int i = 0; i < copy; i++) buf[i] = new System.Numerics.Complex(samples[i], 0);
        for (int i = copy; i < fftSize; i++) buf[i] = System.Numerics.Complex.Zero;

        // naive FFT (Cooley-Tukey)
        FFT(buf);

        int bins = fftSize / 2;
        var mag = new float[bins];
        for (int i = 0; i < bins; i++) mag[i] = (float)buf[i].Magnitude;
        return mag;
    }

    static void FFT(System.Numerics.Complex[] buffer)
    {
        int n = buffer.Length;
        int bits = (int)Math.Log2(n);

        // bit-reverse
        for (int j = 1, i = 0; j < n; j++)
        {
            int bit = n >> 1;
            for (; (i & bit) != 0; bit >>= 1) i &= ~bit;
            i |= bit;
            if (j < i) (buffer[j], buffer[i]) = (buffer[i], buffer[j]);
        }

        for (int len = 2; len <= n; len <<= 1)
        {
            double ang = -2 * Math.PI / len;
            var wlen = new System.Numerics.Complex(Math.Cos(ang), Math.Sin(ang));
            for (int i = 0; i < n; i += len)
            {
                var w = System.Numerics.Complex.One;
                for (int j = 0; j < len / 2; j++)
                {
                    var u = buffer[i + j];
                    var v = buffer[i + j + len / 2] * w;
                    buffer[i + j] = u + v;
                    buffer[i + j + len / 2] = u - v;
                    w *= wlen;
                }
            }
        }
    }

    static float[] MapBands(float[] spectrum, int sampleRate, float[] bandEdges)
    {
        int fftSize = spectrum.Length * 2;
        float binHz = sampleRate / (float)fftSize;
        var bands = new float[bandEdges.Length - 1];
        for (int b = 0; b < bands.Length; b++)
        {
            int startBin = Math.Max(1, (int)(bandEdges[b] / binHz));
            int endBin = Math.Min(spectrum.Length - 1, (int)(bandEdges[b + 1] / binHz));
            float sum = 0f;
            for (int i = startBin; i <= endBin; i++) sum += spectrum[i];
            float avg = sum / Math.Max(1, (endBin - startBin + 1));
            bands[b] = avg;
        }

        // simple dynamic range compression
        for (int b = 0; b < bands.Length; b++)
        {
            float x = bands[b];
            x = (float)Math.Log10(1 + 9 * x); // 0..1-ish
            bands[b] = x;
        }
        return bands;
    }

    static (byte r, byte g, byte b)[] BandsToColors(float[] bands)
    {
        var colors = new (byte r, byte g, byte b)[bands.Length];
        for (int i = 0; i < bands.Length; i++)
        {
            float v = Math.Clamp(bands[i], 0f, 1f);
            var (r, g, b) = HeatColor(v);
            colors[i] = ((byte)r, (byte)g, (byte)b);
        }
        return colors;
    }

    // Maps 0..1 to blue->cyan->green->yellow->red
    static (byte r, byte g, byte b) HeatColor(float v)
    {
        v = Math.Clamp(v, 0f, 1f);
        float h = (1f - v) * 240f; // hue 240 (blue) to 0 (red)
        return HsvToRgb(h, 1f, v * 0.9f + 0.1f);
    }

    static (byte r, byte g, byte b) HsvToRgb(float h, float s, float v)
    {
        float c = v * s;
        float x = c * (1 - Math.Abs((h / 60f) % 2 - 1));
        float m = v - c;
        float r1 = 0, g1 = 0, b1 = 0;
        if (h < 60) { r1 = c; g1 = x; b1 = 0; }
        else if (h < 120) { r1 = x; g1 = c; b1 = 0; }
        else if (h < 180) { r1 = 0; g1 = c; b1 = x; }
        else if (h < 240) { r1 = 0; g1 = x; b1 = c; }
        else if (h < 300) { r1 = x; g1 = 0; b1 = c; }
        else { r1 = c; g1 = 0; b1 = x; }
        byte r = (byte)Math.Round((r1 + m) * 255);
        byte g = (byte)Math.Round((g1 + m) * 255);
        byte b = (byte)Math.Round((b1 + m) * 255);
        return (r, g, b);
    }
}

// NAudio helper to expose WasapiLoopbackCapture as IWaveProvider
sealed class WasapiCaptureProvider : IWaveProvider, IDisposable
{
    private readonly WasapiCapture _capture;
    private readonly BufferedWaveProvider _buffered;

    public WasapiCaptureProvider(WasapiCapture capture)
    {
        _capture = capture;
        _buffered = new BufferedWaveProvider(capture.WaveFormat)
        {
            DiscardOnBufferOverflow = true,
            BufferLength = capture.WaveFormat.AverageBytesPerSecond * 2
        };
        capture.DataAvailable += (s, a) => _buffered.AddSamples(a.Buffer, 0, a.BytesRecorded);
    }

    public WaveFormat WaveFormat => _buffered.WaveFormat;

    public int Read(byte[] buffer, int offset, int count) => _buffered.Read(buffer, offset, count);

    public void Dispose() { _capture?.Dispose(); }
}
