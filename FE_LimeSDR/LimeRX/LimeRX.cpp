#define _CRT_SECURE_NO_WARNINGS

#include "lime/LimeSuite.h"
#include <iostream>
#include <thread>
#include <queue>
#include <vector>
#include <mutex>
#include <condition_variable>
#include <chrono>
#include <cstring>
#include <atomic>
#include <cstdio>
#include <cmath>
#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

using namespace std;

constexpr double CENTER_FREQ = 2492.028e6; // 2492.028 MHz
constexpr double SAMPLE_RATE = 12e6;       // 12 MHz
constexpr float  NORMALIZED_GAIN = (float)0.7;    // RX gain
//constexpr unsigned BUFERSIZE = 16384;      // Complex samples per buffer, 16k
constexpr unsigned BUFERSIZE = 65536;      // Complex samples per buffer, 64k
//constexpr size_t QUEUE_MAX_BUFFERS = 256;  // Max queue size before RX waits
constexpr size_t QUEUE_MAX_BUFFERS = 1024;  // Max queue size before RX waits
const char* DEFAULT_OUTPUT_FILENAME = "iq_samples_2bit.bin";
constexpr double ALPHA = 0.001;        // Smoothing factor for stddev

// Thread-safe queue for IQ buffers
class BufferQueue {
public:
    void push(vector<int16_t>&& buf) {
        unique_lock<mutex> lock(mtx);
        cv_full.wait(lock, [this] { return queue.size() < QUEUE_MAX_BUFFERS; });
        queue.push(move(buf));
        cv_empty.notify_one();
    }
    bool pop(vector<int16_t>& buf) {
        unique_lock<mutex> lock(mtx);
        cv_empty.wait(lock, [this] { return !queue.empty() || finished; });
        if (queue.empty())
            return false;
        buf = move(queue.front());
        queue.pop();
        cv_full.notify_one();
        return true;
    }
    void set_finished() {
        unique_lock<mutex> lock(mtx);
        finished = true;
        cv_empty.notify_all();
    }
    bool is_finished() const {
        return finished;
    }
private:
    queue<vector<int16_t>> queue;
    mutable mutex mtx;
    condition_variable cv_empty, cv_full;
    bool finished = false;
};

// Helper for error handling
void error(lms_device_t* device) {
    if (device) LMS_Close(device);
    exit(-1);
}

// RX thread: receives samples and pushes to queue
void rx_thread_func(lms_device_t* device, lms_stream_t* streamId, BufferQueue* bufq, atomic<bool>* rx_done) {
    while (!rx_done->load()) {
        vector<int16_t> buffer(BUFERSIZE * 2);
        int samplesRead = LMS_RecvStream(streamId, buffer.data(), BUFERSIZE, nullptr, 1000);
        if (samplesRead > 0) {
            buffer.resize(samplesRead * 2);
            bufq->push(move(buffer));
        }
    }
    bufq->set_finished();
}

// Statistics thread: prints RX stats once per second
void stats_thread_func(lms_stream_t* streamId, atomic<bool>* stats_done) {
    using namespace std::chrono_literals;
    while (!stats_done->load()) {
        lms_stream_status_t status;
        if (LMS_GetStreamStatus(streamId, &status) == 0) {
            cout << "RX data rate: " << status.linkRate / 1e6 << " MB/s"
                << " | Dropped: " << status.droppedPackets
                << " | RX FIFO: " << (100 * status.fifoFilledCount) / status.fifoSize << "%\n";
        }
        this_thread::sleep_for(1s);
    }
}

// Quantization function: converts a 16-bit sample to {-3,-1,+1,+3} based on stddev
inline int8_t quantize_2bit(int16_t sample, double stddev) {
    if (sample >= 0) {
        if (abs(sample) < stddev) return +1;
        else return +3;
    }
    else {
        if (abs(sample) < stddev) return -1;
        else return -3;
    }
}

// Writer thread: quantizes in real time using stddev and writes as 2-bit I/Q
void writer_thread_func(BufferQueue* bufq, const char* filename) {
    bool is_stdout = (strcmp(filename, "-") == 0);
    FILE* out = nullptr;
    if (is_stdout) {
#ifdef _WIN32
        _setmode(_fileno(stdout), _O_BINARY);
#endif
        out = stdout;
    }
    else {
        out = fopen(filename, "wb");
        if (!out) {
            cerr << "Error opening output file!" << endl;
            return;
        }
    }

    // 2-bit I/Q
    double mean = 0.0, var = 0.0, stddev = 1.0;
    size_t sample_count = 0;

    vector<int16_t> buffer;
    vector<int8_t> outbuf;

    while (bufq->pop(buffer)) {
        if (buffer.empty()) continue;
        size_t n = buffer.size();
        outbuf.resize(n);

        for (size_t i = 0; i < n; ++i) {
            int16_t x = buffer[i];
            ++sample_count;
            // Smoothing stddev
            mean = (1.0 - ALPHA) * mean + ALPHA * x;
            var = (1.0 - ALPHA) * var + ALPHA * (x - mean) * (x - mean);
            stddev = sqrt(var);
            // Quantize
            outbuf[i] = quantize_2bit(x, stddev);
        }
        fwrite(outbuf.data(), sizeof(int8_t), n, out);
    }

    cout << "Final mean: " << mean << ", stddev: " << stddev << ", total samples: " << sample_count / 2 << " complex\n";
    if (!is_stdout && out) fclose(out);
}

void print_usage(const char* prog) {
    cout << "Usage: " << prog << " [--time <seconds>] [--file <filename>]\n"
        "  --time <seconds>   Set capture duration (default 10 seconds)\n"
        "  --file <filename>  Output file (default: iq_samples_8bit.bin, use '-' for stdout)\n";
}

int main(int argc, char** argv)
{
    int capture_seconds = 10; // default
    string output_filename = DEFAULT_OUTPUT_FILENAME;
    // Parse command line
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--time") == 0 && i + 1 < argc) {
            capture_seconds = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "--file") == 0 && i + 1 < argc) {
            output_filename = argv[++i];
        }
        else if (strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        }
    }
    cout << "Capture duration: " << capture_seconds << " seconds\n";
    cout << "Output file: " << output_filename << "\n";

    lms_device_t* device = nullptr;

    // Find and open device
    int n;
    if ((n = LMS_GetDeviceList(nullptr)) < 0) error(device);
    if (n < 1) {
        cerr << "No devices found.\n";
        return -1;
    }
    lms_info_str_t* list = new lms_info_str_t[n];
    if (LMS_GetDeviceList(list) < 0) error(device);

    if (LMS_Open(&device, list[0], nullptr)) error(device);
    delete[] list;

    if (LMS_Init(device) != 0) error(device);

    // Enable RX channel
    if (LMS_EnableChannel(device, LMS_CH_RX, 0, true) != 0) error(device);

    // Set frequency and sample rate
    if (LMS_SetLOFrequency(device, LMS_CH_RX, 0, CENTER_FREQ) != 0) error(device);
    if (LMS_SetSampleRate(device, SAMPLE_RATE, 8) != 0) error(device);

    // Set LPF, gain, and calibrate
    if (LMS_SetLPFBW(device, LMS_CH_RX, 0, SAMPLE_RATE) != 0) error(device);
    if (LMS_SetNormalizedGain(device, LMS_CH_RX, 0, NORMALIZED_GAIN) != 0) error(device);
    if (LMS_Calibrate(device, LMS_CH_RX, 0, SAMPLE_RATE, 0) != 0) error(device);

    // Stream setup
    lms_stream_t streamId;
    streamId.channel = 0;
    //streamId.fifoSize = 1024 * 1024; // SDR-side buffer
    streamId.fifoSize = 8 * 1024 * 1024; // 8MB SDR-side buffer
    streamId.throughputVsLatency = 1.0;
    streamId.isTx = false;
    streamId.dataFmt = lms_stream_t::LMS_FMT_I12;
    if (LMS_SetupStream(device, &streamId) != 0) error(device);

    LMS_StartStream(&streamId);

    BufferQueue bufq;
    atomic<bool> rx_done(false);
    atomic<bool> stats_done(false);

    // Start threads
    thread rx_thread(rx_thread_func, device, &streamId, &bufq, &rx_done);
    thread writer_thread(writer_thread_func, &bufq, output_filename.c_str());
    thread stats_thread(stats_thread_func, &streamId, &stats_done);

    // Capture for specified seconds
    this_thread::sleep_for(chrono::seconds(capture_seconds));
    rx_done = true;

    // Wait for RX thread to finish and queue to drain
    rx_thread.join();
    writer_thread.join();

    // Stop stats thread after RX is done
    stats_done = true;
    stats_thread.join();

    LMS_StopStream(&streamId);
    LMS_DestroyStream(device, &streamId);
    LMS_Close(device);

    cout << "I/Q data captured and quantized to '" << output_filename << "' for " << capture_seconds << " seconds.\n";
    return 0;
}
