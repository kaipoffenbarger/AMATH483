#include <iostream>
#include <vector>
#include <string>
#include <thread>
#include <chrono>
#include <cstring>
#include <png.h>

// Grayscale conversion for a portion of rows
void grayscaleWorker(png_bytep* image, int startRow, int endRow, int width, int channels) {
    for (int y = startRow; y < endRow; ++y) {
        png_bytep row = image[y];
        for (int x = 0; x < width; ++x) {
            png_bytep px = &(row[x * channels]);
            uint8_t gray = static_cast<uint8_t>(0.3 * px[0] + 0.59 * px[1] + 0.11 * px[2]);
            px[0] = px[1] = px[2] = gray;
            // if RGBA, preserve alpha
        }
    }
}

// Threaded grayscale function
void grayscaleThreaded(
    png_bytep* image,
    int width, int height, int channels,
    int numThreads)
{
    std::vector<std::thread> threads;
    int rowsPerThread = height / numThreads;
    int remainingRows = height % numThreads;
    int currentRow = 0;

    for (int i = 0; i < numThreads; ++i) {
        int rowsToProcess = rowsPerThread + (i < remainingRows ? 1 : 0);
        int startRow = currentRow;
        int endRow = startRow + rowsToProcess;

        threads.emplace_back(grayscaleWorker, image, startRow, endRow, width, channels);

        currentRow = endRow;
    }

    for (auto& t : threads) {
        t.join();
    }
}