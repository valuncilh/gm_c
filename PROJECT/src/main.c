#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <windows.h>
#include <commdlg.h>
#include <direct.h>      // для _mkdir (если понадобится, но мы его не используем)

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "graphics.h"
#include "fft_processor.h"

#define WIDTH 800
#define HEIGHT 600
#define NUM_FRAMES 100
#define MAX_SIGNAL_SIZE 2048

static char szFileFilter[] = "Text Files (*.txt)\0*.txt\0All Files (*.*)\0*.*\0";

static int open_file_dialog(char* out_path, size_t out_size) {
    OPENFILENAMEA ofn = {0};
    char file_name[MAX_PATH] = "";
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = NULL;
    ofn.lpstrFilter = szFileFilter;
    ofn.lpstrFile = file_name;
    ofn.nMaxFile = MAX_PATH;
    ofn.Flags = OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;
    if (GetOpenFileNameA(&ofn)) {
        strncpy(out_path, file_name, out_size);
        out_path[out_size - 1] = '\0';
        return 1;
    }
    return 0;
}

static int read_signal_from_file(const char* filename, double* signal, int max_size) {
    FILE* f = fopen(filename, "r");
    if (!f) return 0;
    int count = 0;
    while (count < max_size && fscanf(f, "%lf", &signal[count]) == 1) {
        count++;
    }
    fclose(f);
    return count;
}

static void generate_test_signal(double* signal, int size) {
    double Fs = size;
    double f1 = 101.3, f2 = 253.7;
    double A1 = 1.0, A2 = 0.7;
    for (int i = 0; i < size; i++) {
        double t = i / Fs;
        signal[i] = A1 * sin(2 * M_PI * f1 * t) + A2 * sin(2 * M_PI * f2 * t);
    }
    printf("Using default test signal (f1=%.1f, f2=%.1f), N=%d\n", f1, f2, size);
}

static void draw_spectrum_colored(uint8_t* pixels, int width, int height,
                                  double* magnitudes, int num_bins,
                                  double max_mag, double factor,
                                  int x_start, int y_base, int bar_width) {
    int max_bar_height = y_base - 30;
    for (int i = 0; i < num_bins; i++) {
        int x = x_start + i * bar_width;
        double amp = magnitudes[i] * factor;
        int bar_height = (int)((amp / max_mag) * max_bar_height);
        if (bar_height < 0) bar_height = 0;
        // Логарифмическая интенсивность для цвета (хорошо показывает малые амплитуды)
        double intensity_linear = amp / max_mag;
        double intensity_log = log10(1.0 + 9.0 * intensity_linear) / log10(10.0);
        if (intensity_log < 0) intensity_log = 0;
        if (intensity_log > 1) intensity_log = 1;
        uint8_t r = (uint8_t)(255 * intensity_log);
        uint8_t g = (uint8_t)(255 * (1 - intensity_log));
        uint8_t b = 0;
        for (int y = 0; y < bar_height; y++) {
            int py = y_base - y;
            for (int dx = 0; dx < bar_width; dx++) {
                set_pixel(pixels, width, height, x + dx, py, r, g, b);
            }
        }
        draw_line(pixels, width, height, x, y_base, x, y_base - bar_height, 0, 0, 0);
    }
}

int main() {
    // Сохраняем кадры прямо в текущую папку (без создания подпапки)
    char cwd[1024];
    if (_getcwd(cwd, sizeof(cwd))) {
        printf("Saving frames to: %s\n", cwd);
    }

    double* signal = (double*)malloc(sizeof(double) * MAX_SIGNAL_SIZE);
    int N = 0;

    printf("Opening file selection dialog...\n");
    char filepath[MAX_PATH];
    if (open_file_dialog(filepath, sizeof(filepath))) {
        N = read_signal_from_file(filepath, signal, MAX_SIGNAL_SIZE);
        if (N > 0) {
            printf("Loaded %d samples from file: %s\n", N, filepath);
        } else {
            printf("File is empty or contains invalid data. Using test signal.\n");
        }
    } else {
        printf("No file selected. Using test signal.\n");
    }

    if (N == 0) {
        N = 1024;
        generate_test_signal(signal, N);
    }

    // Дополнение до степени двойки
    if ((N & (N - 1)) != 0) {
        int newN = 1;
        while (newN < N) newN <<= 1;
        if (newN > MAX_SIGNAL_SIZE) newN = MAX_SIGNAL_SIZE;
        for (int i = N; i < newN; i++) signal[i] = 0.0;
        N = newN;
        printf("Padded to %d samples (power of two).\n", N);
    }

    double* magnitudes = (double*)malloc(sizeof(double) * (N / 2));
    compute_fft(signal, N, magnitudes);

    double max_mag = 0;
    for (int i = 0; i < N / 2; i++) {
        if (magnitudes[i] > max_mag) max_mag = magnitudes[i];
    }
    if (max_mag == 0) max_mag = 1;
    printf("Max magnitude: %f\n", max_mag);

    int x_start = 50;
    int y_base = HEIGHT - 50;
    int bar_width = (WIDTH - x_start - 30) / (N / 2);
    if (bar_width < 1) bar_width = 1;

    for (int frame = 0; frame < NUM_FRAMES; frame++) {
        uint8_t* pixels = (uint8_t*)malloc(WIDTH * HEIGHT * 3);
        memset(pixels, 255, WIDTH * HEIGHT * 3);

        double t = (double)frame / (NUM_FRAMES - 1);
        double factor = t;

        draw_spectrum_colored(pixels, WIDTH, HEIGHT,
                              magnitudes, N / 2, max_mag, factor,
                              x_start, y_base, bar_width);

        char fname[256];
        snprintf(fname, sizeof(fname), "frame_%04d.png", frame);
        if (!stbi_write_png(fname, WIDTH, HEIGHT, 3, pixels, WIDTH * 3)) {
            printf("Error writing file: %s\n", fname);
        }
        free(pixels);
        printf("Frame %d/%d done\n", frame + 1, NUM_FRAMES);
    }

    free(signal);
    free(magnitudes);
    printf("Done! Frames saved in current directory.\n");
    return 0;
}