#include "fft_processor.h"
#include "kiss_fft.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

void generate_signal(double* signal, int size) {
    for (int i = 0; i < size; i++) {
        double t = (double)i / size;
        signal[i] = sin(2 * M_PI * 5 * t) + 0.5 * sin(2 * M_PI * 12 * t) + 0.3 * sin(2 * M_PI * 20 * t);
    }
}

void compute_fft(double* signal, int size, double* magnitudes) {
    kiss_fft_cfg cfg = kiss_fft_alloc(size, 0, NULL, NULL);
    kiss_fft_cpx* in = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * size);
    kiss_fft_cpx* out = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * size);
    for (int i = 0; i < size; i++) {
        in[i].r = (float)signal[i];
        in[i].i = 0.0f;
    }
    kiss_fft(cfg, in, out);
    for (int i = 0; i < size / 2; i++) {
        magnitudes[i] = sqrt(out[i].r * out[i].r + out[i].i * out[i].i);
    }
    free(in);
    free(out);
    free(cfg);
}