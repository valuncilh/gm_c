#ifndef KISS_FFT_H
#define KISS_FFT_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float r;
    float i;
} kiss_fft_cpx;

typedef struct kiss_fft_state* kiss_fft_cfg;

kiss_fft_cfg kiss_fft_alloc(int nfft, int inverse_fft, void* mem, size_t* lenmem);
void kiss_fft(kiss_fft_cfg cfg, const kiss_fft_cpx* fin, kiss_fft_cpx* fout);
void kiss_fft_stride(kiss_fft_cfg cfg, const kiss_fft_cpx* fin, kiss_fft_cpx* fout, int fin_stride);
#define kiss_fft_free free
size_t kiss_fft_cleanup(void);
int kiss_fft_next_fast_size(int n);

#ifdef __cplusplus
}
#endif

#endif