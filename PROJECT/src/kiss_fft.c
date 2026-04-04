#include "_kiss_fft_guts.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

static void kf_bfly2(kiss_fft_cpx* Fout, size_t fstride, kiss_fft_cfg st, int m) {
    kiss_fft_cpx* Fout2 = Fout + m;
    kiss_fft_cpx* tw1 = st->twiddles;
    kiss_fft_cpx t;
    int i;
    for (i = 0; i < m; ++i) {
        C_MUL(t, Fout2[i], tw1[i]);
        C_SUB(Fout2[i], Fout[i], t);
        C_ADDTO(Fout[i], t);
    }
}

static void kf_bfly4(kiss_fft_cpx* Fout, size_t fstride, kiss_fft_cfg st, int m) {
    kiss_fft_cpx* tw1 = st->twiddles;
    kiss_fft_cpx* tw2 = tw1 + m;
    kiss_fft_cpx* tw3 = tw2 + m;
    kiss_fft_cpx scratch[6];
    int
    for (i = 0; i < m; ++i) {
        C_MUL(scratch[0], Fout[i + m], tw1[i]);
        C_MUL(scratch[1], Fout[i + 2*m], tw2[i]);
        C_MUL(scratch[2], Fout[i + 3*m], tw3[i]);
        C_SUB(scratch[5], Fout[i], scratch[1]);
        C_ADDTO(Fout[i], scratch[1]);
        C_ADD(scratch[3], scratch[0], scratch[2]);
        C_SUB(scratch[4], scratch[0], scratch[2]);
        C_SUB(Fout[i + 2*m], Fout[i], scratch[3]);
        C_ADDTO(Fout[i], scratch[3]);
        C_MUL(scratch[0], scratch[5], st->twiddles[fstride * m]);
        C_MUL(scratch[1], scratch[4], st->twiddles[fstride * 2 * m]);
        C_ADD(Fout[i + m], scratch[0], scratch[1]);
        C_SUB(Fout[i + 3*m], scratch[0], scratch[1]);
    }
}

static void kf_work(kiss_fft_cpx* Fout, const kiss_fft_cpx* f, size_t fstride, int in_stride, int* factors, kiss_fft_cfg st) {
    const int p = *factors++;
    const int m = *factors++;
    const kiss_fft_cpx* f_cur = f;
    kiss_fft_cpx* Fout_end = Fout + p * m;
    int i;

    for (; Fout != Fout_end; ++Fout) {
        *Fout = *f_cur;
        f_cur += in_stride;
    }

    if (p == 2)
        kf_bfly2(Fout_end - p * m, fstride, st, m);
    else if (p == 4)
        kf_bfly4(Fout_end - p * m, fstride, st, m);
}

static void kf_factor(int n, int* facbuf) {
    int p = 4;
    double floor_sqrt = floor(sqrt((double)n));
    while (n > 1) {
        while (n % p) {
            switch (p) {
            case 4: p = 2; break;
            case 2: p = 3; break;
            default: p += 2; break;
            }
            if (p > floor_sqrt) p = n;
        }
        n /= p;
        *facbuf++ = p;
        *facbuf++ = n;
    }
}

kiss_fft_cfg kiss_fft_alloc(int nfft, int inverse_fft, void* mem, size_t* lenmem) {
    kiss_fft_cfg st = NULL;
    size_t memneeded = sizeof(struct kiss_fft_state) + sizeof(kiss_fft_cpx) * (nfft - 1);
    if (lenmem == NULL) {
        st = (kiss_fft_cfg)malloc(memneeded);
    } else {
        if (mem != NULL && *lenmem >= memneeded)
            st = (kiss_fft_cfg)mem;
        *lenmem = memneeded;
    }
    if (st) {
        int i;
        st->nfft = nfft;
        st->inverse = inverse_fft;
        for (i = 0; i < nfft; ++i) {
            double phase = -2 * M_PI * i / nfft;
            if (inverse_fft) phase *= -1;
            kf_cexp(st->twiddles + i, phase);
        }
        kf_factor(nfft, st->factors);
    }
    return st;
}

void kiss_fft(kiss_fft_cfg cfg, const kiss_fft_cpx* fin, kiss_fft_cpx* fout) {
    int factors[2 * MAXFACTORS];
    memcpy(factors, cfg->factors, sizeof(int) * (cfg->factors[0] * 2 + 2));
    kf_work(fout, fin, 1, 1, factors, cfg);
    if (cfg->inverse) {
        int i;
        double scaling = 1.0 / cfg->nfft;
        for (i = 0; i < cfg->nfft; ++i) {
            fout[i].r *= (float)scaling;
            fout[i].i *= (float)scaling;
        }
    }
}

int kiss_fft_next_fast_size(int n) {
    while (1) {
        int m = n;
        while ((m % 2) == 0) m /= 2;
        while ((m % 3) == 0) m /= 3;
        while ((m % 5) == 0) m /= 5;
        if (m <= 1) break;
        n++;
    }
    return n;
}

size_t kiss_fft_cleanup(void) { return 0; }