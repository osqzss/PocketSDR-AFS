// 
//  Pocket SDR C Library - GNSS SDR Receiver Channel Functions
//
//  Author:
//  T.TAKASU
//
//  History:
//  2022-07-08  1.0  port sdr_ch.py to C
//  2022-07-16  1.1  modify API sdr_ch_new()
//  2023-12-16  1.2  reduce memory usage
//  2023-12-28  1.3  support type and API changes.
//  2024-01-12  1.4  ch->state: const char * -> int
//  2024-01-16  1.5  add doppler assist for acquisition
//  2024-03-20  1.6  modify API sdr_ch_update()
//  2024-04-28  1.7  modify API sdr_ch_new()
//  2024-06-06  1.8  modify API sdr_ch_new()
//  2024-06-10  1.9  add API sdr_ch_stat_req(), sdr_ch_stat_get()
//
#include <ctype.h>
#include <math.h>
#include "pocket_sdr.h"

// constants and macros --------------------------------------------------------
#define SP_CORR    0.25     // correlator spacing (chip)
#define T_ACQ      0.010    // non-coherent integration time for acquis. (s)
#define T_DLL      0.010    // non-coherent integration time for DLL (s)
#define T_CN0      1.0      // averaging time for C/N0 (s)
#define T_FPULLIN  1.0      // frequency pullin time (s)
#define T_NPULLIN  1.5      // navigation data pullin time (s)
#define B_DLL      0.25     // band-width of DLL filter (Hz)
#define B_PLL      5.0      // band-width of PLL filter (Hz)
#define B_FLL_W    5.0      // band-width of FLL filter (Hz) (wide)
#define B_FLL_N    2.0      // band-width of FLL filter (Hz) (narrow)
#define MAX_DOP    5000.0   // max Doppler for acquisition (Hz)
#define THRES_CN0_L 35.0    // C/N0 threshold (dB-Hz) (lock)
#define THRES_CN0_U 32.0    // C/N0 threshold (dB-Hz) (lost)
#define THRES_SYNC 0.02     // threshold for sec-code sync
#define THRES_LOST 0.002    // threshold for sec-code lost
#define N_CODE     10       // number of resampled code bank
#define ADD_CORR   40       // number of additional correlators

#define DPI        (2.0 * PI)
#define SQR(x)     ((x) * (x)) 
#define MAX(x, y)  ((x) > (y) ? (x) : (y))
#define MIN(x, y)  ((x) < (y) ? (x) : (y))

// global variables ------------------------------------------------------------
double sdr_sp_corr = SP_CORR;
double sdr_t_acq   = T_ACQ;
double sdr_t_dll   = T_DLL;
double sdr_b_dll   = B_DLL;
double sdr_b_pll   = B_PLL;
double sdr_b_fll_w = B_FLL_W;
double sdr_b_fll_n = B_FLL_N;
double sdr_max_dop = MAX_DOP;
double sdr_thres_cn0_l = THRES_CN0_L;
double sdr_thres_cn0_u = THRES_CN0_U;

// upper cases of signal string ------------------------------------------------
static void sig_upper(const char *sig, char *Sig)
{
    int i;
    
    for (i = 0; i < 15; i++) {
        if (!(Sig[i] = (char)toupper(sig[i]))) break;
    }
    Sig[i] = '\0';
}

// new signal acquisition ------------------------------------------------------
static sdr_acq_t *acq_new(const int8_t *code, int len_code, double T, double fs,
    int N)
{
    sdr_acq_t *acq = (sdr_acq_t *)sdr_malloc(sizeof(sdr_acq_t));
    
    acq->code_fft = sdr_cpx_malloc(2 * N);
    sdr_gen_code_fft(code, len_code, T, 0.0, fs, N, N, acq->code_fft);
    acq->fd_ext = 0.0;
    acq->fds = sdr_dop_bins(T, 0.0, sdr_max_dop, &acq->len_fds);
    acq->P_sum = NULL;
    acq->n_sum = 0;
    return acq;
}

// free signal acquisition -----------------------------------------------------
static void acq_free(sdr_acq_t *acq)
{
    if (!acq) return;
    sdr_cpx_free(acq->code_fft);
    sdr_free(acq->fds);
    sdr_free(acq->P_sum);
    sdr_free(acq);
}

// new signal tracking ---------------------------------------------------------
static sdr_trk_t *trk_new(const char *sig, int prn, const int8_t *code,
    int len_code, double T, double fs)
{
    sdr_trk_t *trk = (sdr_trk_t *)sdr_malloc(sizeof(sdr_trk_t));
    int i = 0, npos = (SDR_N_CORR - 5) / 2;
    
    int pos = (int)(sdr_sp_corr * T / len_code * fs) + 1;
    trk->pos[i++] = 0;    // P
    trk->pos[i++] = -pos; // E
    trk->pos[i++] = pos;  // L
    trk->pos[i++] = -80;  // N
    for (pos = -npos; pos <= npos; pos++) { // additional correlators
        trk->pos[i++] = pos;
    }
    trk->npos = 4;
    trk->sec_sync = trk->sec_pol = 0;
    trk->err_phas = 0.0;
    trk->sumP = trk->sumE = trk->sumL = trk->sumN = 0.0;
    int N = (int)(fs * T);
    if (!strcmp(sig, "L6D") || !strcmp(sig, "L6E")) {
        trk->code_fft = sdr_cpx_malloc(N * N_CODE);
        for (int i = 0; i < N_CODE; i++) {
            double coff = -i / fs / N_CODE;
            sdr_gen_code_fft(code, len_code, T, coff, fs, N, 0,
                trk->code_fft + i * N);
        }
    }
    else {
        trk->code = (sdr_cpx16_t *)sdr_malloc(sizeof(sdr_cpx16_t) * N * N_CODE);
        for (int i = 0; i < N_CODE; i++) {
            double coff = -i / fs / N_CODE;
            sdr_res_code(code, len_code, T, coff, fs, N, 0, trk->code + i * N);
        }
    }
    return trk;
}

// free signal tracking --------------------------------------------------------
static void trk_free(sdr_trk_t *trk)
{
    if (!trk) return;
    sdr_free(trk->code);
    sdr_cpx_free(trk->code_fft);
    sdr_free(trk);
}

//------------------------------------------------------------------------------
//  Generate new receiver channel.
//
//  args:
//      sig      (I)  Signal ID as string ('L1CA', 'L1CB', 'L1CP', ....)
//      prn      (I)  PRN number
//      fs       (I)  Sampling frequency (Hz)
//      fi       (I)  IF carrier frequency (Hz)
//
//  return:
//      Receiver channel (NULL: error)
//
sdr_ch_t *sdr_ch_new(const char *sig, int prn, double fs, double fi)
{
    sdr_ch_t *ch = (sdr_ch_t *)sdr_malloc(sizeof(sdr_ch_t));
    
    ch->state = SDR_STATE_IDLE;
    ch->time = 0.0;
    sig_upper(sig, ch->sig);
    ch->prn = prn;
    sdr_sat_id(ch->sig, prn, ch->sat);
    if (!(ch->code = sdr_gen_code(sig, prn, &ch->len_code)) ||
        !(ch->sec_code = sdr_sec_code(sig, prn, &ch->len_sec_code))) {
        sdr_free(ch);
        return NULL;
    }
    ch->fc = sdr_sig_freq(sig);
    ch->fs = fs;
    ch->fi = sdr_shift_freq(sig, prn, fi);
    ch->T = sdr_code_cyc(sig);
    ch->N = (int)(fs * ch->T);
    ch->fd = ch->coff = ch->adr = ch->cn0 = 0.0;
    ch->lock = ch->lost = 0;
    ch->costas = strcmp(ch->sig, "L6D") && strcmp(ch->sig, "L6E");
    ch->acq = acq_new(ch->code, ch->len_code, ch->T, fs, ch->N);
    ch->trk = trk_new(ch->sig, ch->prn, ch->code, ch->len_code, ch->T, fs);
    ch->nav = sdr_nav_new();
    pthread_mutex_init(&ch->mtx, NULL);
    return ch;
}

//------------------------------------------------------------------------------
//  Free receiver channel.
//
//  args:
//      ch       (I) Receiver channel
//
//  return:
//      none
//
void sdr_ch_free(sdr_ch_t *ch)
{
    if (!ch) return;
    acq_free(ch->acq);
    trk_free(ch->trk);
    sdr_nav_free(ch->nav);
    sdr_free(ch);
}

// initialize signal tracking --------------------------------------------------
static void trk_init(sdr_trk_t *trk)
{
    trk->err_phas = 0.0;
    trk->sec_sync = trk->sec_pol = 0;
    trk->sumP = trk->sumE = trk->sumL = trk->sumN = 0.0;
    memset(trk->C, 0, sizeof(sdr_cpx_t) * SDR_N_CORR);
    memset(trk->P, 0, sizeof(sdr_cpx_t) * SDR_N_HIST);
}

// start tracking --------------------------------------------------------------
static void start_track(sdr_ch_t *ch, double time, double fd, double coff,
    double cn0)
{
    ch->state = SDR_STATE_LOCK;
    ch->time = time;
    ch->lock = 0;
    ch->fd = fd;
    ch->coff = coff;
    ch->adr = 0.0;
    ch->cn0 = cn0;
    ch->week = 0;
    ch->tow = -1;
    trk_init(ch->trk);
    sdr_nav_init(ch->nav);
}

// search signal ---------------------------------------------------------------
static void search_sig(sdr_ch_t *ch, double time, const sdr_buff_t *buff,
    int ix)
{
    float *fds = ch->acq->fds, fd_ext[3];
    int n = ch->acq->len_fds;
    
    if (ch->acq->fd_ext != 0.0) { // assist by external Doppler
        for (n = 0; n < 3; n++) {
            fd_ext[n] = ch->acq->fd_ext + (n - 1) * 0.5 / ch->T;
        }
        fds = fd_ext;
    }
    if (!ch->acq->P_sum) {
        ch->acq->P_sum = (float *)sdr_malloc(sizeof(float) * 2 * ch->N * n);
    }
    // parallel code search and non-coherent integration
    sdr_search_code(ch->acq->code_fft, ch->T, buff, ix, 2 * ch->N, ch->fs,
        ch->fi, fds, n, ch->acq->P_sum);
    ch->acq->n_sum++;
    
    if (ch->acq->n_sum * ch->T >= sdr_t_acq) {
        int ix[2];
        
        // search max correlation power
        float cn0 = sdr_corr_max(ch->acq->P_sum, 2 * ch->N, ch->N, n, ch->T, ix);
        
        if (cn0 >= sdr_thres_cn0_l) {
            double fd = sdr_fine_dop(ch->acq->P_sum, 2 * ch->N, fds, n, ix);
            double coff = ix[1] / ch->fs;
            start_track(ch, time, fd, coff, cn0);
            // TEB: Log level 4 -> 3
            sdr_log(3, "$LOG,%.3f,%s,%d,SIGNAL FOUND (%.1f,%.1f,%.7f)", time,
                ch->sig, ch->prn, cn0, fd, coff * 1e3);
        }
        else {
            sdr_sleep_msec(100);
            ch->state = SDR_STATE_IDLE;
            sdr_log(4, "$LOG,%.3f,%s,%d,SIGNAL NOT FOUND (%.1f)", time, ch->sig,
                ch->prn, cn0);
        }
        sdr_free(ch->acq->P_sum);
        ch->acq->P_sum = NULL;
        ch->acq->n_sum = 0;
    }
}

// sync and remove secondary code ----------------------------------------------
static void sync_sec_code(sdr_ch_t *ch, int N)
{
    if (ch->trk->sec_sync == 0) {
        float P = 0.0, R = 0.0;
        for (int i = 0; i < N; i++) {
            P += ch->trk->P[SDR_N_HIST-N+i][0] * ch->sec_code[i] / N;
            R += fabsf(ch->trk->P[SDR_N_HIST-N+i][0]) / N;
        }
        if (fabsf(P) >= R && R >= THRES_SYNC) {
            ch->trk->sec_sync = ch->lock;
            ch->trk->sec_pol = (P > 0.0f) ? 1 : -1;
        }
    }
    else if ((ch->lock - ch->trk->sec_sync) % N == 0) {
        float P = 0.0;
        for (int i = 0; i < N; i++) {
            P += ch->trk->P[SDR_N_HIST-N+i][0] / N;
        }
        if (fabsf(P) < THRES_LOST) {
            ch->trk->sec_sync = ch->trk->sec_pol = 0;
        }
    }
    if (ch->trk->sec_sync > 0) {
        int8_t C = ch->sec_code[(ch->lock - ch->trk->sec_sync - 1) % N] *
            ch->trk->sec_pol;
        for (int i = 0; i < ch->trk->npos; i++) {
            ch->trk->C[i][0] *= C;
            ch->trk->C[i][1] *= C;
        }
        ch->trk->P[SDR_N_HIST-1][0] *= C;
        ch->trk->P[SDR_N_HIST-1][1] *= C;
    }
}

// FLL ------------------------------------------------------------------------
static void FLL(sdr_ch_t *ch)
{
    if (ch->lock >= 2) {
        double IP1 = ch->trk->P[SDR_N_HIST-1][0];
        double QP1 = ch->trk->P[SDR_N_HIST-1][1];
        double IP2 = ch->trk->P[SDR_N_HIST-2][0];
        double QP2 = ch->trk->P[SDR_N_HIST-2][1];
        double dot   = IP1 * IP2 + QP1 * QP2;
        double cross = IP1 * QP2 - QP1 * IP2;
        if (dot != 0.0) {
            double B = ch->lock * ch->T < T_FPULLIN ? sdr_b_fll_w : sdr_b_fll_n;
            double err_freq = ch->costas ? atan(cross / dot) : atan2(cross, dot);
            ch->fd -= B / 0.25 * err_freq / DPI;
        }
    }
}

// PLL -------------------------------------------------------------------------
static void PLL(sdr_ch_t *ch)
{
    double IP = ch->trk->C[0][0];
    double QP = ch->trk->C[0][1];
    if (IP != 0.0) {
        double err_phas = (ch->costas ? atan(QP / IP) : atan2(QP, IP)) / DPI;
        double W = sdr_b_pll / 0.53;
        ch->fd += 1.4 * W * (err_phas - ch->trk->err_phas) +
            W * W * err_phas * ch->T;
        ch->trk->err_phas = err_phas;
    }
}

// DLL -------------------------------------------------------------------------
static void DLL(sdr_ch_t *ch)
{
    int N = MAX(1, (int)(sdr_t_dll / ch->T));
    ch->trk->sumE += sdr_cpx_abs(ch->trk->C[1]); // non-coherent sum 
    ch->trk->sumL += sdr_cpx_abs(ch->trk->C[2]);
    if (ch->lock % N == 0) {
        double E = ch->trk->sumE;
        double L = ch->trk->sumL;
        if (E + L > 0.0) {
            double err_code = (E - L) / (E + L) * 0.5f * ch->T / ch->len_code;
            ch->coff -= sdr_b_dll / 0.25 * err_code * ch->T * N;
            ch->trk->err_code = err_code;
        }
        ch->trk->sumE = ch->trk->sumL = 0.0;
    }
}

// update C/N0 -----------------------------------------------------------------
static void CN0(sdr_ch_t *ch)
{
    ch->trk->sumP += SQR(ch->trk->C[0][0]) + SQR(ch->trk->C[0][1]);
    ch->trk->sumN += SQR(ch->trk->C[3][0]) + SQR(ch->trk->C[3][1]);
    if (ch->lock % (int)(T_CN0 / ch->T) == 0) {
        if (ch->trk->sumN > 0.0) {
            double cn0 = 10.0 * log10(ch->trk->sumP / ch->trk->sumN / ch->T);
            ch->cn0 += 0.5 * (cn0 - ch->cn0);
        }
        ch->trk->sumP = ch->trk->sumN = 0.0;
    }
}

// interpolate correlation -----------------------------------------------------
static void interp_corr(const sdr_cpx_t *C, double x, sdr_cpx_t *c)
{
    int i = (int)x;
    double a1 = x - i, a0 = 1.0 - a1;
    (*c)[0] = a0 * C[i][0] + a1 * C[i+1][0];
    (*c)[1] = a0 * C[i][1] + a1 * C[i+1][1];
}

// decode L6 CSK ---------------------------------------------------------------
static void CSK(sdr_ch_t *ch, const sdr_cpx_t *corr)
{
    double R = (double)ch->N / (ch->len_code / 2); // samples / chips 
    int n = (int)(280 * R);
    sdr_cpx_t *C = sdr_cpx_malloc(2 * n);
    memcpy(C, corr + ch->N - n, sizeof(sdr_cpx_t) * n);
    memcpy(C + n, corr, sizeof(sdr_cpx_t) * n);
    
    // interpolate correlation powers and detect peak
    double P_max = 0.0;
    int ix = 0;
    for (int i = -255; i <= 255; i++) {
        sdr_cpx_t c;
        interp_corr(C, n + i * R, &c);
        double P = sdr_cpx_abs(c);
        if (P <= P_max) continue;
        P_max = P;
        ix = i;
    }
    // add CSK symbol to buffer 
    uint8_t sym = (uint8_t)(255 - ix % 256);
    sdr_add_buff(ch->nav->syms, SDR_MAX_NSYM, &sym, sizeof(sym));
    
    // generate correlator outputs
    for (int i = 0; i < ch->trk->npos; i++) {
        interp_corr(C, n + ix * R + ch->trk->pos[i], ch->trk->C + i);
    }
    sdr_cpx_free(C);
}

// update TOW ------------------------------------------------------------------
static void update_tow(sdr_ch_t *ch, double sec)
{
    if (ch->tow < 0) return;
    ch->tow = (ch->tow + (int)(sec / 1e-3)) % (86400 * 7 * 1000);
}

// track signal ----------------------------------------------------------------
static void track_sig(sdr_ch_t *ch, double time, const sdr_buff_t *buff, int ix)
{
    double tau = time - ch->time;   // time interval (s) 
    double fc = ch->fi + ch->fd;    // IF carrier frequency with Doppler (Hz) 
    ch->adr += ch->fd * tau;        // accumulated Doppler (cyc)
    ch->coff -= ch->fd / ch->fc * tau; // carrier-aided code offset (s) 
    ch->time = time;
    
    // adjust code offset within 1 cycle range (0 <= ch->coff < ch->T)
    if (ch->coff >= ch->T) {
        ch->coff -= ch->T;
        update_tow(ch, -ch->T);
        ch->lock--;
        memmove(ch->trk->P + 1, ch->trk->P, sizeof(sdr_cpx_t) *
            (SDR_N_HIST - 1));
    }
    else if (ch->coff < 0.0) {
        ch->coff += ch->T;
        update_tow(ch, ch->T);
        ch->lock++;
        memmove(ch->trk->P, ch->trk->P + 1, sizeof(sdr_cpx_t) *
            (SDR_N_HIST - 1));
    }
    // code position (samples) and carrier phase (cyc) 
    int i = (int)(ch->coff * ch->fs);
    int j = (int)((ch->coff * ch->fs - i) * N_CODE); // code bank index
    double phi = ch->fi * tau + ch->adr + fc * i / ch->fs;
    
    if (!strcmp(ch->sig, "L6D") || !strcmp(ch->sig, "L6E")) {
        sdr_cpx_t *corr = sdr_cpx_malloc(ch->N);
        
        // FFT correlator 
        sdr_corr_fft(buff, ix + i, ch->N, ch->fs, fc, phi,
            ch->trk->code_fft + j * ch->N, corr);
        
        // decode L6 CSK 
        CSK(ch, corr);
        
        sdr_cpx_free(corr);
    }
    else {
        // standard correlator
        sdr_corr_std(buff, ix + i, ch->N, ch->fs, fc, phi,
            ch->trk->code + j * ch->N, ch->trk->pos, ch->trk->npos, ch->trk->C);
    }
    // add P correlator outputs to history 
    sdr_add_buff(ch->trk->P, SDR_N_HIST, ch->trk->C[0], sizeof(sdr_cpx_t));
    update_tow(ch, ch->T);
    ch->lock++;
    
    // sync and remove secondary code 
    if (ch->len_sec_code >= 2 && ch->lock * ch->T >= T_NPULLIN) {
        sync_sec_code(ch, ch->len_sec_code);
    }
    // FLL/PLL, DLL and update C/N0 
    if (ch->lock * ch->T <= T_FPULLIN) {
        FLL(ch);
    }
    else {
        PLL(ch);
    }
    DLL(ch);
    CN0(ch);
    
    // decode navigation data 
    if (ch->lock * ch->T >= T_NPULLIN) {
        sdr_nav_decode(ch);
    }
    if (ch->cn0 < sdr_thres_cn0_u) { // signal lost 
        ch->state = SDR_STATE_IDLE;
        ch->lock = 0;
        ch->trk->sec_sync = ch->trk->sec_pol = 0;
        ch->nav->ssync = ch->nav->fsync = ch->nav->rev = 0;
        ch->lost++;
        sdr_log(4, "$LOG,%.3f,%s,%d,SIGNAL LOST (%s, %.1f)", ch->time, ch->sig,
            ch->prn, ch->sig, ch->cn0);
    }
}

//------------------------------------------------------------------------------
//  Update a receiver channel. A receiver channel is a state machine which has
//  the following internal states indicated as ch.state. By calling the function,
//  the receiver channel search and track GNSS signals and decode navigation
//  data in the signals. The results of the signal acquisition, tracking and
//  navigation data decoding are output as log messages. The internal status are
//  also accessed as object instance variables of the receiver channel after
//  calling the function. The function should be called in the cycle of GNSS
//  signal code with 2-cycle samples of digitized IF data (which are overlapped
//  between previous and current). 
//
//    SDR_STATE_SRCH : signal acquisition state
//    SDR_STATE_LOCK : signal tracking state
//    SDR_STATE_IDLE : waiting for a next signal acquisition cycle
//
//  args:
//      ch       (I)  Receiver channel
//      time     (I)  Sampling time of the end of digitized IF data (s)
//      buff     (I)  IF data buffer
//      ix       (I)  index of IF data buffer
//
//  return:
//      none
//
void sdr_ch_update(sdr_ch_t *ch, double time, const sdr_buff_t *buff, int ix)
{
    if (ch->state == SDR_STATE_SRCH) {
        search_sig(ch, time, buff, ix);
    }
    else if (ch->state == SDR_STATE_LOCK) {
        pthread_mutex_lock(&ch->mtx);
        track_sig(ch, time, buff, ix);
        pthread_mutex_unlock(&ch->mtx);
    }
}

// set receiver channel correlator ---------------------------------------------
void sdr_ch_set_corr(sdr_ch_t *ch, int npos)
{
    pthread_mutex_lock(&ch->mtx);
    ch->trk->npos = npos;
    pthread_mutex_unlock(&ch->mtx);
}

// get receiver correlator status ----------------------------------------------
int sdr_ch_corr_stat(sdr_ch_t *ch, double *stat, int *pos, sdr_cpx_t *C)
{
    pthread_mutex_lock(&ch->mtx);
    int n = SDR_N_CORR;
    stat[0] = ch->state;
    stat[1] = ch->fs;
    stat[2] = ch->lock * ch->T;
    stat[3] = ch->cn0;
    stat[4] = ch->coff * 1e3;
    stat[5] = ch->fd;
    memcpy(pos, ch->trk->pos, sizeof(int) * n);
    memcpy(C, ch->trk->C, sizeof(sdr_cpx_t) * n);
    pthread_mutex_unlock(&ch->mtx);
    return n;
}

// get receiver correlator history ---------------------------------------------
int sdr_ch_corr_hist(sdr_ch_t *ch, double tspan, double *stat, sdr_cpx_t *P)
{
    pthread_mutex_lock(&ch->mtx);
    int n = MIN((int)(tspan / ch->T), SDR_N_HIST);
    stat[0] = ch->time;
    stat[1] = ch->T;
    memcpy(P, ch->trk->P + SDR_N_HIST - n, sizeof(sdr_cpx_t) * n);
    pthread_mutex_unlock(&ch->mtx);
    return n;
}
