#include "pocket_sdr.h"

#define MAX_CHAR (512)
#define MAX_SAT (12)
#define MAX_OBS (120)

// Defined in rtklib.h
//#define PI 3.1415926535898
//#define R2D 57.2957795131

#define POW2_M19 1.907348632812500e-6
#define POW2_M31 4.656612873077393e-10
#define POW2_M32 2.328306436538696e-10
#define POW2_M43 1.136868377216160e-13

#define GM_MOON 4.9028e12
#define R_MOON 1737.4e3

#define SECONDS_IN_WEEK 604800.0
#define SECONDS_IN_HALF_WEEK 302400.0

#define SPEED_OF_LIGHT 2.99792458e8

#define TRUE  1
#define FALSE 0
#define SINGULAR    0 
#define NONSINGULAR 1
#define MAXCHANNELS 12

typedef struct
{
    int week;
    double sec;
} gpstime_t;

typedef struct
{
    int vflg;

    gpstime_t toe;
    gpstime_t toc;

    int prn;

    double ecc;   // Eccentricity
    double sqrta; // SQRT(A) 
    double m0;    // Mean Anom
    double omg0;  // Right Ascen at Week
    double inc0;  // Orbital Inclination
    double aop;   // Argument of Perigee

    double af0;
    double af1;

    // Working variables
    double n; // Mean motion
    double A; // Semi-major axis
    double sq1e2; // sqrt(1-e^2)

} ephem_t;

typedef struct
{
    gpstime_t g;
    int nsat;
    int prn[MAX_SAT];
    double range[MAX_SAT];
} obsrv_t;

// Global variables
static double gdop, pdop;
static ephem_t eph[MAX_SAT] = { 0 };

void xyz2llh(const double *xyz, double *llh)
{
	double a, x, y, z, rho2;

	a = R_MOON;

	x = xyz[0];
	y = xyz[1];
	z = xyz[2];

	rho2 = x * x + y * y;

	llh[0] = atan2(z, sqrt(rho2));
	llh[1] = atan2(y, x);
	llh[2] = sqrt(rho2 + z*z) - a;

    return;
}
static void satpos(ephem_t eph, gpstime_t g, double* pos, double* vel, double* clk)
{
    double tk;
    double mk;
    double ek;
    double ekold;
    double OneMinusecosE;
    double cek, sek;
    double ekdot;
    double uk;
    double cuk, suk;
    double ukdot;
    double rk;
    double rkdot;
    double ik;
    double cik, sik;
    double xpk, ypk;
    double xpkdot, ypkdot;
    double ok;
    double cok, sok;

    tk = g.sec - eph.toe.sec;

    if (tk > SECONDS_IN_HALF_WEEK)
        tk -= SECONDS_IN_WEEK;
    else if (tk < -SECONDS_IN_HALF_WEEK)
        tk += SECONDS_IN_WEEK;

    // Mean anomaly
    mk = eph.m0 + eph.n * tk;

    // Eccentric anomaly
    ek = mk;
    ekold = ek + 1.0;

    OneMinusecosE = 1.0;

    while (fabs(ek - ekold) > 1.0E-14)
    {
        ekold = ek;
        OneMinusecosE = 1.0 - eph.ecc * cos(ekold);
        ek = ek + (mk - ekold + eph.ecc * sin(ekold)) / OneMinusecosE;
    }

    sek = sin(ek);
    cek = cos(ek);

    ekdot = eph.n / OneMinusecosE;

    // True anomaly + Argument of perigee
    uk = atan2(eph.sq1e2 * sek, cek - eph.ecc) + eph.aop;
    suk = sin(uk);
    cuk = cos(uk);
    ukdot = eph.sq1e2 * ekdot / OneMinusecosE;

    // Range and range rate
    rk = eph.A * OneMinusecosE;
    rkdot = eph.A * eph.ecc * sek * ekdot;

    xpk = rk * cuk;
    ypk = rk * suk;
    xpkdot = rkdot * cuk - ypk * ukdot;
    ypkdot = rkdot * suk + xpk * ukdot;

    // Inclination
    ik = eph.inc0;

    sik = sin(ik);
    cik = cos(ik);

    // RAAN
    ok = eph.omg0;
    sok = sin(ok);
    cok = cos(ok);

    // Moon-centered inertial coordinates
    pos[0] = xpk * cok - ypk * cik * sok;
    pos[1] = xpk * sok + ypk * cik * cok;
    pos[2] = ypk * sik;

    vel[0] = xpkdot * cok - ypkdot * cik * sok;
    vel[1] = xpkdot * sok + ypkdot * cik * cok;
    vel[2] = ypkdot * sik;

    // Satellite clock correction
    tk = g.sec - eph.toc.sec;

    if (tk > SECONDS_IN_HALF_WEEK)
        tk -= SECONDS_IN_WEEK;
    else if (tk < -SECONDS_IN_HALF_WEEK)
        tk += SECONDS_IN_WEEK;

    clk[0] = eph.af0 + tk * eph.af1;
    clk[1] = eph.af1;
}

static void subVect(double* y, double* x1, double* x2)
{
    y[0] = x1[0] - x2[0];
    y[1] = x1[1] - x2[1];
    y[2] = x1[2] - x2[2];

    return;
}

static double normVect(double* x)
{
    return(sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]));
}

static double computeRange(ephem_t eph, gpstime_t t_rcv, double xyz[], double dt_rcv, double drho_dr[])
{
    gpstime_t g;
    double pos[3], vel[3], clk[2];
    double los[3];
    double tau;
    double r, rho;

    g.week = t_rcv.week;
    g.sec = t_rcv.sec - dt_rcv;
    if (g.sec < 0.0) {
        g.sec += SECONDS_IN_WEEK;
        g.week--;
    }

    // SV position at time of the pseudorange observation.
    satpos(eph, g, pos, vel, clk);

    // Receiver to satellite vector and light-time.
    subVect(los, pos, xyz);
    tau = normVect(los) / SPEED_OF_LIGHT;

    // Extrapolate the satellite position backwards to the transmission time.
    pos[0] -= vel[0] * tau;
    pos[1] -= vel[1] * tau;
    pos[2] -= vel[2] * tau;

    // New observer to satellite vector and satellite range.
    subVect(los, pos, xyz);
    r = normVect(los);
    rho = r + dt_rcv * SPEED_OF_LIGHT; // - dt_sat * SPEED_OF_LIGHT

    // Partials
    drho_dr[0] = -los[0] / r;
    drho_dr[1] = -los[1] / r;
    drho_dr[2] = -los[2] / r;

    return(rho);
}

#define MAX_M 4
#define MAX_N 12

static int invMat(double a[MAX_M][MAX_M], int m)
{
    int i,j,k;
    double b[MAX_M][MAX_M+MAX_M];

    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            b[i][j] = a[i][j];
            if (i == j) {
                b[i][j + m] = 1.0;
            } else {
                b[i][j + m] = 0.0;
            }
        }
    }
    
    // Gaussian elimination
    for (i = 0; i < m; i++) {
        if (fabs(b[i][i]) <= 1.0E-10) {
            return(SINGULAR);
        }
        for (j = m + m - 1; j >= i; j--) {
            b[i][j] /= b[i][i];
        }

        for (k = 0; k < m; k++) {
            if (k != i) {
                for (j = m + m -1; j >= i; j--) {
                    b[k][j] -= b[k][i] * b[i][j];
                }
            }
        }
    }

    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            a[i][j] = b[i][j + m];
        }
    }

    return(NONSINGULAR);
}

static int invCovMat(double sm[MAX_N][MAX_M], int n, double icm[MAX_M][MAX_M])
{
    int i,j,k;
    int m = MAX_M;

    // First determine the product sm'*sm.
    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            icm[i][j] = 0.0;
            for (k = 0; k < n; k++) {
                icm[i][j] += sm[k][i] * sm[k][j];
            }
        }
    }

    // Now invert to produce the inverse covariance matrix.
    if (invMat(icm, m) == SINGULAR)
    {
        gdop = 0.0;
        pdop = 0.0;
        return(FALSE);
    }
    else
    {
        gdop = sqrt(icm[0][0] + icm[1][1] + icm[2][2] + icm[3][3]);
        pdop = sqrt(icm[0][0] + icm[1][1] + icm[2][2]);
        return(TRUE);
    }
}

static void out_log_pos(double time, const sol_t *sol, int nsat) 
{
    double ep[6], pos[3];

    time2epoch(timeadd(sol->time, sol->dtr[0]), ep);
    xyz2llh(sol->rr, pos);

    sdr_log(3, "$POS,%.3f,%.0f,%.0f,%.0f,%.0f,%.0f,%.3f,%.9f,%.9f,%.3f,%d,%d",
        time, ep[0], ep[1], ep[2], ep[3], ep[4], ep[5],
        pos[0] * R2D, pos[1] * R2D, pos[2], 5, nsat);
}

void update_sol_afs(sdr_pvt_t *pvt)
{
    double time = pvt->ix * SDR_CYC;
    
    int i, j, k, itr;
    obsrv_t obs;
    int week;
    double tow;

    int prn, sv, nsat;
    double rho_c;
    double drho_dr[3];

    int m = MAX_M;
    double omp[MAX_N];
    double nm[MAX_M][MAX_N];
    double sm[MAX_N][MAX_M];
    double icm[MAX_M][MAX_M];

    double r[3], dt_rcv = 0.0;
    double position;
    double position_update[4];

    sol_t *sol = pvt->sol;

    //sdr_log(3, "$DBG,%.3f,UPDATE AFS SOLUTION,obs->n=%d", time, pvt->obs->n);

    sol->stat = SOLQ_NONE;
    sol->time = pvt->obs->data[0].time;

    pvt->nsat = pvt->obs->n;

    // Read observations
    nsat = 0;
    obs.nsat = pvt->obs->n;
    tow = time2gpst(pvt->obs->data[0].time, &week);
    obs.g.week = week;
    obs.g.sec = tow;

    for (i = 0; i < obs.nsat; i++) {
        const obsd_t *data = pvt->obs->data + i;

        satsys(data->sat, &prn);

        if (eph[prn - 1].vflg != 1) continue;

        for (int j = 0; j < NFREQ + NEXOBS; j++) {
            if (!data->code[j]) continue;

            obs.prn[nsat] = prn;
            obs.range[nsat] = data->P[j];
            break;
        }

        nsat++; // Number of satellites with valid ephemeris and observation

        //sdr_log(3, "$DBG,%.3f,%d,%d,%.3f,%02d,%.3f", time, nsat, obs.g.week, obs.g.sec, obs.prn[i], obs.range[i]);
    }

    if (nsat < 4) {
        //sdr_log(3, "$DBG,%.3f,nsat=%d", time, nsat);
        return;
    }

    // Initial position
    r[0] = 0.0; r[1] = 0.0; r[2] = -R_MOON; // South pole
    dt_rcv = 0.0;

    for (itr = 0; itr < 4; itr++) {
        for (i = 0; i < nsat; i++) {

            sv = obs.prn[i] - 1;
            rho_c = computeRange(eph[sv], obs.g, r, dt_rcv, drho_dr);

            omp[i] = obs.range[i] - rho_c; // observed minus predicted range

            sm[i][0] = drho_dr[0];
            sm[i][1] = drho_dr[1];
            sm[i][2] = drho_dr[2];
            sm[i][3] = 1.0;
        }

        if (invCovMat(sm, nsat, icm) == TRUE) {

            for (i = 0; i < m; i++) {
                for (j = 0; j < nsat; j++) {
                    nm[i][j] = 0.0;
                    for (k = 0; k < m; k++) {
                        nm[i][j] += icm[i][k] * sm[j][k];
                    }
                }
            }

            for (i = 0; i < m; i++) {
                position = 0.0;
                for (k = 0; k < nsat; k++) {
                    position += nm[i][k] * omp[k];
                }
                position_update[i] = position;
            }

            r[0] += position_update[0];
            r[1] += position_update[1];
            r[2] += position_update[2];
            dt_rcv += position_update[3] / SPEED_OF_LIGHT;
        }
        else {
            sdr_log(3, "$PVT,%.3f,ERROR: Solution matrix is singular.", time);
            sol->ns = 0;
            return;
        }
    }

    // Save PVT solutions
    sol->type = 0;
    sol->time = timeadd(pvt->obs->data[0].time, -dt_rcv);
    sol->dtr[0] = dt_rcv;
    for (j=0;j<6;j++) sol->rr[j] = j<3?r[j]:0.0;
    sol->age = sol->ratio = 0.0;
    sol->stat = SOLQ_SINGLE + MAXSOLQ; // for LANS AFS navigation solutions
    sol->ns = nsat;

    out_log_pos(time, pvt->sol, nsat);
}

int decode_afs_frame(const uint8_t *buff, int prn)
{
    int sv = prn - 1;
    int i = 0;
    int wn;

    if (sv >= MAX_SAT) return 0;

    wn = getbitu(buff, i, 13); i += 22;

    eph[sv].toe.week = wn;
    eph[sv].toe.sec  = getbitu(buff, i, 16) * 16.0;          i += 16;
    eph[sv].ecc      = getbitu(buff, i, 32) * POW2_M32;      i += 32;
    eph[sv].sqrta    = getbitu(buff, i, 32) * POW2_M19;      i += 32;
    eph[sv].inc0     = getbits(buff, i, 32) * POW2_M31 * PI; i += 32;
    eph[sv].omg0     = getbits(buff, i, 32) * POW2_M31 * PI; i += 32;
    eph[sv].aop      = getbits(buff, i, 32) * POW2_M31 * PI; i += 32;
    eph[sv].m0       = getbits(buff, i, 32) * POW2_M31 * PI; i += 32;

    eph[sv].toc.week = wn;
    eph[sv].toc.sec  = getbitu(buff, i, 16) * 16.0;          i += 16;
    eph[sv].af0      = getbits(buff, i, 22) * POW2_M31;      i += 22;
    eph[sv].af1      = getbits(buff, i, 16) * POW2_M43;

    eph[sv].prn = prn; // not used

    // Update the working variables
    eph[sv].A = eph[sv].sqrta * eph[sv].sqrta;
    eph[sv].n = sqrt(GM_MOON / (eph[sv].A * eph[sv].A * eph[sv].A));
    eph[sv].sq1e2 = sqrt(1.0 - eph[sv].ecc * eph[sv].ecc);

    // Valid ephemeris
    eph[sv].vflg = 1;

    return 1;
}