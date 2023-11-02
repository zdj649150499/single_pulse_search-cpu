#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "fftw3.h"

#include "errno.h"
#include "limits.h"

#include "gsl/gsl_statistics.h"
#include "gsl/gsl_fit.h"
#include "gsl/gsl_multifit.h"





#define NUMSCOPES  8
#define NUMBANDS   6
#define MAXNUMONOFF 40


#ifdef __cplusplus
extern "C" {
#endif

    void help();
    void help_H();

    // int useffts;
    // int dosearch;

    // int fftlen;   // Should be a power-of-two for best speed
    // int chunklen;  // Must be at least max_downfact less than fftlen

    // float opt_maxwidth;
    // float opt_threshold;
    // float opt_start;
    // float opt_end;
    // char opt_glob[128];
    // int opt_fast;
    // int opt_nobadblocks;
    // int opt_detrendlen;

    // int detrendlen;
    // int blocks_per_chunk;
    // int overlap;
    // int worklen;  // currently it is fftlen...
    // int max_downfact;
    // int default_downfacts[14];

    long long next2_to_n(long long x);
    void cutExtName(char *name, char *marker);

    typedef struct INFODATA {
    double ra_s;		/* Right ascension seconds (J2000)       */
    double dec_s;		/* Declination seconds (J2000)           */ 
    double N;		        /* Number of bins in the time series     */
    double dt;	 	        /* Width of each time series bin (sec)   */
    double fov;			/* Diameter of Beam or FOV in arcsec     */
    double mjd_f;		/* Epoch of observation (MJD) frac part  */
    double dm;			/* Radio -- Dispersion Measure (cm-3 pc) */
    double freq;		/* Radio -- Low chan central freq (Mhz)  */
    double freqband;		/* Radio -- Total Bandwidth (Mhz)        */
    double chan_wid;		/* Radio -- Channel Bandwidth (Mhz)      */
    double wavelen;		/* IR,Opt,UV -- central wavelength (nm)  */
    double waveband;		/* IR,Opt,UV -- bandpass (nm)            */
    double energy;		/* x-ray,gamma -- central energy (kev)   */
    double energyband;		/* x-ray,gamma -- energy bandpass (kev)  */
    double onoff[MAXNUMONOFF*2];/* Bin number pairs where obs is "on"    */
    int num_chan;		/* Radio -- Number Channels              */
    int mjd_i;			/* Epoch of observation (MJD) int part   */
    int ra_h;			/* Right ascension hours (J2000)         */
    int ra_m;			/* Right ascension minutes (J2000)       */
    int dec_d;			/* Declination degrees (J2000)           */
    int dec_m;			/* Declination minutes (J2000)           */  
    int bary;			/* Barycentered?  1=yes, 0=no            */
    int numonoff;		/* The number of onoff pairs in the data */ 
    int breaks;         /*  Any breaks in the data? (1 yes, 0 no) */
    char notes[500];		/* Any additional notes                  */
    char name[200];		/* Data file name without suffix         */
    char object[100];		/* Object being observed                 */ 
    char instrument[100];	/* Instrument used                       */
    char observer[100];		/* Observer[s] for the data set          */
    char analyzer[100];		/* Who analyzed the data                 */
    char telescope[40];		/* Telescope used                        */
    char band[40];		/* Type of observation (EM band)         */
    char filt[7];		/* IR,Opt,UV -- Photometric Filter       */
} infodata;

typedef float rawtype_part;

typedef struct FCOMPLEX {
  rawtype_part r, i;
} fcomplex;


typedef struct {
    double DM;
    double time;
    float val;
    int bin;
    int downfact;
} Candidate;

typedef struct {
    Candidate *candidates;
    int count;
} CandidateList;


    FILE *chkfopen(char *path, const char *mode);
    double chk_str2double(char *instr, char *desc);
    long chk_str2long(char *instr, char *desc);
    char *rmtrail(char *str);
    char *rmlead(char *str);
    char *remove_whitespace(char *str);
    void read_inf_line_valstr(FILE * infofile, char *valstr, char *errdesc);
    void ra_dec_from_string(char *radec, int *h_or_d, int *m, double *s);
    void readinf(infodata * data, char *filenm);

    fcomplex *gen_cvect(long length);
    void vect_free(void *vect);
    void realfft(float idata[], long n, int isign);
    void fftwcall(fcomplex * indata, long nn, int isign);
    void tablesixstepfft(fcomplex * indata, long nn, int isign);
    void read_wisdom(void);
    long long good_factor(long long nn);
    // static int TOMS_gcd(int a, int b);
    short transpose_fcomplex(fcomplex * a, int nx, int ny, unsigned char *move,
                         int move_size);


    float median(float arr[], int n);

    void swap(float* a, float* b);
    void swap_d(double* a, double* b) ;
    int partition(float arr[], int low, int high);
    int partition_d(double arr[], int low, int high) ;
    void quickSort(float arr[], int low, int high);
    void quickSort_d(double arr[], int low, int high) ;
    void quickSortWrapper(float arr[], int size);
    void quickSortWrapper_d(double arr[], int size) ;

    int dofit(const gsl_multifit_robust_type *T,
      const gsl_matrix *X, const gsl_vector *y,
      gsl_vector *c, gsl_matrix *cov);

    void fit_line_robust(float *indata, gsl_vector *x, gsl_vector *y, gsl_matrix *X, gsl_matrix *cov, gsl_vector *c, int N);
    void fit_line_samp(double *indata, double *xarray, double *c0, double *c1, int N);


    void fft_convolve(float *fftd_data, float *fftd_kern, float *prod, int N);


    void prune_related1(int *hibins, float *hivals, int N, int downfact, int *off);
    void prune_related2(Candidate *dm_candlist, int dm_candlist_count, int downfacts[], int N, int *off);

    int mycmp(float a, float b) ;
    int cmp(const void *a, const void *b);

    void addCandidate(CandidateList *candlist, float DM, float val, double time, int bin, int valid);

    void prune_border_cases(Candidate *dm_candlist, int dm_candlist_count, long offregions[], int *off);


#ifdef __cplusplus
}
#endif



#endif