#include "functions.h"

#define strMove(d,s) memmove(d,s,strlen(s)+1)
#define _FILE_OFFSET_BITS 64
#define TWOPI         6.2831853071795864769252867665590057683943387987502
#define COMPLEXFFT    fftwcall
#define BIGFFTWSIZE           200000000
#define is_aligned(POINTER, BYTE_COUNT) \
    ((uintptr_t)(const void *)(POINTER)) % (BYTE_COUNT)
#define ALIGNSIZE 64

char bands[NUMBANDS][40] = { "Radio", "IR", "Optical", "UV", "X-ray", "Gamma" };

char scopes[NUMSCOPES][40] =
    { "None (Artificial Data Set)", "Arecibo", "Parkes", "VLA",
    "MMT", "Las Campanas 2.5m", "Mt. Hopkins 48in", "Other"
};


void help()
{
    printf("    Usage: single_pulse_search.py [options] .dat files\n");
    printf("Options:");
    printf("    -h, --help              show this help message and exit\n");
    // printf("    -x, --xwin              Don't make a postscript plot, just use an X-window\n");
    // printf("    -p, --noplot            Look for pulses but do not generate a plot\n");
    printf("    -m  MAXWIDTH, --maxwidth=MAXWIDTH [float]\n");
    printf("                            Set the max downsampling in sec (see below for default)\n");
    printf("    -t THRESHOLD, --threshold=THRESHOLD, [float]\n");
    printf("                            Set a different threshold SNR (default=5.0)\n");
    // printf("    -s T_START, --start=T_START [float]\n");
    // printf("                            Only plot events occuring after this time (s)\n");
    // printf("    -e T_END, --end=T_END [float]\n");
    // printf("                            Only plot events occuring before this time (s)\n");
    // printf("    -g GLOBEXP, --glob=GLOBEXP [string]\n");
    // printf("                            Process the files from this glob expression\n");
    printf("    -f, --fast [0/1]        Use a faster method of de-trending (2x speedup)\n");
    printf("    -b, --nobadblocks [0/1] Don't check for bad-blocks (may save strong pulses)");
    printf("    -d DETRENDFACT, --detrendlen=DETRENDFACT [int]\n");
    printf("                            Chunksize for detrending (pow-of-2 in 1000s)\n");
    printf("\n");
}
void help_H()
{
    printf("***********************************************************************\n");
    printf("usage:  single_pulse_search.py [options] .dat files\n");
    printf("[-h, --help]        : Display this help\n");
    printf("[-m, --maxwidth]    : Set the max downsampling in sec (see below for default)\n");
    // printf("[-p, --noplot]      : Look for pulses but do not generate a plot\n");
    printf("[-t, --threshold]   : Set a different threshold SNR (default=5.0)\n");
    // printf("[-x, --xwin]        : Don't make a postscript plot, just use an X-window\n");
    // printf("[-s, --start]       : Only plot events occuring after this time (s)\n");
    // printf("[-e, --end]         : Only plot events occuring before this time (s)\n");
    // printf("[-g, --glob]        : Use the files from these glob expressions (in quotes)\n");
    printf("[-f, --fast]        : Use a less-accurate but much faster method of detrending\n");
    printf("[-b, --nobadblocks] : Don't check for bad-blocks (may save strong pulses)\n");
    printf("[-d, --detrendlen]  : Chunksize for detrending (pow-of-2 in 1000s, default=1)\n\n");

    printf("Perform a single-pulse search (or simply re-plot the results of a\n");
    printf("single-pulse search) on a set of de-dispersed time series (.dat\n");
    printf("files).\n\n");

    printf("The search attempts to find pulses by matched-filtering the data with\n");
    printf("a series of different width boxcar functions.  The possible boxcar\n");
    printf("sizes are [1, 2, 3, 4, 6, 9, 14, 20, 30, 45, 70, 100, 150, 220, 300]\n");
    printf("bins.  By default the boxcars <= 30 are used.  You can specify\n");
    printf("that the larger boxcars are used with the -m (or --maxwidth) option.\n\n");

    printf("The matched filtering (and accounting for all the possible 'phase'\n");
    printf("offsets of each boxcar) is accomplished by convolving the boxcars\n");
    printf("with the full resolution data.  'Duplicate' candidates from this\n");
    printf("process are filtered, leaving only the most significant.  The time\n");
    printf("series are initially smoothed (by default) using a piecewise linear\n");
    printf("fit to the data where each piece is 1000 data points long.\n\n");

    printf("If the input files are .singlepulse files, we won't actually perform\n");
    printf("a search, we'll only read in the output .singlepulse files and make\n");
    printf("a plot using the information they contain (along with the\n");
    printf("corresponding .inf files).\n\n");

    printf("Notes on usage and performance:\n\n");

    printf("-- single_pulse_search.py is tuned for finding *narrow* pulses\n");
    printf("   (i.e. those of only a few bins width).  Because of this, you\n");
    printf("   should always search appropriately downsampled data (as\n");
    printf("   recommended by DDplan.py, for instance) where dispersion\n");
    printf("   smearing is <~ 1 time series bin.\n\n");

    printf("-- the linear-piecewise detrending is very useful in long\n");
    printf("   observations with modern instrumentation where you can see\n");
    printf("   long timescale power fluctuations.  Strong pulses can skew the\n");
    printf("   statistics of the 1000-bin chunks, though, and caused some\n");
    printf("   suppression in the detection levels of bright pulses (weak\n");
    printf("   pulses are mostly unaffected since they don't strongly change\n");
    printf("   the statistics).  If your data have no long-timescale\n");
    printf("   fluctuations (for instance, if you are processing old 1-bit\n");
    printf("   analog filterbank data which is AC-coupled or if you remove\n");
    printf("   rednoise via realfft/rednoise/(inverse-)realfft), I recommend\n");
    printf("   using the -f/--fast flag.  And if you want to find wide\n");
    printf("   pulses, it might be worth making the chunksize bigger (i.e.\n\n");
    printf("   4000 or 8000).\n\n");

    printf("-- The bad-block detection and removal code can and does remove\n");
    printf("   blocks that have very strong, and particularly strong and broad,\n");
    printf("   pulses in them.  It can also quite effectively remove RFI-\n");
    printf("   infused portions of the data.  Whether to turn it on or off\n");
    printf("   depends on your data.  Note that if there are multiple pulses,\n");
    printf("   only the brightest will usually be \"bad-blocked\" and removed.\n\n");

    printf("-- The fourier-domain matched filtering used here has no phase-\n");
    printf("   dependent effects.  So a 15-bin pulse can be found with equal\n");
    printf("   significance no matter which bin it starts in in the time series.\n\n");

    printf("-- The definition of \"sigma\" used is possibly slightly different\n");
    printf("   from that used in other codes for S/N:\n");
    printf("       sigma = sum(signal-bkgd_level)/RMS/sqrt(boxcar_width)\n");
    printf("   where the bkgd_level is typically 0 after detrending and RMS=1\n");
    printf("   after normalization.  This definition has the advantage that\n");
    printf("   you will get (basically) the same sigma for any pulse no\n");
    printf("   matter how much the input time series has been downsampled as\n");
    printf("   long as the pulse is still resolved.\n\n");

    printf("Copyright Scott Ransom <sransom@nrao.edu>, 2015\n");
}


long long next2_to_n(long long x)
/* Return the first value of 2^n >= x */
{
    long long i = 1;

    while (i < x)
        i <<= 1;
    return i;
}

void cutExtName(char *name, char *marker)
{
  int i,len;
  len=strlen(name);
  for(i=len-1;i>=0;i--)
    if(name[i]==*marker)
      break;
  if(i>0)
    name[i]='\0';
}

FILE *chkfopen(char *path, const char *mode)
{
    FILE *file;

    if ((file = fopen(path, mode)) == NULL) {
        perror("\nError in chkfopen()");
        printf("   path = '%s'\n", path);
        exit(-1);
    }
    return (file);
}

double chk_str2double(char *instr, char *desc)
{
    char tmp[100], *sptr = instr, *endptr;
    double retval;

    retval = strtod(sptr, &endptr);
    if (retval == 0.0 && endptr == instr) {
        sprintf(tmp,
                "Error:  can not convert '%s' to a double (%s) in chk_str2double()\n",
                instr, desc);
        perror(tmp);
        exit(EXIT_FAILURE);
    }
    return retval;
}

long chk_str2long(char *instr, char *desc)
{
    char tmp[100], *sptr = instr, *endptr;
    long retval;

    errno = 0;
    retval = strtol(sptr, &endptr, 10);
    if ((errno == ERANGE && (retval == LONG_MAX || retval == LONG_MIN))
        || (errno != 0 && retval == 0)) {
        sprintf(tmp,
                "Error:  can not convert '%s' to an int/long (%s) in chk_str2long()\n",
                instr, desc);
        perror(tmp);
        exit(EXIT_FAILURE);
    }
    if (endptr == instr) {
        sprintf(tmp,
                "Error:  No digits were found in '%s' for %s in chk_str2long()\n",
                instr, desc);
        perror(tmp);
        exit(EXIT_FAILURE);
    }
    return retval;
}

char *rmtrail(char *str)
/* Removes trailing space from a string */
{
    int i;

    if (str && 0 != (i = strlen(str))) {
        while (--i >= 0) {
            if (!isspace(str[i]))
                break;
        }
        str[++i] = '\0';
    }
    return str;
}

char *rmlead(char *str)
/* Removes leading space from a string */
{
    char *obuf;

    if (str) {
        for (obuf = str; *obuf && isspace(*obuf); ++obuf);
        if (str != obuf)
            strMove(str, obuf);
    }
    return str;
}


char *remove_whitespace(char *str)
/* Remove leading and trailing space from a string */
{
    return rmlead(rmtrail(str));
}

void read_inf_line_valstr(FILE * infofile, char *valstr, char *errdesc)
{
    char line[250], *sptr = NULL;
    int ii, slen;

    sptr = fgets(line, 250, infofile);
    if (sptr != NULL && sptr[0] != '\n' && 0 != (ii = strlen(sptr))) {
        // Check to see if this is a "standard" .inf line
        // which has an '=' in character 40
        if (ii >= 40 && line[40] == '=') {
            sptr = line + 41;
        } else {
            // Else, look for the first '=' back from the end of the line
            while (--ii >= 0) {
                if (sptr[ii] == '=')
                    break;
            }
            if (ii + 1 == 0) {
                sprintf(line,
                        "Error:  no '=' to separate key/val while looking for '%s' in readinf()\n",
                        errdesc);
                perror(line);
                exit(EXIT_FAILURE);
            }
            sptr = line + ii + 1;
        }
        sptr = remove_whitespace(sptr);
        slen = strlen(sptr);
        if (slen) {
            if ((strcmp(errdesc, "data->name") == 0 && slen > 199) ||
                (strcmp(errdesc, "data->telescope") == 0 && slen > 39) ||
                (strcmp(errdesc, "data->band") == 0 && slen > 39) ||
                (strcmp(errdesc, "data->name") != 0 && slen > 99)) {
                sprintf(line,
                        "Error:  value string is too long (%d char) while looking for '%s' in readinf()\n",
                        slen, errdesc);
                perror(line);
                exit(EXIT_FAILURE);
            }
            strcpy(valstr, sptr);
        } else {
            strcpy(valstr, "Unknown");
        }
        return;
    } else {
        if (feof(infofile)) {
            sprintf(line,
                    "Error:  end-of-file while looking for '%s' in readinf()\n",
                    errdesc);
        } else {
            sprintf(line,
                    "Error:  found blank line while looking for '%s' in readinf()\n",
                    errdesc);
        }
        perror(line);
        exit(EXIT_FAILURE);
    }
    // Should never get here....
}

void ra_dec_from_string(char *radec, int *h_or_d, int *m, double *s)
/* Return a values for hours or degrees, minutes and seconds        */
/* given a properly formatted RA or DEC string.                     */
/*   radec is a string with J2000 RA  in the format 'hh:mm:ss.ssss' */
/*   or a string with J2000 DEC in the format 'dd:mm:ss.ssss'       */
{
    int retval;

    radec = remove_whitespace(radec);
    retval = sscanf(radec, "%d:%d:%lf\n", h_or_d, m, s);
    if (retval != 3) {
        char tmp[100];
        sprintf(tmp,
                "Error:  can not convert '%s' to RA or DEC in ra_dec_from_string()\n",
                radec);
        perror(tmp);
        exit(1);
    }
    if (radec[0] == '-' && *h_or_d == 0) {
        *m = -*m;
        *s = -*s;
    }
}



void readinf(infodata * data, char *filenm)
{
    char tmp1[100], tmp2[100], tmp3[100], *infofilenm, *sptr;
    int ii, retval, noteslen = 0;
    FILE *infofile;

    infofilenm = malloc(strlen(filenm) + 5);
    sprintf(infofilenm, "%s.inf", filenm);
    infofile = chkfopen(infofilenm, "r");
    

    read_inf_line_valstr(infofile, data->name, "data->name");
    read_inf_line_valstr(infofile, data->telescope, "data->telescope");
    /* If not using makedata */
    if (strcmp(data->telescope, scopes[0]) != 0) {
        read_inf_line_valstr(infofile, data->instrument, "data->instrument");
        read_inf_line_valstr(infofile, data->object, "data->object");
        read_inf_line_valstr(infofile, tmp1, "RA string");
        ra_dec_from_string(tmp1, &data->ra_h, &data->ra_m, &data->ra_s);
        read_inf_line_valstr(infofile, tmp1, "DEC string");
        ra_dec_from_string(tmp1, &data->dec_d, &data->dec_m, &data->dec_s);
        read_inf_line_valstr(infofile, data->observer, "data->observer");
        read_inf_line_valstr(infofile, tmp1, "MJD string");
        retval = sscanf(tmp1, "%d.%s", &data->mjd_i, tmp2);
        if (retval != 2) {
            sprintf(tmp3, "Error:  can not parse MJD string '%s' in readinf()'\n",
                    tmp1);
            perror(tmp3);
            exit(EXIT_FAILURE);
        }
        sprintf(tmp3, "0.%s", tmp2);
        data->mjd_f = chk_str2double(tmp3, "data->mjd_f");
        read_inf_line_valstr(infofile, tmp1, "data->bary");
        data->bary = chk_str2long(tmp1, "data->bary");
    } else {
        data->mjd_i = -1;
        strcpy(data->object, "fake pulsar");
    }
    read_inf_line_valstr(infofile, tmp1, "data->N");
    data->N = chk_str2double(tmp1, "data->N");
    
    read_inf_line_valstr(infofile, tmp1, "data->dt");
    data->dt = chk_str2double(tmp1, "data->dt");

    read_inf_line_valstr(infofile, tmp1, "data->numonoff");
    data->numonoff = chk_str2long(tmp1, "data->numonoff");
    if (data->numonoff) 
    {
        ii = 0;
        do {
            read_inf_line_valstr(infofile, tmp1, "on-off pairs");
            retval = sscanf(tmp1, "%lf %*[ ,] %lf",
                            &data->onoff[ii], &data->onoff[ii + 1]);
            if (retval != 2) {
                sprintf(tmp3,
                        "Error:  can not parse on-off pair (%d) in readinf()\n",
                        ii / 2);
                perror(tmp3);
                exit(EXIT_FAILURE);
            }
            ii += 2;
        } while (data->onoff[ii - 1] < data->N - 1 && ii < 2 * MAXNUMONOFF);
        data->numonoff = ii / 2;
        if (data->numonoff >= MAXNUMONOFF) {
            sprintf(tmp3,
                    "Error:  number of onoff pairs (%d) >= MAXNUMONOFF (%d) in readinf().\n",
                    data->numonoff, MAXNUMONOFF);
            perror(tmp3);
            exit(EXIT_FAILURE);
        }
    } 
    else 
    {
        data->numonoff = 1;
        data->onoff[0] = 0;
        data->onoff[1] = data->N - 1;
    }
    /* If not using makedata */
    if (strcmp(data->telescope, scopes[0]) != 0) {
        read_inf_line_valstr(infofile, data->band, "data->band");
        if (strcmp(data->band, bands[0]) == 0) {
            read_inf_line_valstr(infofile, tmp1, "data->fov");
            data->fov = chk_str2double(tmp1, "data->fov");
            read_inf_line_valstr(infofile, tmp1, "data->dm");
            data->dm = chk_str2double(tmp1, "data->dm");
            read_inf_line_valstr(infofile, tmp1, "data->freq");
            data->freq = chk_str2double(tmp1, "data->freq");
            read_inf_line_valstr(infofile, tmp1, "data->freqband");
            data->freqband = chk_str2double(tmp1, "data->freqband");
            read_inf_line_valstr(infofile, tmp1, "data->num_chan");
            data->num_chan = chk_str2long(tmp1, "data->num_chan");
            read_inf_line_valstr(infofile, tmp1, "data->chan_wid");
            data->chan_wid = chk_str2double(tmp1, "data->chan_wid");
        } else if ((strcmp(data->band, bands[4]) == 0) ||
                   (strcmp(data->band, bands[5]) == 0)) {
            read_inf_line_valstr(infofile, tmp1, "data->fov");
            data->fov = chk_str2double(tmp1, "data->fov");
            read_inf_line_valstr(infofile, tmp1, "data->energy");
            data->energy = chk_str2double(tmp1, "data->energy");
            read_inf_line_valstr(infofile, tmp1, "data->energyband");
            data->energyband = chk_str2double(tmp1, "data->energyband");
        } else {
            read_inf_line_valstr(infofile, data->filt, "data->filt");
            read_inf_line_valstr(infofile, tmp1, "data->fov");
            data->fov = chk_str2double(tmp1, "data->fov");
            read_inf_line_valstr(infofile, tmp1, "data->wavelen");
            data->wavelen = chk_str2double(tmp1, "data->wavelen");
            read_inf_line_valstr(infofile, tmp1, "data->waveband");
            data->waveband = chk_str2double(tmp1, "data->waveband");
        }
    }
    read_inf_line_valstr(infofile, data->analyzer, "data->analyzer");
    // The following is description line for the 'Notes' part
    sptr = fgets(tmp1, 100, infofile);
    // Now read all the notes lines
    while (1) {
        sptr = fgets(tmp1, 100, infofile);
        if (noteslen + strlen(tmp1) > 500)
            break;
        if (sptr) {
            if (noteslen == 0)
                strcpy(data->notes + noteslen, rmlead(tmp1));
            else
                strcpy(data->notes + noteslen, tmp1);
            noteslen += strlen(data->notes);
        } else {
            if (feof(infofile))
                break;
        }
    }
    fclose(infofile);

    char shell[128];
    sprintf(shell, "cat %s | grep \"Any breaks in the data\" | awk '{print $11}'", infofilenm);
    infofile = popen(shell, "r");
    fscanf(infofile, "%d", &data->breaks);
    pclose(infofile);

    // read_inf_line_valstr(infofile, tmp1, "data->breaks");
    // data->breaks = chk_str2long(tmp1, "data->breaks");


    free(infofilenm);
    
}


void realfft(float idata[], long n, int isign)
/*  This is a modified version of the NR routine with correct (-)  */
/*  exponent.  It uses the above tablesixstepfft making it very    */
/*  fast.  The forward transform (i.e. normal FFT) is isign=-1     */
{
    long nby2, il, ih;
    double cc, h1r, h1i, h2r, h2i, h2rwr, h2iwr, h2rwi, h2iwi;
    double wr, wi, wpr, wpi, tmp1, theta;
    fcomplex *data;

    if (n % 2) {
        printf("\nrealfft() arrays lengths must be evenly divisible by 2.\n\n");
        exit(-1);
    }
    nby2 = n >> 1;
    data = (fcomplex *) idata;
    if (isign == -1) {
        cc = -0.5;
        theta = -TWOPI / (double) n;
        COMPLEXFFT(data, nby2, -1);
    } else {
        cc = 0.5;
        theta = TWOPI / (double) n;
        /* Numerical Recipes gives a sign error for */
        /* the imaginary part of frequency n/2.     */
        if ((n + 2) % 4)
            data[(n >> 2)].i = -data[(n >> 2)].i;
    }
    /* Prep the trig recursion */
    wr = cos(theta);
    wi = sin(theta);
    tmp1 = sin(0.5 * theta);
    wpr = -2.0 * tmp1 * tmp1;
    wpi = wi;
    il = 1;                     /* n     */
    ih = nby2 - il;             /* N/2-n */
    for (; il <= (n >> 2); il++, ih--) {
        h1r = 0.5 * (data[il].r + data[ih].r);
        h1i = 0.5 * (data[il].i - data[ih].i);
        h2r = -cc * (data[il].i + data[ih].i);
        h2i = cc * (data[il].r - data[ih].r);
        h2rwr = h2r * wr;
        h2rwi = h2r * wi;
        h2iwr = h2i * wr;
        h2iwi = h2i * wi;
        data[il].r = h1r + h2rwr - h2iwi;
        data[il].i = h1i + h2iwr + h2rwi;
        data[ih].r = h1r - h2rwr + h2iwi;
        data[ih].i = -h1i + h2iwr + h2rwi;
        tmp1 = wr;
        wr = tmp1 * wpr - wi * wpi + wr;
        wi = wi * wpr + tmp1 * wpi + wi;
    }
    if (isign == -1) {
        /* Set data[0].r to Freq 0 value  */
        /* Set data[0].i to Nyquist value */
        tmp1 = data[0].r;
        data[0].r = tmp1 + data[0].i;
        data[0].i = tmp1 - data[0].i;
    } else {
        /* Numerical Recipes gives a sign error for */
        /* the imaginary part of frequency n/2.     */
        if ((n + 2) % 4)
            data[(n >> 2)].i = -data[(n >> 2)].i;
        tmp1 = data[0].r;
        data[0].r = 0.5 * (tmp1 + data[0].i);
        data[0].i = 0.5 * (tmp1 - data[0].i);
        COMPLEXFFT(data, nby2, 1);
        tmp1 = 2.0 / (double) n;
        for (il = 0; il < n; il++)
            idata[il] *= tmp1;
    }
}

fcomplex *gen_cvect(long length)
{
    fcomplex *v;

#ifdef USE_FFTW_MALLOC
    v = (fcomplex *) fftwf_malloc((size_t) (sizeof(fcomplex) * length));
#else
    v = (fcomplex *) aligned_alloc(ALIGNSIZE, (size_t) (sizeof(fcomplex) * length));
#endif
    if (!v) {
        perror("\nError in gen_cvect()");
        printf("\n");
        exit(-1);
    }
    return v;
}

void vect_free(void *vect)
{
#ifdef USE_FFTW_MALLOC
    fftwf_free(vect);
#else
    free(vect);
#endif
}

void fftwcall(fcomplex * indata, long nn, int isign)
/* This routine calls the FFTW complex-complex FFT using stored wisdom */
/* files.  It is VERY fast.  nn does _not_ have to be a power of two   */
/* size.  indata is a complex array but stored as floats.              */
{
    fftwf_plan *plan_forward, *plan_inverse;
    fftwf_complex *dataptr = (fftwf_complex *) indata;
    int ii, indata_align, slot, incache = 0, oldestplan = 0;
    //static int goodct = 0, badct = 0;
    static fftwf_plan plancache_forward[4] = { NULL, NULL, NULL, NULL };
    static fftwf_plan plancache_inverse[4] = { NULL, NULL, NULL, NULL };
    static int aligncache[4] = { -99, -99, -99, -99 };
    static int firsttime = 1;
    static int lastslot = 0, lastused[4] = { 0, 0, 0, 0 };
    static long nncache[4] = { 0, 0, 0, 0 };

    // This determines the alignment of the input array.  Allows
    // more flexible calling of FFTW using its plans.
    // A return value of 0 is "properly" aligned.
    indata_align = is_aligned(indata, 16);
    //if (indata_align)
    //    printf("Data not properly aligned (%d)!\n", indata_align);

    // Call the six-step algorithm if the FFT is too big to be
    // efficiently handled by FFTW
    if (nn > BIGFFTWSIZE) {
        tablesixstepfft(indata, nn, isign);
        return;
    }
    // If calling for the first time, read the wisdom file
    if (firsttime)
        read_wisdom();

    // If we used the same plan during the last few calls, use it
    // again.  We keep, in effect, a stack of the 4 most recent plans.
    ii = 0;
    slot = lastslot;
    while (ii < 4) {
        if (nn == nncache[slot] && indata_align == aligncache[slot]) {
            plan_forward = &plancache_forward[slot];
            plan_inverse = &plancache_inverse[slot];
            lastused[slot] = 0;
            lastused[(slot + 1) % 4]++;
            lastused[(slot + 2) % 4]++;
            lastused[(slot + 3) % 4]++;
            //printf("Found plan in slot %d (iter = %d):  nn=%ld  align=%d  number=%d\n",
            //       slot, ii, nn, aligncache[slot], goodct++);
            lastslot = slot;
            incache = 1;
            break;
        }
        slot = (slot + 1) % 4;
        ii++;
    }
    if (!incache) {
        unsigned int planflag;
        fcomplex *datacopy;
        if (!firsttime) {
            for (ii = 3; ii >= 0; ii--)
                if (lastused[ii] >= oldestplan)
                    oldestplan = ii;
            // Delete the old plans to prevent memory leaks
            if (plancache_forward[oldestplan])
                fftwf_destroy_plan(plancache_forward[oldestplan]);
            if (plancache_inverse[oldestplan])
                fftwf_destroy_plan(plancache_inverse[oldestplan]);
        }
        //printf("Making a new plan for nn=%ld, align=%d (dropping nn=%ld) %d\n",
        //       nn, indata_align, nncache[oldestplan], badct++);
        // We don't want to wait around to measure huge transforms
        planflag = (nn > 16384) ? FFTW_ESTIMATE : FFTW_MEASURE;
        // FFTW_MEASURE will destroy the input/output data, so copy it
        datacopy = gen_cvect(nn);
        memcpy(datacopy, dataptr, nn * sizeof(fcomplex));
        // Actually make the plans
        plancache_forward[oldestplan] =
            fftwf_plan_dft_1d(nn, dataptr, dataptr, -1, planflag);
        plancache_inverse[oldestplan] =
            fftwf_plan_dft_1d(nn, dataptr, dataptr, +1, planflag);
        // Now copy the input data back
        memcpy(dataptr, datacopy, nn * sizeof(fcomplex));
        vect_free(datacopy);
        nncache[oldestplan] = nn;
        aligncache[oldestplan] = indata_align;
        plan_forward = &plancache_forward[oldestplan];
        plan_inverse = &plancache_inverse[oldestplan];
        lastused[oldestplan] = 0;
        lastused[(oldestplan + 1) % 4]++;
        lastused[(oldestplan + 2) % 4]++;
        lastused[(oldestplan + 3) % 4]++;
        lastslot = oldestplan;
    }
    // Call the transform using the "new-array" functionality of FFTW
    if (isign == -1) {
        fftwf_execute_dft(*plan_forward, dataptr, dataptr);
    } else {
        fftwf_execute_dft(*plan_inverse, dataptr, dataptr);
    }
    firsttime = 0;
}

void tablesixstepfft(fcomplex * indata, long nn, int isign)
/*  This is a modified version of a six-step-FFT routine from the    */
/*  apfloat() package.  It is a standard complex FFT.                */
/*  It uses a split-radix, table-look-up, six-step FFT.              */
/*  It is very fast for huge transforms due to high memory locality. */
/*  The forward transform (i.e. normal FFT) is isign=-1              */
{
    long n1, n2, jj, kk, kind;
    double wpr, wpi, wr, wi, theta, tmp1, tmp2;
    int move_size;
    unsigned char *move;
    fftwf_plan plan;
    static fftwf_plan last_plan = NULL;
    static int firsttime = 1, lastn = 0, lastisign = 0;

    /* If calling for the first time, read the wisdom file */
    if (firsttime)
        read_wisdom();

    if (nn < 2)
        return;

    /* Treat the input data as a n1 (rows) x n2 (cols) */
    /* matrix.  Make sure that n2 >= n1.               */

    n1 = good_factor(nn);
    if (n1 == 0) {
        printf("\nLength of FFT in tablesixstepfft() must be factorable\n\n");
        exit(0);
    }
    n2 = nn / n1;

    /* transpose scratch space */

    move_size = (n1 + n2) / 2;
    move = (unsigned char *) malloc(move_size);

    /* first transpose the matrix */

    transpose_fcomplex(indata, n1, n2, move, move_size);

    /* then do n2 transforms of length n1 */

    /* Use FFTW for the small transforms if available. */

    if (n1 == lastn && isign == lastisign) {
        plan = last_plan;
    } else {
        const int N1[1] = { n1 };
        if (firsttime)
            firsttime = 0;
        else
            fftwf_destroy_plan(last_plan);
        plan = fftwf_plan_many_dft(1, N1, n2,
                                   (fftwf_complex *) indata, NULL, 1, n1,
                                   (fftwf_complex *) indata, NULL, 1, n1,
                                   isign, FFTW_ESTIMATE);
        last_plan = plan;
        lastn = n1;
        lastisign = isign;
    }
    fftwf_execute(plan);

    /* transpose the matrix */

    transpose_fcomplex(indata, n2, n1, move, move_size);

    /* then multiply the matrix A_jk by exp(isign * 2 pi i j k / nn) */
    /* Use recursion formulas from NR */

    for (jj = 1; jj < n1; jj++) {
        theta = isign * jj * TWOPI / nn;
        wr = cos(theta);
        wi = sin(theta);
        tmp1 = sin(0.5 * theta);
        wpr = -2.0 * tmp1 * tmp1;
        wpi = wi;
        kind = jj * n2 + 1;
        for (kk = 1; kk < n2; kk++, kind++) {
            tmp1 = indata[kind].r;
            tmp2 = indata[kind].i;
            indata[kind].r = tmp1 * wr - tmp2 * wi;
            indata[kind].i = tmp2 * wr + tmp1 * wi;
            tmp1 = wr;
            wr = tmp1 * wpr - wi * wpi + wr;
            wi = wi * wpr + tmp1 * wpi + wi;
        }
    }

    /* then do n1 transforms of length n2 */

    /* Use FFTW for the small transforms if available. */

    if (n2 == lastn && isign == lastisign) {
        plan = last_plan;
    } else {
        const int N2[1] = { n2 };
        fftwf_destroy_plan(last_plan);
        plan = fftwf_plan_many_dft(1, N2, n1,
                                   (fftwf_complex *) indata, NULL, 1, n2,
                                   (fftwf_complex *) indata, NULL, 1, n2,
                                   isign, FFTW_ESTIMATE);
        last_plan = plan;
        lastn = n2;
        lastisign = isign;
    }
    fftwf_execute(plan);

    // Comment this out so it matches FFTW
    // Scale the FFT if it is an inverse FFT
    //if (isign == 1) {
    //   tmp1 = 1.0 / (double) nn;
    //   for (jj = 0; jj < n1 * n2; jj++) {
    //      indata[jj].r *= tmp1;
    //      indata[jj].i *= tmp1;
    //   }
    //}

    /* last transpose the matrix */

    transpose_fcomplex(indata, n1, n2, move, move_size);
    free(move);
}


void read_wisdom(void)
{
    FILE *wisdomfile;
    static char wisdomfilenm[120];

    /* First try to import the system wisdom if available */
    fftwf_import_system_wisdom();
    sprintf(wisdomfilenm, "%s/lib/fftw_wisdom.txt", getenv("PRESTO"));
    wisdomfile = fopen(wisdomfilenm, "r");
    if (wisdomfile == NULL) {
        printf("Warning:  Couldn't open '%s'\n"
               "          You should run 'makewisdom'.  See $PRESTO/INSTALL.\n",
               wisdomfilenm);
    } else {
        fftwf_import_wisdom_from_file(wisdomfile);
        // if (!fftwf_import_wisdom_from_file(wisdomfile))
        //     printf("Warning:  '%s' is not up-to-date.\n"
        //            "          You should run 'makewisdom'.  See $PRESTO/INSTALL.\n",
        //            wisdomfilenm);
        fclose(wisdomfile);
    }
    // The following resets errno if one of the wisdom files was not found
    errno = 0;
}


long long good_factor(long long nn)
/* Return the factor of a number that is closest to its sqrt. */
/* If the number is prime, return 0.                          */
{
    long long pp;

    /* Optimal factoring is one factor twice the size of the other */
    /* Try this brute force...                                     */

    pp = (long long) sqrt((double) (nn / 2));
    if (2 * pp * pp == nn)
        return pp;

    /* Calculate the best (closest to each other) factors */
    /* This is certainly not the best way to do this...   */

    pp = (long long) sqrt((double) nn);
    while (pp > 1) {
        if (nn % pp == 0)
            return pp;
        pp--;
    }
    return 0;
}

static int TOMS_gcd(int a, int b)
/* Return the greatest common denominator of 'a' and 'b' */
{
    int r;
    do {
        r = a % b;
        a = b;
        b = r;
    } while (r != 0);

    return a;
}

short transpose_fcomplex(fcomplex * a, int nx, int ny, unsigned char *move,
                         int move_size)
/*
 * TOMS Transpose.  Revised version of algorithm 380.
 * 
 * These routines do in-place transposes of arrays.
 * 
 * [ Cate, E.G. and Twigg, D.W., ACM Transactions on Mathematical Software, 
 *   vol. 3, no. 1, 104-110 (1977) ]
 * 
 * C version by Steven G. Johnson. February 1997.
 *
 * "a" is a 1D array of length ny*nx which contains the nx x ny matrix to be
 * transposed.  "a" is stored in C order (last index varies fastest).  move
 * is a 1D array of length move_size used to store information to speed up
 * the process.  The value move_size=(ny+nx)/2 is recommended.
 * 
 * The return value indicates the success or failure of the routine. Returns 0
 * if okay, -1 if ny or nx < 0, and -2 if move_size < 1. The return value
 * should never be positive, but it it is, it is set to the final position in
 * a when the search is completed but some elements have not been moved.
 * 
 * Note: move[i] will stay zero for fixed points.
 */
{
    int i, j, im, mn;
    fcomplex b, c, d;
    int ncount;
    int k;

    /* check arguments and initialize: */
    if (ny < 0 || nx < 0)
        return -1;
    if (ny < 2 || nx < 2)
        return 0;
    if (move_size < 1)
        return -2;

    if (ny == nx) {
        /*
         * if matrix is square, exchange elements a(i,j) and a(j,i):
         */
        for (i = 0; i < nx; ++i)
            for (j = i + 1; j < nx; ++j) {
                b = a[i + j * nx];
                a[i + j * nx] = a[j + i * nx];
                a[j + i * nx] = b;
            }
        return 0;
    }
    ncount = 2;                 /* always at least 2 fixed points */
    k = (mn = ny * nx) - 1;

    for (i = 0; i < move_size; ++i)
        move[i] = 0;

    if (ny >= 3 && nx >= 3)
        ncount += TOMS_gcd(ny - 1, nx - 1) - 1; /* # fixed points */

    i = 1;
    im = ny;

    while (1) {
        int i1, i2, i1c, i2c;
        int kmi;

    /** Rearrange the elements of a loop
	and its companion loop: **/

        i1 = i;
        kmi = k - i;
        b = a[i1];
        i1c = kmi;
        c = a[i1c];

        while (1) {
            i2 = ny * i1 - k * (i1 / nx);
            i2c = k - i2;
            if (i1 < move_size)
                move[i1] = 1;
            if (i1c < move_size)
                move[i1c] = 1;
            ncount += 2;
            if (i2 == i)
                break;
            if (i2 == kmi) {
                d = b;
                b = c;
                c = d;
                break;
            }
            a[i1] = a[i2];
            a[i1c] = a[i2c];
            i1 = i2;
            i1c = i2c;
        }
        a[i1] = b;
        a[i1c] = c;

        if (ncount >= mn)
            break;              /* we've moved all elements */

    /** Search for loops to rearrange: **/

        while (1) {
            int max;

            max = k - i;
            ++i;
            if (i > max)
                return i;
            im += ny;
            if (im > k)
                im -= k;
            i2 = im;
            if (i == i2)
                continue;
            if (i >= move_size) {
                while (i2 > i && i2 < max) {
                    i1 = i2;
                    i2 = ny * i1 - k * (i1 / nx);
                }
                if (i2 == i)
                    break;
            } else if (!move[i])
                break;
        }
    }

    return 0;
}


#define ELEM_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }

float median(float arr[], int n)
{
    int low, high;
    int median;
    int middle, ll, hh;

    low = 0;
    high = n - 1;
    median = (low + high) / 2;
    for (;;) {
        if (high <= low)        /* One element only */
            return arr[median];
        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]);
            return arr[median];
        }
        /* Find median of low, middle and high items; swap into position low */
        middle = (low + high) / 2;
        if (arr[middle] > arr[high])
            ELEM_SWAP(arr[middle], arr[high]);
        if (arr[low] > arr[high])
            ELEM_SWAP(arr[low], arr[high]);
        if (arr[middle] > arr[low])
            ELEM_SWAP(arr[middle], arr[low]);
        /* Swap low item (now in position middle) into position (low+1) */
        ELEM_SWAP(arr[middle], arr[low + 1]);
        /* Nibble from each end towards middle, swapping items when stuck */
        ll = low + 1;
        hh = high;
        for (;;) {
            do
                ll++;
            while (arr[low] > arr[ll]);
            do
                hh--;
            while (arr[hh] > arr[low]);
            if (hh < ll)
                break;
            ELEM_SWAP(arr[ll], arr[hh]);
        }
        /* Swap middle item (in position low) back into correct position */
        ELEM_SWAP(arr[low], arr[hh]);
        /* Re-set active partition */
        if (hh <= median)
            low = ll;
        if (hh >= median)
            high = hh - 1;
    }
}

#undef ELEM_SWAP

// 交换两个元素的值
void swap(float* a, float* b) 
{
    float temp = *a;
    *a = *b;
    *b = temp;
}

void swap_d(double* a, double* b) 
{
    double temp = *a;
    *a = *b;
    *b = temp;
}

// 快速排序的分割函数
int partition(float arr[], int low, int high) 
{
    float pivot = arr[low]; // 选择第一个元素作为基准
    while (low < high) {
        while (low < high && arr[high] >= pivot)
            high--;
        arr[low] = arr[high]; // 将比基准小的元素移到低端
        while (low < high && arr[low] <= pivot)
            low++;
        arr[high] = arr[low]; // 将比基准大的元素移到高端
    }
    arr[low] = pivot; // 基准元素归位
    return low;
}

int partition_d(double arr[], int low, int high) 
{
    double pivot = arr[low]; // 选择第一个元素作为基准
    while (low < high) {
        while (low < high && arr[high] >= pivot)
            high--;
        arr[low] = arr[high]; // 将比基准小的元素移到低端
        while (low < high && arr[low] <= pivot)
            low++;
        arr[high] = arr[low]; // 将比基准大的元素移到高端
    }
    arr[low] = pivot; // 基准元素归位
    return low;
}

// 快速排序递归函数
void quickSort(float arr[], int low, int high) 
{
    if (low < high) {
        int pivotPos = partition(arr, low, high);
        quickSort(arr, low, pivotPos - 1); // 对左半部分进行快速排序
        quickSort(arr, pivotPos + 1, high); // 对右半部分进行快速排序
    }
}

void quickSort_d(double arr[], int low, int high) 
{
    if (low < high) {
        int pivotPos = partition_d(arr, low, high);
        quickSort_d(arr, low, pivotPos - 1); // 对左半部分进行快速排序
        quickSort_d(arr, pivotPos + 1, high); // 对右半部分进行快速排序
    }
}

// 快速排序接口函数
void quickSortWrapper(float arr[], int size) 
{
    quickSort(arr, 0, size - 1);
}

void quickSortWrapper_d(double arr[], int size) 
{
    quickSort_d(arr, 0, size - 1);
}



int dofit(const gsl_multifit_robust_type *T,
      const gsl_matrix *X, const gsl_vector *y,
      gsl_vector *c, gsl_matrix *cov)
{
  int s;
  gsl_multifit_robust_workspace * work
    = gsl_multifit_robust_alloc (T, X->size1, X->size2);

  s = gsl_multifit_robust (X, y, c, cov, work);
  gsl_multifit_robust_free (work);

  return s;
}

void fit_line_robust(float *indata, gsl_vector *x, gsl_vector *y, gsl_matrix *X, gsl_matrix *cov, gsl_vector *c, int N)
{
    int i;

    /* generate linear dataset */
    for (i = 0; i < N; i++)
    {
        gsl_vector_set (x, i, 1.0*i);
        gsl_vector_set (y, i, (double)(indata[i]));

        gsl_matrix_set (X, i, 0, 1.0);
        gsl_matrix_set (X, i, 1, 1.0*i);
    }
    /* construct design matrix X for linear fit */
    // for (i = 0; i < N; ++i)
    // {
    //     double xi = gsl_vector_get(x, i);

    //     gsl_matrix_set (X, i, 0, 1.0);
    //     gsl_matrix_set (X, i, 1, xi);
    // }

    dofit(gsl_multifit_robust_bisquare, X, y, c, cov);

}

void fit_line_samp(double *indata, double *xarray, double *c0, double *c1, int N)
{
    double cov00, cov01, cov11, chisq;
    gsl_fit_linear(xarray, 1, indata, 1, N, c0, c1, &cov00, &cov01, &cov11, &chisq);
}


void fft_convolve(float *fftd_data, float *fftd_kern, float *prod, int N)
{
    // """
    // fft_convolve(fftd_data, fftd_kern, lo, hi):
    // Perform a convolution with the complex floating point vectors
    //        'fftd_data' and 'fftd_kern'.  The returned vector will start at
    //        at bin 'lo' (must be an integer), and go up to but not
    //        include bin 'hi' (also an integer).
    // """

    int i;
    float a, b;
    float c, d;
    for(i=1; i<N/2; i++)
    {
        a = fftd_data[i*2];
        b = fftd_data[i*2+1];
        c = fftd_kern[i*2];
        d = fftd_kern[i*2+1];
        prod[i*2] = a*c - b*d;
        prod[i*2+1] = a*d + b*c;
    }
    prod[0] = fftd_data[0]*fftd_kern[0];
    prod[1] = fftd_data[1]*fftd_kern[1];
    realfft(prod, N, 1);
}

void prune_related1(int *hibins, float *hivals, int N, int downfact, int *off)
{
    // Remove candidates that are close to other candidates
    // but less significant.  This one works on the raw 
    // candidate arrays and uses the single downfact
    // that they were selected with.

    int i, j;
    int *toremove = malloc(sizeof(int)*N);

    memset(toremove, 0, N * sizeof(int));

    int xbin;
    float xsigma;
    int ybin;
    float ysigma;

    for(i=0; i<N-1; i++)
    {
        if (toremove[i] == 1) continue;

        xsigma = hivals[i];
        xbin = hibins[i];
        for(j=i+1; j<N; j++)
        {
            ysigma = hivals[j];
            ybin = hibins[j];
            if (abs(ybin-xbin) > (int)(downfact/2))  break;
            else 
            {
                if (toremove[j] == 1) continue;
                if(xsigma > ysigma)
                    toremove[j] = 1;
                else
                    toremove[i] = 1;
            }
        }
    }

    int offset = 0;
    for(i=0; i<N; i++)
    {
        if (toremove[i] == 1) 
        {
            offset++;
        } else {
            hibins[i - offset] = hibins[i];
            hivals[i - offset] = hivals[i];
        }
    }

    *off = offset;
    free(toremove);
}

static int max(int a, int b, int c) 
{
    int max_value = a;

    if (b > max_value) {
        max_value = b;
    }
    if (c > max_value) {
        max_value = c;
    }
    return max_value;
}


void prune_related2(Candidate *dm_candlist, int dm_candlist_count, int downfacts[], int N, int *off)
{
    int* toremove = (int*)malloc(dm_candlist_count * sizeof(int));
    int ii, jj;
    int xbin, ybin ;
    int prox;
    float xsigma, ysigma;
    Candidate xx, yy;

    memset(toremove, 0, dm_candlist_count * sizeof(int));

    for ( ii = 0; ii < dm_candlist_count - 1; ii++)
    {
        if (toremove[ii] == 1) continue;
        
        xx = dm_candlist[ii];
        xsigma = xx.val;
        xbin = xx.bin;
    
        for ( jj = ii + 1; jj < dm_candlist_count; jj++)
        {
            if (toremove[jj] == 1) continue;

            yy = dm_candlist[jj];
            ysigma = yy.val;
            ybin = yy.bin;
            
            if (abs(ybin - xbin) > downfacts[N-1] / 2)  break;
            else
            {
                prox = max(xx.downfact / 2, yy.downfact / 2, 1);

                if (abs(ybin - xbin) <= prox) 
                {
                    if (xsigma > ysigma) {
                        toremove[jj] = 1;
                    } else {
                        toremove[ii] = 1;
                    }
                }
            }
        }
        
    }

    int offset = 0;
    for (ii = 0; ii < dm_candlist_count; ii++) {
        if (toremove[ii] == 1) {
            offset++;
        } else {
            dm_candlist[ii - offset] = dm_candlist[ii];
        }
    }

    *off = offset;

    free(toremove);
}


// Define comparison function for sorting candidate objects
int mycmp(float a, float b) 
{
    return ((a > b) - (a < b));
}

// Define comparison functions for sorting candidate objects
int cmp(const void *a, const void *b) 
{
    Candidate *ca = (Candidate *)a;
    Candidate *cb = (Candidate *)b;
    return mycmp(ca->bin, cb->bin);
}


// Function to add a candidate to the list
void addCandidate(CandidateList *candlist, float DM, float val, double time, int bin, int valid) 
{
    candlist->candidates = realloc(candlist->candidates, (candlist->count + 1) * sizeof(Candidate));
    candlist->candidates[candlist->count].DM = DM;
    candlist->candidates[candlist->count].val = val;
    candlist->candidates[candlist->count].time = time;
    candlist->candidates[candlist->count].bin = bin;
    candlist->candidates[candlist->count].downfact = valid;
    candlist->count++;
}


void prune_border_cases(Candidate *dm_candlist, int dm_candlist_count, long offregions[], int *off)
{
    int* toremove = (int*)malloc(dm_candlist_count * sizeof(int));
    memset(toremove, 0, dm_candlist_count * sizeof(int));

    int ii;

    Candidate cand;
    int loside, hiside;
    
    for(ii=dm_candlist_count-1; ii>=0; ii++)
    {
        cand = dm_candlist[ii];
        loside = cand.bin - cand.downfact/2;
        hiside = cand.bin + cand.downfact/2;
        if(hiside < offregions[0])   break;
        if(hiside > offregions[0] && loside < offregions[1])
            toremove[ii] = 1;
    }

    int offset = 0;
    for(ii=0; ii<dm_candlist_count; ii++)
    {
        if(toremove[ii] == 1)
            offset++;
        else
            dm_candlist[ii-offset] = dm_candlist[ii];
    }
    *off = offset;

    free(toremove);
}