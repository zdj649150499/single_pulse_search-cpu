#include <unistd.h>
#include <getopt.h>

#include "functions.h"

int main(int argc,char **argv)
{
    int i, j, k, ii;


    int useffts = 1;

    float opt_maxwidth = 0.0;
    float opt_threshold = 5.0;
    float opt_start=0.0;
    float opt_end=1e9;
    int opt_fast = 0;
    int opt_badblocks=1;
    int opt_detrendlen = 1;

    char opt_glob[128];

    int fftlen = 8192;  // Should be a power-of-two for best speed
    int chunklen = 8000;  // Must be at least max_downfact less than fftlen

    char filename[512];
    char filenmbase[512];
    char infoname[512];
    char outfile[512];

    if(argc < 2)
    {
        help();
        exit(0);
    }
    for(i=0;i<argc;i++)
    {
        if (strcmp(argv[i],"-h")==0)
        {
            help();
            exit(0);
        }
        else if (strcmp(argv[i],"-H")==0)
        {
            help();
            help_H();
            exit(0);
        }
        else if (strcmp(argv[i],"-m")==0)
        {
            sscanf(argv[++i],"%f",&opt_maxwidth);
        }
        else if (strcmp(argv[i],"-t")==0)
        {
            sscanf(argv[++i],"%f",&opt_threshold);
        }
        else if (strcmp(argv[i],"-s")==0)
        {
            sscanf(argv[++i],"%f",&opt_start);
        }
        else if (strcmp(argv[i],"-e")==0)
        {
            sscanf(argv[++i],"%f",&opt_end);
        }
        else if (strcmp(argv[i],"-g")==0)
        {
            sscanf(argv[++i],"%s",opt_glob);
        }
        else if (strcmp(argv[i],"-f")==0)
        {
            opt_fast = 1;
        }
        else if (strcmp(argv[i],"-b")==0)
        {
            opt_badblocks = 0;
        }
        else if (strcmp(argv[i],"-d")==0)
        {
            sscanf(argv[++i],"%d",&opt_detrendlen);
        }
    }
    strcpy(filename, argv[argc-1]);
    strcpy(filenmbase, filename);
    cutExtName(filenmbase,".");
    strcpy(infoname, filenmbase);
    strcat(infoname, ".inf");
    strcpy(outfile, filenmbase);
    strcat(outfile, ".singlepulse");

    infodata info;
    readinf(&info, filenmbase);


    if(opt_detrendlen != 1 &&  opt_detrendlen != 2 && opt_detrendlen != 4 && opt_detrendlen != 8 && opt_detrendlen != 16 && opt_detrendlen != 32)
    {
        printf(" -d   detrendlen mast be 1,2,4,8,16,32\n");
        exit(0);
    }

    int detrendlen = opt_detrendlen*1000;
    if(detrendlen > chunklen)
    {
        chunklen = detrendlen;
        fftlen = (int)(next2_to_n(chunklen));
    }

    int blocks_per_chunk = chunklen/detrendlen;
    int overlap = (fftlen - chunklen)/2;
    int worklen = chunklen + 2*overlap;  // currently it is fftlen...

    int max_downfact = 30;
    int downfacts_n = 14;
    int default_downfacts[14] = {2, 3, 4, 6, 9, 14, 20, 30, 45, 70, 100, 150, 220, 300};


    
    // float candlist;
    // float num_v_DMstr;



    // float DMstr = info.dm;
    long N = info.N;
    double dt = info.dt;
    // double obstime = N*dt;

    // Choose the maximum width to search based on time instead
    // of bins.  This helps prevent increased S/N when the downsampling
    // changes as the DM gets larger.

    int downfacts_N = 0;
    if(opt_maxwidth > 0.0f)
    {
        for(i=0; i<downfacts_n; i++)
        {
            if(default_downfacts[i]*dt <= opt_maxwidth)
                downfacts_N++;
        }
    }
    else
    {
        for(i=0; i<downfacts_n; i++)
        {
            if(default_downfacts[i] <= max_downfact)
                downfacts_N++;
        }
    }
    if(downfacts_N==0)
        downfacts_N = 1;

    int downfacts[downfacts_N];
    
    if(downfacts_N ==1 )
        downfacts[0] = default_downfacts[0];
    else 
    {
        downfacts_N = 0;
        if(opt_maxwidth > 0.0f)
        {
            for(i=0; i<downfacts_n; i++)
            {
                if(default_downfacts[i]*dt <= opt_maxwidth)
                {
                    downfacts[downfacts_N] = default_downfacts[i];
                    downfacts_N++;
                }
            }
        }
        else
        {
            for(i=0; i<downfacts_n; i++)
            {
                if(default_downfacts[i] <= max_downfact)
                {
                    downfacts[downfacts_N] = default_downfacts[i];
                    downfacts_N++;
                }
            }
        }
    }

    // long orig_N = N;
    // double orig_dt = dt;

    float *fftd_kerns_f = malloc(sizeof(float)*fftlen*downfacts_N);
    for(j=0; j< fftlen*downfacts_N; j++)
    {
        fftd_kerns_f[j] = 0.0f;
    }
    int downfact;
    for(i=0; i<downfacts_N; i++)
    {
        downfact = default_downfacts[i];
        if(downfact%2)
        {
            for(j=0; j<downfact/2+1; j++)
                fftd_kerns_f[i*fftlen+j] = 1.0/sqrt(1.0*downfact);
            for(j=fftlen-(downfact/2); j<fftlen; j++)
                fftd_kerns_f[i*fftlen+j] = 1.0f/sqrt(1.0*downfact);
        }
        else
        {
            for(j=0; j<downfact/2+1; j++)
                fftd_kerns_f[i*fftlen+j] = 1.0f/sqrt(1.0*downfact);
            if(downfact>2)
                for(j=fftlen-(downfact/2-1); j<fftlen; j++)
                    fftd_kerns_f[i*fftlen+j] = 1.0f/sqrt(1.0*downfact);
        }
        realfft(fftd_kerns_f+i*fftlen, fftlen, -1);
    }

    long offregions[2];
    offregions[0] = (long) info.onoff[1];
    offregions[1] = (long) info.onoff[3];
    
    if(info.breaks)
    {
        // If last break spans to end of file, don't read it in (its just padding)
        if(offregions[1] == N-1)
            N = offregions[0] + 1;
    }
    
    // Compute the file length in detrendlens
    int roundN = ((int)(N/detrendlen))*detrendlen;
    int numchunks = roundN/chunklen;


    // Read in the file
    // printf("Reading %s ...\n", filename);
    float *timeseries = malloc(sizeof(float)*roundN);
    FILE *FP;
    FP = fopen(filename, "r");
    fread(timeseries, sizeof(float), roundN, FP);
    fclose(FP);

    // Split the timeseries into chunks for detrending
    int numblocks = roundN/detrendlen;
    double *stds = malloc(sizeof(double)*numblocks);
    int  *bad_blocks = malloc(sizeof(int)*numblocks);
    memset(bad_blocks, -1, sizeof(int)*numblocks);
    
    int bad_blocks_size = 0;
    double *sort_stds = malloc(sizeof(double)*numblocks);

    // for(i=0; i<numblocks; i++)
    //     stds[i] = 0.0;

    float *tmpchunk = malloc(sizeof(float)*detrendlen);
    double sum;
    int edge = detrendlen/40;
    if(opt_fast)
    {
        // use median removal instead of detrending (2x speedup)
        for(i=0; i<numblocks; i++)
        {
            memcpy(tmpchunk, timeseries+i*detrendlen, sizeof(float)*detrendlen);
            quickSortWrapper(tmpchunk, detrendlen);
            float med = tmpchunk[detrendlen/2];
            
            sum = 0.0;
            for(j=0; j<detrendlen; j++)
            {
                timeseries[j+i*detrendlen] -= med;
                tmpchunk[j] -= med;

                if(j >= edge && j< detrendlen-edge)
                    sum += (tmpchunk[j]*tmpchunk[j]);
            }
            sort_stds[i] = stds[i] = sqrt(sum/(0.95*detrendlen))*1.148;
        }
    }
    else
    {
        double c0, c1;
        double sum;
        /* sample linear fit , good */
        double *xtime_d = malloc(sizeof(double)*detrendlen);
        double *timeseries_d = malloc(sizeof(double)*detrendlen);
        for(i=0;i<detrendlen;i++)
        {
            xtime_d[i] = i*1.0f;
        }
        // The detrend calls are the most expensive in the program
        for(i=0; i<numblocks; i++)
        {
            for(j=0; j<detrendlen; j++)
            {
                timeseries_d[j] = timeseries[j+i*detrendlen];
            }
            fit_line_samp(timeseries_d, xtime_d, &c0, &c1, detrendlen);
            
            for(j=0; j<detrendlen; j++)
            {
                tmpchunk[j] = timeseries[j+i*detrendlen] = timeseries[j+i*detrendlen] - (c0+c1*j*1.0f);
                // tmpchunk[j] = timeseries[j+i*detrendlen];
            }
            quickSortWrapper(tmpchunk, detrendlen);
            sum = 0.0;
            for(j=edge; j<detrendlen-edge; j++)
            {
                sum += (tmpchunk[j]*tmpchunk[j]);
            }
            // The following gets rid of (hopefully) most of the 
            // outlying values (i.e. power dropouts and single pulses)
            // If you throw out 5% (2.5% at bottom and 2.5% at top)
            // of random gaussian deviates, the measured stdev is ~0.871
            // of the true stdev.  Thus the 1.0/0.871=1.148 correction below.
            // The following is roughly .std() since we already removed the median
            sort_stds[i] = stds[i] = sqrt(sum/(0.95*detrendlen))*1.148;
        }
        free(xtime_d);
        free(timeseries_d);

        /* robust fit, slow*/
        // gsl_matrix *X, *cov;
        // gsl_vector *c;
        // gsl_vector *x, *y;

        // X = gsl_matrix_alloc (detrendlen, 2);
        // c = gsl_vector_alloc (2);
        // cov = gsl_matrix_alloc (2, 2);
        // x = gsl_vector_alloc (detrendlen);
        // y = gsl_vector_alloc (detrendlen);

        // for(i=0; i<numblocks; i++)
        // {
        //     fit_line_robust(timeseries+i*detrendlen, x, y, X, cov, c, detrendlen);
        //     c0=gsl_vector_get(c,0);
        //     c1=gsl_vector_get(c,1);
        //     for(j=0; j<detrendlen; j++)
        //     {
        //         timeseries[j+i*detrendlen] -= (c0+c1*j*1.0f);
        //         tmpchunk[j] = timeseries[j+i*detrendlen];
        //     }
        // }
        // gsl_matrix_free (X);
        // gsl_vector_free (c);
        // gsl_matrix_free (cov);
        // gsl_vector_free (x);
        // gsl_vector_free (y);
    }
    free(tmpchunk);
    
    // sort the standard deviations and separate those with
    // very low or very high values
    quickSortWrapper_d(sort_stds, numblocks);
    // identify the differences with the larges values (this
    // will split off the chunks with very low and very high stds
    int locut, hicut;
    double sort_stds_max = sort_stds[1] - sort_stds[0];
    locut = 1;
    for(i=2; i<numblocks/2+1; i++)
    {
        if(sort_stds_max < (sort_stds[i] - sort_stds[i-1]))
        {
            sort_stds_max = sort_stds[i] - sort_stds[i-1];
            locut = i;
        }
    }
    sort_stds_max = sort_stds[numblocks/2+1] - sort_stds[numblocks/2];
    hicut = numblocks/2+1;
    for(i=numblocks/2+2; i<numblocks; i++)
    {
        if(sort_stds_max < (sort_stds[i] - sort_stds[i-1]))
        {
            sort_stds_max = sort_stds[i] - sort_stds[i-1];
            hicut = i;
        }
    }
    hicut = hicut - 3;
    double std_stds = sqrt(gsl_stats_variance(sort_stds+locut, 1, hicut-locut+1));
    double median_stds = sort_stds[(locut+hicut)/2];
    // printf("    pseudo-median block standard deviation = %.2f\n", median_stds);

    if(opt_badblocks)
    {
        double lo_std = median_stds - 4.0*std_stds;
        double hi_std = median_stds + 4.0*std_stds;
        // Determine a list of "bad" chunks.  We will not search these.
        bad_blocks_size=0;
        for(i=0; i<numblocks; i++)
        {
            if(stds[i]<lo_std || stds[i]>hi_std)
            {
                bad_blocks[bad_blocks_size] = i;
                bad_blocks_size++;
                stds[i] = median_stds;
                int loind = i*detrendlen;
                int hiind = (i+1)*detrendlen;
                for(j=loind; j<hiind; j++)
                {
                    timeseries[j] = 0.0f;
                }
            }
            else
            {
                // Now normalize all of the data and reshape it to 1-D
                for(j=0; j<detrendlen; j++)
                    timeseries[j+i*detrendlen] /= stds[i];
            }
        }

        // printf("    identified %d bad blocks out of %d (i.e. %.2f%%)\n", bad_blocks_size, numblocks, 100.0*bad_blocks_size/numblocks);
    }

    // printf("  Now searching...\n");

    int chunknum;
    int chunknumPchunklen;
    int loind;
    // int hiind;
    int lowblock;
    int localgoodblocks_len;

    float *fftd_chunk = malloc(sizeof(float)*worklen);
    float *chunk = malloc(sizeof(float)*worklen);
    float *prod = malloc(sizeof(float)*worklen);
    int *hibins = malloc(sizeof(int)*chunklen);
    int *hiblocks = malloc(sizeof(int)*chunklen);
    float *hivals = malloc(sizeof(float)*chunklen);
    int *currentblocks = malloc(sizeof(int)*blocks_per_chunk);

    CandidateList  dm_candlist;
    dm_candlist.candidates = NULL;
    dm_candlist.count = 0;

    int off;

    for(chunknum = 0; chunknum < numchunks; chunknum++)
    {
        loind = chunknum*chunklen-overlap;
        // hiind = (chunknum+1)*chunklen+overlap;
        chunknumPchunklen = chunknum*chunklen;


        // Take care of beginning and end of file overlap issues
        if(chunknum == 0)  // Beginning of file
        {
            memcpy(chunk+overlap, timeseries+loind+overlap, sizeof(float)*(chunklen+overlap));
            for(i=0; i<overlap; i++)
                chunk[i] = 0.0f;
        }
        else if(chunknum == numchunks-1)  // end of the timeseries
        {
            memcpy(chunk, timeseries+loind, sizeof(float)*(chunklen+overlap));
            for(i=worklen-overlap; i<worklen; i++)
                chunk[i] = 0.0f;
        }
        else
        {
            memcpy(chunk, timeseries+loind, sizeof(float)*worklen);
        }

        // Make a set with the current block numbers
        lowblock = blocks_per_chunk * chunknum;
        for(i=0; i<blocks_per_chunk; i++)
            currentblocks[i] = i+lowblock;
        int found ;
        localgoodblocks_len = 0;
        for(i=0; i<blocks_per_chunk; i++)
        {
            found = 0;
            for(j=0; j<bad_blocks_size; j++)
                if(currentblocks[i] == bad_blocks[j])
                {
                    found = 1;
                    break;
                }
            if(!found)
                localgoodblocks_len++;
        }
        

        // Search this chunk if it is not all bad
        if(localgoodblocks_len)
        {
            // memcpy(goodchunk, chunk+overlap, sizeof(float)*chunklen);
            
            // need to pass blocks/chunklen, localgoodblocks
            // dm_candlist, dt, opts.threshold to cython routine

            // Search non-downsampled data first
            // NOTE:  these nonzero() calls are some of the most
            //        expensive calls in the program.  Best bet would 
            //        probably be to simply iterate over the goodchunk
            //        in C and append to the candlist there.
            int hibins_counts = 0;
            int hibins_counts_thisDown;
            
            for(i=0; i<chunklen; i++)
            {
                if(chunk[i+overlap] > opt_threshold)
                {
                    hibins[hibins_counts] = i + chunknumPchunklen;
                    hiblocks[hibins_counts] = hibins[hibins_counts]/detrendlen;
                    hivals[hibins_counts] = chunk[i+overlap];
                    hibins_counts++;
                }
            }

            // Add the candidates (which are sorted by bin)
            int bin, block;
            int blockIsgood;
            float val;
            
            for(i=0; i<hibins_counts; i++)
            {
                bin = hibins[i];
                block = hiblocks[i];
                val = hivals[i];
                blockIsgood = 1;

                for(j=0; j<bad_blocks_size; j++)
                {
                    if(block == bad_blocks[j])
                    {
                        blockIsgood = 0;
                        break;
                    }
                }
                if(blockIsgood)
                    addCandidate(&dm_candlist, info.dm, val, bin*dt, bin, 1);
            }

            // Prepare our data for the convolution
            if(useffts)
            {
                memcpy(fftd_chunk, chunk, sizeof(float)*worklen);
                realfft(fftd_chunk, worklen, -1);
            }


            // Now do the downsampling...
            for(ii=0; ii< downfacts_N; ii++)
            {
                downfact = downfacts[ii];
                if(useffts)
                {
                    fft_convolve(fftd_chunk, fftd_kerns_f+ii*fftlen, prod, fftlen);
                    // memcpy(goodchunk, prod+overlap, sizeof(float)*chunklen);
                }
                // else
                // { 

                // }

                hibins_counts_thisDown=0;
                for(j=0; j<chunklen; j++)
                {
                    if(prod[j+overlap] > opt_threshold)
                    {
                        hibins[hibins_counts_thisDown] = j + chunknumPchunklen;
                        hiblocks[hibins_counts_thisDown] = hibins[hibins_counts_thisDown]/detrendlen;
                        hivals[hibins_counts_thisDown] = prod[j+overlap];
                        hibins_counts_thisDown++;
                    }
                }
                
                // Now walk through the new candidates and remove those
                // that are not the highest but are within downfact/2
                // bins of a higher signal pulse
                prune_related1(hibins, hivals, hibins_counts_thisDown, downfact, &off);
                hibins_counts_thisDown -= off;
                // Insert the new candidates into the candlist, but
                // keep it sorted...
                for(j=0; j<hibins_counts_thisDown; j++)
                {
                    bin = hibins[j];
                    block = hiblocks[j];
                    val = hivals[j];
                    blockIsgood = 1;

                    for(k=0; k<bad_blocks_size; k++)
                    {
                        if(block == bad_blocks[k])
                        {
                            blockIsgood = 0;
                            break;
                        }
                    }
                    if(blockIsgood)
                        addCandidate(&dm_candlist, info.dm, val, bin*dt, bin, downfact);
                }
            }
        }
    }

    

    // Now walk through the dm_candlist and remove the ones that
    // are within the downsample proximity of a higher
    // signal-to-noise pulse
    qsort(dm_candlist.candidates, dm_candlist.count, sizeof(Candidate), cmp);
    prune_related2(dm_candlist.candidates, dm_candlist.count, downfacts, downfacts_N, &off );
    dm_candlist.count -= off;
    // printf("  Found %d pulse candidates\n",dm_candlist.count);

    // Get rid of those near padding regions
    if(info.breaks)
    {
        prune_border_cases(dm_candlist.candidates, dm_candlist.count, offregions, &off);
        dm_candlist.count -= off;
    }


    //  Write the pulses to an ASCII output file
    FP = fopen(outfile, "w");
    fprintf(FP, "# DM      Sigma      Time (s)     Sample    Downfact\n");
    Candidate cand;
    for(i=0; i<dm_candlist.count; i++)
    {
        cand = dm_candlist.candidates[i];
        // dm_candlist.sort(cmp_sigma)
        fprintf(FP, "  %.2f    %.2f    %.6f    %d    %d\n", cand.DM, cand.val, cand.time, cand.bin, cand.downfact);
    }
    fclose(FP);





    free(fftd_kerns_f);
    free(timeseries);
    free(stds);
    free(bad_blocks);
    free(sort_stds);

    free(fftd_chunk);
    free(chunk);
    free(prod);
    free(hibins);
    free(hiblocks);
    free(hivals);
    free(currentblocks);
    free(dm_candlist.candidates);
    

    return 0;
}

