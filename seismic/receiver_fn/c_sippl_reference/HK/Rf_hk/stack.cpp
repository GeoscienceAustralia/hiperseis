#include "recordrf.h"
#include "nrutil.h"

/******************************************************************************/
/* Calculate travel time, pick up amplitude
   Stacking all stations                                                      */
/******************************************************************************/
void nth_stk(RECORD *rfrec, HKS_INFO *hks)
{
    int izr;
    int irec;
    int i;
    double tmpx1, tmpx2, tmpx3;
                
/*----------------------------------------------------------------------------*/
/* linear stacking                                                            */
/*----------------------------------------------------------------------------*/
    if (hks->nroot == 1) {  
        for (irec = 0; irec < hks->nrfrec; irec++) {
            rfrec[irec].LoadData();
            printf("Processing rfrec[%d/%d]...\n", irec, hks->nrfrec);

            for (izr = 0; izr < hks->nzr; izr++) {
                rfrec[irec].Ttps(hks, izr);  //Note how the mh and rv sort into 1-D array
                    
                rfrec[irec].PeakAmp(hks);

                //if (rfrec[irec].a_0p1s == 0.0) printf("0p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
                //if (rfrec[irec].a_2p1s == 0.0) printf("2p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
                //if (rfrec[irec].a_1p2s == 0.0) printf("1p2s: irec = %d, izr = %d, == 0.0\n", irec, izr);
                    
                hks->a_0p1s[izr] += rfrec[irec].a_0p1s;
                hks->a_2p1s[izr] += rfrec[irec].a_2p1s;
                hks->a_1p2s[izr] += rfrec[irec].a_1p2s;
            }
            rfrec[irec].FreeData();
        }
        // Sign on a_1p2s here might be incorrect, it is supposed to be negated, e.g. see Chen et al. (2010)
        for (izr = 0; izr < hks->nzr; izr++) {
            hks->a_0p1s[izr] /= hks->nrfrec;
            hks->a_2p1s[izr] /= hks->nrfrec;
            hks->a_1p2s[izr] /= hks->nrfrec;
        }
    }
/*----------------------------------------------------------------------------*/
/* nth-root stacking n == 2                                                   */
/*----------------------------------------------------------------------------*/
    else if (hks->nroot == 2) {
        // N-th root stacking (N = 2)
        for (irec = 0; irec < hks->nrfrec; irec++) {
            rfrec[irec].LoadData();
            printf("Processing rfrec[%d/%d]...\n", irec, hks->nrfrec);

            for (izr = 0; izr < hks->nzr; izr++) {
                rfrec[irec].Ttps(hks, izr);  //Note how the mh and rv sort into 1-D array

                rfrec[irec].PeakAmp(hks);

                //if (rfrec[irec].a_0p1s == 0.0) printf("0p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
                //if (rfrec[irec].a_2p1s == 0.0) printf("2p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
                //if (rfrec[irec].a_1p2s == 0.0) printf("1p2s: irec = %d, izr = %d, == 0.0\n", irec, izr);
                    
                tmpx1 = SIGN1(rfrec[irec].a_0p1s); 
                tmpx2 = SIGN1(rfrec[irec].a_2p1s); 
                tmpx3 = SIGN1(rfrec[irec].a_1p2s); 
    
                hks->a_0p1s[izr] += tmpx1 * sqrt(fabs(rfrec[irec].a_0p1s));
                hks->a_2p1s[izr] += tmpx2 * sqrt(fabs(rfrec[irec].a_2p1s));
                hks->a_1p2s[izr] += tmpx3 * sqrt(fabs(rfrec[irec].a_1p2s));
            }
            rfrec[irec].FreeData();
        }
        // This loop computes the mean, then raises each phase component to the power of nroot
        // as per Muirhead (1968).
        for (izr = 0; izr < hks->nzr; izr++) {
                    
            hks->a_0p1s[izr] /= hks->nrfrec;
            hks->a_2p1s[izr] /= hks->nrfrec;
            hks->a_1p2s[izr] /= hks->nrfrec;
    
            tmpx1 = 1.0;
            tmpx2 = 1.0;
            tmpx3 = 1.0;
            for (i = 0; i < hks->nroot; i++)    {
                tmpx1 *= fabs(hks->a_0p1s[izr]);
                tmpx2 *= fabs(hks->a_2p1s[izr]);
                tmpx3 *= fabs(hks->a_1p2s[izr]);
            }   

            hks->a_0p1s[izr] =  SIGN1(hks->a_0p1s[izr]) * tmpx1;
            hks->a_2p1s[izr] =  SIGN1(hks->a_2p1s[izr]) * tmpx2;
            hks->a_1p2s[izr] = -SIGN1(hks->a_1p2s[izr]) * tmpx3;
        }
    }
/*----------------------------------------------------------------------------*/
/* nth-root stacking n == 4                                                   */
/*----------------------------------------------------------------------------*/
    else if (hks->nroot == 4) {  
        for (irec = 0; irec < hks->nrfrec; irec++) {
            rfrec[irec].LoadData();
            printf("Processing rfrec[%d/%d]...\n", irec, hks->nrfrec);

            for (izr = 0; izr < hks->nzr; izr++) {
                //if (izr % (hks->nzr/10) == 0) {
                //  printf("%d-root Stacking: %d, iz = %d, ir = %d\n", 
                //      hks->nroot, izr, izr / hks->nrv, izr % hks->nrv);
                //}
                            
                rfrec[irec].Ttps(hks, izr);  //Note how the mh and rv sort into 1-D array
                rfrec[irec].PeakAmp(hks);

                //if (rfrec[irec].a_0p1s == 0.0) printf("0p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
                //if (rfrec[irec].a_2p1s == 0.0) printf("2p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
                //if (rfrec[irec].a_1p2s == 0.0) printf("1p2s: irec = %d, izr = %d, == 0.0\n", irec, izr);
                    
                tmpx1 = SIGN1(rfrec[irec].a_0p1s); 
                tmpx2 = SIGN1(rfrec[irec].a_2p1s); 
                tmpx3 = SIGN1(rfrec[irec].a_1p2s); 
    
                hks->a_0p1s[izr] += tmpx1 * sqrt(sqrt(fabs(rfrec[irec].a_0p1s)));
                hks->a_2p1s[izr] += tmpx2 * sqrt(sqrt(fabs(rfrec[irec].a_2p1s)));
                hks->a_1p2s[izr] += tmpx3 * sqrt(sqrt(fabs(rfrec[irec].a_1p2s)));
            }
            rfrec[irec].FreeData();
        }
            
        for (izr = 0; izr < hks->nzr; izr++) {
                    
            hks->a_0p1s[izr] /= hks->nrfrec;
            hks->a_2p1s[izr] /= hks->nrfrec;
            hks->a_1p2s[izr] /= hks->nrfrec;
    
            tmpx1 = 1.0;
            tmpx2 = 1.0;
            tmpx3 = 1.0;
            for (i = 0; i < hks->nroot; i++)    {
                tmpx1 *= fabs(hks->a_0p1s[izr]);
                tmpx2 *= fabs(hks->a_2p1s[izr]);
                tmpx3 *= fabs(hks->a_1p2s[izr]);
            }   

            hks->a_0p1s[izr] =  SIGN1(hks->a_0p1s[izr]) * tmpx1;
            hks->a_2p1s[izr] =  SIGN1(hks->a_2p1s[izr]) * tmpx2;
            hks->a_1p2s[izr] = -SIGN1(hks->a_1p2s[izr]) * tmpx3;
        }
    }
/*----------------------------------------------------------------------------*/
/* nth-root stacking                                                          */
/*----------------------------------------------------------------------------*/
    else {  
        for (irec = 0; irec < hks->nrfrec; irec++) {
            rfrec[irec].LoadData();
            printf("Processing rfrec[%d/%d]...\n", irec, hks->nrfrec);

            for (izr = 0; izr < hks->nzr; izr++) {
                rfrec[irec].Ttps(hks, izr);  //Note how the mh and rv sort into 1-D array

                rfrec[irec].PeakAmp(hks);

                //if (rfrec[irec].a_0p1s == 0.0) printf("0p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
                //if (rfrec[irec].a_2p1s == 0.0) printf("2p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
                //if (rfrec[irec].a_1p2s == 0.0) printf("1p2s: irec = %d, izr = %d, == 0.0\n", irec, izr);
                    
                tmpx1 = SIGN1(rfrec[irec].a_0p1s); 
                tmpx2 = SIGN1(rfrec[irec].a_2p1s); 
                tmpx3 = SIGN1(rfrec[irec].a_1p2s); 
    
                hks->a_0p1s[izr] += tmpx1 * pow(fabs(rfrec[irec].a_0p1s),  hks->nrootpw);
                hks->a_2p1s[izr] += tmpx2 * pow(fabs(rfrec[irec].a_2p1s),  hks->nrootpw);
                hks->a_1p2s[izr] += tmpx3 * pow(fabs(rfrec[irec].a_1p2s),  hks->nrootpw);
            }
            rfrec[irec].FreeData();
        }
            
        for (izr = 0; izr < hks->nzr; izr++) {
                    
            hks->a_0p1s[izr] /= hks->nrfrec;
            hks->a_2p1s[izr] /= hks->nrfrec;
            hks->a_1p2s[izr] /= hks->nrfrec;
    
            tmpx1 = 1.0;
            tmpx2 = 1.0;
            tmpx3 = 1.0;
            for (i = 0; i < hks->nroot; i++)    {
                tmpx1 *= fabs(hks->a_0p1s[izr]);
                tmpx2 *= fabs(hks->a_2p1s[izr]);
                tmpx3 *= fabs(hks->a_1p2s[izr]);
            }   

            hks->a_0p1s[izr] =  SIGN1(hks->a_0p1s[izr]) * tmpx1;
            hks->a_2p1s[izr] =  SIGN1(hks->a_2p1s[izr]) * tmpx2;
            hks->a_1p2s[izr] = -SIGN1(hks->a_1p2s[izr]) * tmpx3;
        }
    }

    return;
}   

