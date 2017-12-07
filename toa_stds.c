#include <stdio.h>
#include <math.h>
    
void toa_stds(FILE *toaptr_stds,char *fname,double *stdamps,int stdnbin,double *prfamps,double prfnbin,int imjdmid,double fmjdmid,double pfold,double freq, double *terr, double *maxbin2)
{
  
  double phase,shift,eshift,ephase,snrfft,esnrfft,maxbin[1],shift2,rms,b;
  int idmid, j, i, ishift, smax;
  double dmjd;
  
  fflush(stdout);
  
  fftconv_(&stdnbin,prfamps,stdamps,&shift,&eshift,&snrfft,&esnrfft,&b,&rms); 
  printf("shift %f\n",shift);
  if (shift < 0.0) (double)(shift += 1.0*stdnbin);
  phase = (double)(shift/(1.0*stdnbin));
  printf("shift %f\n",shift);
  ephase = eshift / (double)stdnbin;
  
  /* Referring the mjd to the position of the standard (which should have its peak at 0) */
  idmid = imjdmid;
  dmjd  = (double) fmjdmid + (double)(phase * pfold/(86400.0));
  *terr  = ephase * pfold * 1000.0 * 1000.0;
  if (dmjd > 1.0) {idmid++; dmjd-=1.0;}
  if (dmjd < 0.0) {idmid--; dmjd+=1.0;}
  
  return;
  fflush(stdout);
}
