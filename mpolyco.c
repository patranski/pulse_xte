#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#define MPROG "MPOLYCO: "

int verb1,verb2,verb3,verb4;

void mpolyco(char *dfname,char *psrname,long int imjd_m,double fmjd_m,char *nsite,int *nspan,int *ncoeff)

{
  FILE     *tzfile,*datefile;
  
  char     unique[32],tzin[16],date[16],mvpoly[64],genpoly[64];
  
  int      maxha; 
  double   freq;
  
  /* Get these from catalogue at a later date */
  *ncoeff = 15;     //#  Number of coefficient in the polynomial expansion
  *nspan = 120;     //#  Number of minutes for tempo polynomial
  maxha = 12;       //#  Maximum hour angle for observations (hr)
  freq = 1408.0;    //#  Default observing frequency
  
  /* First we need to make the tz.in file */
  
  strcpy(tzin,psrname);
  strcat(tzin,".tz");
  tzfile = fopen(tzin,"w");

  fprintf(tzfile,"    %s   12 180  15 1408.      (Defaults for nsite, maxha, nspan, ncoeff, freq)\n",nsite);
  fprintf(tzfile,"  Name     Nspan  Ncoeffs Maxha Freq (For each PSR you can override the defaults)\n");
  fprintf(tzfile,"--------------------------------------------------------------------------------\n");
  fprintf(tzfile,"%s %d %d %d %f\n",psrname,*nspan,*ncoeff,maxha,freq); 
  fclose(tzfile);
  //  printf("ncoeff = %d\n",*ncoeff);
  /* Okay now we need a file to put the date in */
  
 
  strcpy(date,psrname);
  strcat(date,".dt");
  datefile = fopen(date,"w");
  fprintf(datefile,"%lf %lf\n",imjd_m+fmjd_m-0.5,imjd_m+fmjd_m+0.5);
  printf("Mpolyco range: %lf %lf\n",imjd_m+fmjd_m-0.5,imjd_m+fmjd_m+0.5);
  fclose(datefile);

  /* Actually generate the polyco */
  
  sprintf(genpoly,"tempo -z %s < %s >> /tmp/%s",tzin,date,dfname);
  system(genpoly);

  /* Move the polyco to a unique name */

  sprintf(mvpoly,"mv polyco.dat %s.polyco\n",dfname);
  system(mvpoly);
}








