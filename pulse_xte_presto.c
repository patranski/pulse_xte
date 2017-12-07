#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <malloc.h>
#include <time.h>
#include <stddef.h>
#include "matrix.h"
#include <stdlib.h>
#include <ctype.h>
#include "fitsio.h"
#include "get_channels.h"
#include "choose_channels.h"
    
#define PROG "FXRAY: "
    
#define nfiles 1
#define nmetafiles 1
#define max_metafiles 1
#define max_channels 256
#define max_files 1
#define NUM_FILES 5000
#define pi 3.141592654
    
#define error "Memory Overflow\n"
    
    
int ioffset[max_files][max_metafiles];
int ibit_chan[max_files][max_metafiles];
int channels[max_channels][max_files][max_metafiles];
int channels1[max_channels][max_files][max_metafiles];
int channels2[max_channels][max_files][max_metafiles];
int nchannels[max_files][max_metafiles];
int nhistos[max_files][max_metafiles];
    
int column,firstcall,tformat;          //Dump subroutine global variables
char metafiles[80],*infilename;        // global variable for input fits filenames    

void printerror_C( int status);
double *selectrows(double tsamp, int column, int incl_excl,int numfiles);  //this is the subroutine to read and filter TOAs in FITS files 
void profiles(long int imjd_m, double fmjd_m, int nout,int counter);
int *toa(FILE *toaptr,FILE *toatempo2,char *dfname,double *count1,int nbin,double *prfamps, int simjd_m, double sfmjd_m, double spobs_m ,double freq, double *terr,char *site,double EMAX, double maxbin, double *deltph,double *refph1);
void toa_stds(FILE *toaptr_stds,char *fname,double *stdamps,int stdnbin,double *prfamps,double prfnbin,int imjdmid,double fmjdmid,double pfold,double freq, double *terr, double *maxbin2);
int dump(int column, int tformat);
int pure_sin(void);
void timelags(void);
float gammq(float aa, float xx);
int check_pcus(double counts_pcu[]);
void prepare_info(int count_info);
void mpolyco(char *dfname,char *psrname,long int imjd_m,double fmjd_m,char *nsite,int *nspan,int *ncoeff);
void ppolyco(char *dfname,long int imjd_m,double fmjd_m,double *pobs,double *refph);

void usage(int status)   
{
  printf("ERROR IN INPUT COMMAND LINE !\n");
  printf("%sCommand line must be: pulse_xte_presto filename\n",PROG);
  printf("%sRemember to put the names of your FITS file in metafile (metafiles are defined as @filename)\n",PROG);
}
    
// Global variables to store the baricentric data from fits files.
double *barytime, counts_pcu[6], *bary_f, var,rms;
int *pcubary, prof_act_pcus, first_toa_call;
double *bary_init_i,*bary_init_f, *bary_fin_i,*bary_fin_f, *time_col;
double timeoffset,*bary_i,EMAX;      
double prfmean,prob,offsetph;
    
double *ph,*deltaph,*refph1,refph2;
int nh, profi, num_pcu, good_pcu[6], pcu_switch, intdph, counterprofile;
int switch_bary, presto_switch, even_bin;
int count_info;
char pcu_choice[20], root[6];
int lcount,nbits,num,chunk,incl_excl,finalpcu,pcuoffset;                    
int noutrows,numfiles, NN,refchan,maxbinprf;
long int imjd_m;
char mjd_s[32];
double *stdamps,*standprof,*prfamps,*terr,*count1,*count0,*totp, gti_i[100000],gti_f[100000],toat[100000],totcs;
double gti_mjd_i,gti_mjd_f,bkgcount[1000],chan[257],dchan[257], bkgshift; 
int numbr,i,j,numb,offon,nbin,ndmp_o,counter,gticount,lines,init, counter_ph, counter_ph_offset;
long int simjd_m;
double sfmjd_m,spobs_m,freq,refph,dphase,phs,fmjd_m,tdmp_o,tdmp_2,store,tsamp[1],toatonoff, delta_pixel, sampling;
double *count2,*maxbin2,standard[200],smax,smaxprf,maxbin;
char site[2],dfname[32], buf[2000], gtiyesno, bkgyesno,toasyn;
FILE *toaptr, *toatempo2,*ndmp, *indivpulse,*totalprof, *indivphase, *indivnorm;
FILE *noshift, *phase,*stdprof,*stdprof2,*timeint,*gtifile, *bkgfile,*toas, *indiv_std, *lc;
FILE *gti_mjd,*test,*intervals,*pcus,*totcounts,*amplitudes, *lightcurves;
/*  New section with fortran definitions to avoid  memory leaking */
struct complex {float real, imag; };
double time_r_i, time_r_f, bintime, rimjd_f, spobs, pobs, fmjd, sph, dph,frequency;
float *opr, ict;
int indx, ib, binct, time_ct;
double time_bell;
long int imjd, rimjd_i;
int    itoa, chkdct, chkct, switch_chan, totnum;
extern void shiftbyfft_(double *,int *,double *);
extern void fftconv_(int *,double *,double *,double *,double *, double *, double *, double *, double *);
extern void cprof_(double *,int *,int *, struct complex *, double *, double *);
extern void ffft_(struct complex *, int *, int *, int *);
extern void fftfit_(double *, double *, double *, int *, double *, double *, double *, double *);
extern void fccf_(double *, double *, double *);
/*                      End of Fortran definitions                               */

/*REBINNING VARIABLES*/
int ibin, ibin_old;
double np;  
int switch_counter, n;
double tbin, tbin_i, tbin_f;
float check_cts;
int  check_bins;
char *prestofile,strfile[3], strfile2[3];
FILE *presto_out, *standa;

int main(int argc, char **argv)
{
    
  fitsfile *infptr;
  int status,hdutype,nfound, nkeys, chansel, foldvar;
  int c,k,kk,nreads,nsamp_o,operation;
  int ncts_r,iopt,nspan_o,ncoeff_o, l;
  long nn;
  
  double time_initial_i, time_initial_f, time_final_i, time_final_f, fix_time_i, fix_time_f; 
  float *rdata;
  int *phpos;

  double limit, step_limit, minbin,timemjd;

  double time,period_s;
  char   prfname[40],*psrname,*string;
  double *x, *p, array[200];
  int    nout,selec,std,ll;
  char toaname[40], toaname2[40],*s;
  FILE *prfptr,*toaptr_stds, *fp,*data, *dataunbary;
  char metafilename[]="meta";
  char yesno;
  
    
  /* Initialization of variables for channels */
  int nnfiles =1;
  int nnmetafiles= 1;
  int nmax_metafiles= 1;
  int nmax_channels= 256;
  int nmax_files= 1;
  int istart,iend;
    

  /* Inizialization of output files*/

  amplitudes=fopen("amplitudes.dat","w");
  fclose(amplitudes);
  test=fopen("test.dat","w");
  fclose(test);
  intervals=fopen("intervals.dat","w");
  fclose(intervals);
  psrname = (char *) calloc(16,1);
  printf("Allocating psrname\n");
  data=fopen("data.dat","w");
  dataunbary=fopen("dataunbary.dat","w");
  fclose(dataunbary);
  fclose(data);
  indiv_std=fopen("indiv_std.dat","w");
  fclose(indiv_std);
  pcus=fopen("pcus.dat","w");
  fclose(pcus);
  lightcurves=fopen("lightcurves.dat","w");
  fclose(lightcurves);    
  timeint=fopen("timeint.dat","w");
  ndmp=fopen("ndmp.dat","w");
  indivpulse=fopen("indivpulse.dat","w");
  indivphase=fopen("indivphase.dat","w");
  indivnorm=fopen("indivnorm.dat","w");
  totalprof =fopen("totalprofile.dat", "w");
  noshift = fopen("noshift.dat","w");
  stdprof=fopen("stdprof.dat","w");
  stdprof2=fopen("stdprof2.dat","w");
  fclose(indivpulse);
  fclose(indivphase);
  fclose(indivnorm);
  fclose(totalprof);
  fclose(noshift);
  fclose(stdprof);
  fclose(stdprof2);
  phase=fopen("phase.dat","w");
  fclose(phase);
  gti_mjd=fopen("gti_mjd.dat","w");
  totcounts=fopen("totcounts.dat","w");
  fclose(totcounts);
    
  /* Initialization of some other variables */
  for(i=0; i<6; i++)
    root[i]=0;
    
  first_toa_call = 1;
  std=0;
  itoa = 0;
  freq = 9999.999;
  status=0;
  nfound=0;
  chansel=1;
  lcount=1;
  foldvar=0;
  offon=0;
  numb=0;
  num=0;
  gticount=0;
  bary_init_i = (double *) malloc(sizeof(double));
  bary_init_f = (double *) malloc(sizeof(double));
  bary_fin_i = (double *) malloc(sizeof(double));
  bary_fin_f = (double *) malloc(sizeof(double));
    

  /* Define usage of shell command line */
    
  if (argc <=1)
    {
      usage(0);
    }
    
  for (iopt = optind; iopt < argc; iopt++)
    {
    
    
      printf("    BINARY PULSAR SELECTION\n\n");
      printf("    [1] SAX J1808.4-3658\n");
      printf("    [2] IGR J0291+5934\n");
      printf("    [3] XTE J0929-314\n");
      printf("    [4] XTE J1751-305\n");
      printf("    [5] XTE J1807-294\n");
      printf("    [6] HETE J1900-2455\n");
      printf("    [7] XTE J1814-338\n");
      printf("    [8] Dump ToAs in an ASCII file (.dmp)\n");
      printf("    [9] SAX J1748-2021\n");
      printf("    [10] User defined (1000-0000.par)\n");
      printf("    [11] Swift 1756.9-2508\n");
      
      scanf("%d", &selec);
      strcpy(dfname,argv[iopt]);
      
      switch(selec)
        {
        case 1:
          strcpy(psrname,"J1808-3658\0");
          break;
        case 2:
          strcpy(psrname,"J0291+5934\0");
          break;
        case 3:
          strcpy(psrname,"J0929-314\0");
          break;
        case 4:
          strcpy(psrname,"J1751-305\0");
          break;
        case 5:
          strcpy(psrname,"J1807-294\0");
          break;
        case 6:
          strcpy(psrname,"J1900-2455\0");
          break;
        case 7:
          strcpy(psrname,"J1814-338\0");
          break;
        case 8:
          firstcall=1;
          break;
        case 9:
          strcpy(psrname,"J1748-2021\0");
          break;
        case 10:
          strcpy(psrname,"J1000-0000\0");
          break;
        case 11:
          strcpy(psrname,"J1756-2508\0");
          break;
        }
      if(selec!=8)
        {
          printf("Binary selected:  %s\n\n\n",psrname);
          
         
          /* Stripping the pulsar name */
          if (strlen(psrname) > 10 )
            strncpy(psrname,psrname+5,strlen(psrname)-4);
          else
            {
              if ( strncmp("J",psrname,1) == 0 || strncmp("B",psrname,1) == 0 )
                strncpy(psrname,psrname+1,strlen(psrname));
              else
                strncpy(psrname,psrname,strlen(psrname)+1);
            }
          
          printf("---> %s\n",psrname);
          strcpy(site,"@");   // @ for barycentering, 0 for geocentering
          
          /* Opening a few ascii files */
          strcpy(toaname,dfname);
          strcat(toaname,".tim");
          toaptr = fopen(toaname,"a");
          strcpy(toaname2,"tempo2");
          strcat(toaname2,".tim");
          toaptr = fopen(toaname,"a");
          toatempo2 = fopen(toaname2,"a");
          strcpy(prfname,dfname);
          strcat(prfname,".prf");
          prfptr = fopen(prfname,"a");
          amplitudes=fopen("amplitudes.dat","a");
          lightcurves= fopen("lightcurves.dat","a");


      /*      NEW PROCEDURE TO READ FILES       */
      /* USES A GLOBAL VARIABLE infilename[80]. */
         

	  printf("Use Barycentered times (0. No | 1. Yes) ?\n");
	  scanf("%d",&switch_bary);
	  printf("Use PRESTO mode ? (create time series and info files for PRESTO)\n");
	  printf("0. No  | 1. Yes\n");
	  scanf("%d",&presto_switch);

          nbin = 32;                           // num. of bin in each chunk-profile.
          nh=nbin/2;

	  /* count1 is the array where the counts of the standard profile */
	  /* (standard.st) are stored                                     */

          count1 = (double *) calloc(nbin,sizeof(double));
	  totp = (double *) calloc(nbin,sizeof(double));
	  for(i=0; i<nbin; i++)
	    {count1[i]=0; totp[i]=0;}

	  printf("Allocating count1 and totp\n");

          /* Selecting/Generating  a new standard profile: */
          printf("Standard profile: if you need a new standard profile\n");
          printf("the program will build one using all\n");
          printf("the fits files provided in your metafile\n\n\n");
          printf("If you are analyzing only a subpart of data, please build the standard profile\n");
          printf("using all the available data, then launch again the program and use your new high S/N profile\n\n");
	  printf("Do you need to build a new standard profile ? (y/n)\n");
          scanf("%s",&yesno);
	  yesno=tolower(yesno);
	  printf("You said: %c\n",yesno);
          
          switch(yesno)
            {
            case 'y': 
              offon=1;
              printf ("OFFON =%d\n",offon);
              printf("Would you like to build:\n");
              printf("1. A sequence of profiles (one for each FITS file) \n");
              printf("2. Use a pure sinusoid at a fixed frequency\n");
              scanf("%d",&std);
              printf("You said: %d\n",std);
              if(std==2)
                {
                  standa=fopen("standard.st","w");
                  pure_sin();
                  fclose(standa);
                  exit(1);
                }
  

	      /* Select GTI's filter file (gti.dat in ASCII format)*/
	      /* First column: Initial GT; Second column: Final GT;*/

	      printf("Do you want to use a GTI file (gti.dat) ?\n");
              scanf("%s",&gtiyesno);
	      gtiyesno=tolower(gtiyesno);
              printf("You said: %c\n",gtiyesno);
              printf("Do you want to exclude (0) or include (1) those gtis ?\n");
              scanf("%d",&incl_excl);
              if(gtiyesno=='y')
                {
                  gtifile=fopen("gti.dat","r");
                  
                  for(;;)
                    {   
                      if(gticount>1 && gti_i[gticount]==gti_i[gticount-1])
                        break;
                      fscanf(gtifile,"%lf %lf",&gti_i[gticount],&gti_f[gticount]);
                      printf("GTI INTERVALS: %lf %lf\n",gti_i[gticount],gti_f[gticount]);
                      fprintf(gti_mjd,"%16.11lf %16.11lf\n",gti_i[gticount]/86400.+49353.00073567628,gti_f[gticount]/86400.+49353.00073567628);
		      fflush(gti_mjd);
                      gticount++;
                    }
                }
              break;
            case 'n':
              offon=0;
              printf ("OFFON =%d\n",offon);  
              printf("Do you want to use a GTI file (gti.dat) ?\n");
              scanf("%s",&gtiyesno);
	      gtiyesno=tolower(gtiyesno);
              printf("You said: %c\n",gtiyesno);
              printf("Do you want to exclude (0) or include (1) those gtis ?\n");
              scanf("%d",&incl_excl);
              if(gtiyesno=='y')
                {
                  gtifile=fopen("gti.dat","r");
                  for(;;)
                    { 
                      if(gticount>1 && gti_i[gticount]==gti_i[gticount-1])
                        break;
                      fscanf(gtifile,"%lf %lf",&gti_i[gticount],&gti_f[gticount]);
                      printf("GTI INTERVALS: %lf %lf\n",gti_i[gticount],gti_f[gticount]);
                       fprintf(gti_mjd,"%16.11lf %16.11lf\n",gti_i[gticount]/86400.+49353.00073567628,gti_f[gticount]/86400.+49353.00073567628);
		       fflush(gti_mjd);
                      gticount++;
                    }
                }
    
              break;
            }
	  /* Select BKG file (bkg.dat) with format Ct/S/PCU  */
	  /* IMPORTANT: the BKG values need to be sorted according */
	  /* to the INPUT FITS file names */
	  /* E.G., if you have file1.fits and file2.fits in your @meta file*/
	  /* then you need to have bkg1 and bkg2 in the bkg.dat file */
          printf("Do you want to use a BKG file (bkg.dat) ?\n");
          scanf("%s",&bkgyesno);
	  bkgyesno=tolower(bkgyesno);
          printf("You said: %c\n",bkgyesno);
          lines=0;
          if(bkgyesno=='y')
            {
              bkgfile=fopen("bkg.dat","r");
	      if(bkgfile==NULL)
		{
		  printf("File bkg.dat does not exist !!!\n");
		  printf("Error. Create background file and restart.\n");
		  exit(1);
		}
              for(;;)
                {
                  string=fgets(buf,1500,bkgfile);
                  if(string==NULL)
                    {fclose(bkgfile);break;}
                  lines++;
                }
    
              bkgfile=fopen("bkg.dat","r");
              for(i=0;i<lines;i++)
                {
                  fscanf(bkgfile,"%lf",&bkgcount[i]);
                  printf("BKG rate (ct/s) for all active PCUs: %lf\n",bkgcount[i]);
                  printf("Number of BKG files=%d\n",i);
                  printf("Number of lines    =%d\n",lines);
                }
            }
    
    
	  /*Check if the standard.st file exists otherwise exit*/
          if(offon==0 || offon==1 && std==2)
            {
              printf("Using the standard profile stored in 'standard.st'\n\n\n ");
              standa = fopen("standard.st","r");
              if(standa == NULL)
                {
                  printf("NO VALID 'standard.st' file present !!! \n\n\n");
                  exit(1);
                }
              for(j=0;j<nbin;j++)
                {
                  fscanf(standa,"%d %lf\n",&j,&standard[j]);
                }
            }
        }


      /* Select TOAs file option: give an arbitrary list of time intervals*/
      /* that you want to use in the folding process, which are subintervals*/
      /* of the GTIs. */
      printf("Do you want to use a TOAs file ?\n");
      scanf("%s",&toasyn);
      toasyn=tolower(toasyn);
      if(toasyn=='y')
        {
          toatonoff=1.;
          lines=0;
          if((toas = fopen("toas.st", "r")) == NULL)
            { 
              printf("No toas.st file found\n");
              exit(1);
            }
          for(;;)
            {
              string=fgets(buf,1500,toas);
              if(string==NULL)
                {fclose(toas);break;}
              lines++;
            }
    
              toas=fopen("toas.st","r");
              for(i=0;i<lines;i++)
                {
                  fscanf(toas,"%lf",&toat[i]);
                  printf("TOAs: %lf\n",toat[i]);
                  printf("Number of TOA files=%d\n",i);
                  printf("Number of lines    =%d\n",lines);
                }            
        }
      else
        {
          toatonoff=0.;
        }


      /*Select input metafile name with lists of FITS files*/
      printf("Insert metafilename (DEFAULT = @meta)\n");
      scanf("%s",metafilename);
      
      if((fp = fopen(metafilename, "r")) == NULL)
        { 
          printf("Using default metafile\n");
          fp = fopen("@meta","r");
        }      
      else
        {
          printf("Your selected metafile: %s\n",metafilename);
          
        }
            
      /* Select folding option: variable chunk, fixed chunk or two different chunk lengths*/
      /* (Option 2 needs to be further tested) */
      if(selec!=8)
        {
          if(!offon)
            {
              printf("Do you want ti use:\n");
              printf("1. A variable time interval for each chunk\n");
              printf("2. Two different variable time interval for each chunk\n");
              printf("3. A fixed time interval for each chunck\n");
              scanf("%d",&chunk);
              if(chunk==1)
                {
                  printf("What is the time interval for each chunk of data ? (e.g. 128 s, 200 s, etc...)\n");
                  scanf("%lf",&tdmp_o);
                  printf("Time interval =%e\n",tdmp_o);
                  store=tdmp_o;
                  timemjd=1.e+9;
                }
              else if(chunk==2)
                {
                  printf("Where is the MJD where you want to change the time interval ?\n");
                  scanf("%lf",&timemjd);
                  printf("What is the first time interval for each chunk of data ? (e.g. 128 s, 200 s, etc...)\n");
                  scanf("%lf",&tdmp_o);
                  printf("Time interval =%e\n",tdmp_o);
                  store=tdmp_o;
                  printf("What is the second time interval for each chunk of data ? (e.g. 128 s, 200 s, etc...)\n");
                  scanf("%lf",&tdmp_2);
                  
                }
              else if(chunk==3)
                {
                  printf("What is the time interval for each chunk of data ? (e.g. 128 s, 200 s, etc...)\n");
                  scanf("%lf",&tdmp_o);
                  ////      scanf("%lf",&tdmp_o);
                  printf("Time interval =%e\n",tdmp_o);
                  store=tdmp_o;
                  timemjd=1.e+9;
                }
    
            }
          else         // In this case we need a value of tdmp_o just to make the folding and binning. 
            {          // It's always the same, whatever the value you choose (not the same for timelags() of course).
              tdmp_o=3000.;
              store=tdmp_o; 
            }
        }

      /* Selects the MAXIMUM allowed error (in microsec) on the pulses ToAs*/
      /* Retains only profiles with error smaller or equal than EMAX*/
      /* If you don't want to use this option give some very large value in input*/
      /* (something like 1E+9 or so...)*/
    
      printf("Which EMAX you want to use ?\n");
      scanf("%lf",&EMAX);
      printf("EMAX = %lf\n",EMAX);
    

      	  
      /*SELECT PCU NUMBER*/
      pcu_switch=0;     
      printf("Which PCU would you like to use ?\n");
      printf("0 1 2 3 4\n");
      printf("Please sort them in increasing order (e.g., 1,3,4, not 3,1, 4)\n");
      scanf("%s",pcu_choice);
      
      i=0;
      num_pcu=0;
      if(pcu_switch==0)
	for(i=0;i<10;i=i+2)
	  {
		if(pcu_choice[i]!=','&& pcu_choice[i]!='\0')
		  num_pcu++;
	  }
      printf("NUMBER OF SELECTED PCUs=%d\n",num_pcu);
      i=0;
      for(i=0; i<num_pcu; i++)
	{
	  good_pcu[i]=atoi(&pcu_choice[2*i]);
	}
	  
	 
      /////////////////////////////////////////
      /*                                     */
      /*                                     */
      /* HERE IT STARTS TO READ THE METAFILE */
      /* AND TO ANALYZE THE SEQUENTIAL DATA  */
      /*                                     */
      /*                                     */
      /////////////////////////////////////////
    
      printf("STARTING NUMFILES CICLE\n");
      init=0; //switch to read the standard profile in case numfile==0 --> skip (less than 10 rows)
      totcounts=fopen("totcounts.dat","a"); //stores total counts in each profile
    
      for(numfiles=0;numfiles<NUM_FILES;numfiles++) /* Loop for 5000 (potential) fits files. */
                                                    /* If you need more, change the variable NUM_FILES */
	                                            /* and recompile */
        {
          time=0.0;
          tdmp_o=store;
          fflush(fp);
          s=fgets(buf,199,fp);
    
          
         /////////////////////////////////////////////////////////
         //  CHECK END OF ANALYSIS AND OTHER PARs for each cycle//
         /////////////////////////////////////////////////////////
          if(s==NULL && lcount != 0)
            {
              nn=ftell(fp);
              fseek(fp, -nn, 1);
              lcount=0;
              numfiles=0;             
              continue;
            }
          else if(lcount == 0 && foldvar == 1)
            {
	      totalprof=fopen("totalprof.dat","w");
	      for(i=0; i<nbin; i++)
		fprintf(totalprof, "%d %d\n",i, (int)(totp[i]));
	      fclose(totalprof);
	      
              printf("DATA ANALYSIS SUCCESSFUL !\n");
              if(offon==1 && std==1)
                {
		  free(standprof);
                  printf("Now start again pulse_xte_presto with the new standard profile (NOT IDEAL, FIX LATER)\n");
                }
              exit(1);
            }
          else if(lcount == 0 && foldvar == 0)
            {
              foldvar++;
              continue;
            }
          //////////////////////////////////   
    
          
          /* READS FITS FILENAMES AND SELECTS CHANNELS */
    
          nn=ftell(fp);
          chkdct = -1;                        
          i=0;
          
          while(buf[i]!='\n')
	    i++;
	  
	  printf("Allocating *infilename (i=%d)\n",i);
          infilename = (char *) calloc(i,4);
          i=0;
          while(buf[i]!='\n')
            {
              infilename[i]=buf[i];
              i++;
            }
          
          
          printf("Infilename=%s\n",infilename);
          

	  /* CFITSIO SUBROUTINES to manage the I/O process*/
	  /* open the existing FITS files */
          if ( fits_open_file(&infptr, infilename,  READONLY,  &status))
            printerror_C( status );
          /* move to the 2nd HDU in the input file (a binary table in this case) */
          if ( fits_movabs_hdu(infptr, 2, &hdutype, &status))
            printerror_C( status );
          /* read the  time sample keyword */
          if(fits_read_key(infptr, TDOUBLE, "TIMEDEL", tsamp, 0, &status))
            printerror_C( status );
          printf("tsamp=%22.21f\n",tsamp[0]);
          
         
	  /* get_channels.f and choose_channels.f are old FORTRAN77 subrountines that however, work very well.*/	 
          /* A:  Select channels  */
          if(numfiles==0 && chansel == 1)
            {
	      get_channels_(infilename,&channels1[0][0][0],&channels2[0][0][0],&nchannels[0][0],&nhistos[0][0],&ioffset[0][0],&ibit_chan[0][0],&nnfiles,&nnmetafiles,&nmax_files,&nmax_metafiles,&nmax_channels);
    
	      choose_channels_(metafilename,infilename,&channels[0][0][0],&channels1[0][0][0],&channels2[0][0][0],&nchannels[0][0],&nnfiles,&nnmetafiles,&nmax_files,&nmax_metafiles,&nmax_channels,&istart,&iend);
    
              chansel=0;
              switch_chan=0;
              totnum=0;
              for (k=0;k<nchannels[0][0];k++)
                {
                  if(istart==channels1[k][0][0] && switch_chan==0)
                    {
                      refchan=k;
                      switch_chan=1;
                    }
                  if(channels[k][0][0]==1)
                    {
                      totnum++;
                    }
                }
              
              i=0;
            }
        

	  /*This software can handle Events_125us and GoodXenon data*/
	  /* So it checks the time resolution and infers the total */
	  /* number of bits to filter on energies, PCUs and time-markers*/
          printf ("The processing infilename is: %s\n",buf);
          if(tsamp[0]==0.0001220703125)
            {
              nbits=16;
              printf("THIS IS AN EVENT_125us FITS FILE  nbits=%d\n\n",nbits);
	      printf("Sampling time: %e\n",tsamp[0]);
	      sampling=tsamp[0];
	      sprintf(root,"event_FILE");
	    }
          else
            {
              nbits=24;
              printf("THIS IS A GoodXenon FITS FILE, nbits=%d\n\n",nbits);
	      printf("Sampling time: %e\n",tsamp[0]);
	      sampling=tsamp[0];
	      sprintf(root,"goodX_FILE");			
            }
          
          if (fits_close_file(infptr, &status))
            printerror_C( status );
	  //////////////////////////////////////////////////////////

	  /*Here you can define your customer resolution*/
	  /*It will be used later to rebin*/
	  //tsamp[0]=0.0001220703125;
	  

	  /*Uncomment the three lines below and recomplile if you use an */
	  /*old version of Xenon2fits with the known half time-pixel bug. */
	  //	  if(nbits==24)
	  //	    delta_pixel=-0.0001220703125/2.; 
	  //	  else
	  delta_pixel=0.;
	  
	  /*Define Nyquist Frequency*/
          frequency = 1/(2*tsamp[0]);          
          chkct = frequency;      
          

	  /*Filter with GTIs if you want to dump ToAs into an ASCII file*/
	  /*This chunk of code is in a weird place, better to move it*/
	  /*somewhere else to improve readability*/
	  if(selec==8)
            {
              if(numfiles==0)
                {
                  printf("Do you want to use a GTI file (gti.dat) ?\n");
                  scanf("%s",&gtiyesno);
		  gtiyesno=tolower(gtiyesno);
                  printf("You said: %c\n",gtiyesno);
                  printf("Do you want to exclude (0) or include (1) those gtis ?\n");
                  scanf("%d",&incl_excl);
                  if(gtiyesno=='y')
                    {
                      gtifile=fopen("gti.dat","r");

                      for(;;)
                        {
                          if(gticount>1 && gti_i[gticount]==gti_i[gticount-1])
                            break;
                          fscanf(gtifile,"%lf %lf",&gti_i[gticount],&gti_f[gticount]);
                          printf("GTI INTERVALS: %lf %lf\n",gti_i[gticount],gti_f[gticount]);
                          fprintf(gti_mjd,"%16.11lf %16.11lf\n",gti_i[gticount]/86400.+49353.00073567628,gti_f[gticount]/86400.+49353.00073567628);
			  fflush(gti_mjd);
                          gticount++;
                
                        }
                    }
                }
              tformat=dump(column,tformat);
              firstcall++;
              free(bary_i);
              free(bary_f);
              free(infilename);
              continue;
            }
          else
            {
              column=1;
            }

          /* B: Read barytimes */
          /* Calls the subroutine selectrows that reads filters and */
	  /* stores barycentered data from the BARYTIME column*/
          selectrows(tsamp[0],column,incl_excl,numfiles);        //read all the rows in the fits file. 
          nout=noutrows;               // counts the total number of rows.
          printf("NOUTROWS = %d\n",noutrows); 
          if(noutrows<=10)
            {
              printf("This file contains less than 10 rows. Skipped\n");
	      if(numfiles==0)
		numfiles--;
              continue;
            }
    
          free(infilename);
          printf("Infilename pointer freed\n");
    


          /* imjd_m+fmjd_m = middle of the obsevation  (refers to the WHOLE FITS file data)*/
    
          imjd_m = (long int) ((*bary_init_i+*bary_fin_i)/(86400.*2.)+(*bary_init_f+*bary_fin_f)/(86400.*2.));// in units of days
          fmjd_m = (double) (*bary_fin_f/(2.*86400.)) ;
          fmjd_m+= (double) (*bary_init_f/(2.*86400.));
          fmjd_m+= (double) (*bary_fin_i/(2.*86400.));
          fmjd_m-= (double) imjd_m; 
          fmjd_m+= (double) (*bary_init_i/(86400.*2.));
	  printf("imjd_m = %ld fmjd_m = %16.15f\n",imjd_m,fmjd_m);

	  /*Defines the total number of profiles (ndmp_o) and the EFFECTIVE time*/
	  /*interval for each chunk (tdmp_o)*/
          if(!offon && imjd_m+fmjd_m>timemjd)
            {
              tdmp_o=tdmp_2;
              printf("We are using the second time interval of %lf\n",tdmp_o);
            }
          if(chunk!=3)
            {
              ndmp_o = (int) (var/tdmp_o);  //number of chunks 
              tdmp_o = (double) ((bary_f[nout-1]-bary_f[0])/ndmp_o+(bary_i[nout-1]-bary_i[0])/ndmp_o); // seconds for each chunk (e.g. 128 sec.)
            }
          else
            {
              ndmp_o = (int) (var/tdmp_o);
            }

          /* If ndmp_o is zero then adjust the tdmp_o time so that ndmp_o = 1 */
	  /* In other words: use the whole data in the file for a single profile*/
          if(ndmp_o==0)
            {
              printf("barytime[nout-1]= %e  barytime[0] = %e \n",bary_f[nout-1],bary_f[0]);
              printf("ndmp_o=0 changing tdmp_o so we have only one profile\n");
              tdmp_o=(double) ((bary_f[nout-1]-bary_f[0]) + (bary_i[nout-1]-bary_i[0]));
              ndmp_o = (int) (var/tdmp_o);//(bary_i[nout-1]-bary_i[0])/tdmp_o);  //number of chunks  
              tdmp_o = (double) ((bary_f[nout-1]-bary_f[0])/ndmp_o+(bary_i[nout-1]-bary_i[0])/ndmp_o); 
            }

          if(offon==1 && std==1)  //i.e. build a standard (offon=1) and local profiles (std==1)
            {
	      tdmp_o=(double) ((bary_f[nout-1]-bary_f[0]));
	      tdmp_o+=(bary_i[nout-1]-bary_i[0]);
	      ndmp_o = (int) (var/tdmp_o); 
            }
    
          printf("Total time in this file: %lf\n",-bary_i[0]-bary_f[0]+bary_i[nout-1]+bary_f[nout-1]);
          printf("dfname=%s\n",dfname);
    
          
	  /* Calculate now the coefficients of the polynomial expansion */
	  /* It uses the reference file in ~tempo/tzpar/namefile.par    */
	  /* The polynomial expansion refers to the middle              */
	  /* of the observations (imjd_m+fmjd_m)                        */
	  mpolyco(dfname,psrname,imjd_m,fmjd_m,site,&nspan_o,&ncoeff_o);
	  printf("ncoeff_o=%d\n",ncoeff_o);
	  /* determine the period and phase at the middle of the observation */
	  /* from the polynomial expansion created with mpolyco*/
	  ppolyco(dfname,imjd_m,fmjd_m,&pobs,&refph); 
	
	  
	  if(numfiles==0)
	    {
	      refph1=(double *) malloc(sizeof(double));
	      *refph1=refph;
	      refph2=refph;
	    }
          else
	    {
	      refph2=refph;
	    }
          

	  /*Re-adjust number of profiles*/
          if(chunk!=3)
            ndmp_o = (int) (var/tdmp_o + 0.5); 
          
	  /*Write output file ndmp.dat with filenames and number of profiles created*/
          printf("nbin %d tdmp_o %f ndmp_o %d\n\n\n\n",nbin,tdmp_o,ndmp_o);
          fprintf(ndmp,"Filename: %s  ndmp_o %d\n",buf,ndmp_o);
	  fflush(ndmp);
          fprintf(timeint,"Initial time:%lf    Final time=%lf  Number of chunks:%d  Time interval for each chunk:%lf  Infilename:%s\n\n",(bary_i[0]+bary_f[0])/86400.0,(bary_i[nout-1]+bary_f[nout-1])/86400.0,ndmp_o,tdmp_o,buf);
    

	  /*Allocates an array to contain the whole binned unfolded lightcurve */
          if(chunk!=3)
            {
              opr = (float *) calloc(nbin*ndmp_o+1,sizeof(float)); 
            }
          else
            {
              opr = (float *) calloc(nbin*(ndmp_o+1)+1,sizeof(float));
            }
	  printf("Allocating *opr\n");	  
          
	  /*Initialize the arrays to store the lightcurve for each PCU*/
          indx = 0;
          for(k=0; k<6; k++)
	    counts_pcu[k]=0.;

	  k=0;
          time_r_i = bary_i[0];
          time_r_f = bary_f[0]; // Round off errors ?!? We can lose the precision required. Change it ?!?
	  counterprofile=0;
	  bintime=time_col[0];

	  /////////////////////////////////////
          /* Main Rebinning & Folding Engine */
	  /////////////////////////////////////

	  /* REBNNING */
	  j=0;
	  np=(double)(tdmp_o/tsamp[0]);
	  check_bins=(int)(ndmp_o*np)+1;
	  rdata= (float *) calloc((int)(ndmp_o*np)+1,sizeof(float));
	  phpos = (int *) calloc((int)(ndmp_o*np)+1,sizeof(int));
	  tbin_i=0.;
	  tbin_f=0.;
	  tbin=0.;
	  ibin=0;
	  ibin_old=0;
	  counter_ph_offset=0;
	  counter_ph=0;
	  check_cts=0.;
	  counter=0;
	  time_initial_i=bary_i[0];
	  time_initial_f=bary_f[0];
	  fix_time_i=time_initial_i;
	  fix_time_f=time_initial_f;
	  time_final_i=bary_i[0]+(double)((int)(tdmp_o));
	  time_final_f=bary_f[0]+tdmp_o-(double)((int)(tdmp_o));
	  while(time_final_f > 1.0) {time_final_f--; time_final_i++;}

	  for(j=0;j<(int)(np*ndmp_o)+1;j++)
	    {
	      rdata[j]=0.;
	      phpos[j]=0;
	    }
	  
	  for(j=0; j<ndmp_o; j++)
	    {
	      if(presto_switch)
		{
		  prestofile = (char *) calloc(22, sizeof(char));
		  strcat(prestofile,root);
		  sprintf(strfile2,"%d",numfiles+1);
		  strcat(prestofile,strfile2);
		  sprintf(strfile,"%d",j+1);
		  strcat(prestofile,"_TS");
		  strcat(prestofile,strfile);
		  strcat(prestofile,".dat");
		  presto_out=fopen(prestofile,"wb");
		  printf("Rebinning time series %d\n",j+1);
		}

	      switch_counter=0;
	      counter_ph_offset+=counter_ph;
	      counter_ph=0;
	      while(switch_counter==0)
		{
		  if((bary_i[counter]-time_final_i>0)	||(bary_i[counter]-time_final_i>=0 && bary_f[counter]-time_final_f>-tsamp[0]/2.))
		    {
		      switch_counter=1;
		    }
		  else
		    switch_counter=0;
		    
		  if(!switch_counter)
		    {
		      tbin= bary_i[counter]-fix_time_i;
		      tbin=tbin+bary_f[counter];
		      tbin= tbin-fix_time_f;
		      ibin_old=ibin;
		      ibin=(int)(tbin/tsamp[0]);
		      rdata[ibin] = rdata[ibin] + 1.0;
		      if(ibin!=ibin_old)
			{
			  phpos[counter_ph+counter_ph_offset] = ibin;
			  counter_ph++;
			}
		      counter++;
		    }
		}
	      if(presto_switch)
		{
		  if((int)(np) % 2 ==0)
		    {
		      printf("Time series is even, writing on file;\n");  
		      even_bin=0;
		    }
		  else
		    {
		      printf("Time series is odd, removing last bin and writing on file;\n");
		      even_bin=1;
		    }
		  for(n=(int)(np*j);n<(int)(np*(j+1))-even_bin; n++)
		    fwrite(&rdata[n], sizeof(float),1, presto_out);
		}
	      printf("Final ibin index = %d\n",ibin);
	      printf("Maximum index in *phpos array --> %d\n",counter_ph+counter_ph_offset);
	      printf("Counter returned %d array elements in barytimes scan\n", counter);

	      time_initial_i= time_final_i;
	      time_initial_f= time_final_f;
	      time_final_i=bary_i[0]+(int)((tdmp_o)*(j+2));
	      time_final_f=bary_f[0]+tdmp_o*(j+2)-(int)(tdmp_o*(j+2));
	      while(time_final_f > 1.0) {time_final_f--; time_final_i++;}
	      
	      if(presto_switch)
		{
		  free(prestofile);
		  fclose(presto_out);
		}
	    }
	  
	  for(n=0; n<check_bins; n++)
	    {
	      check_cts+=rdata[n];
	    }
	  printf("TOTAL COUNTS IN RDATA = %f    TOTAL NUMBER OF BINS=%d\n",check_cts, check_bins);
	  i=0;
	  j=0;
	  counter=0;
	  ibin=0;
	  check_cts=0;
	 

	  /* A: Binning the data of the FITS file   */ 
	  /*    according to the number of pulse profiles */	  
	  for(n=0; n<check_bins; n++)
	    {
	      if(rdata[n])
		{
		  ibin=phpos[indx];
		  ict=rdata[ibin];
		  check_cts+=rdata[ibin];
		  

		  counts_pcu[pcubary[indx]]++;
		  if(indx==0)printf("Counter Treshold=%19.16lf    tdmp_o*(counterprofile+1)=%16.15lf\n",bary_i[0]+bary_f[0] + tdmp_o*(counterprofile+1), tdmp_o*(counterprofile+1));		 		 
		  if((bary_i[(int)(check_cts)] > bary_i[0]+ (int)(tdmp_o*(counterprofile+1))) ||  (bary_i[(int)(check_cts)] == bary_i[0]+ floor(tdmp_o*(counterprofile+1)) && bary_f[(int)(check_cts)]>=tdmp_o*(counterprofile+1)-floor(tdmp_o*(counterprofile+1))))
		  {
		      check_pcus(counts_pcu);
		      bkgshift= bkgcount[numfiles]*tdmp_o;
		      fprintf(lightcurves,"%14.13lf %14.13lf %lf %lf %lf %lf %lf %d %lf %lf\n",time_col[0]/86400.+timeoffset+(tdmp_o*counterprofile)/86400., bary_i[0]/86400.+bary_f[0]/86400.+(tdmp_o*counterprofile)/86400.,counts_pcu[0],counts_pcu[1],counts_pcu[2],counts_pcu[3],counts_pcu[4], prof_act_pcus, bkgshift, tdmp_o);		       
		      fflush(lightcurves);
		      counterprofile++;
		      for(kk=0; kk<6; kk++)
			counts_pcu[kk]=0.;
		    }

		  rimjd_i = (long int) (bary_i[(int)(check_cts)]/86400.); 
		  rimjd_f = (double) (bary_i[(int)(check_cts)]/86400.) ;
		  rimjd_f-= (double) rimjd_i;
		  rimjd_f+= (double) (bary_f[(int)(check_cts)]/86400.) ;
		  
		  if(indx==0)
		    {
		      printf("rimjd_i = %10.3ld rimjd_f=%18.17lf  indx =%f   bary_i[%f]=%lf    bary_f[%d]=%17.16lf\n",rimjd_i,rimjd_f, check_cts, check_cts, bary_i[(int)(check_cts)],indx,bary_f[(int)(check_cts)]);
		    }		 
		  indx++;
		}
	      else
		{
		  rimjd_i = (long int)   (time_r_i/86400.0);
		  rimjd_f = (double) (time_r_i/86400.);
		  rimjd_f-= (double) rimjd_i;
		  rimjd_f+= time_r_f/86400. ; 
		  ict = 0.0;
		}
	      
	      if (chkct > frequency - 1)  //in the first loop of the while cycle chkct=frequency, then it's zero incresing by one per time
		{
		  imjd = rimjd_i;
		  fmjd = rimjd_f;
		  
		  /* Polycode for the photon index-th, calculates the phase and period at the moment of its arrival*/
		  ppolyco(dfname,imjd,fmjd,&spobs,&phs);
		  
		  sph = fmod(phs - refph,1.0); //fractional part of the difference between the photon phase (phs) and the middle of the obs. phase (refph)
		  
		  if(sph < 0.0) sph += 1.0;                   
		  dph = (double) (nbin * tsamp[0]/spobs);
		  sph = sph*nbin;  
		  sph = sph + dph/2.0;
		  chkct = 0;
		  chkdct++;
		}
	      
	      time = (double) (chkdct*frequency*tsamp[0] + chkct*tsamp[0]);  //chkdct here is zero and incereases only when chkct > frequency - 1

	      time_ct = (int) (time/tdmp_o);                                 
	      if(sph > (1.0*(nbin-1))+0.5) sph = sph - 1.0*nbin;             	     	     
	      ib = (int) (sph+0.5);   /* The integer bin number */           
	      binct = ib+time_ct*nbin;                                       
	      opr[binct] = opr[binct] + 1.0*ict;                             
	      sph = sph+dph;
	      chkct++;
	      time_r_f += tsamp[0];
	      if(time_r_f > 1.0) {time_r_f--; time_r_i++;}          


	    }
	  printf("TOTAL COUNTS in while cycle=%f\n",check_cts);

	  check_cts=0;
	  for(counter=0; counter<nbin*ndmp_o+1; counter++)
	    check_cts+=opr[counter];
	  printf("TOTAL OPR COUNTS=%f\n",check_cts);

	  free(rdata);
	  free(phpos);
	  deltaph = (double *) malloc(sizeof(double));

	  /* Global or Local standard profile (local means a profile for each FITS file)*/          
          switch(std)
            {
            case 1:
              printf("Building local profiles\n");
              counter++;
              profiles(imjd_m,fmjd_m,nout,counter);
              timelags();
              continue; //stops the numciles and brings you back to the beginning
            }
          count2 =  (double *) calloc(nbin,sizeof(double));
          stdamps = (double *) calloc(nbin,sizeof(double));
          prfamps = (double *) calloc(nbin+1,sizeof(double));
          maxbin2 = (double *) malloc(sizeof(double));
          terr    = (double *) malloc(sizeof(double));
          printf ("Allocating count2, stdamps, prfamps,maxbin2, terr\n");
          
          /* Builds the standard profile */          
          smax=0.0;
	  if( s != NULL && lcount != 0)
	    {
              if(offon==1 && lcount == 1)
                {
                  for(j=0;j<nbin;j++)
                    {
                      for (i=0;i<ndmp_o;i++)
                        {
                          count2[j]  += opr[j+i*nbin];
                          stdamps[j] += opr[i*nbin+j];
                          if (smax < stdamps[j]) {smax = stdamps[j]; maxbin = 1.0*j;} 
                        }
                      if(bkgyesno=='y')
                        count2[j]-=bkgcount[numfiles]*tdmp_o/nbin;
                    }
                }  
              else if(offon == 0 && init==0)
                {
                  init=1;
                  printf("Using the previous standard profile\n");
                  for(j=0;j<nbin;j++) 
                    {
                      count1[j]=standard[j];          
    
                    }
                  for(j=0;j<nbin;j++)
                    {
                      if (smax < count1[j]) {smax = count1[j]; maxbin = 1.0*j;}
                    }
		  standa=fopen("standard.st","w");
                  for(j=0;j<nbin;j++)
                    fprintf(standa,"%d %lf\n",j, count1[j]);
                  
                  fclose(standa);
                }
              
              if(numfiles==0 && offon == 1)
                { 
		  for(j=0;j<nbin;j++)
                    {
                      count1[j]=stdamps[j];
                    } 
                  lcount++; 
                  continue; 
                } 
              if(offon==1  && lcount >=2)
                    {  
                      for(i=0;i<ndmp_o;i++)
                        {
                          for(j=0;j<nbin;j++)
                            {
                              count2[j]+=opr[j+i*nbin];
                            }
                        }
                      toa_stds(toaptr_stds,dfname,count1,nbin,count2,nbin,imjd_m,fmjd_m,pobs,freq,terr,maxbin2);
                      smax=0.0;
                      for(j=0;j<nbin;j++)
                        {
                          count1[j]+=count2[j];
                          if (smax < count1[j]) {smax = count1[j]; maxbin = 1.0*j;}
                        }
		      free(maxbin2);
                      free(terr);
		      free(deltaph);
                      free(bary_f);
                      free(bary_i);
                      free(opr);  
                      free(count2);
                      free(stdamps);
                      free(prfamps);
                      printf("Pointers *maxbin2, *terr, *bary_i, *bary_f, *opr, *count2, *stdamps and *prfamps freed\n");
                  continue;
                    }  
            }  /*close the main IF statement*/
          
          
          free(count2);
          free(stdamps);
          printf("Pointers *stdamps and *count2 freed\n");
	  
          /* Now we start the analysis for each chunck of data */
	  /* i.e., we calculate the ToAs for each profile*/
	  /* via cross-correlation with a standard template*/        
          numbr=0;
          intervals=fopen("intervals.dat","a");
          for (i=0;i<ndmp_o;i++)
            {
              prfmean=0.;
              smaxprf=0.;
              maxbinprf=0.;
	      if(ndmp_o>1)
		{
		  /* Need to calculate the mjd for the pulse arrival phase for each subint */
		  simjd_m = (long int) ((bary_i[0]/86400.+bary_f[0]/86400.)+(double)(((tdmp_o * i)+tdmp_o/2.)/86400.0));
		  sfmjd_m = bary_i[0]/86400. + bary_f[0]/86400.+ (double)(((tdmp_o * i)+tdmp_o/2.)/86400.0)-simjd_m;     
		}
	      else
		{
		  simjd_m = imjd_m;
		  sfmjd_m = fmjd_m;
		}
	      
	      if (sfmjd_m > 1.0) {simjd_m += 1; sfmjd_m -= 1.0;}



	      /* Print middle of each chunk and chunk length */
	      fprintf(intervals,"%ld %16.15lf %ld %16.15lf %lf\n",simjd_m,sfmjd_m,imjd_m, fmjd_m, tdmp_o);		  
	      ppolyco(dfname,simjd_m,sfmjd_m,&spobs_m,&phs); 
              dphase = fmod(phs - refph,1.0);
              sfmjd_m = sfmjd_m - dphase*spobs_m/86400.0;
              if (sfmjd_m > 1.0) {simjd_m += 1; sfmjd_m -= 1.0;}
              for (j=0;j<nbin;j++)
                {
                  prfamps[j] = opr[i*nbin+j];
		  if(bkgyesno=='y')
                    {
                      prfamps[j] =prfamps[j] -bkgcount[numfiles]*tdmp_o/nbin;
		    }
                  
                  if (smaxprf < prfamps[j]) {smaxprf = prfamps[j]; maxbinprf = j;}
		  prfmean+=prfamps[j];
                }
              
	      fflush(intervals);
              if(offon==1)
                {
                  continue;
                }
              else
                {
                  numb++;
                  numbr++;
                  num++;
		  for(profi=0; profi<nbin; profi++)
		    if(prfamps[i]<0.)
		      {
			fprintf(amplitudes, "Amplitude < 0. Skipping file nr. %d\n", numfiles);
			fflush(amplitudes);		       	    
			continue;
		      }
		  *deltaph=fmod(refph2-*refph1,1.0)*nbin;
                  count_info=i;
		  toa(toaptr,toatempo2,dfname,count1,nbin,prfamps,simjd_m,sfmjd_m,spobs_m,freq,terr,site,EMAX,maxbin,deltaph,refph1);
                }
              if(offon==0 && i==ndmp_o-1)
                { 
                  printf("Pointer *maxbin2 *deltaph and *terr freed\n");
                  free(maxbin2); 
                  free(terr); 
		  free(deltaph);
                }            
              totcs=0.;
              for(k=0;k<nbin;k++)
                totcs+=prfamps[k];
              fprintf(totcounts,"%lf %ld %17.16lf %d %lf %lf %lf %lf %d\n",totcs,simjd_m,sfmjd_m,finalpcu,totcs+bkgcount[numfiles]*tdmp_o/nbin,bkgcount[numfiles]*tdmp_o/nbin,bkgcount[numfiles],tdmp_o, num_pcu);
              fflush(totcounts);
            } //End of ndmp_o loop
          
          fclose(intervals);
          free(bary_f);
          free(bary_i);
	  free(time_col);
	  free(opr);
          free(prfamps);
	  free(pcubary);
	  printf("Pointer *bary_i, *bary_f, *pcubary, *prfamps, *opr and *time_col freed\n");
         
          printf("File %s correctly processed\n",buf);
        } // End "numfile" loop 
    
    }  // End of "optind" loop
    

 
  free(refph1);
  free(psrname);
  printf("Pointer *psrname & refph1 freed\n");
  free(count1);
  free(totp);
  printf("Pointer *count1  and *totp freed\n");
  fclose(toaptr);  
  printf("toaptr file closed\n");
    
  fclose(standa);
  printf("standa file closed\n");
    
  fclose(fp);
  printf("fp file closed\n");
    
  fclose(ndmp);
  printf("ndmp file closed\n");
    
  fclose(prfptr);
  printf("prfptr file closed\n");
    
  fclose(amplitudes);
  printf("amplitudes file closed");

  fclose(lightcurves);
  printf("lightcurves file closed");

  
  fclose(totcounts);
  fclose(gti_mjd);
    
  free(bary_init_i);
  free(bary_init_f);
  free(bary_fin_i);
  free(bary_fin_f);
    
    
  return(0);
}       
    
    
    
/*-------------------------------------------------------------------*/
/************************************************************/
/* Reads rows in FITS files and selects the good Event data */
/************************************************************/
    
double *selectrows(double tsamp, int column, int incl_excl, int numfiles)
{
  fitsfile *infptr;  /* pointer to input and output FITS files*/
  int status, hdutype, nkeys, keypos, nfound, colnum,colnum_time, anynulls;
  long naxes[2], frow, felem, irow, firstrow,finrow;
  float nullval;
  double  mjdrefi[1], timezero[1], mjdreff[1]; //Change here later...
  int  exit_pcu,i,ndet,ii,l,m,nr,n,ioff,j,k,ichan, ibit, npcus,pcu[6];
  int active_pcus;
  double flux_pcu[6];
  char array[nbits];
  double intbary,fbary;
  //"faxbFS4f_10915ff0.fits";  /* name for existing FITS file   */ 
  char *s;
  FILE *data, *dataunbary;
    
  data=fopen("data.dat","a");
  dataunbary=fopen("dataunbary.dat","a");

  pcus=fopen("pcus.dat","a");
  status = 0;
  anynulls=0;
  npcus=0;
  i=0;
    
  /* open the existing FITS files */
  if ( fits_open_file(&infptr, infilename,  READONLY,  &status))
    printerror_C( status );
  /* move to the 2nd HDU in the input file (a binary table in this case) */
  if ( fits_movabs_hdu(infptr, 2, &hdutype, &status))
    printerror_C( status );
    
  if (hdutype != BINARY_TBL)  {
    printf("Error: expected to find a binary table in this HDU\n");
    return(0);
  }
    
  /* read the NAXIS1 and NAXIS2 keyword to get table size */
  if (fits_read_keys_lng(infptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
    printerror_C( status );
    
    
  /* read the  TIMEZERO keyword */
  if(fits_read_key(infptr, TDOUBLE, "TIMEZERO", timezero, 0, &status))
    printerror_C( status );
    
  /* read the  MJDREFF keyword */
  if(fits_read_key(infptr, TDOUBLE, "MJDREFF", mjdreff, 0, &status))
        printerror_C( status );
    
  /* read the MJDREFI keyword */
  if(fits_read_key(infptr, TDOUBLE, "MJDREFI", mjdrefi, 0, &status))
    printerror_C( status );
    
  /* find which column contains the BARYTIME values */
  if(column==1 && switch_bary)
    {
      if ( fits_get_colnum(infptr, CASEINSEN, "BARYTIME", &colnum, &status) )
	printerror_C( status );
      if ( fits_get_colnum(infptr, CASEINSEN, "Time", &colnum_time, &status) )
	printerror_C( status );
    }
  else if(column==1 && !switch_bary)
    {
      if ( fits_get_colnum(infptr, CASEINSEN, "Time", &colnum, &status) )
        printerror_C( status );
      colnum_time=colnum;
    }
  else
    {
      if ( fits_get_colnum(infptr, CASEINSEN, "Time", &colnum, &status) )
        printerror_C( status );
    }
    
  printf("TOAs are in column number %d\n",colnum);
  printf("MJDREFI = %5.1f  MJDREFF = %20.19f\n",mjdrefi[0],mjdreff[0]);
  timeoffset= (double) ((mjdrefi[0]+mjdreff[0])+timezero[0]/86400.);
  printf("Time offset=%16.11lf\n",timeoffset);
    
  if(numfiles==0 && toasyn=='y')
    {
      for(i=0;i<lines;i++)
        {
          toat[i]-=timeoffset;
          toat[i]=toat[i]*86400.;
        }
    }
    
    
  /* read the BARYTIME column values */
  frow = 1;
  felem = 1;
  nullval = -99.;
  nr=1;
  printf("frow=%ld, numrows=%ld\n",frow, naxes[1]);
  printf("Allocating barytime\n");
  barytime   = (double *) calloc(naxes[1], sizeof(double));
  time_col       = (double *) calloc(naxes[1], sizeof(double));
  pcubary     = (int *) calloc(naxes[1], sizeof(int));
    
  // Reads Barytimes. (fits_read_col = ffgcv CFITSIO subroutine)
  if ( fits_read_col(infptr, TDOUBLE, colnum, frow, felem, naxes[1], 
                     &nullval, barytime, &anynulls, &status) )
    printerror_C(status);
  if ( fits_read_col(infptr, TDOUBLE, colnum_time, frow, felem, naxes[1], 
                     &nullval, time_col, &anynulls, &status) )
    printerror_C(status);
  
  
  

  printf("Ok, debugging\n");
  printf("UNCLEANED rows=%ld\n",naxes[1]);
  printf("Removing Delta Pixel\n");
  printf("Barytime[0] = %18.6lf\n",barytime[0]);
  printf("Barytime[naxes[1]-1] = %18.6lf\n",barytime[naxes[1]-1]);
  printf("Delta Pixel = %e\n",delta_pixel); 
  ii=0;
  firstrow=1;
  
  for(i=0;i<5;i++)
    pcu[i]=10;
  ndet=0;
  for(k=0;k<6;k++)
    flux_pcu[k]=0.;
  k=0;
  
  
  for(frow=1;frow<naxes[1];frow++)
    
    {
      // read each Event row in array[i] once per time
      fits_read_col_bit(infptr, 2, frow, 1,nbits, array, &status);
      //NOTE THAT IF YOU ARE USING GTIFILTER() THE FIRST ROW CAN 
      // BE DIFFERENT FROM THE FIRST IN THE FITS FILE. 
      // AND THIS IS NOT AN ERROR: IT IS CORRECT!
      //  array[0] = 0 means FLAG TIME MARKER
      //  array[0] = 1 means PHOTON TIME OF ARRIVAL 
      // (array[0] is the most significant bit)
      // If there's a channels range selection, array[0] is modified
      // and we use array[0]-->0 if it's not in the selected channel range
      
      if(array[0]==1)   
        {
          //Channel selection:
          ichan = 0 ;
          ibit  = 0 ;
    
          for(ibit=0;ibit<ibit_chan[0][0];ibit++) //e.g. for Event mode --> ibit_chan = 6 --> 2**6=64 chan.
            {
              if(nbits==16)
                ioffset[0][0]=5;
              else if(nbits==24)
                ioffset[0][0]=17;
              else
                {printf("Strange data format, nbits does not match. Check offset.\n"); exit(1);}
	      
              if(array[ibit+ioffset[0][0]-1])  //e.g. ioffset = 5, but array starts from 0 and not 1, so ioffset[0][0]]-1
                ichan = ichan + pow(2,(ibit_chan[0][0]-1-ibit));
           }
          //ichan++;  //e.g. for Event mode converts 0-63 into 1-64 channels (1-2**n)

    
          //PCU counting         
	  // I check which PCU is hitted by each photon. (PCU Nr. = [0:4]. "a" is the PCU nr. here)
	  // Then I see at the end how many different detectors were hitted.
	  // THIS IS STILL WRONG, SINCE THE NUMBER OF active PCUs
	  // CAN CHANGE WITHIN ONE FITS FILE, SO THIS CAN BE FAR
	  // FROM BEING CORRECT
	  exit_pcu=0;
          npcus=0;
          if(nbits==16)
	    {
	      pcuoffset=1;
	    }
          else if(nbits==24)
	    {
	      pcuoffset=7;
	    }
          
          for(i=0; i<3; i++)
            {
              if(array[i+pcuoffset])  // i = bit, +1 or +7= offest
                {
                  npcus= npcus + pow(2,(2-i));
                }
            }
          
	  for(i=0;i<num_pcu;i++)
	    { 
	      if(good_pcu[i]!=npcus)
		{
		  exit_pcu=1;
		  continue;
		}
	      else if(good_pcu[i]==npcus)
		{
		  exit_pcu=0;
		  break;
		}
	    }
	  if(exit_pcu==1)
	    continue;
	  

          if(ndet==0)
            {pcu[ndet]=npcus;ndet++;}

          for(i=0;i<ndet+1;i++)
            {
              if(npcus==pcu[i])
                break;
              else if(npcus!=pcu[i] && pcu[i]==10)
                {pcu[i]=npcus; ndet++; break;}
            }
          
	

          *bary_init_i=barytime[0];
          finrow=naxes[1];
          fbary=barytime[finrow-1];
          *bary_fin_i=barytime[finrow-1];
          
          if(toasyn=='y')
            {
            for(k=0;k<lines;k++)
              {
                if(barytime[frow-1]>(toat[k]-tdmp_o/2.) && barytime[frow-1] < (toat[k] + tdmp_o/2.))
                  {
                    if(channels[0][0][ichan] && gtiyesno=='y')
                      {  
                        for(i=0;i<gticount-1;i++)
                          {
                            if(incl_excl==1)
                              {
                                if( (time_col[frow-1] > gti_i[i] && time_col[frow-1]<gti_f[i]) || (i==gticount-2 && time_col[frow-1] >  gti_f[gticount-2]))
                                  {
                                    barytime[ii]=barytime[frow-1];
				    time_col[ii]=time_col[frow-1];
				    pcubary[ii]=npcus;
				    ii++;
				    flux_pcu[npcus]=flux_pcu[npcus]+1;
                                  }
                              }
                            else if(incl_excl==0)
                              {
                                if( (time_col[frow-1] < gti_i[0] && i==0) || (time_col[frow-1] > gti_f[i] && time_col[frow-1] < gti_i[i+1])|| (i==gticount-2 && time_col[frow-1] >  gti_f[gticount-2]))
                                  {
                                    barytime[ii]=barytime[frow-1];
				    time_col[ii]=time_col[frow-1];
				    pcubary[ii]=npcus;
				    ii++;
				    flux_pcu[npcus]=flux_pcu[npcus]+1;
                                  } 
                              }
                          }
                      }
                    else if (channels[0][0][ichan] && gtiyesno=='n')
                      {
                        barytime[ii]=barytime[frow-1];
			time_col[ii]=time_col[frow-1];
			pcubary[ii]=npcus;
                        ii++;
			flux_pcu[npcus]=flux_pcu[npcus]+1;
                      }
                  }
              }
            }
          else
            {
	      if(channels[0][0][ichan] && gtiyesno=='y')
		{  
		  for(i=0;i<gticount-1;i++)
		    {
		      if(incl_excl==1)
			{
			  if((time_col[frow-1] > gti_i[i] && time_col[frow-1]<gti_f[i]))
			    {
			      barytime[ii]=barytime[frow-1];
			      time_col[ii]=time_col[frow-1];
			      pcubary[ii]=npcus;
			      ii++;
			      flux_pcu[npcus]=flux_pcu[npcus]+1;
			    }
			}
		      else if(incl_excl==0)
			{
			  if((time_col[frow-1] < gti_i[0] && i==0) || (time_col[frow-1] > gti_f[i] && time_col[frow-1] < gti_i[i+1]) || (i==gticount-2 && time_col[frow-1] >  gti_f[gticount-2])) 
			    {
			      barytime[ii]=barytime[frow-1];
			      time_col[ii]=time_col[frow-1];
			      pcubary[ii]=npcus;
			      ii++;
			      flux_pcu[npcus]=flux_pcu[npcus]+1;
			    } 
			}
		    }
		}
	      else if (channels[0][0][ichan] && gtiyesno=='n')
		{
		  barytime[ii]=barytime[frow-1];
		  time_col[ii]=time_col[frow-1];
		  pcubary[ii]=npcus;
		  ii++;
		  flux_pcu[npcus]=flux_pcu[npcus]+1;
		}
            }            
        }
      firstrow++;
    }
  
    
  noutrows=ii;
  finalpcu=ndet;
  printf("NUMBER OF ACTIVE PCUs=%d\n",finalpcu);
  printf("cleaned rows=%d\n",noutrows);
  printf("Barytime[0] = %18.6lf\n",barytime[0]);
  printf("Barytime[%d] = %18.6lf\n",noutrows-1,barytime[noutrows-1]);
  printf("PCUs=%d\n",k+1);
  printf("COUNTS PCU0=%lf  PCU1=%lf  PCU2=%lf  PCU3=%lf PCU4=%lf\n",flux_pcu[0],flux_pcu[1],flux_pcu[2],flux_pcu[3],flux_pcu[4]);
 
  if(noutrows<10)
    {
      free(barytime);
      free(pcubary);
      free(time_col);
      return(0);
    }
    
  printf("Re-allocating barytime\n");
  if(!realloc(barytime, noutrows*sizeof(double)))
    printf("Error in Re-alloc barytime\n");;
   printf("Re-allocating time_col\n");
  if(!realloc(time_col, noutrows*sizeof(double)))
    printf("Error in Re-alloc time_col\n");;
  if(!realloc(pcubary,  noutrows*sizeof(int)))
    printf("Error in Re-alloc pcubary\n");
  
  if(barytime==NULL)
    {
      printf(error);
      exit(1);
    }
    
   printf("Allocating bary_i and bary_f\n");   
   bary_f = (double *) calloc(noutrows, sizeof(double));
   if(bary_f==NULL)
     {
       printf(error);
       exit(1);
     }
   bary_i = (double *) calloc(noutrows, sizeof(double));
   if(bary_i==NULL)
     {
       printf(error);
       exit(1);
     }
    
  for(i=0;i<noutrows;i++)
    {
      bary_i[i]= (long int) (barytime[i]); 
      bary_f[i]= (double) barytime[i] - bary_i[i];// +0.016725271940231; 
      /* delta_pixel accounts for the GoodXenon vs. Event pixel difference of 0.5 bins*/
    }
      
  intbary=*bary_init_i;
  fbary  =*bary_init_i;
  *bary_init_i = (long int) intbary;
  *bary_init_f = (double)   fbary-*bary_init_i;
  intbary=*bary_fin_i;
  fbary  =*bary_fin_i;
  *bary_fin_i = (long int) intbary;
  *bary_fin_f = (double)   (fbary - (long int)intbary);
    
  i=0;
    
    
  while(i<noutrows)
    {
      if(i==0) 
        {
          printf("Barytime[%d]=  %20.15lf\n",i,barytime[i]);
          printf("Mjdref seconds  =  %20.15lf\n",mjdreff[0]*86400.);
          printf("Timezero  =  %20.15lf\n",timezero[0]);  
	  printf("Delta Pixel       =  %20.15lf\n",delta_pixel);
	  timezero[0]-=delta_pixel;
	  printf("Timezero (Delta Pixel removed) =  %20.15lf\n",timezero[0]);  
	}
      bary_i[i]+= (long int) mjdrefi[0]*86400.0;
      bary_i[i]+= (long int)(timezero[0]+mjdreff[0]*86400.0);
      if(i==0)
        {
          *bary_init_i+=(long int) mjdrefi[0]*86400.0;
          *bary_init_i+=(long int)(timezero[0]+mjdreff[0]*86400.0);
          *bary_fin_i+=(long int) mjdrefi[0]*86400.0;
          *bary_fin_i+=(long int)(timezero[0]+mjdreff[0]*86400.0);
        }
      bary_f[i]+= fmod(mjdreff[0]*86400.0,1.0)+ fmod(timezero[0],1.0);
      if(i==0)
        {
          *bary_init_f+=fmod(mjdreff[0]*86400.0,1.0)+ fmod(timezero[0],1.0);
          *bary_fin_f+=fmod(mjdreff[0]*86400.0,1.0)+ fmod(timezero[0],1.0);
        }
    
      if(bary_f[i] > 1.0) {bary_f[i]--; bary_i[i]++;}
      if(*bary_init_f > 1.0) {*bary_init_f=*bary_init_f-1; *bary_init_i=*bary_init_i+1;}
      if(*bary_fin_f > 1.0) {*bary_fin_f=*bary_fin_f-1; *bary_fin_i=*bary_fin_i+1;}
    
      
      var=barytime[noutrows-1]-barytime[0];
    
      if(i==noutrows-1) {printf("Bary_f[0]= %20.15lf    Bary_f[%d]= %20.15lf\n",bary_f[0],noutrows,bary_f[noutrows-1]);}
      if(i==noutrows-1) {printf("Bary_i[0]= %20.1f     Bary_i[%d]= %20.1lf\n",bary_i[0],noutrows,bary_i[noutrows-1]);}
      i++;
    }
    
  if(numfiles==0) 
    {
      fprintf(data,"FILE NAME                                                                                       START(MJD)                FINISH(MJD)           Total time (sec)    Num Bins\n"); 
      fprintf(dataunbary,"FILE NAME                                                                                 START(MJD)                FINISH(MJD)           Total time (sec)    Num Bins\n"); 
    }
  fprintf(data,"%s  %15.14lf %15.14lf %15.14lf %lf\n",infilename,bary_i[0]/86400.+bary_f[0]/86400.,bary_i[noutrows-1]/86400.+bary_f[noutrows-1]/86400.,var,var/sampling);
  fprintf(dataunbary,"%s  %15.14lf %15.14lf %15.14lf %lf\n",infilename,time_col[0]/86400.+timeoffset,time_col[noutrows-1]/86400.+timeoffset,time_col[noutrows-1]-time_col[0],(time_col[noutrows-1]-time_col[0])/sampling);
  fflush(data);
  fflush(dataunbary);
  free(barytime);
  printf("Pointer *barytime freed\n");
    
  if (fits_close_file(infptr, &status))
    printerror_C( status );
  
  active_pcus=0;
  for(i=0; i<6;i++)
    if(flux_pcu[i])
      active_pcus++;
  
  fprintf(pcus,"%lf %lf %lf %lf %lf %lf %d\n",bary_i[0]/86400.+bary_f[0]/86400.,flux_pcu[0],flux_pcu[1],flux_pcu[2],flux_pcu[3],flux_pcu[4], active_pcus);
  fclose(data);
  fclose(dataunbary);
  fclose(pcus);
  return(0);
    
}
    
/*--------------------------------------------------------------------------*/
void printerror_C( int status)
{
  /*****************************************************/
  /* Print out cfitsio error messages and exit program */
  /*****************************************************/
    
    
  if (status)
    {
      fits_report_error(stderr, status); /* print error report */
      
      exit( status );    /* terminate the program, returning error status */
    }
  return;
}
    
    
void profiles(long int imjd_m, double fmjd_m, int nout, int counter)
{
  int j,i;
  char mjd_s[32];
  double *stdamps,smax,maxbin,step,tint;
  FILE *stdprof,*stdprof2;
    
  stdprof=fopen("stdprof.dat","a");
  stdprof2=fopen("stdprof2.dat","a");
  smax=0.0;
  stdamps = (double *) calloc(nbin,8);  
  count0  = (double *) calloc(nbin,8);  
  printf("Allocating stdamps\n");
  for(j=0;j<nbin;j++)
    {
      for (i=0;i<ndmp_o;i++)
        {
          if(bkgyesno=='y')
            stdamps[j] += opr[i*nbin+j]-bkgcount[numfiles]*tdmp_o/nbin;
          else
            stdamps[j] += opr[i*nbin+j];
          if (smax < stdamps[j]) {smax = stdamps[j]; maxbin = 1.0*j;} 
        }
      count0[j]=stdamps[j];
    }
    
    
  shiftbyfft_(&count0[0],&nbin,&maxbin); 
  tint=-bary_i[0]-bary_f[0]+bary_i[nout-1]+bary_f[nout-1];
  sprintf(mjd_s,"%lf",imjd_m+fmjd_m);
    
    
  for(j=0;j<nbin;j++)
    {
      if(j==0) 
        { 
          step=imjd_m+fmjd_m;
          step+=(bary_i[0]+bary_f[0]-bary_i[nout-1]-bary_f[nout-1])/(2.*86400.);
          
          fprintf(stdprof,"%2d %10.9f %10.9f %12.6f %12.6f %12.6f %s",j,stdamps[j]/smax,sqrt(stdamps[j])/smax,step,imjd_m+fmjd_m,tint,buf);
          
          fprintf(stdprof2,"%3d %10.3f %19s %5d\n ",j,stdamps[j],mjd_s,counter);
	}
      else
        {
          step=imjd_m+fmjd_m;
          step+=(bary_i[0]+bary_f[0]-bary_i[nout-1]-bary_f[nout-1])/(2.*86400.);
          step+=(-bary_i[0]-bary_f[0]+bary_i[nout-1]+bary_f[nout-1])*j/(nbin*86400.);
          fprintf(stdprof,"%2d %10.9f %10.9f %12.6f %12.6f\n",j,stdamps[j]/smax,sqrt(stdamps[j])/smax,step,imjd_m+fmjd_m);
          fprintf(stdprof2,"%2d %10.3f %19s %5d\n ",j,stdamps[j],mjd_s,counter);
	} 
    }
    
    
  fprintf(stdprof,"\n\n");
  fprintf(stdprof2,"\n\n");
  free(stdamps);
  free(count0);
  printf("Pointer *stdamps and count0 freed\n");
  fflush(stdprof);
  fclose(stdprof);
  fflush(stdprof2);
  fclose(stdprof2);
}
        
int pure_sin(void)
{
  int i,t;
  double y;
  FILE *sinus;
    
    
  sinus=fopen("standard.st","w");
    
  for(i=0;i<nbin;i++)
    {
      y=cos(i*pi/(nbin/2)-4.*pi/nbin/2);
      fprintf(sinus,"%d %e\n",i,y);
      printf("y=%e\n",y);
    }
  fclose(sinus);
  return(0);
}
    
    
int *toa(FILE *toaptr,FILE *toatempo2, char *dfname,double *count1,int nbin,double *prfamps, int simjd_m, double sfmjd_m, double spobs_m ,double freq, double *terr,char *site,double EMAX,double maxbin, double *deltaph,double *refph1)
{
  double phase_toa,shift,eshift,ephase,snrfft,esnrfft,b,smax,*glbprf,*glbprf2, offsetph, ishift, *shifted;
  int j,k,intdph, intdph_ref;
  int refpulseph; 
  double dmjd,offset;
  char dmjd_s[16];
  
  indivpulse = fopen("indivpulse.dat","a");
  indivphase = fopen("indivphase.dat","a");
  indivnorm = fopen("indivnorm.dat","a");

  indiv_std=  fopen("indiv_std.dat","a");
  noshift= fopen("noshift.dat","a");
  phase=fopen("phase.dat","a");
  test=fopen("test.dat","a");

  refpulseph = (int)(fmod(refph, 1.0)*nbin)+1;
  fftconv_(&nbin,prfamps,count1,&shift,&eshift,&snrfft,&esnrfft,&b,&rms); 
  if (shift < 0.0) (double)(shift += 1.0*nbin);
  phase_toa= (double)(shift/(1.0*nbin));
  ephase = eshift / (double)nbin;
  /* Referring the mjd to the position of the standard (which should have its peak at 0) */
  dmjd  = (double) sfmjd_m + (double)(phase_toa * spobs_m/(86400.0));
  *terr  = ephase * spobs_m * 1000.0 * 1000.0;
    
  if (dmjd > 1.0) {simjd_m++; dmjd-=1.0;}
  if (dmjd < 0.0) {simjd_m--; dmjd+=1.0;}
  
  sprintf(dmjd_s,"%15.13lf",dmjd);
  sprintf(mjd_s,"%d",simjd_m);
  strcat(mjd_s,dmjd_s+1);
  
  smax=0;
  for(j=0;j<nbin;j++)
    if (smax < prfamps[j]) {smax = prfamps[j]; maxbin = 1.0*j;} 
  
  //printf("*terr=%8.1f\n\n",*terr);
  if(*terr<EMAX) //empirical selection criterion for good pulses (Jake)
    {
      
     
      fprintf(test,"%19s %16.15lf %15.14lf %lf %lf %lf %lf %lf %lf\n",mjd_s,bary_f[0],spobs_m,shift,eshift,snrfft,esnrfft,b,rms);
       if(presto_switch)
	prepare_info(count_info);
       
      fflush(test);
      fprintf(toaptr," %17s %3d %2d %8.3f  %19s    %4.2f%8.1f        %1s\n",dfname, 0, 0, freq, mjd_s, 0.0,*terr,site);
      fflush(toaptr);
      /*TEMPO2 FORMAT*/
      fprintf(toatempo2,"200100503.tt    %8.3f   %19s  %7.2f i -tel wsrt\n",freq, mjd_s,*terr);
      fflush(toatempo2);

      offsetph=0.;
      offsetph=fmod(*refph1,1.0);
      
      intdph=0;
      intdph=(int)(*deltaph+offsetph*nbin);
      if(intdph<0)
        intdph+=nbin;
      if(intdph>=nbin)
        intdph-=nbin;
    

      if(numfiles==0 && first_toa_call == 1 )
      	{intdph_ref=0; first_toa_call=0;}

      glbprf = (double *) calloc(nbin,sizeof(double));
      glbprf2 = (double *) calloc(nbin,sizeof(double));
      shifted = (double *) calloc(nbin+1,sizeof(double));
      
      for(j=0;j<nbin;j++)
	{
	  if(j-intdph<0)
	    {
	      glbprf[j]+=prfamps[j-intdph+nbin];
	    }
	  else
	    glbprf[j]+=prfamps[j-intdph];
      	}
      
      for(j=0; j<nbin; j++)
	glbprf2[j]=prfamps[j];           

      ishift = (double) ((int)(shift+intdph_ref));
      shiftbyfft_(&glbprf2[0],&nbin,&maxbin);

      for(j=0;j<nbin;j++)
        {	  
	  k=j;
          if(tdmp_o>10)
	    fprintf(indivpulse,"%2d %10.3lf %19s %5d\n",j,prfamps[k],mjd_s,finalpcu); 
	  totp[j]+=glbprf2[j];
	  fprintf(noshift,"%2d %10.3lf %19s %5d\n",j,glbprf[j],mjd_s,finalpcu);
        }
      fprintf(noshift,"\n\n");
      if(tdmp_o>10)
	fprintf(indivpulse,"\n\n");          

      fprintf(phase,"%lf %16.15lf %16.15lf %16.15lf %19s %16.15lf %16.15lf %16.15lf\n",shift, phase_toa,ephase,spobs_m,mjd_s,sfmjd_m, refph, dphase);
      k=0;
      for(j=0; j<nbin; j++)
	{
	  if(j==0)
	    {
	    k=ceil(shift);
	    }	  
	  if(k>nbin-1)
	    k-=nbin;
	  shifted[j]=prfamps[k];
	  k++;
	}
      for(j=0; j<nbin; j++)
	{
	  fprintf(indivphase,"%2d %10.3lf %19s %lf\n",j,shifted[j],mjd_s, tdmp_o);
	  fprintf(indivnorm,"%2d %10.3lf %19s %lf\n",j,shifted[j]/smax,mjd_s, tdmp_o);
	}
      fprintf(indivphase, "\n\n");
      fprintf(indivnorm, "\n\n");
      fflush(indivphase);
      fflush(indivnorm);
    }
  else
    {
      num--;
    }
    
  free(glbprf);
  free(glbprf2);
  free(shifted);
  fclose(indivpulse);
  fclose(noshift);
  fclose(indiv_std);
  fclose(indivphase);
  fclose(indivnorm);
  fclose(phase);
  fclose(test);
  return(0);
}
    
    
    
int dump(int column, int tformat)
    
{
  FILE *output, *output2;
  int i,binyesno, imjd;
  char mjd[100],dmjd[50];
  double fmjd;

  column=1;    
  selectrows(tsamp[0],column, incl_excl,numfiles);  
  strcat(infilename,".dmp");
  output=fopen(infilename,"w");
  output2=fopen("barytimes.dat","w");
    
  tformat=1;
  if(tformat==1)
    {
      i=0;
      while(i<noutrows)
        {
          imjd = (long int) (bary_i[i]/86400.); 
          fmjd = (double) (bary_i[i]/86400.) ;
          fmjd-= (double) imjd;
          fmjd+= (double) (bary_f[i]/86400.) ;
    
    
          sprintf(dmjd,"%18.17lf",fmjd);
          sprintf(mjd,"%d",imjd);
          strcat(mjd,dmjd+1);
          fprintf(output,"%24s\n",mjd);
          fprintf(output2,"%lf\n",bary_i[i]+bary_f[i]); 
          i++;
        }
      
    }

  fclose(output);
  fclose(output2);
    
  return(tformat);
}
    
    
    
    
    
void timelags(void)    
{
  //  FILE *tlags;
  count2 =  (double *) calloc(nbin,8);
  stdamps = (double *) calloc(nbin,8);
  prfamps = (double *) calloc(nbin,8);
  maxbin2 = (double *) malloc(8);
  terr    = (double *) malloc(8);
    
  //  tlags=fopen("tlags.dat","a");
  printf ("Allocating count2, stdamps, prfamps, maxbin2, terr\n");
  printf("Using the previous standard profile\n");
  for(j=0;j<nbin;j++) 
    {
      count1[j]=standard[j];
      
    }
  for(j=0;j<nbin;j++)
    {
      if (smax < count1[j]) {smax = count1[j]; maxbin = 1.0*j;}
    }
    
  free(count2);
  free(stdamps);
  printf("Pointers *stdamps and *count2 freed\n");
    
  numbr=0;
  for (i=0;i<ndmp_o;i++)
    {
      /* Need to calculate the mjd for the pulse arrival phase for each subint */
    
    
      simjd_m = (long int) ((bary_i[0]/86400.)+bary_f[0]/86400.)+(double)(((tdmp_o * i)+tdmp_o/2.)/86400.0);
      //              var_int = (int) (((tdmp_o * i)+tdmp_o/2)/86400.0);
    
      sfmjd_m = (double) bary_i[0]/86400.;
      sfmjd_m+= bary_f[0]/86400.;
      sfmjd_m+= (double)(((tdmp_o * i)+tdmp_o/2.)/86400.0);  
      sfmjd_m-=simjd_m;
      
      printf("simjd_m = %ld sfmjd_m = %16.15f\n",simjd_m,sfmjd_m);
      if (sfmjd_m > 1.0) {simjd_m += 1; sfmjd_m -= 1.0;}
      ppolyco(dfname,simjd_m,sfmjd_m,&spobs_m,&phs); 
      
      dphase = fmod(phs - refph,1.0);
      sfmjd_m = sfmjd_m - dphase*spobs_m/86400.0;
      if (sfmjd_m > 1.0) {simjd_m += 1; sfmjd_m -= 1.0;}
      for (j=0;j<nbin;j++)
        {
          if(bkgyesno=='y')
            prfamps[j] = opr[i*nbin+j]- bkgcount[numfiles]*tdmp_o/nbin;
          else 
            prfamps[j] = opr[i*nbin+j];
        }
      
      numb++;
      numbr++;
      num++;
      toa(toaptr,toatempo2, dfname,count1,nbin,prfamps,simjd_m,sfmjd_m,spobs_m,freq,terr,site,EMAX,maxbin,deltaph,refph1);
      
      /* SHAPE of each individual pulse of tdmp_o seconds  */         
      if(offon==0 && i==ndmp_o-1)
        { 
          free(maxbin2); 
          free(terr); 
          printf("Pointer *maxbin2 and *terr freed\n");
        }            
    } // End of "ndmp_o" loop
  free(bary_f);
  free(bary_i);
  free(opr);
  free(prfamps);
  free(time_col);
  return;
    
}
    
float gammq(float aa, float xx)
{
        void gcf(float *gammcf, float a, float xx, float *gln);
        void gser(float *gamser, float a, float xx, float *gln);
        void nrerror(char error_text[]);
        float gamser,gammcf,gln;
    
        if (xx < 0.0 || aa <= 0.0) nrerror("Invalid arguments in routine gammq");
        if (xx < (aa+1.0)) {
                gser(&gamser,aa,xx,&gln);
                return 1.0-gamser;
        } else {
                gcf(&gammcf,aa,xx,&gln);
                return gammcf;
        }
}
    
int check_pcus(double counts_pcu[])
{
  
  int counter_pcus; 

  prof_act_pcus=0;
  for(counter_pcus=0; counter_pcus<5; counter_pcus++)
    {
      if(counts_pcu[counter_pcus])
	prof_act_pcus++;
    }
  return(0);
}
     
void prepare_info(int i)
{
  FILE *info;
  char geninfo[100];


  info=fopen("input_info", "w");
  
  fprintf(info,"%s%d_TS%d\n",root,numfiles+1,i+1);
  fprintf(info,"8\n");
  fprintf(info,"X\n");
  fprintf(info,"HETE J1900\n");
  fprintf(info,"19:11:16.05\n");
  fprintf(info,"00:35:05.80\n");
  fprintf(info,"A. Patruno\n");
  fprintf(info,"%15.14lf\n",bary_i[0]/86400.+bary_f[0]/86400.);
  fprintf(info,"1\n");
  fprintf(info,"%lf\n",floor(tdmp_o/tsamp[0])-even_bin);
  fprintf(info,"%14.13lf\n",tsamp[0]);
  fprintf(info,"0\n");
  fprintf(info,"5\n");
  fprintf(info,"3600\n");
  fprintf(info,"7\n");
  fprintf(info,"14\n");
  fprintf(info,"A.Patruno\n");
  fprintf(info,"\n");
  fflush(info);
  fclose(info);
  sprintf(geninfo,"makeinf < input_info");
  system(geninfo);

  return;

}
