c 
c $Log: fftfit.f,v $
c Revision 1.2  1996/03/15 01:14:00  jsandhu
c removed Id string (Log string is sufficient).
c
c Revision 1.1  1994/11/25 03:22:44  jsandhu
c Initial revision
c

*       fftfit.f,v 1.2 94/01/28 14:11:51 will 
C       @(#)fftfit.f	1.1 9/7/90
	subroutine fftfit(prfamps,amp,phase,nbin,shift
     *,eshift,snrfft,esnrfft,b,rms)
	
C       Fourier transform domain routine for determining pulse TOAs.
C  Input data:
C	prfamps(nbin)	profile
C	s(nh)		standard profile amplitude
C	phase(nh)	standard profile phase
C	nbin		length of prof

C  Outputs:
C	shift		shift required to align std with prof, in bins
C	eshift      	uncertainty in shift
C	snrfft		signal to noise ratio of prof
C	esnrfft		uncertainty in snrfft

C  Method:
C  It is assumed that prfamps(j)=a + b*std(j-tau) + noise(j), where the
C  (j-tau) subscript is supposed to mean that the time-shift between the
C  observed and standard profiles need not be an integer number of bins.
C  The algorithm is a straightforward Chi-squared minimization with respect
C  to a, b, and tau.  (The result for a, being of no physical or instrumental
C  interest, is not actually evaluated -- though it easily could be.)
C  First and second partial derivatives of Chisqr with respect to b and tau
C  are computed, and used to evaluate s, snrfft, and their "standard errors,"
C  or one-sigma uncertainties.  The only special trick is that the expression
C  for the best-fit value of tau defines it implicitly, and cannot be solved
C  analytically.  It is solved numerically instead, finding the minimum of
C  Chisqr near a best guess from a CCF at 32 lags done in the Fourier domain.

C  Also note that it may
C  be desirable to return the scale factor b relating prof to std, instead
C  of snrfft.  In that case you could also return the noise estimate rms.

*	MODIFICATIONS:
*		real*4 -> real*8.
*		MAXSAM -> nbin in some declarations.
*		Check that nbin <= MAXSAM.  WD 14aug92.

***
	implicit none 
	real*8 twopi, fa,fb,fc,s1,s2,s3,shift,eshift,snrfft,esnrfft
	real*8 errb,errtau,rms,sq,cosfac,b,a,ftau,edtau,dtau,tau,fac
	real*8 dchisqr,zbrent
	integer MAXSAM,nbin,ia,ib,ic,isum,nsum,ntries,k,nh
***

	parameter (twopi=6.2831853,MAXSAM=8192)
	real*8 prfamps(nbin),p(nbin/2),theta(nbin/2)
	real*8 amp(nbin/2),phase(nbin/2),r(nbin/2),tmp(nbin/2)
	complex cp(0:nbin/2)
	logical low,high

CCCCC	if (phase(1).ne.0) then
CCCCC	  write(6,*) 'fftfit: Phase of fundamental not zero.'
CCCCC	  call exit(1)
CCCCC	end if
	if (nbin .gt. MAXSAM) then
	  write(6,*) 'fftfit: nsam > MAXSAM = ', MAXSAM
	  write(6,*) 'Recompile fftfit.f with larger MAXSAM'
	  call exit(1)
	endif


	nh=nbin/2
	call cprof(prfamps,nbin,nh,cp,p,theta)
	do 10 k=1,nh
	 tmp(k)=p(k)*amp(k)
10	 r(k)=theta(k)-phase(k)

	fac=nbin/twopi

	call fccf(tmp,r,shift)

C  The "DO 60" loop solves the transcendental equation yielding the best-fit
C  value of tau.  Here the number starts at 16 (number used in CCF)
C  and increases by factors of 2, at each step finding the function
C  zero closest to the starting point from the previous iteration.

	tau=shift
c	print *, 'Uncorr tau = ',tau
	do 60 isum=1,99
	nsum=2.0**isum
	if(nsum.gt.nh) go to 70
	dtau=twopi/(nsum*5)
	edtau=1./(2.*nsum+1.)
	if (nsum.gt.(nh/2.+.5)) edtau=1.e-4

	ntries=0
	low=.false.
	high=.false.
50	ftau=dchisqr(tau,tmp,r,nsum)
	ntries=ntries+1
	if(ftau.lt.0.0) then
	  a=tau
	  fa=ftau
	  tau=tau+dtau
	  low=.true.
	else
	  b=tau
	  fb=ftau
	  tau=tau-dtau
	  high=.true.
	end if
	if (ntries.gt.10) then
	  shift=0.
	  eshift=999.
	  snrfft=0.
	  esnrfft=0.
	  return
	end if
	if (low.neqv.high) go to 50
	tau=zbrent(a,b,fa,fb,edtau,tmp,r,nsum)
60	continue
	

70	s1=0.
c	print *, 'Corr tau = ',tau
	s2=0.
	s3=0.
	do 80 k=1,nh
	cosfac=cos(-r(k)+k*tau)
	s1=s1 + tmp(k)*cosfac
	s2=s2 + amp(k)**2
	s3=s3 + k**2 *tmp(k)*cosfac
c       print *, 'k',k,'tmp= ',tmp(k),'s3 = ',s3,'cosfac = ',cosfac
c   	print *, 'sum = ',k**2 *tmp(k)*cosfac
 80	continue
	b=s1/s2
	s1=0.
	do 90 k=1,nh
	sq=p(k)**2-2.*b*p(k)*amp(k)*cos(r(k)-k*tau)+(b*amp(k))**2
90	s1=s1+sq
	rms=sqrt(s1/nh)


	errb=rms/sqrt(2.0*s2)
	errtau=rms/sqrt(2.0*b*s3)

	snrfft=2.0*sqrt(2.0*nh)*b/rms
c	write(*,*) 'b',b,'s1',s1,'s2',s2

c Ale: transforms shift and eshift in units of bins instead of radians
c Ale: fac=nbin/2pi

	shift=fac*tau   
	eshift=fac*errtau

	esnrfft=snrfft*errb/b
c	write(*,*) 'eshift=',eshift,'errtau=',errtau,'rms=',rms
c	write(*,*) 'b=',b
c	write(*,*) 's1=',s1
c	write(*,*) 's2=',s2
c	write(*,*) 's3=',s3
c	write(*,*) 'cosfac=',cosfac,'k=',k
c	do k=1,16
c	   write(*,*) 'tmp*cosfac', tmp(k)*cos(-r(k)+k*tau)
c	enddo
c        write(*,*) 's(k)=',s
	return
	end
