c $Id: fccf.f,v 1.1 1994/11/25 03:17:24 jsandhu Exp $
c $Log: fccf.f,v $
c Revision 1.1  1994/11/25  03:17:24  jsandhu
c Initial revision
c
c Revision 1.1  1994/11/25  03:17:24  jsandhu
c Initial revision
c

* fccf.f,v 1.1 93/10/29 22:51:12 will 
C @(#)fccf.f	1.1 9/7/90
	subroutine fccf(tamp,pha,shift)

C  Input: Fourier-domain CCF of two profiles: tamp = A1*A2, pha = phi1 - phi2.
C  Output: Using 16 harmonics of tamp and pha array, finds maximum of CCF
C  at 64 lags over the pulsar period, and returns value of shift in radians.
*	MODIFICATIONS:
*		real*4 -> real*8.

***
	implicit none
	
	real*8   twopi,imax,ib,shift
	integer nbin,maxsam,rc,cmax,fa,fb,fc, nh,i,ia,ic

***

	parameter (nbin=32)
	parameter (MAXSAM=8192,twopi=6.2831853)
	real*8 tamp(nbin/2),pha(nbin/2)
	complex ccf(0:127)

	nh=nbin/2
	ccf(0)=(0.,0.)
	do 10 i=1,nh/2
	 ccf(i)=cmplx(tamp(i)*cos(pha(i)),tamp(i)*sin(pha(i)))
10	 ccf(nbin-i)=conjg(ccf(i))
	do 20 i=nh/2+1,nh
	 ccf(i)=(0.,0.)
20	 ccf(nbin-i)=(0.,0.)
	call ffft(ccf,nbin,-1,0)
	cmax=-1.e30
	do 30 i=0,127
	 rc=real(ccf(i))
c	 print *, 'rc=',rc
	 if (rc.gt.cmax) then
	   cmax=rc
	   imax=i
	 end if
30	continue

	fb=cmax
	ia=imax-1
	if(ia.eq.-1) ia=nbin-1
	fa=real(ccf(ia))
	ic=imax+1
	if(ic.eq.nbin) ic=0
	fc=real(ccf(ic))
	if ((2*fb-fc-fa).ne.0) then
	  shift=imax+0.5*(fa-fc)/(2*fb-fc-fa)
	else
	  shift=imax
	end if
	if(shift.gt.nh) shift=shift-nbin
	shift=shift*twopi/nbin
c	print *, 'Shift =',shift
	return
	end
