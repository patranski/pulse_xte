c $Id: cprof.f,v 1.1 1994/11/25 03:22:13 jsandhu Exp $
c $Log: cprof.f,v $
c Revision 1.1  1994/11/25  03:22:13  jsandhu
c Initial revision
c
c Revision 1.1  1994/11/25  03:22:13  jsandhu
c Initial revision
c

* cprof.f,v 1.1 93/10/29 22:51:06 will 
C @(#)cprof.f	1.1 9/7/90
	subroutine cprof2(count1,nbin,nh,amp,phase)

C  Compute FFT of profile in array y(nbin), and return amplitude and phase
C  in arrays amp(nh) and phase(nh).  Note that nh=nbin/2, and that the DC term
C  is returned in c(0), fundamental in c(1), ..., Nyquist freq in c(nh).
*	MODIFICATIONS:
*		real*4 -> real*8.
*		Checks the MAXSAM limit. WD 14aug92.

	implicit none

	integer MAXSAM,i,nbin,nh
	parameter (MAXSAM=8192)
	real*8 count1(nbin),amp(nh),phase(nh)
	complex space(0:nh),temp(MAXSAM)

	if (nh .ge. MAXSAM) then
	    write(6,*) 'cprof(): nh > MAXSAM: ', nh, MAXSAM
	    call exit(1)
	endif
	do 10 i=1,nh
10	 temp(i)=cmplx(count1(2*i-1),count1(2*i))
	call ffft(temp,nbin,1,1)
	space(0)=temp(1)
	do 20 i=1,nh
	 space(i)=temp(i+1)
	 amp(i)=cabs(space(i))
	 phase(i)=0.
20	if(amp(i).gt.0.) phase(i)=aimag(clog(space(i)))
	return
	end
