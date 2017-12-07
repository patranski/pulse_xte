c $Id: shiftbyfft.f,v 1.1 1994/11/25 03:11:31 jsandhu Exp $
c $Log: shiftbyfft.f,v $
c Revision 1.1  1994/11/25  03:11:31  jsandhu
c Initial revision
c
c Revision 1.1  1994/11/25  03:11:31  jsandhu
c Initial revision
c

*       shiftbyfft.f,v 1.1 93/10/29 22:56:53 will 
*	****************************************************************
	subroutine shiftbyfft(prf, nphase, shift)
*	****************************************************************

	implicit none 

c	include '../include/implicit.h'

*	Phase-shifts a profile.

*	INPUT: Length of profile
	integer nphase

*	INPUT/OUTPUT: Profile to shift
	real*8 prf(0:nphase-1)

*	INPUT: Amount to shift (in units of phase bins)
	real*8 shift

*	* * * * * 

	integer MAXNPRF
	parameter (MAXNPRF=128)
	complex x(0:MAXNPRF-1)

	real*8 ampl, phase, shiftrad, temp(2)
	integer i

	double precision twopi
	parameter(twopi = 6.2831853071795865d0)

	if (nphase .gt. MAXNPRF) then
	    write(6,*) 'lll shiftout: nphase (', nphase,
     :			') > MAXNPRF (', MAXNPRF, ')'
	    call exit(1)
	endif

	temp(1) = 0
	temp(2) = nphase - 1
CCCCC	call display(nphase, 1, temp, prf, 1, 1,
CCCCC     :		'Sample', 'Amplitude', 'Before Phase-Shifting')

	do i = 0, nphase-1
	    x(i) = cmplx(prf(i), 0.)
	enddo

	call ffft(x, nphase, 1, 0)

	shiftrad = twopi * shift / nphase
*	X(0) doesn't change phase.  Only need to do 1..nphase/2-1 and
*	the corresponding negative frequencies.
C X(0) non cambia perche' questa e' la FFT di un numero reale. 
C Inoltre si somma fino a nphase/2-1 perche' le componenti
C A_n=A_(nphase-n) sono appunto uguali.
C Questo perche' per una funzione reale A si ha A(n)* = A(-n)
C e quindi nn c'e' bisogno di calcolare l'intero spettro, ma solo meta'.
 
 
	do i = 1, nphase/2-1
	    ampl = abs(x(i))
	    phase = 0.
	    if (ampl .gt. 0) phase = aimag(log(x(i)))
	    phase = mod(phase - i*shiftrad, real(twopi))
	    x(i) = ampl * cmplx(cos(phase), sin(phase))
	    x(nphase-i) = conjg(x(i))
	enddo

	call ffft(x, nphase, -1, 0)

CCCCC	call display(nphase, 2, temp, x, 1, 2,
CCCCC     :		'Sample', 'Amplitude', 'Phase-Shifted Complex Profile')

*	Transfer back to profile array and normalize the result
	do i = 0, nphase-1
	    prf(i) = Real(x(i)) / nphase
	enddo

CCCCC	call display(nphase, 1, temp, prf, 1, 1,
CCCCC     :		'Sample', 'Amplitude', 'Phase-Shifted Profile')

	return
	end
