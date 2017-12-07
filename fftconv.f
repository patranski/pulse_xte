c $Id: fftconv.f,v 1.1 1994/11/25 03:15:16 jsandhu Exp mbailes $
c $Log: fftconv.f,v $
c Revision 1.1  1994/11/25 03:15:16  jsandhu
c Initial revision
c
c------------------------------------------------------------
	subroutine fftconv(nbin,prfamps,count1,shift,eshift,snrfft
     *       ,esnrfft,b,rms)
c
c       MXB 18-2-1994
c


c       count1= standard amplitudes
	implicit none

	integer nbin              ! Nbins in both profile + standard
	real*8 prfamps(nbin)             ! profile to be convolved
	real*8 count1(nbin)            ! standard profile
	real*8 shift                ! shift in bins between pr+count1
	real*8 eshift               ! error in that shift
	real*8 snrfft                  ! signal to noise ratio
	real*8 esnrfft                 ! error in signal to noise ratio
	real*8 b                      ! Ale: amplitude scale factor
	real*8 rms                    ! Ale: noise contribution in amplitude

c	integer nbinmax
	integer nh
c	parameter (nbinmax=4*8192)
	parameter (nh=64)

	complex space(0:nh)
	real*8 amp(nh),phase(nh)

c	if (nbin.gt.nbinmax) then
c	   write(*,*) 'Increase nbinmax in fftconv and recompile'
c	   write(*,*) 'nbinmax = ',nbinmax
c	   stop
c	end if

	call cprof(count1,nbin,nh,space,amp,phase)
	call fftfit(prfamps,amp,phase,nbin,shift,eshift,snrfft
     *,esnrfft,b,rms)
c	write(*,*) 'eshift = ',eshift
	return
	end












