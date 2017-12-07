
C                          $Author: martin $
C                          $Date: 1998/02/16 17:06:44 $
C                          $Id: des2chan.f,v 1.8 1998/02/16 17:06:44 martin Exp $
C     $Source: /home/xray/sources/subroutines/des2chan.f,v $
      
      subroutine des2chan(descriptor,ichans,ichane,nchan,ioffset,ibit_ch
     *an)
C     ** subroutine returns the start and end channels as described in the
C     ** string descriptor in the array ichans and ichane
C     ** for a description of how the descriptor should look see the XFF documentation
C     
C     nchan = number of channels (output)
C     ioffset = "start" bit of event which deals with channel information (output)
C     ibit_chan = number of bit that deals with channel information (output)
      
      common /debug/sw_debug
      logical sw_debug
      
      character*2048 descriptor, temp*16
      integer ichans(*),ichane(*)

C--------------------------------------------------------------------------
C These comments have been added by Mariano on Feb 24 1996.
C The changes were necessary to be able to parse Good_Xenon files.
C After these changes fft_event works for this mode too.
C--------------------------------------------------------------------------
C
C This routine deals with the "Event descriptor" of the FITS file.
C This descriptor is stored in the FITS Header in the Keyword
C TEVTB. FFT program reads this keyword, and pass it here.
C Hopefully this routine can deal with it.
C
C Here are couple of examples of the "Event descriptor" as read from 2
C typical FITS files.
C
C This is for the Event Mode 8C_0 (see Appendix F, p. 191; I've splitted
C it into several lines, although the program reads it in long string):
C
C (M[1]{1},D[0:4]{3},C[0~12,13~16,17~20,21~25,26~31,32~39,40~55,56~249]{3}
C ,T[0.0:0.00390625;7.62939453125e-06]{9})^(M[127]{8},S[One5]{5}
C ,S[FirstFlag]{1},S[SpillFlag]{1},S[Zero1]{1})^(M[1]{8},S[LostEvents0]
C {8})^(M[2]{8},S[LostEvents1]{8})^(M[4]{8},S[Zero]{7},S[Spillage]
C {1})^(M[8]{8},S[Zero]{5},S[ModeSpecific]{3})'                         
C
C
C This is for a Good_Xenon Mode (Good_Xenon should always have 256
C channels, 0 to 255; see Appendix F, p. 194):
C
C (M[1]{1},S[Zero]{6},D[0:4]{3},E[0:63]{6},C[0:255]{8})
C
C In both cases this routine deals with the part between C[]{}, to
C figure out how the energy channels are binned.
C
C A short (and certainly incomplete) explanation of what you may get
C here, and its meaning:
C
C A descriptor like C[0~4,5,6,7~9,10:12,13~255]{3} tells you that:
C
C
C  1. Your observing mode has 8 Channels (the {3} means 2**3 channels).
C  2. The channels that you have are formed as follows:
C
C               Original PCA Channles   ---->  Rebinned to channel
C
C                  0 to   4             ---->          1
C                         5             ---->          2
C                         6             ---->          3
C                  7 to   9             ---->          4
C                        10             ---->          5
C                        11             ---->          6
C                        12             ---->          7
C                 13 to 255             ---->          8
C
C For more information on the Data Description Language (DDL), see:
C
C http://heasarc.gsfc.nasa.gov/docs/xte/abc/data_files.html#event
C
C   or
C
C http://heasarc.gsfc.nasa.gov/docs/xte/abc/ddl.html
C
C--------------------------------------------------------------------------

C
C First, let's find the beginning and end of the channel description
C within the whole descriptor.
C
C This is the beginning:
C (this comment was written by Mariano; the code was written by Tim).
C
      ipos1 = istrpos(descriptor,'C[')
      ii = 1
      do while (ipos1 + ii .lt. len_used(descriptor))
         ii = ii +1
C
C And this is the end
C (this comment was written by Mariano; the code was written by Tim).
C
         ipos2 = istrpos(descriptor(ipos1:ipos1+ii),']')
         if (ipos2 .ne. 0) goto 111
      enddo
      stop ' End of channel description not found'
  111 ipos2 = ipos2 + ipos1
C
C Now read what's between `{}'. This gives the number of (binned) channels
C for the mode.
C (this comment was written by Mariano; the code was written by Tim).
C                      
      iii = 0
      do while (ipos2 + iii .lt. len_used(descriptor))
         iii = iii + 1
         ipos3 = istrpos(descriptor(ipos2:ipos2+iii),'{')
         if (ipos3 .ne. 0) goto 112
      enddo
  112 ipos3 = ipos2 + ipos3
c
      read(descriptor(ipos3:ipos3),'(i1)') ibit_chan
      nchan = 2**ibit_chan
      if (sw_debug) 
     &    write (*,'('' This data file has: '',i3,'' channels'')') nchan

C     *** FIND OUT WHICH CHANNELS ARE USED.
C     This piece of unreadable code is used to find out
C     which original channels correspond to the rebinned channels.
C     Note that rebinned channels are numbered 1 to 64 (if there
C     are 64 channels of course), while the original channels
C     are numbered 0 to 255.
C this comment was written by Tim.
C
       i = 1
      ii = ipos1+3
      ilast = ipos1+2
      ic = nchan
           
      do while (ii .le. ipos2)
         ii = ii+1
C
C Let's check if there is a `,' or a `]'.
C (this comment was written by Mariano; the code was written by Tim).
C             
         if(istrpos(descriptor(ii:ii),',').ne.0
     &        .or.istrpos(descriptor(ii:ii),']') .ne.0) then
C     *** A comma is found (or a ] at the end)
C this comment was written by Tim.
            read(descriptor(ilast:ii-1),'(a)') temp
            
C     
C     No `~' or `:', and as it was no `,', then it's a single channel, so
C     no rebinning applyied.
C     (this comment was written by Mariano; the code was written by Tim).
C     
            
            if( istrpos(temp,'~') .eq. 0) then 
               if (istrpos(temp,':') .eq. 0) then
                  read(temp,*) ichans(i)
                  ichane(i) = ichans(i)
               else
                  
C     
C     Here we deal with the `:'. Channels before and after `:' are not binned, so
C     you keep the instrument resolution there.
C     
C     In fact this is the only thing I've changed here.
C     Mariano.
C
                  icolon= istrpos(temp,':')
                  read(temp(1:icolon-1),*) istart
                  read(temp(icolon+1:16),*) iend
                  do j=istart,iend
                     ichans(i)=j
                     ichane(i)=j
                     i=i+1
                  enddo
                  i=i-1
               endif
            else
C     
C     Now we deal with the `~'. Channels before and after the `~' are binned
C     together.
C     (this comment was written by Mariano; the code was written by Tim).
C     
               itilde= istrpos(temp,'~')
c               print *, 'itilde ',itilde
c               print *, 'temp ',temp
c               print *, 'temp(1:itilde-1) ', temp(1:itilde-1)
c               print *, 'temp(itilde+1:16) ', temp(itilde+1:16)
c               print *, 'i ',i
               read(temp(1:itilde-1),*) ichans(i)
               read(temp(itilde+1:16),*) ichane(i)

            endif
            i = i+1
            ilast=ii+1
         endif
      enddo
      
     

C
C I didn't get into this (Mariano).
C

      ioffset = 1
      iposl=ipos1
      iposr=istrpos(descriptor(1:iposl-1),'}')
      iposl=istrpos(descriptor(1:iposl-1),'{')
            do while (iposl .ne. 0)
         read(descriptor(iposl+1:iposr-1),*) i
         ioffset = ioffset +i
         iposr=istrpos(descriptor(1:iposl-1),'}')
         iposl=istrpos(descriptor(1:iposl-1),'{')
      enddo
C
C If in debug mode, then write to the screen the rebinning.
C (this I added too, but it's harmless. Mariano)
C
      if (sw_debug) then
         write(*,*) 'Channel Binning is:'
         do i=1,nchan
            write(*,*) ichans(i),ichane(i)
         enddo
      endif
C
      return
      end
