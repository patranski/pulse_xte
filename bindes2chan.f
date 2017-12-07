**
c** Subroutine bindes2chan
c**      
      subroutine bindes2chan(descriptor,ichans,ichane,nchan)

      character*2048 descriptor, temp*16
      integer ichans(*),ichane(*)


      print *, 'Ale test'
      ipos1 = istrpos(descriptor,'C[')

      ii = 1
      do while (ipos1 + ii .le. len_used(descriptor) )
         ii = ii +1
         ipos2 = istrpos(descriptor(ipos1:ipos1+ii),']')
         if (ipos2 .ne. 0) goto 111
      enddo
      stop ' End of channel description not found'
 111  ipos2 = ipos2 + ipos1

c      read(descriptor(ipos2+1:ipos2+1),'(i1)') ibit_chan

c      nchan = 2**ibit_chan
c      write (*,'('' This data file has: '',i3,'' channels'')') nchan

C     *** FIND OUT WHICH CHANNELS ARE USED.
C     This piece of unreadable code is used to find out
C     which original channels correspond to the rebinned channels.
C     Note that rebinned channels are numbered 1 to 64 (if there
C     are 64 channels of course), while the original channels
C     are numbered 0 to 255.
      i = 1
      ii = ipos1+3
      ilast = ipos1+2
      ic = nchan
      do while (ii .le. ipos2)
         ii = ii+1
         if (istrpos(descriptor(ii:ii),',') .ne. 0
     * .or. istrpos(descriptor(ii:ii),']') .ne. 0) then
C     *** A comma is found (or a ] at the end)
            read(descriptor(ilast:ii-1),'(a)') temp
            if (istrpos(temp,'~') .ne. 0) then
               itilde= istrpos(temp,'~')
               if (istrpos(temp,'(') .eq. 0) then
               read(temp(1:itilde-1),'(I1)') ichans(i)
               read(temp(itilde+1:16),'(I1)') ichane(i)
               i = i+1
            else
               isemi=istrpos(temp,';')
               read(temp(2:itilde-1),'(I1)') ia
               read(temp(itilde+1:isemi-1),'(I1)') ib
               read(temp(isemi+1:len_used(temp)-1),'(I1)') ic
               do ichan =ia,ib,ic
                  ichans(i) = ichan
                  ichane(i) = ichan+ic-1
                  i = i +1
               enddo
            endif
            elseif (istrpos(temp,':') .ne. 0) then
               icolon= istrpos(temp,':')
               read(temp(1:icolon-1),'(I1)') imin
               read(temp(icolon+1:16),'(I1)') imax
               do ichan=imin,imax
                  ichans(i) = ichan
                  ichane(i) = ichan
                  i =i+1
               enddo
            else
               read(temp,'(I1)') ichans(i)
               ichane(i) = ichans(i)
               i = i+1
            endif
            ilast=ii+1
         endif
      enddo

c**
c     *** SANITY CHECK: see if we have found as many channels as is in other
c     *** keyword. If not this can be caused strange data selection etc.
c     *** at the moment TDDES2 is used. Might be wrong.
c** 
      if (i-1 .ne. nchan) then
         print *, 'Error in determination of channels from descriptor'
         print *, 'Check that your data selection is consistent'
         print *, 'If not: contact Tim to debug unreadable code ;-)'
         stop
      endif

      return
      end
