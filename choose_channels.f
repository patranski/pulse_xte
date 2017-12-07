      subroutine choose_channels(metafiles,infilename,channels,channels1
     *,channels2,nchannels,nfiles,nmetafiles,max_files,max_metafiles,
     *max_channels,istart,iend)
c     
      implicit none
c
      integer       nfiles,nmetafiles,max_files
      integer       max_metafiles,max_channels
      character*180 metafiles
      character*380  infilename
      integer       channels(max_channels,max_files,max_metafiles)
      integer       channels1(max_channels,max_files,max_metafiles)
      integer       channels2(max_channels,max_files,max_metafiles)
      integer       temp1(8),temp2(8)
      integer       nchannels(max_files,max_metafiles)
      integer       nhistos(max_files,max_metafiles),len_used
c
      logical       sw_debug,swok
      common /debug/sw_debug
c
      integer   max_chan
      parameter(max_chan=256)
c
      integer       i, j, k, ichan, nrang, nranges, istart, iend
      integer       ichanstart(max_chan), ichanend(max_chan)
      character*16  cpar
      character*80  ask_string
      data swok/.true./
      character     option
c
  112 do j = 1, nmetafiles
         do i = 1, nfiles
            do k = 1, nchannels(i,j)
               channels(k,i,j) = 0
            enddo
         enddo
      enddo
c
      write(*,*)
      print *,'Overview of available channels'
      do j = 1, nmetafiles
         write(*,*)
         write(*,'(8(I3,A1,I3,A1))') (channels1(k,1,j),'-'
     &        ,channels2(k,1,j),' ',k=1,nchannels(1,j))
      enddo
      write(*,*)
      write(*,*)
      print *,'             Give ABSOLUTE channel numbers !!!'
c
      cpar = 'ALL'
      write (*,*) 'Channel selection menu'
      write (*,*) 'All channels summed (A) '
      write (*,*) 'Continuous range (R) '
      write (*,*) 'Multiple ranges (M)'
      write (*,*)
      read  (*,*)  cpar
c     
      print *, cpar 
            
      if (cpar .eq. 'A') then
         do j = 1, nmetafiles
            do i = 1, nfiles
               do k = 1, nchannels(i,j)
                  channels(k,i,j) = 1
               enddo
            enddo
         enddo
      elseif (cpar .eq. 'R') then
 113     write (*,*) 'Start channel for FFT'
         read (*,*) istart
         write (*,*) 'End channel for FFT'
         read (*,*) iend
         if (iend .lt. istart) goto 113
         do ichan = istart, iend
            do j = 1, nmetafiles
               do i = 1, nfiles
                  do k = 1, nchannels(i,j)
                     if ((ichan .ge. channels1(k,i,j)) .and.
     &                    (ichan .le. channels2(k,i,j)))  then
                     channels(k,i,j) = 1   
                  endif
                  enddo
               enddo
            enddo
         enddo
      elseif (cpar .eq. 'M') then
         write(*,*) 'How many ranges ? '
         read (*,*) nranges
 115     do nrang= 1, nranges
            write(*,*) 'Start channel for range',nrang,' ='
            read (*,*) ichanstart(nrang)
            write(*,*) 'End channel for range',nrang,' ='
            read (*,*) ichanend(nrang)
            if (ichanend(nrang) .lt. ichanstart(nrang)) goto 115            
         enddo

       do nrang = 1, nranges
            do ichan = ichanstart(nrang), ichanend(nrang)
               do j = 1, nmetafiles
                  do i = 1, nfiles
                     do k = 1, nchannels(i,j)
                        if ((ichan .ge. channels1(k,i,j)) .and.
     &                      (ichan .le. channels2(k,i,j)))
     *                       channels(k,i,j) = 1
                     enddo
                  enddo
               enddo
            enddo
         enddo         
      endif
      
c     Place the cutted part here 
      
      i=0
      j=1
      do k=1, nchannels(1,j)
         if (channels(k,1,j) .gt. 0) then
            i = i + 1
            temp1(i) = channels1(k,1,j)
            temp2(i) = channels2(k,1,j)
            if (i.eq.8) then
               write(*,'(8(I3,A1,I3,A1))')
     *              (temp1(i),'-',temp2(i),' ',i=1,8)
               i = 0
            endif
         endif
      enddo
      if (i.gt.0) write(*,'(8(I3,A1,I3,A1))')
     *     (temp1(k),'-',temp2(k),' ',
     &     k=1,i)
     
      write(*,*)
c     
      
      return
      end
