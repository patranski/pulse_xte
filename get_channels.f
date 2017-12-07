      subroutine get_channels(infilename,channels1,channels2,nchannels,n
     *histos,ioffset,ibit_chan,nfiles,nmetafiles,max_files,max_metafi
     *les,max_channels)
c     
      implicit none
c     
      integer       nfiles,nmetafiles,max_files
      integer       max_metafiles,max_channels
C filenames is not necessary anymore here
c      character*180 filenames(max_files,max_metafiles)
      integer       channels(max_channels,max_files,max_metafiles)
      integer       channels1(max_channels,max_files,max_metafiles)
      integer       channels2(max_channels,max_files,max_metafiles)
      integer       nchannels(max_files,max_metafiles)
      integer       nhistos(max_files,max_metafiles)
      integer       ioffset(max_files,max_metafiles)
      integer       ibit_chan(max_files,max_metafiles)
c     
      logical       sw_debug
      common /debug/sw_debug
      character*280     infilename;
c
      integer   max_dim,   max_chan,     max_col
      parameter(max_dim=4, max_chan=256, max_col=256)
c
      integer       i, j, k, nhisto, nchans, ipos, istrpos
      integer       len_used, ntemp
      integer       status, unit, readwrite, blocksize, nkeys
      integer       nspace, hdutype,
     &              nfound
      integer       nrows, tfields, varidat, naxes(max_dim)
      real*8        timedel
      character*6   type
      character*16  instrument
      character*10  tform(max_col)
      character*32  ttype(max_col)
      character*80  tdim(max_col), tunit(max_col), temp_string, comment, 
     &              datamode, extname
c
      character*2048 event_description, tevtb(16), tddes(max_col)
      integer        ichans(max_chan), ichane(max_chan)
c


      print *, 'ciao'
      do 10 j = 1, nmetafiles
         do 20 i = 1, nfiles
            status = 0
            readwrite = 0
c
        
c Get an unused Logical Unit Number to use to open the FITS file
c
            call ftgiou(unit,status)
c
c Open the FITS file, with read-only access
c
            call ftopen(unit,infilename,readwrite,blocksize,status)
            if (status .ne. 0) call printerror(status)
            if (status .ne. 0) stop 'No FITS file in get_channels'
c
c Determine the number of keywords in the header of the FITS file'
c
            call ftghsp(unit,nkeys,nspace,status)
            if (status .ne. 0) stop 'Cannot determine nr. keywords'
c
c Try moving to the next extension in the FITS file, if it exists
c
            call ftmrhd(unit,1,hdutype,status)
            if (status .ne. 0) stop 'Cannot move to the next extension'
c
c Determine the dimension of the next extension in the FITS file
c
            call ftgknj(unit,'NAXIS',1,4,naxes,nfound,status)
            if (status .ne. 0) stop 'Cannot determine dimensions NAXIS'
c
            call ftghbn(unit,max_col, nrows, tfields, ttype, tform,
     &                  tunit, extname, varidat, status)
            if (status .ne. 0) 
     &          stop 'Cannot determine data type in the FITS file'
c
            call ftgkys(unit,'INSTRUME',instrument,comment,status)
            if (status .ne. 0) stop 'Cannot determine INSTUMENT'
c     
c Determine the datamode of the observation'
c
            call ftgkys(unit,'DATAMODE',datamode,comment,status)
            if (status .ne. 0) stop 'Cannot determine the obs. datamode'
    

            if (extname.eq.'XTE_HK') then
               type = 'xte_hk'
               write(*,*)
               print *,'Not a proper data file'
               print *,'Only know SA, SE, GX, HEXTE data files'
               print *,'XTE_HK will not be supported'
               print *,'Data file is: ', extname
               write(*,*)
               stop 'Check your metafiles!'
            endif
c     
            if ((extname.eq.'XTE_SA').and.
     *           (instrument(1:3).eq.'PCA')) then
               type = 'pca_sa'
            endif
            if ((extname.eq.'XTE_SE').and.
     *           (instrument(1:3).eq.'PCA')) then
               type = 'pca_se'
               if (datamode(1:9).eq.'GoodXenon') type = 'pca_gx'
            endif
            if ((extname.eq.'XTE_SA').and.
     *           (instrument(1:5).eq.'HEXTE')) then
               type = 'hxt_sa'
            endif
            if ((extname.eq.'XTE_SE').and.
     *           (instrument(1:5).eq.'HEXTE')) then
               type = 'hxt_se'
            endif
c     
            if (type.eq.'') then
               write(*,*)
               print *,'Not a proper data file'
               print *,'Only know SA, SE, GX, HEXTE data files'
               print *,'Data file is: ', extname
               write(*,*)
               stop 'Check your metafiles!'
            endif
c     
c Determine the time resolution of the observation'
c
            call ftgkyd(unit,'TIMEDEL',timedel,comment,status)
            if (status .ne. 0) stop 'Cannot determine time res'
c
c Determine the TDIM of the observation'
c
c
            call local_ftgkns(unit,'TDIM',2,tfields,tdim,nfound,status)
            if (status .ne. 0) stop 'Cannot determine TDIM'
           
            nhisto = 1
            nchans = 1
c
            if (type(5:6) .eq. 'sa') then
               do 30 k = 2, tfields
                  if (tdim(k) .ne. '') then
                     if (tdim(k)(1:1) .eq. '(') then
                        temp_string = tdim(k)(2:len_used(tdim(k))-1)
                        ipos = istrpos(temp_string,',')
                        nchans = 1
                        if (ipos .ne. 0 ) then                       
c
c KEYWORD TDIMn CONTAINS A "," SO THIS MUST BE A MULTI-DIMENSIONAL ARRAY
c
                           read(temp_string,*) nhisto, nchans
                        else
c
c ONLY ONE VALUE IN KEYWORD
c
                           read(temp_string,*) nhisto
                        endif
                        if (datamode(1:9) .eq. 'Standard2') then
c
c This is really a kludge! Should read NCPIXM keyword to find out
c whether we are dailing with CHANNELS or TIME.
c
                           ntemp  = nhisto
                           nhisto = nchans
                           nchans = ntemp
c
c ONLY for nhisto = 1 (ALWAYS in standard2) the indexing will go right
c
                           if (nhisto .ne. 1) then
                              stop 'NHISTO not 1, wrong results'
                           endif
                        endif
                        goto 35
                     endif
                  endif
 30            continue
            endif
c     
 35         if (type(5:6) .eq. 'sa') then
               call local_ftgkns(unit,'TDDES',2,tfields,tddes,
     *              nfound,status)
               if (status .ne. 0) stop ' Cannot read TDDES field '
               if (type(1:3).eq.'pca') event_description = tddes(2)
               if (type(1:3).eq.'hxt') event_description = tddes(4)
c     
c     INTERPRET THE MEANING OF THE BITS IN THE EVENTS (BINNED DATA)
c     
               print *, 'Calling bindes2chan'
               call bindes2chan(event_description,ichans,ichane,nchans)
               print *, 'THERE MIGHT BE A PROBLEM HERE !!!'
            else
               if ((type(5:6).eq.'se').or.(type(5:6).eq.'gx')) then
                  call local_ftgkns(unit,'TEVTB',2,tfields
     *                 ,tevtb,nfound,status)
                  if (status .ne. 0) stop ' Cannot read TEVTB field '
                  if (type(1:3).eq.'pca') event_description = tevtb(2)
                  if (type(1:3).eq.'hxt') event_description = tevtb(4)
c     
c     INTERPRET THE MEANING OF THE BITS IN THE EVENTS (EVENT DATA)
c     
                  if (sw_debug) print *,infilename
                  if (sw_debug) print *,event_description
C                  print *, 'eventdescr=',event_description, 'ichans',ichans, 'ichane',ichane, 'nchans',nchans, 'ioffs',ioffset(i,j),'ibit', ibit_chan(i,j)
                  call des2chan(event_description,ichans,ichane,nchans,
     &                 ioffset(i,j),ibit_chan(i,j))
               endif
            endif
c     
            nhistos(i,j)   = nhisto
            nchannels(i,j) = nchans
c     
            do k = 1, nchans
               channels1(k,i,j) = ichans(k)
               channels2(k,i,j) = ichane(k)
            enddo
c     
c     Close the file
c     
 66         call ftclos(unit,status)
            call ftfiou(unit,status)
            if (status .ne. 0) stop 'Cannot close the FITS file'
c     
 20      continue
 10   continue
c     
      return
      end
