c**
c** Subroutine local_ftgkns
c**

C
C This routine should replace all the calls to the FITSIO routine ftgkns.
C This is because in the new release of FITSIO ftgkns cannot longer read
C LONG keywords, so we had to change things a bit to still be able to get
C the LONG keywords (for instance, the channel description in a FITS file
C is usually longer than 80 characters, and therefore is one of those).
C This routine should make things transparent to the user and the
C programmer, as long as he/she uses this one instead of the old one.
C One of the advantages of the new FITSIO is that you can work with
C compressed (either .Z or .gz) files.
C
C                             Mariano (7/8/98).
C-----------------
C
C The calling sequence was kept identical to that of the
C FITSIO routine ftgkns.
C From the FITSIO Manual:
C
C Get a sequence of numbered keyword values 
C FTGKN[EDJLS](unit,keyroot,startno,max_keys, > keyvals,nfound,status)
C
C This routine (which only replaces FTGKNS, i.e., for strings) does
C the same, but it can deal with long keywords (using the continuation
C & and the keyword CONTINUE) while FTGKNS (in the new FITSIO releases)
C cannot.
C
C <<< Input
C unit:    Unit number of the FITS file.
C keyroot: Root part of the name of the desired keyword.
C seq_no:  First keyword to get (usually 1).
C tfields: Maximum number of keywords to look for.
C >>> Output
C data:    Contents of the keywords found in the file.
C nfound:  Number of keywords found matching the given root name.
C status:  Status (0 = Ok, > 0 error).
C
C -- Calls FTGKNS, which tells how many keywords are there that
C    match the given root name, although it cannot longer provide
C    the values, if they are LONG ... yes, I know I said this
C    several times. It should be clear now :^).
C -- Calls FTKEYN which appends numbers to the root name of
C    the given keyword (e.g., the user provides TDIM, and
C    this routine builds up TDIM1, TDIM2, etc...).
C -- Calls FTGKYS which actually gets the value for the given
C    keyword. This is the routine that can read LONG keywords.
C -- Calls PRINTERROR to print error messages (if status > 0).
C
C
      subroutine local_ftgkns(unit,keyroot,seq_no,
     *tfields,data,nfound,status)
c


      parameter(max_col=256)
      character*(*) data(max_col),keyroot
      character*8   keyword
      character     comment*80
      integer       unit,seq_no,nfound,status
      
      nfound=0
      status=0
      
C-- Use FTGKNS to get the number of keywords that match with
C   the (user) given keyroot.
      call ftgkns(unit,keyroot,seq_no,tfields,data,nfound,status)
      call printerror(status)
C-- And now loop through the whole header looking for those
C   keywords that matched.
C   Start on the one given by the user (seq_no), and loop until
C   we looked at all of them (nfound).
    
      do jj=seq_no,seq_no+nfound-1
         status=0
C-- Use the keyroot to build up the keyword name as follows:
C   keyword = keyroot // jj, e.g., TDIM + 1 -> TDIM1.
C MM -- Here there seems to be a problem I do not see in Amsterdam.
C       The first appearence of a keyword seems not to be numbered.
C       So, e.g., TDIM goes as TDIM, TDIM2, TDIM3, etc. (no TDIM1).
C       This causes a "Keyword not found" error which apparently is
C       harmless, but nevertheless annoying. The problem seems to be
C       with TEVTB, which starts from 2. Other keywords may have the
C       same problem.
C
C       The reason they don't have this in Amsterdam seems to be that
C       they use FTGKNS instead of this routine, but most probably they
C       are still linking the old FITSIO (given that in the new one 
C       FTGKNS no longer supports long strings).
C
        
         call FTKEYN(keyroot,jj,keyword,status)
         call printerror(status)
         status=0
C-- Get the keyword value (even if it is LONG).
         call FTGKYS(unit,keyword,data(jj),comment,status)
         call printerror(status)
      enddo
      return
      end
      

