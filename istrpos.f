
C                          $Author: tim $
C                          $Date: 1996/06/21 13:55:25 $
C                          $Id: istrpos.f,v 1.4 1996/06/21 13:55:25 tim Exp $
C                          $Source: /home/xray/sources/subroutines/istrpos.f,v $

	integer function istrpos(string,search_string)

c*****  This function returns the position of the first occurence
c*****  of first character of the search_string in string as seen from the back.
c*****  If the search_string is not found in string istrpos is returned
c*****  as zero.
c*****  
c*****  History: 
c*****          created 23-01-1992 : Tim Oosterbroek


	character*(*) string
	character*(*)   search_string
	length_search=len(search_string)
	length=len(string) -length_search		

        if (length .eq. 0) then
           istrpos = 0
           if (search_string .eq. string) istrpos = 1
           return
        endif
        do while (.true.)
           if (string(length:length+length_search-1).eq.
     *       search_string) goto 36
           length = length - 1
           if (length .le. 0) goto 36
        enddo
 36	istrpos=length
	return
	end
