
C                          $Author: tim $
C                          $Date: 1996/06/27 12:03:03 $
C                          $Id: len_used.f,v 1.2 1996/06/27 12:03:03 tim Exp $
C                          $Source: /home/xray/sources/subroutines/len_used.f,v $

  	integer function len_used(string)
C ** Gives the number of used characters (excluding the right spaces padding)
C ** in a string. If the string is just spaces, then a value of 1 is returned.
  	character string*(*)
  	do i = len(string),1,-1
  	  if (string(i:i).ne.' ') goto 1
  	enddo
        len_used = 1
        return

1  	len_used = i
  	return
  	end
