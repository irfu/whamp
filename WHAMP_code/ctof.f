*
************************************************************************
* 
       SUBROUTINE CTOF_CONVERT_STRING(CSTRING,FSTRING)
* 
       PARAMETER(MAX_LENGTH = 80)
       CHARACTER C*1 
       CHARACTER*(*) CSTRING,FSTRING
*
*       CALL CEOS(C)
       DO 10 I=1,MAX_LENGTH
         IF(CSTRING(I:I).EQ.C)THEN
           GOTO 15
         ENDIF
   10  CONTINUE
       RETURN 

   15  CONTINUE

       FSTRING=CSTRING(1:I-1)

       RETURN
       END

