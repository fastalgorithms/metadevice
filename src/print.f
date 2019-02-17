c
c
c
        subroutine prin3(msg,a,n)
        implicit double precision (a-h,o-z)
        dimension a(1)
        character*1 msg(1),AST
        data AST / '*' /
        data iw1  / 13 /, iw2 / 6 /
c
c       Output an array of double words as a formatted column; show the 
c       most significant 15 digits of each value.
c
 0100 format(1X,80A)
 0200 format(4(2X,E22.15))
c
        len = 0
        do 1000 j=1,10 000
        if (msg(j) .eq. AST) goto 1100
        len=len+1
 1000 continue
 1100 continue
        write (iw1,0100) (msg(j),j=1,len)
        write (iw1,0200) (a(j),j=1,n)
c
        write (iw2,0100) (msg(j),j=1,len)
        write (iw2,0200) (a(j),j=1,n)

c
        end

c
c
c
c


        subroutine elapsed(t)
        implicit double precision (a-h,o-z)
        integer*8 i,irate
c
c       Return the elapsed time from the beginning of program execution
c       in seconds as a double word.  
c
c       NOTE: the resolution of the timer which is queried is operating
c       system and compiler dependent.  This routines seems to produce
c       reasonable results on a wide range of modern systems, however.
c
        call system_clock(i,irate)
c
        dd = i
        dd = dd/irate
        t = dd
c
c
        return
        end
