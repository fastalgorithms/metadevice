
C
C
C
C
        SUBROUTINE XPRINI(IP1,IQ1)
        CHARACTER *1 MES(1), AA(1)
         save
        REAL *4 A(1)
        REAL *8 A2(1)
        COMPLEX *8 AC(1)
        COMPLEX *16 AC2(1)
ccc        INTEGER *4 IA(1)
        INTEGER IA(1)
        INTEGER *4 IA1(1)
        INTEGER *2 IA2(1)
        IP=IP1
        IQ=IQ1
        RETURN

C
C
C
 1200   FORMAT(1(2X,E24.8))
 1300   FORMAT(1(3X,'(',E22.16,',',3X,E22.16,')'))
 1400   FORMAT(1(2X,E24.16))
 1600   FORMAT(1(1X,I15))
C
C
C
        ENTRY XPRIN(MES,A,N)
        ENTRY XPRIN1(MES,A,N)
        CALL  XMESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1200)(A(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1200)(A(J),J=1,N)
        RETURN
C
C
C
        ENTRY XPRIN2(MES,A2,N)
        CALL XMESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1400)(A2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1400)(A2(J),J=1,N)
        RETURN
C
C
C
C
        ENTRY XPRINC(MES,AC,N)
        ENTRY XPRINC1(MES,AC,N)
        CALL  XMESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1200)(AC(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1200)(AC(J),J=1,N)
        RETURN

        ENTRY XPRINC2(MES,AC2,N)
        CALL  XMESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1400)(AC2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1400)(AC2(J),J=1,N)
        RETURN
c
        ENTRY XPRINC2B(MES,AC2,N)
        CALL  XMESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1300)(AC2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1300)(AC2(J),J=1,N)
        RETURN

C
        ENTRY XPRINF(MES,IA,N)
        CALL XMESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA(J),J=1,N)
        RETURN
C
        ENTRY XPRINF1(MES,IA1,N)
        CALL XMESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA1(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA1(J),J=1,N)
        RETURN
C
        ENTRY XPRINF2(MES,IA2,N)
        CALL XMESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA2(J),J=1,N)
        RETURN
C
C
C
C
        ENTRY XPRINA(MES,AA,N)
        CALL XMESSPR(MES,IP,IQ)
 2000 FORMAT(1X,80A1)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,2000)(AA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,2000)(AA(J),J=1,N)
        RETURN
        END
c
c
c
        SUBROUTINE XMESSPR(MES,IP,IQ)
        CHARACTER *1 MES(1),AST
        DATA AST/'*'/
C
C         DETERMINE THE LENGTH OF THE MESSAGE
C
        I=0
        DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
         IF ( (I1.NE.0) .AND. (IP.NE.0) )
     1     WRITE(IP,1800) (MES(I),I=1,I1)
         IF ( (I1.NE.0) .AND. (IQ.NE.0) )
     1     WRITE(IQ,1800) (MES(I),I=1,I1)
 1800 FORMAT(1X,80A1)
        RETURN
        END
C
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        SUBROUTINE XREADI(IP1)
        CHARACTER *1 MES(1), AA(1)
         save
        REAL *4 A(1)
        REAL *8 A2(1)
        COMPLEX *8 AC(1)
        COMPLEX *16 AC2(1)
ccc        INTEGER *4 IA(1)
        INTEGER IA(1)
        INTEGER *4 IA1(1)
        INTEGER *2 IA2(1)
        IP=IP1
        RETURN
C
C
C
 1200   FORMAT(1(2X,E24.8))
 1300   FORMAT(1(3X,'(',E22.16,',',3X,E22.16,')'))
 1400   FORMAT(1(2X,E24.16))
C
C
C
        ENTRY XREAD(MES,A,N)
        ENTRY XREAD1(MES,A,N)
        CALL XMESSRE(IP,MES,LEN)
        IF(IP.NE.0 .AND. N.NE.0) READ(IP,*)(A(J),J=1,N)
        RETURN
C
C
C
        ENTRY XREAD2(MES,A2,N)
        CALL XMESSRE(IP,MES,LEN)
        IF(IP.NE.0 .AND. N.NE.0) READ(IP,*)(A2(J),J=1,N)
        RETURN
C
C
        ENTRY XREADC(MES,AC,N)
        ENTRY XREADC1(MES,AC,N)
        CALL XMESSRE(IP,MES,LEN)
        IF(IP.NE.0 .AND. N.NE.0) READ(IP,*)(AC(J),J=1,N)
        RETURN
C
C
        ENTRY XREADC2(MES,AC2,N)
        CALL XMESSRE(IP,MES,LEN)
        IF(IP.NE.0 .AND. N.NE.0) READ(IP,*)(AC2(J),J=1,N)
        RETURN
c
c
        ENTRY XREADC2B(MES,AC2,N)
        CALL XMESSRE(IP,MES,LEN)
        IF(IP.NE.0 .AND. N.NE.0) READ(IP,*)(AC2(J),J=1,N)
        RETURN

c
        ENTRY XREADF(MES,IA,N)
        CALL XMESSRE(IP,MES,LEN)
        IF(IP.NE.0 .AND. N.NE.0) READ(IP,*)(IA(J),J=1,N)
        RETURN
c
        ENTRY XREADF1(MES,IA1,N)
        CALL XMESSRE(IP,MES,LEN)
        IF(IP.NE.0 .AND. N.NE.0) READ(IP,*)(IA1(J),J=1,N)
 1600 FORMAT(10(1X,I15))
C
        ENTRY XREADF2(MES,IA2,N)
        CALL XMESSRE(IP,MES,LEN)
        IF(IP.NE.0 .AND. N.NE.0) READ(IP,*)(IA2(J),J=1,N)
        RETURN
C
C
C
C
        ENTRY XREADA(MES,AA,N)
 2000 FORMAT(1X,80A1)
        CALL XMESSRE(IP,MES,LEN)
        IF(IP.NE.0 .AND. N.NE.0) READ(IP,*)(AA(J),J=1,N)
        RETURN
        END
c
c
c
C
c
        SUBROUTINE XMESSRE(IP,MES,LEN)
        CHARACTER *1 MES(1),AST
        CHARACTER *1 MES1(1)
        DATA AST/'*'/
C
C         DETERMINE THE LENGTH OF THE MESSAGE
C
        READ(IP,1800) MES1
C
        I=0
        DO 1400 I=1,10000
        IF(MES1(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE

        N=I1

 1800 FORMAT(1X,80A1)
        RETURN
        END
C
C
