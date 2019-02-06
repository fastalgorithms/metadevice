
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging routines and the beginning 
c       of the complex GMRES routines
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine cgmres(ier,n,a,multa,par1,par2,y,eps,numit,
     1     x,niter,errs,ngmrec,w)
        implicit real *8 (a-h,o-z)
        dimension errs(1)
        complex *16 a(n,n),x(1),y(1),w(1)
c
        external multa
c
c     This subroutine solves a complex linear system Ax=y by means
c     of GMRES algorithm. This is a memory management routine for 
c     cgmres1 which performs the actual work.
c
c       Cyclical GMRES algorithm
c
c     Note #1: This routine will use a cyclical buffer to store 
c     cgmres vectors. This improves convergence if cgmres needs to be
c     restarted (e.g, due to memory limitation reasons)
c
c
c
c     Input parameters:
c
c     n - the dimensionality of the linear system
c     y - the right hand side
c     a - the matrix of the system (or whatever other parameter)
c     eps - the required accuracy
c     multa - the user-defined matrix-vector multiplication subroutine
c
c     the calling sequence for multa must be
c
c     multa(a,par1,par2,x,y,n)               (1)
c
c       in (1), a is a matrix the system or whatever other parameter,
c       par1, par2 are whatever other parameters, x is an input vector,
c       y is a product Ax, and n is the dimensionality of a, x, and y.
c
c     numit - the maximum number of iteration permitted
c     ngmrec - the maximum number of iteration 
c            which GMRES algorithm needs to be restarted
c
c     w - must be at least (ngmrec*2+4)*n complex *16 elements long
c
c     Output parameters:
c
c     ier - error return code
c        ier=0 normal execution of the subroutine
c        ier=4 means that the maximum number iterations numit
c           has been reached without achieving the required accuracy eps
c        ier=8 means that the errors failed to decrease before the maximum
c           number of iterations has been reached or the required accuracy
c           eps has benn reached. 
c
c     x - the solution of the system
c     niter - the number of iterations performed 
c     errs - the array of errors produced by the algorithm. 
c        errs(i)=||y-Ax_i||, where x_i is a solution obtained on the i-th
c        GMRES iteration.
c
        ie=1
        lie=n*ngmrec
c
        iae=ie+lie
        liae=n*ngmrec
c
        iz=iae+liae
        lz=n
c
        iw=iz+lz
        lw=n
c
        ixr=iw+lw
        lxr=n
c
        izr=ixr+lxr
        lzr=n
c
        call cgmres1(ier,n,a,multa,par1,par2,y,eps,numit,x,niter,errs,
     1     ngmrec,w(ie),w(iae),w(iz),w(iw),w(ixr),w(izr))

        return
        end
c
c
c
c
c
        subroutine cgmres1(ier,n,a,multa,par1,par2,y,eps,numit,
     1     x,niter,errs,ngmrec,e,ae,z,w,xr,zr)
        implicit real *8 (a-h,o-z)
c
        dimension errs(1)
        complex *16 a(n,n),x(1),y(1),e(n,1),ae(n,1),
     1       z(1),w(1),xr(1),zr(1),d
c
        external multa
c
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        ier=0
        niter=0
c
        rhsnorm=0
        do 1000 i=1,n
c
        rhsnorm=rhsnorm+abs(conjg(y(i))*y(i))
        ifrec=0
        x(i)=0
        z(i)=y(i)
c
c     set x=y for cgmrespp  
c
ccc        ifrec=1
ccc        x(i)=y(i)
 1000 continue
        rhsnorm=sqrt(rhsnorm/n)
c
c
        dold=1d+100
c
        m=0
        ifcycle=0
c
        do 6000 iter=1,numit
c     
c       ... restart GMRES
c
        if( m .eq. ngmrec ) then
           m=0
           ifcycle=1
        endif
c
        m=m+1
        niter=niter+1
c
cccc        call prinf('niter=*',niter,1)
        write(*,*) 'in cgmres, niter=', niter
cccc        call prinf('in cgmres, m=*',m,1)
c
        nsteps=10
c
ccc        if( m .eq. 1 ) ifrec=1
        if( mod(niter,nsteps) .eq. 1 ) ifrec=1
        if( ifrec .eq. 0 ) goto 1600
c
cccc        call prinf('in cgmres, recalculate the residual at iter=*'
cccc     1     ,iter,1)
c
c       ... this is the first iteration, calculate the residual z=y-Ax
c       also, recalculate the residual every 10 steps
c
        call multa(a,par1,par2,x,z,n)
c
        do 1200 i=1,n
        z(i)=y(i)-z(i)
 1200 continue
c
 1600 continue
c
c
c       ... compute the error
c
        d=0
        do 2000 i=1,n
        d=d+dconjg(z(i))*z(i)
 2000 continue
c
        d=sqrt(d)
        errs(niter)=dble(d)
c
ccc        if( abs(d) .lt. eps ) return
        if( abs(d) .lt. eps*errs(1) ) return
c
        if( abs(d) .gt. dold ) then
c
c       ... the errors stopped decreasing, abort
c
        niter=niter-1
        do 2200 i=1,n
        x(i)=xr(i)
 2200 continue 
c
        ier=8
        return
        endif
c
 2400 continue
c
        dold=d
c
c       ... compute the new direction w=Az
c
        call multa(a,par1,par2,z,w,n)
c
c       ... orthogonalize w to all preceeding ae
c
        do 3000 i=1,n
        zr(i)=z(i)
 3000 continue
c
        if( ifcycle .eq. 0 ) then
        do 3600 j=1,m-1
        d=0
        do 3200 i=1,n
        d=d+dconjg(ae(i,j))*w(i)
 3200 continue
        do 3400 i=1,n
        w(i)=w(i)-d*ae(i,j)
        zr(i)=zr(i)-d*e(i,j)
 3400 continue
 3600 continue
        endif
c       
        if( ifcycle .eq. 1 ) then
        do 3650 j=1,ngmrec
        d=0
        do 3250 i=1,n
        d=d+dconjg(ae(i,j))*w(i)
 3250 continue
        do 3450 i=1,n
        w(i)=w(i)-d*ae(i,j)
        zr(i)=zr(i)-d*e(i,j)
 3450 continue
 3650 continue
        endif
c       
c       
c       ... normalize the current direction w
c
        d=0
        do 3800 i=1,n
        d=d+dconjg(w(i))*w(i)
 3800 continue
        d=1/sqrt(d)
c
c     ... store e and ae 
c
        do 4000 i=1,n
        ae(i,m)=w(i)*d
        e(i,m)=zr(i)*d
 4000 continue
c
ccc        call prinf('m=*',m,1)
ccc        call prin2('d=*',d,1)
ccc        call prin2('e=*',e(1,m),2*n)
ccc        call prin2('Ae=*',Ae(1,m),2*n)
c
c       ... double orthogonalize the current ae to all preceeding ae
c
        if( ifcycle .eq. 0 ) then
        do 4030 j=1,m-1
        d=0
        do 4010 i=1,n
        d=d+dconjg(ae(i,j))*ae(i,m)
 4010 continue
        do 4020 i=1,n
        ae(i,m)=ae(i,m)-d*ae(i,j)
        e(i,m)=e(i,m)-d*e(i,j)
 4020 continue
 4030 continue
        endif
c
        if( ifcycle .eq. 1 ) then
        do 4035 j=1,ngmrec
        if( j.eq.m ) goto 4035
        d=0
        do 4015 i=1,n
        d=d+dconjg(ae(i,j))*ae(i,m)
 4015 continue
        do 4025 i=1,n
        ae(i,m)=ae(i,m)-d*ae(i,j)
        e(i,m)=e(i,m)-d*e(i,j)
 4025 continue
 4035 continue
        endif
c
        d=0
        do 4040 i=1,n
        d=d+dconjg(ae(i,m))*ae(i,m)
 4040 continue
        d=1/sqrt(d)
c
cccc        call prin2('d=*',d-1,1)
c
c     ... store e and ae 
c
        do 4050 i=1,n
        ae(i,m)=ae(i,m)*d
        e(i,m)=e(i,m)*d
 4050 continue
c
c
c       ... update the solution x and the residual z
c
        d=0
        do 4200 i=1,n
        d=d+dconjg(ae(i,m))*z(i)
 4200 continue
c
        do 4400 i=1,n
        z(i)=z(i)-d*ae(i,m)
 4400 continue
c
        do 4600 i=1,n
        xr(i)=x(i)
        x(i)=x(i)+d*e(i,m)
 4600 continue
c
c
c       ... update the error
c
        d=0
        do 4800 i=1,n
        d=d+dconjg(z(i))*z(i)
 4800 continue
c
        d=sqrt(d/n)
        errs(niter)=dble(d)
        write(*,*) 'in cgmres, rms error=', errs(niter)
        write(*,*) 'in cgmres, rel error=', errs(niter)/rhsnorm
c
ccc        if( abs(d) .lt. eps ) return
        if( abs(d) .lt. eps*errs(1) ) return
c
 6000 continue
c
c       ... the maximum number of iterations has been reached, abort
c
        ier=4
c
        return
        end
c
c
c
c
c
        subroutine cgmrespp(w,ngmrec,nvec,n,ainv,a,multb,par3,par4,aw)
        implicit real *8 (a-h,o-z)
        complex *16 w(n,ngmrec,2),ainv(n,n),a(n,n),cd,aw(1),alpha
c
        external multb
c
c       This subroutine finds the inverse (or an approximation of the
c       inverse) of the matrix a by using the rank-1 updates of the UNIT
c       matrix, as contructed via a preceding call to cgmres.
c
c       This subroutine assumes that we will start with the UNIT matrix,
c       which is to say that, the matrix a is assumed to be a low rank
c       preturbation of a unit matrix. In cgmres, the initial
c       approximation of the solution x MUST be set to the right hand
c       side y.
c
c         Input parameters:
c
c   a - the matrix to be inverted
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   multb - the user-supplied adjoint matrix-vector multiplication routine.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      description of the calling system of subroutine multb
c
c       the calling sequence of multb is as follows:
c
c          multb(a,par3,par4,x,y,n)
c
c             with the following parameters:
c      a - the matrix of the system (or whatever other parameter) (input)
c      par3,par4 - other parameters (input)
c      x - the vector to which conjugate transpose of a
c             is to be applied (input)
c      y -= b(x) = conjugate transpose of a(x) (output)
c      n - the dimensionality of a,x,y (input)
c
c      end of description of the calling system of subroutine multb
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   w - the rank-1 update vectors from cgmres routine
c   ngmrec - the maximum number of iteration 
c       before GMRES algorithm needs to be restarted
c   nvec - the number of rank-1 update vectors (must less or equal to niter-1,
c       where niter - the actual number iteration performed 
c       by the GMRES algorithm)
c   n - the dimensionality of the problem
c     
c   aw - work array. must be at least n complex *16 elements long
c
c         Output parameters:
c
c   ainv - the inverse of a
c     
c
        alpha=1
c        
        do 1400 i=1,n
        do 1200 j=1,n
        ainv(i,j)=0
 1200 continue
        ainv(i,i)=alpha
 1400 continue
c
c
        do 2600 k=1,nvec
c       
        call multb(a,par3,par4,w(1,k,2),aw,n)
c       
        do 2400 i=1,n
        do 2200 j=1,n
cccc        ainv(i,j)=ainv(i,j)+dconjg(w(i,k,1))*(w(j,k,2)-aw(j)*alpha) 
        ainv(i,j)=ainv(i,j)+(w(i,k,1))*dconjg(w(j,k,2)-aw(j)*alpha) 
 2200 continue
 2400 continue
c       
 2600 continue
c       
        return
        end
c
c
c
c
c
        subroutine cgmrespo(w,ngmrec,nvec,n,ainv,a)
        implicit real *8 (a-h,o-z)
        complex *16 w(n,ngmrec,2),ainv(n,n),a(n,n),cd,alpha
c
c       This subroutine finds the inverse (or an approximation of the
c       inverse) of the matrix a by using the rank-1 updates of the UNIT
c       matrix, as contructed via a preceding call to cgmres.
c
c       This subroutine assumes that we will start with the UNIT matrix,
c       which is to say that, the matrix a is assumed to be a low rank
c       preturbation of a unit matrix. In cgmres, the initial
c       approximation of the solution x MUST be set to the right hand
c       side y.
c
c         Input parameters:
c
c   a - the matrix to be inverted
c   w - the rank-1 update vectors from cgmres routine
c   ngmrec - the maximum number of iteration 
c       before GMRES algorithm needs to be restarted
c   nvec - the number of rank-1 update vectors (must less or equal to niter-1,
c       where niter - the actual number iteration performed 
c       by the GMRES algorithm)
c   n - the dimensionality of the problem
c     
c         Output parameters:
c
c   ainv - the inverse of a
c     
c
        alpha=1
c        
        do 1400 i=1,n
        do 1200 j=1,n
        ainv(i,j)=0
 1200 continue
        ainv(i,i)=alpha
 1400 continue
c
c
        do 2600 k=1,nvec
c       
        do 2400 i=1,n
        do 2200 j=1,n
        ainv(i,j)=ainv(i,j)+(w(i,k,1)-w(i,k,2)*alpha)*dconjg(w(j,k,2))
 2200 continue
 2400 continue
c       
 2600 continue
c       
        return
        end
c
c
c
c
c
        subroutine cgmrespv(w,ngmrec,nvec,n,a,multb,par3,par4,aw)
        implicit real *8 (a-h,o-z)
        complex *16 w(n,ngmrec,2),a(n,n),cd,aw(1),alpha
c
        external multb
c
c       This subroutine finds the rank-1 update vectors of the inverse
c       (or an approximation of the inverse) of the matrix a by using
c       the rank-1 updates of the UNIT matrix, as contructed via a
c       preceding call to cgmres. The inverse matrix itself is not
c       formed.
c
c       This subroutine assumes that we will start with the UNIT matrix,
c       which is to say that, the matrix a is assumed to be a low rank
c       preturbation of a unit matrix. In cgmres, the initial
c       approximation of the solution x MUST be set to the right hand
c       side y.
c
c         Input parameters:
c
c   a - the matrix to be inverted
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   multb - the user-supplied adjoint matrix-vector multiplication routine.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      description of the calling system of subroutine multb
c
c       the calling sequence of multb is as follows:
c
c          multb(a,par3,par4,x,y,n)
c
c             with the following parameters:
c      a - the matrix of the system (or whatever other parameter) (input)
c      par3,par4 - other parameters (input)
c      x - the vector to which conjugate transpose of a
c             is to be applied (input)
c      y -= b(x) = conjugate transpose of a(x) (output)
c      n - the dimensionality of a,x,y (input)
c
c      end of description of the calling system of subroutine multb
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   w - the rank-1 update vectors from cgmres routine
c   ngmrec - the maximum number of iteration 
c       before GMRES algorithm needs to be restarted
c   nvec - the number of rank-1 update vectors (must less or equal to niter-1,
c       where niter - the actual number iteration performed 
c       by the GMRES algorithm)
c   n - the dimensionality of the problem
c     
c   aw - work array. must be at least n complex *16 elements long
c
c         Output parameters:
c
c   w - the rank-1 update vectors of the inverse matrix
c
c     
c
        alpha=1        
c
        do 2600 k=1,nvec
c       
        call multb(a,par3,par4,w(1,k,2),aw,n)
c
        do 2200 j=1,n
        w(j,k,2)=w(j,k,2)-aw(j)*alpha
 2200 continue
c       
 2600 continue
c       
        return
        end
c
c
c
c
c
        subroutine cgmrespz(w,ngmrec,nvec,n)
        implicit real *8 (a-h,o-z)
        complex *16 w(n,ngmrec,2),cd,alpha
c
c       This subroutine finds the rank-1 update vectors of the inverse
c       (or an approximation of the inverse) of the matrix a by using
c       the rank-1 updates of the UNIT matrix, as contructed via a
c       preceding call to cgmres. The inverse matrix itself is not
c       formed.
c
c       This subroutine assumes that we will start with the UNIT matrix,
c       which is to say that, the matrix a is assumed to be a low rank
c       preturbation of a unit matrix. In cgmres, the initial
c       approximation of the solution x MUST be set to the right hand
c       side y.
c
c         Input parameters:
c
c   w - the rank-1 update vectors from cgmres routine
c   ngmrec - the maximum number of iteration 
c       before GMRES algorithm needs to be restarted
c   nvec - the number of rank-1 update vectors (must less or equal to niter-1,
c       where niter - the actual number iteration performed 
c       by the GMRES algorithm)
c   n - the dimensionality of the problem
c     
c         Output parameters:
c
c   w - the rank-1 update vectors of the inverse matrix
c
c     
c       
        alpha=1
c
        do 2600 k=1,nvec
c       
        do 2200 j=1,n
        w(j,k,1)=w(j,k,1)-w(j,k,2)*alpha
 2200 continue
c       
 2600 continue
c       
        return
        end
c
c
c
c
c
        subroutine cgmrespr(w,ngmrec,nvec,x,y,n)
        implicit real *8 (a-h,o-z)
        complex *16 w(n,ngmrec,2),x(1),y(1),d,alpha
c
        alpha=1
c       
        do 1200 i=1,n
        y(i)=x(i)*alpha
 1200 continue
c
        do 2800 j=1,nvec
        d=0
        do 2200 i=1,n
cccc        d=d+w(i,j,2)*x(i)
        d=d+dconjg(w(i,j,2))*x(i)
 2200 continue        
        do 2400 i=1,n
cccc        y(i)=y(i)+dconjg(w(i,j,1))*d
        y(i)=y(i)+w(i,j,1)*d
 2400 continue
 2800 continue        
c
        return
        end
c
c
c
c
c
        subroutine cgmrespr1(w,ngmrec,nvec,x,y,n)
        implicit real *8 (a-h,o-z)
        complex *16 w(n,ngmrec,2),x(1),y(1),d,alpha
c
        alpha=1
c
        do 1200 i=1,n
        y(i)=x(i)*alpha
 1200 continue
c
        do 2800 j=1,nvec
        d=0
        do 2200 i=1,n
cccc        d=d+w(i,j,2)*x(i)
        d=d+dconjg(w(i,j,2))*x(i)
 2200 continue        
        do 2400 i=1,n
cccc        y(i)=y(i)+dconjg(w(i,j,1))*d
        y(i)=y(i)+w(i,j,1)*d
 2400 continue
 2800 continue        
c
        return
        end
c
c
c
c
c
        subroutine cgmrespr2(w,ngmrec,nvec,x,y,n)
        implicit real *8 (a-h,o-z)
        complex *16 w(n,ngmrec,2),x(1),y(1),d,alpha
c
        alpha=1
c
        do 1200 i=1,n
        y(i)=x(i)*alpha
 1200 continue
c
        do 2800 j=1,nvec
        d=0
        do 2200 i=1,n
cccc        d=d+w(i,j,1)*x(i)
        d=d+dconjg(w(i,j,1))*x(i)
 2200 continue        
        do 2400 i=1,n
cccc        y(i)=y(i)+dconjg(w(i,j,2))*d
        y(i)=y(i)+w(i,j,2)*d
 2400 continue
 2800 continue        
c
        return
        end
c
c
c
c
c
