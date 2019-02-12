c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       this is the end of the debugging code and the beginning of the
c       actual Bi-conjugate(stab) gradient algorithm routines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine cbicgstab(ier,n,a,multa,par1,par2,
     $     y,eps,numit,x,niter,errs,w)
        implicit real *8 (a-h,o-z)
        complex *16 a(1),y(1),x(1),w(1)
        dimension errs(1)
c
        external multa
c
c
c     This subroutine solves a complex linear system Ax=y by means
c     of BiCG(stab) algorithm. This is a memory management routine for 
c     cbicgstab1 which performs the actual work.
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
c
c     w - must be at least 11*n complex *16 elements long
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
c        iteration.
c
c
        iaxk=1
        iek=iaxk+n
        iaek=iek+n
        iekm1=iaek+n
        iaekm1=iekm1+n
c
        ibek=iaekm1+n
        ibaek=ibek+n
        ibekm1=ibaek+n
        ibaekm1=ibekm1+n
        iesk=ibaekm1+n
        ietk=iesk+n
c
        call cbicgstab1(ier,n,a,y,numit,x,w(iaxk),
     1     w(iek),w(iaek),w(iekm1),w(iaekm1), w(iesk),w(ietk),
     2     w(ibek),w(ibaek),w(ibekm1),w(ibaekm1),errs,nrec,
     2     multa,par1,par2,
     4     eps,niter)
c
        return
        end
c
c
c
c
c
        subroutine cbicgstab1(ier,n,a,y,numit,xk,axk,
     1     ek,aek,ekm1,aekm1,esk,etk,
     2     bek,baek,bekm1,baekm1,errs,nrec,
     2     multa,par1,par2,
     4     eps,niter)
        implicit real *8 (a-h,o-z)
        complex *16 a(1),y(1),xk(1),axk(1),
     1     ek(1),aek(1),ekm1(1),aekm1(1),
     2     bek(1),baek(1),bekm1(1),baekm1(1),
     3     er1,er2,cd1,cd2,alpha,beta,gamma,omega,esk(1),etk(1)
        dimension errs(1)
c
        external multa
c
c       initialize the Bi-conjugate-STAB gradient iterations
c
        done=1
        ifrec=0
        ier=0
        dold=1.0d60
cccc        call prin2('y=*',y,n*2)
c        
        do 1400 i=1,n
ccc        xk(i)=0
ccc        axk(i)=0
        xk(i)=y(i)
ccc        xk(i)=i
 1400 continue
c
        if( n .eq. 1 ) xk(1)=1
c
        niter=1
        write(*,*) 'in cbicgstab, niter=',niter
        call multa(a,par1,par2,xk,axk,n)
c
        if( n .eq. 1 ) xk(1)=y(1)/axk(1)
        if( n .eq. 1 ) axk(1)=y(1)
c
c       ... find the first direction
        er=0
        er1=0
        do 1600 i=1,n
        ek(i)=y(i)-axk(i)
        ekm1(i)=ek(i)
c
        bek(i)=ek(i)
        bekm1(i)=bek(i)
c
        er=er+dconjg(ek(i))*ek(i)
        er1=er1+dconjg(bek(i))*ek(i)
 1600 continue
c
        er=dsqrt(abs(er)/n)
c
cccc        er=dznrm2(n,ek,1)
c
        errs(1)=er
c
        write(*,*) 'in cbicgstab, error=',er
ccc        if( er .lt. eps ) return
        if( er .eq. 0 ) return
c
c
c      conduct Bi-conjugate-STAB gradient iterations
c
        do 4000 k=2,numit
        niter=k
ccc        call prinf('k=*',k,1)
        write(*,*) 'in cbicgstab, niter=',niter
c
        call multa(a,par1,par2,ekm1,aekm1,n)
c
        cd1=0
        cd2=0
        do 2200 i=1,n
        cd1=cd1+dconjg(bek(i))*ek(i)
        cd2=cd2+dconjg(bek(i))*aekm1(i)
 2200 continue
c
ccc        cd1=dznrm2(n,ek,1)**2
c
        alpha=cd1/cd2
c
ccc        call prin2('cd1=*',cd1,2)
ccc        call prin2('cd2=*',cd2,2)
ccc        call prin2('alpha=*',gamma,2)
c
        do i=1,n
        esk(i)=ek(i)-alpha*aekm1(i)
        enddo
        call multa(a,par1,par2,esk,etk,n)
c
        cd1=0
        cd2=0
        do i=1,n
        cd1=cd1+dconjg(etk(i))*esk(i)
        cd2=cd2+dconjg(etk(i))*etk(i)
        enddo
c
        omega=cd1/cd2
c
ccc        call prin2('cd1=*',cd1,2)
ccc        call prin2('cd2=*',cd2,2)
ccc        call prin2('omega=*',omega,2)
c
c       ... update the solution
c
        do 2400 i=1,n
        xk(i)=xk(i)+alpha*ekm1(i)+omega*esk(i)
 2400 continue
c
c       ... update the residual
c
        er=0
        er2=0
        do 2600 i=1,n
        ek(i)=esk(i)-omega*etk(i)
        er=er+dconjg(ek(i))*ek(i)
        er2=er2+dconjg(bek(i))*ek(i)
 2600 continue
c
        er=er/n
        er2=er2/n
c
        er=dsqrt(abs(er))
c
        errs(niter)=er
c
        write(*,*) 'in cbicgstab, error=',er
c
ccc        if( er .lt. eps ) return
        if( er .lt. eps*errs(1) ) return
        if( er .eq. 0 ) return
c
c       ... update the search direction
c
        beta=(er2/er1) * (alpha/omega)
        er1=er2
c
cccc        call prin2('beta=*',beta,2)
c
        do 3200 i=1,n
        ekm1(i)=ek(i)+beta*(ekm1(i)-omega*aekm1(i))
 3200 continue
c
 4000 continue
c        the maximum permitted number of iterations has been
c        performed without the desired accuracy being accuracy
c        achieved. inform the user and exit
c
        ier=4
        return
        end
c
