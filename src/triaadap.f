c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the 
c        start of the actual quadrature routines.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine triaadap(ier,vert1,vert2,vert3,fun,
     1      par1,par2,m,eps,rint,maxrec,numfunev,w)
        implicit real *8 (a-h,o-z)
        dimension w(1),par1(1),par2(1),vert1(2),vert2(2),vert3(2)
c
c       this subroutine uses the adaptive tensor product gaussian
c       quadrature to integrate a user-supplied function on a 
c       triangle in the plane (also user-supplied)
c
c                       input parameters:
c
c  vert1,vert2,vert3 - the vertices in the plane of the triangle
c       over which the function is to be integrated
c  fun - the user-supplied function to be integrated. the calling
c       sequence of fun must be 
c
c        call fun(x,y,par1,par2,rint).                            (1)
c
c        in (1), (x,y) is a point in the plane where 
c        the function is to be evaluated, and par1, par2
c        are two parameters to be used by fun; they can be 
c        variables or arrays, real or integer, as desired. 
c        rint is assumed to be real *8.
c  par1, par2 - parameters to be used by the user-supplied 
c       subroutine fun (see above)
c  m - the order of the quadrature to me used on each subinterval
c  eps - the accuracy (absolute) to which the integral will be 
c       evaluated
c
c                       output parameters:
c
c  ier - error return code. 
c          ier=0 means normal conclusion
c          ier=8 means that at some point, the depth of recursion
c                reached 200. this is a fatal error.
c          ier=16 means that the total number of subtriangles in the
c                adaptive subdivision of [a,b] turned out to be greater 
c                than nnmax*4.  this is a fatal error.
c                
c  rint - the integral as evaluated
c  maxrec - the maximum depth to which the recursion went at its 
c         deepest point. can not be greater than 200, since at that
c         point ier is set to 8 and the execution of the subroutine
c         terminated.
c  numfunev - the total number of function evaluations (calls to the
c         user-supplied subroutine fun) that has occured in the 
c         calculation
c         
c                         work arrays:
c
c  w - must be at least 3*m**2+4*m+1500 real *8 elements long
c
c        . . . integrate the user-supplied function using the 
c              adaptive gaussian quadratures
c
        nnmax=100000
        maxdepth=200
        iiiw=21
c
c        allocate memory for the subroutine trianrec
c
        istack=1
        lstack=1207
c
        iw=istack+lstack
        lw=3*m**2+50
c
        ivals=iw+lw
        lvals=207
c
        iww=ivals+lvals
        lww=4*m+10
c
c       . . . integrate
c
        call trianrec(ier,w(istack),vert1,vert2,vert3,fun,
     1      par1,par2,w(iw),m,w(ivals),nnmax,eps,
     2      rint,maxdepth,maxrec,numint,iiiw,numfunev,w(iww))
c
        return
        end
c
c
c
c
c
        subroutine trianrec(ier,stack,z1,z2,z3,fun,
     1      par1,par2,w,m,vals,nnmax,eps,
     2      rint,maxdepth,maxrec,numint,iw,numfunev,ww)
        implicit real *8 (a-h,o-z)
        dimension stack(2,3,1),w(1),vals(1),par1(1),par2(1)
        dimension z1(2),z2(2),z3(2),zout(2,3,10),ww(1)
c
c       this subroutine uses the adaptive tensor product gaussian
c       quadrature to integrate a user-supplied function on a 
c       triangle in the plane (also user-supplied)
c
c                       input parameters:
c
c  z1,z2,z3 - the vertices in the plane of the triangle
c       over which the function is to be integrated
c  fun - the user-supplied function to be integrated. the calling
c       sequence of fun must be 
c
c        call fun(x,y,par1,par2,rint).                            (1)
c
c        in (1), (x,y) is a point in the plane where 
c        the function is to be evaluated, and par1, par2
c        are two parameters to be used by fun; they can be 
c        variables or arrays, real or integer, as desired. 
c        rint is assumed to be real *8.
c  par1, par2 - parameters to be used by the user-supplied 
c       subroutine fun (see above)
c  m - the order of the quadrature to be used on each subinterval
c  nnmax - maximum permitted number of subdivisions
c  eps - the accuracy (absolute) to which the integral will be 
c       evaluated
c  maxrec - maximum permitted depth of recursion
c  iw - the FORTRAN unit number on which the subroutine will
c       plot the subdivision of the triangle used for the integration.
c       setting ls .leq. 0 suppresses the plotting
c
c                       output parameters:
c
c  ier - error return code. 
c          ier=0 means normal conclusion
c          ier=8 means that at some point, the depth of recursion
c                reached 200. this is a fatal error.
c          ier=16 means that the total number of subtriangles in the
c                adaptive subdivision of the user-supplied one turned 
c                out to be greater than nnmax*4.  this is a fatal error.
c  rint - the integral as evaluated
c  maxrec - the maximum depth to which the recursion went at its 
c         deepest point. can not be greater than 200, since at that
c         point ier is set to 8 and the execution of the subroutine
c         terminated.
c  numint - the total number of triangles in the subdivision divided by 
c         four. can not be greater than nnmax,  since at that
c         point ier is set to 16 and the execution of the subroutine
c         terminated.
c  numfunev - the total number of calls to the subroutine fun that
c         have been performed ("the number of function evaluations")
c
c                         work arrays:
c
c  stack - must be at least 1206 real *8 elements long
c  w - must be at least 3*m**2+50 real *8 elements long
c  vals - must be at least 200 real *8 elements long
c  ww - must be at leats m*4+10 real *8 elements long
c
c       . . . start the recursion
c
        numfunev=0
c
        call triaarrm(z1,stack(1,1,1),2)
        call triaarrm(z2,stack(1,2,1),2)
        call triaarrm(z3,stack(1,3,1),2)
c
        numfunev=numfunev+m**2
c
        call triaint1(m,z1,z2,z3,fun,
     1      par1,par2,vals(1),w,ww)
c
c       recursively integrate the thing
c
        j=1
        rint=0
        ier=0
        maxrec=0
        do 3000 i=1,nnmax
ccc        call prinf('i=*',i,1)
        numint=i
        if(j .gt. maxrec) maxrec=j
cccc        call prinf('j=*',j,1)
c
c       subdivide the current triangle
c
        call triadiv(stack(1,1,j),stack(1,2,j),stack(1,3,j),zout)
c
        numfunev=numfunev+m**2
        call triaint1(m,zout(1,1,1),zout(1,2,1),zout(1,3,1),fun,
     1      par1,par2,value1,w,ww)
c
        numfunev=numfunev+m**2
        call triaint1(m,zout(1,1,2),zout(1,2,2),zout(1,3,2),fun,
     1      par1,par2,value2,w,ww)
c
        numfunev=numfunev+m**2
        call triaint1(m,zout(1,1,3),zout(1,2,3),zout(1,3,3),fun,
     1      par1,par2,value3,w,ww)
c
        numfunev=numfunev+m**2
        call triaint1(m,zout(1,1,4),zout(1,2,4),zout(1,3,4),fun,
     1      par1,par2,value4,w,ww)
c
        dd=dabs(value1+value2+value3+value4-vals(j))
c
cccc         call prin2('in trianrec, dd=*',dd,1)
c
        ifdone=0
        if(dd .le. eps) ifdone=1
c
c       if the function on this subinterval has been 
c       integrated with sufficient accuracy - add the 
c       value to that of the global integral and move up
c       in the stack
c
        if(ifdone  .eq. 0) goto 2000
c
        rint=rint+value1+value2+value3+value4
        j=j-1
c
c        if the whole thing has been integrated - return
c
        if(j .eq. 0) return
        goto 3000
 2000 continue
c        
c       if the function on this subinterval has not been 
c       integrated with sufficient accuracy - move 
c       down the stack
c
        call triaarrm(zout(1,1,1),stack(1,1,j),6)
        call triaarrm(zout(1,1,2),stack(1,1,j+1),6)
        call triaarrm(zout(1,1,3),stack(1,1,j+2),6)
        call triaarrm(zout(1,1,4),stack(1,1,j+3),6)
c
        vals(j)=value1
        vals(j+1)=value2
        vals(j+2)=value3
        vals(j+3)=value4
c
        j=j+3
c
c        plot the newly created triangles
c
        if(iw .le. 0) goto 2200  
c        call triaplot(iw,zout(1,1,1))
c        call triaplot(iw,zout(1,1,2))
c        call triaplot(iw,zout(1,1,3))
c        call triaplot(iw,zout(1,1,4))
c
 2200 continue
c     
c       if the depth of the recursion has become excessive - bomb
c
        if(j .le. maxdepth) goto 3000
        ier=8
        return
 3000 continue
        ier=16
        return
        end
c
c
c
c  
        subroutine triaint1(n,vert1,vert2,vert3,fun,
     1      par1,par2,rint,w,work)
        implicit real *8 (a-h,o-z)
        dimension vert1(1),vert2(1),vert3(1),w(1),work(450)
        data nold/-10/
c
c        construct the quadrature formula on this triangle
c
        ifinit=1
        if(n .eq. nold) ifinit=0
c
        irnodes=1
        lrnodes=n**2 *2 +10
c
        iweights=irnodes+lrnodes
        lweights=n**2+5
c
        if( n .le. 50 ) then
        call triasymq(n,vert1,vert2,vert3,w(irnodes),
     1      w(iweights),numnodes)
        endif
c
        if( n .gt. 50 ) then
        call triagauc(n,vert1,vert2,vert3,w(irnodes),
     1      w(iweights),ifinit,work)
        numnodes=n*n
        endif
c
        call triagauc(n,vert1,vert2,vert3,w(irnodes),
     1      w(iweights),ifinit,work)
c
c       integrate the user-specified function 
c
        call triaint0(n,w(irnodes),w(iweights),fun,par1,par2,rint)
c
        return
        end
c
c
c
c
c
        subroutine triaint0(n,rnodes,weights,fun,par1,
     1    par2,rint)
        implicit real *8 (a-h,o-z)
        dimension rnodes(2,n,n),weights(n,n),par1(1),par2(1)
c
c       this subroutine integrates the user-specified function 
c       fun over a triangle in the plane; it uses the quadrature 
c       formula with n**2 nodes rnodes and weights weights; 
c       presumably, these have been created by a prior call to
c       the subroutine triagauc (see)
c
        rint=0
        do 1400 i=1,n
        do 1200 j=1,n
c
        call fun(rnodes(1,j,i),rnodes(2,j,i),par1,par2,d)

cccc        call prin2('in triinte, d=*',d,1)
c
        rint=rint+d*weights(j,i)
 1200 continue
 1400 continue
c
        return
        end
c
c
c
c
c
        
