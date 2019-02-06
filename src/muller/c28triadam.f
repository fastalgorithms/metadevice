c
c
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
        subroutine c28triaadam(ier,vert1,vert2,vert3,fun,nm,
     1      p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,m,eps,cints,
     $      maxrec,numfunev,w,nrec,ipatch,targinfo,info)
        implicit real *8 (a-h,o-z)
        dimension vert1(2),vert2(2),vert3(2)
        complex *16 w(1),cints(1)
c
c       This subroutine uses the the adaptive tensor product gaussian
c       quadrature to integrate a collection of complex *16
c       user-supplied functions R^2 \to R^1 on a triangle in the plane
c       (also user-supplied). In fact, this is simply a memory
c       management routine; all actual work is done by the subroutine
c       c28trianrem (see below).
c
c       Optionally, this subroutine can the adaptive symmetric
c       quadrature (for n.le.50) or the adaptive tensor product gaussian
c       quadrature (for n.ge.51)
c
c                       input parameters:
c
c  vert1,vert2,vert3 - the vertices in the plane of the triangle
c       over which the function is to be integrated
c  fun - the complex *16 user-supplied function to be integrated. the calling
c       sequence of fun must be 
c
c        call fun(x,y,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,f).     (1)
c
c        in (1), (x,y) is a point in the plane where 
c        the function is to be evaluated, and 
c        p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11
c        are twelve parameters to be used by fun; they can be 
c        variables or arrays, real or integer, as desired. 
c        f is assumed to be a complex *16 vector of length nm
c  nm - the length of the vectors cints in the calling sequence of the
c        subroutine c28triaadam (see above), and of the vector f in (1) 
c        above
c  par1, par2 - parameters to be used by the user-supplied 
c       subroutine fun (see above)
c  m - the order of the quadrature to me used on each subinterval
c  eps - the accuracy (absolute) to which the integral will be 
c       evaluated
c  nrec - the depth of recursion
c
c
c                       output parameters:
c
c  ier - error return code. 
c          ier=0 means normal conclusion
c          ier=4 means that at some point, the depth of recursion
c                reached nrec. 
c          ier=16 means that the total number of subtriangles in the
c                adaptive subdivision of [a,b] turned out to be greater 
c                than nnmax*4.  this is a fatal error.
c                
c  cints - the integrals as evaluated (nm of them)
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
c        . . . integrate the complex *16 user-supplied function using the 
c              adaptive gaussian quadratures
c
        nnmax=100000
        maxdepth=200
c
c        allocate memory for the subroutine c28trianrem
c
        istack=1
        lstack=1207
c
        iw=istack+lstack
        lw=3*m**2+50
c
        ivals=iw+lw
        lvals=207*nm+1000
c
cccc        call prinf('lvals=*',lvals,1)
c
        iww=ivals+lvals
        lww=4*m+10
c
        ivalue1=iww+lww
        lvalue1=nm+4
c
        ivalue2=ivalue1+lvalue1
        lvalue2=nm+4
c
        ivalue3=ivalue2+lvalue2
        lvalue3=nm+4
c
        ivalue4=ivalue3+lvalue3
        lvalue4=nm+4
c
        ivals2=ivalue4+lvalue4
        lvals2=nm+4
c
        iiiw=21
c
c       . . . integrate
c
        call c28trianrem(ier,w(istack),vert1,vert2,vert3,fun,nm,
     1      p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,
     $      w(iw),m,w(ivals),nnmax,eps,
     2      cints,maxdepth,maxrec,numint,iiiw,numfunev,w(iww),
     3      w(ivalue1),w(ivalue2),w(ivalue3),w(ivalue4),
     4      w(ivals2), nrec,ipatch,targinfo,info)
c
        return
        end
c
c
c
c
c
        subroutine c28trianrem(ier,stack,z1,z2,z3,fun,nm,
     1      p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,w,m,vals,nnmax,eps,
     2      cints,maxdepth,maxrec,numint,iw,numfunev,ww,
     3      value1,value2,value3,value4,vals2,nrec,ipatch,targinfo,info)
c
        implicit real *8 (a-h,o-z)
        dimension stack(2,3,1)
        dimension z1(2),z2(2),z3(2),zout(2,3,10)
        complex *16 ww(1),cints(1),w(1),vals(nm,1),
     1     value1(1),value2(1),value3(1),value4(1),vals2(1)
c
        dimension targinfo(12),info(*)
        dimension xyz0(3),dxyzduv0(3,2)
        dimension xyz1(3),dxyzduv1(3,2)
        dimension xyz2(3),dxyzduv2(3,2)
        dimension xyz3(3),dxyzduv3(3,2)
c
        dimension irec(207)
c
        external p0
c
c
c       This subroutine uses the adaptive symmetric quadrature (for
c       n.le.50) or the adaptive tensor product gaussian quadrature (for
c       n.ge.51) to integrate a complex *16 user-supplied function R^2 \to R^1
c       on a triangle in the plane (also user-supplied).
c
c                       input parameters:
c
c  z1,z2,z3 - the vertices in the plane of the triangle
c       over which the function is to be integrated
c  fun - the complex *16 user-supplied function to be integrated. the calling
c       sequence of fun must be 
c
c        call fun(x,y,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,f).     (1)
c
c        in (1), (x,y) is a point in the plane where 
c        the function is to be evaluated, and 
c        p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11
c        are twelve parameters to be used by fun; they can be 
c        variables or arrays, real or integer, as desired. 
c        f is assumed to be a complex *16 vector of length nm
c  nm - the length of the vectors cints in the calling sequence of the
c        subroutine c28triaadam (see above), and of the vector f in (1) 
c        above
c  par1, par2 - parameters to be used by the user-supplied 
c       subroutine fun (see above)
c  m - the order of the quadrature to be used on each subinterval
c  nnmax - maximum permitted number of subdivisions
c  eps - the accuracy (absolute) to which the integral will be 
c       evaluated
c  maxrec - maximum permitted depth of recursion
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
c  cints - the integrals as evaluated (nm of them)
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
        info(1)=0
        info(2)=0
c
cccc        call prinf('in c28trianrem, nm=*',nm,1)
c
        call triaarrm(z1,stack(1,1,1),2)
        call triaarrm(z2,stack(1,2,1),2)
        call triaarrm(z3,stack(1,3,1),2)
c
c        
c        plot the top level triangle
c
        if( iw .le. 0) goto 1000
cccc        call triaplot(iw,stack(1,1,1) )
c
 1000   continue
c
        call c28triainm1(m,z1,z2,z3,fun,
     1      p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,
     $      vals(1,1),w,ww,nm,vals2,numnodes)
        numfunev=numfunev+numnodes
c
c       recursively integrate the thing
c
        j=1
        do 1200 ij=1,nm
        cints(ij)=0
 1200 continue
c
        irec(j)=1
c
cccc        cints(1)=0
c
        ier=0
        maxrec=0
        do 3000 i=1,nnmax
        numint=i
        if(j .gt. maxrec) maxrec=j
c
c
c       first, check if the target point is in the far field 
c       and act accordingly
c
ccc        iffar=0
ccc        call p0(ipatch,z1(1),z1(2),xyz2,dxyzduv2,p1,p2,p3,p4)
ccc        call p0(ipatch,z2(1),z2(2),xyz2,dxyzduv2,p1,p2,p3,p4)
ccc        call p0(ipatch,z3(1),z3(2),xyz3,dxyzduv3,p1,p2,p3,p4)
ccc        xyz0(1)=(xyz1(1)+xyz2(1)+xyz3(1))/3
ccc        xyz0(2)=(xyz1(2)+xyz2(2)+xyz3(2))/3
ccc        xyz0(3)=(xyz1(3)+xyz2(3)+xyz3(3))/3
        iffar=0
        u1=stack(1,1,j)
        v1=stack(2,1,j)
        call p0(ipatch,u1,v1,xyz1,dxyzduv1,p1,p2,p3,p4)
        u2=stack(1,2,j)
        v2=stack(2,2,j)
        call p0(ipatch,u2,v2,xyz2,dxyzduv2,p1,p2,p3,p4)
        u3=stack(1,3,j)
        v3=stack(2,3,j)
        call p0(ipatch,u3,v3,xyz3,dxyzduv3,p1,p2,p3,p4)
        u0=(stack(1,1,j)+stack(1,2,j)+stack(1,3,j))/3
        v0=(stack(2,1,j)+stack(2,2,j)+stack(2,3,j))/3
        call p0(ipatch,u0,v0,xyz0,dxyzduv0,p1,p2,p3,p4)
        dx=targinfo(1)-xyz0(1)
        dy=targinfo(2)-xyz0(2)
        dz=targinfo(3)-xyz0(3)
        r0=sqrt(dx*dx+dy*dy+dz*dz)
        dx=xyz1(1)-xyz0(1)
        dy=xyz1(2)-xyz0(2)
        dz=xyz1(3)-xyz0(3)
        r1=sqrt(dx*dx+dy*dy+dz*dz)
        dx=xyz2(1)-xyz0(1)
        dy=xyz2(2)-xyz0(2)
        dz=xyz2(3)-xyz0(3)
        r2=sqrt(dx*dx+dy*dy+dz*dz)
        dx=xyz3(1)-xyz0(1)
        dy=xyz3(2)-xyz0(2)
        dz=xyz3(3)-xyz0(3)
        r3=sqrt(dx*dx+dy*dy+dz*dz)
c
        rmax=r1
        if( r2 .gt. rmax ) rmax=r2
        if( r3 .gt. rmax ) rmax=r3

ccc        if( r0 .gt. 1.4*rmax) iffar=1 
ccc        if( r0 .gt. 2*rmax) iffar=1 
        if( r0 .gt. 3*rmax) iffar=1 
c
        if( iffar .eq. 1 ) then        
ccc        call prinf('iffar=*',iffar,1)
        endif
c
c
        if( irec(j) .ge. nrec ) then 
        ier=4 
        iffar=1
        endif
c
c       
        if( iffar .eq. 0 ) then
c
c       subdivide the current triangle
c
        call triadiv(stack(1,1,j),stack(1,2,j),stack(1,3,j),zout)
c
c
        call c28triainm1(m,zout(1,1,1),zout(1,2,1),zout(1,3,1),fun,
     1     p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,
     $     value1,w,ww,nm,vals2,numnodes)
        numfunev=numfunev+numnodes
c
        call c28triainm1(m,zout(1,1,2),zout(1,2,2),zout(1,3,2),fun,
     1     p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,
     $     value2,w,ww,nm,vals2,numnodes)
        numfunev=numfunev+numnodes
c
        call c28triainm1(m,zout(1,1,3),zout(1,2,3),zout(1,3,3),fun,
     1     p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,
     $     value3,w,ww,nm,vals2,numnodes)
        numfunev=numfunev+numnodes
c
        call c28triainm1(m,zout(1,1,4),zout(1,2,4),zout(1,3,4),fun,
     1     p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,
     $     value4,w,ww,nm,vals2,numnodes)
        numfunev=numfunev+numnodes
c
        ifdone=1
c
        do 1600 ij=1,nm
c
        dd=abs(value1(ij)+value2(ij)+value3(ij)
     1      +value4(ij)-vals(ij,j))
c
        if(dd .gt. eps) ifdone=0
c
 1600 continue
c
c
c       if the function on this subinterval has been 
c       integrated with sufficient accuracy - add the 
c       value to that of the global integral and move up
c       in the stack
c
        if(ifdone  .eq. 0) goto 2000
c
        do 1800 ij=1,nm
        cints(ij)=cints(ij)+value1(ij)+value2(ij)+value3(ij)+value4(ij)
 1800 continue
c
cccc        cints(1)=cints(1)+value1(1)+value2(1)+value3(1)+value4(1)
        j=j-1
c
        info(1)=info(1)+1
c
        endif
c
        if( iffar .eq. 1 ) then
c
c       ... target is in the far field, 
c       add the contribution due to the smooth quadrature
c
c        call c28triainm1(m,stack(1,1,j),stack(1,2,j),stack(1,3,j),fun,
c     1      p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,vals(1,j),
c     $      w,ww,nm,vals2,numnodes)
c        numfunev=numfunev+numnodes
c
        do 1900 ij=1,nm
        cints(ij)=cints(ij)+vals(ij,j)
 1900 continue
c
        j=j-1
        ifdone=1
c
        info(2)=info(2)+1
c
        endif
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
        do 2200 ij=1,nm
c
        vals(ij,j)=value1(ij)
        vals(ij,j+1)=value2(ij)
        vals(ij,j+2)=value3(ij)
        vals(ij,j+3)=value4(ij)
 2200 continue
c
cccc        vals(j)=value1(1)
cccc        vals(j+1)=value2(1)
cccc        vals(j+2)=value3(1)
cccc        vals(j+3)=value4(1)
c       
        irectmp=irec(j)+1
        irec(j)=irectmp
        irec(j+1)=irectmp
        irec(j+2)=irectmp
        irec(j+3)=irectmp
c
        j=j+3
c

c        plot the newly created triangles
c
        if(iw .le. 0) goto 2300  
cccc        call triaplot(iw,zout(1,1,1) )
cccc        call triaplot(iw,zout(1,1,2) )
cccc        call triaplot(iw,zout(1,1,3) )
cccc        call triaplot(iw,zout(1,1,4) )
c
 2300 continue
c     
c       if the depth of the recursion has become excessive - bomb
c
        if((j+3) .le. maxdepth) goto 3000
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
c
        subroutine c28triainm1(n,vert1,vert2,vert3,fun,
     1      p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,
     $      cints,w,work,nm,vals,numnodes)
        implicit real *8 (a-h,o-z)
        dimension vert1(1),vert2(1),vert3(1)
        complex *16 cints(1),vals(1000),w(1),work(450)
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
        itype=1
        if( n .le. 50 ) itype = 1
        if( n .gt. 50 ) itype = 2
c
ccc        itype=2

        if( itype .eq. 1 ) then
        call triasymq(n,vert1,vert2,vert3,w(irnodes),
     1      w(iweights),numnodes)
        endif
c
        if( itype .eq. 2 ) then
        call triagauc(n,vert1,vert2,vert3,w(irnodes),
     1      w(iweights),ifinit,work)
        numnodes=n*n
        endif
c
c       integrate the user-specified function 
c
        call c28triainm0(numnodes,w(irnodes),w(iweights),fun,
     $     p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,
     1     cints,nm,vals)
c
        return
        end
c
c
c
c
c
        subroutine c28triainm0(numnodes,rnodes,weights,fun,
     $     p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,cints,nm,vals)
        implicit real *8 (a-h,o-z)
        dimension rnodes(2,numnodes),weights(numnodes)
        complex *16 vals(1),cints(1)
c
c       this subroutine integrates the user-specified function 
c       fun over a triangle in the plane; it uses the quadrature 
c       formula with numnodes nodes rnodes and weights weights; 
c
        do 1050 ij=1,nm
        cints(ij)=0
 1050 continue
c
        do 1400 i=1,numnodes
c
        call fun(rnodes(1,i),rnodes(2,i),
     $     p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,vals)
c
        do 1100 ij=1,nm
        cints(ij)=cints(ij)+vals(ij)*weights(i)
 1100 continue
c
 1400 continue
c
        return
        end
c
c
c
c
c
