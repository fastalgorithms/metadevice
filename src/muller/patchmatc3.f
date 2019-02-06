c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       the triangulation refinement routines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine genrefineinfo(noversamp,npatches,
     $     npatchesout,ipatchinfo,refineinfo)
        implicit real *8 (a-h,o-z)
        dimension ipatchinfo(1),refineinfo(4,1)
c
c       This subroutine generates a refinement table, used by patch
c       descriptor routines.
c
c       i--+
c       |\ |
c       |_\|
c
        npatchesout=npatches*noversamp**2
c
        done=1
        dx=done/noversamp
        dy=done/noversamp
        kk=0
c
        do 1400 i=1,npatches
c
        do 1300 iy=1,noversamp
        do 1200 ix=1,noversamp
c
        if( ix+iy .le. noversamp+1 ) then
c
        kk=kk+1
c
c       ... lower sub-triangles
c
        ipatchinfo(kk)=i
        refineinfo(1,kk)=dx
        refineinfo(2,kk)=dx*(ix-1)
        refineinfo(3,kk)=dy
        refineinfo(4,kk)=dy*(iy-1)
c
        endif
c
        if( ix+iy .le. noversamp ) then
c
        kk=kk+1
c
c       ... upper sub-triangles
c
        ipatchinfo(kk)=i
        refineinfo(1,kk)=-dx
        refineinfo(2,kk)=dx*(ix)
        refineinfo(3,kk)=-dy
        refineinfo(4,kk)=dy*(iy)
c
        endif
c
 1200   continue
 1300   continue
c
 1400   continue
c
ccc        call prinf('ipatchinfo=*',ipatchinfo,kk)
ccc        call prin2('refineinfo=*',refineinfo,4*kk)
ccc        pause
c
        return
        end
c
c
c
c
c
        subroutine genrefineinfo1(noversamp,npatches,
     $     npatchesout,w,lused)
        implicit real *8 (a-h,o-z)
        dimension ipatchinfo(1),refineinfo(4,1)
        dimension w(1)
c
        npatchesout=npatches*noversamp**2
c
        iipatchinfo=1
        lipatchinfo=npatchesout
c
        irefineinfo=iipatchinfo+lipatchinfo
        lrefineinfo=npatchesout*4
c
        lused=irefineinfo+lrefineinfo
c
        call genrefineinfo(noversamp,npatches,
     $     npatchesout,w(iipatchinfo),w(irefineinfo))
c
        return
        end
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       the patch descriptor routines 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine fpatchpnt
     $     (ipatch,u,v,xyz,dxyzduv,triainfo,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension xyz(3),dxyzduv(3,2)        
        dimension triainfo(3,3,1)
c
c       This subroutine maps a standard simplex triangle (0,0), (1,0), (0,1)
c       into R^3
c
c       ... setup a flat triangle in R^3
c
c
c            2     
c          .   . 
c         .     .
c        .       .
c       0 .. . .. 1
c
c
        x0=triainfo(1,1,ipatch)
        y0=triainfo(2,1,ipatch)
        z0=triainfo(3,1,ipatch)
c
        x1=triainfo(1,2,ipatch)
        y1=triainfo(2,2,ipatch)
        z1=triainfo(3,2,ipatch)
c
        x2=triainfo(1,3,ipatch)
        y2=triainfo(2,3,ipatch)
        z2=triainfo(3,3,ipatch)
c
c
c       ... process the geometry, return the point location and the
c       derivatives with respect to u and v
c
        x=x0+u*(x1-x0)+v*(x2-x0)
        y=y0+u*(y1-y0)+v*(y2-y0)
        z=z0+u*(z1-z0)+v*(z2-z0)
c
        xyz(1)=x
        xyz(2)=y
        xyz(3)=z
c
        dxyzduv(1,1)=x1-x0
        dxyzduv(1,2)=x2-x0
c
        dxyzduv(2,1)=y1-y0
        dxyzduv(2,2)=y2-y0
c
        dxyzduv(3,1)=z1-z0
        dxyzduv(3,2)=z2-z0
c
        return
        end
c
c
c
c
c
        subroutine qpatchpnt
     $     (ipatch,u,v,xyz,dxyzduv,qtriainfo,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension xyz(3),dxyzduv(3,2)        
        dimension qtriainfo(3,6,1)
c
c       This subroutine maps a standard simplex triangle (0,0), (1,0), (0,1)
c       into R^3
c
c
c       ... setup a quadratic triangle in R^3
c
c
c            2     
c          .   . 
c         C     B
c        .       .
c       0 .. A .. 1
c
c
        x0=qtriainfo(1,1,ipatch)
        y0=qtriainfo(2,1,ipatch)
        z0=qtriainfo(3,1,ipatch)
c
        x1=qtriainfo(1,2,ipatch)
        y1=qtriainfo(2,2,ipatch)
        z1=qtriainfo(3,2,ipatch)
c
        x2=qtriainfo(1,3,ipatch)
        y2=qtriainfo(2,3,ipatch)
        z2=qtriainfo(3,3,ipatch)
c
        xa=qtriainfo(1,4,ipatch)
        ya=qtriainfo(2,4,ipatch)
        za=qtriainfo(3,4,ipatch)
c
        xb=qtriainfo(1,5,ipatch)
        yb=qtriainfo(2,5,ipatch)
        zb=qtriainfo(3,5,ipatch)
c
        xc=qtriainfo(1,6,ipatch)
        yc=qtriainfo(2,6,ipatch)
        zc=qtriainfo(3,6,ipatch)
c
c
c       ... process the geometry, return the point location and the
c       derivatives with respect to u and v
c
c
        xu=-3*x0+4*xa-x1
        xv=-3*x0+4*xc-x2
        xuu=x0-2*xa+x1
        xuv=x0-xa+xb-xc
        xvv=x0-2*xc+x2
c
        yu=-3*y0+4*ya-y1
        yv=-3*y0+4*yc-y2
        yuu=y0-2*ya+y1
        yuv=y0-ya+yb-yc
        yvv=y0-2*yc+y2
c
        zu=-3*z0+4*za-z1
        zv=-3*z0+4*zc-z2
        zuu=z0-2*za+z1
        zuv=z0-za+zb-zc
        zvv=z0-2*zc+z2
c
        x=x0+u*xu+v*xv+2*(u*(u*xuu+v*xuv)+v*(u*xuv+v*xvv))
        y=y0+u*yu+v*yv+2*(u*(u*yuu+v*yuv)+v*(u*yuv+v*yvv))
        z=z0+u*zu+v*zv+2*(u*(u*zuu+v*zuv)+v*(u*zuv+v*zvv))
c
        xyz(1)=x
        xyz(2)=y
        xyz(3)=z
c
        dxyzduv(1,1)=xu+4*(u*xuu+v*xuv)
        dxyzduv(1,2)=xv+4*(u*xuv+v*xvv)
c
        dxyzduv(2,1)=yu+4*(u*yuu+v*yuv)
        dxyzduv(2,2)=yv+4*(u*yuv+v*yvv)
c
        dxyzduv(3,1)=zu+4*(u*zuu+v*zuv)
        dxyzduv(3,2)=zv+4*(u*zuv+v*zvv)
c
        return
        end
c
c
c
c
c
        subroutine cpatchpnt
     $     (ipatch,u,v,xyz,dxyzduv,ctriainfo,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension xyz(3),dxyzduv(3,2)        
        dimension ctriainfo(3,10,1)
c
c       This subroutine maps a standard simplex triangle (0,0), (1,0), (0,1)
c       into R^3
c
c
c       ... setup a cubic triangle in R^3
c
c
c              2
c             . .     
c           C1   B2 
c          .       .
c         C2   M   B1
c        .           .
c       0 . A1 . A2 . 1
c
c
        x0=ctriainfo(1,1,ipatch)
        y0=ctriainfo(2,1,ipatch)
        z0=ctriainfo(3,1,ipatch)
c
        x1=ctriainfo(1,2,ipatch)
        y1=ctriainfo(2,2,ipatch)
        z1=ctriainfo(3,2,ipatch)
c
        x2=ctriainfo(1,3,ipatch)
        y2=ctriainfo(2,3,ipatch)
        z2=ctriainfo(3,3,ipatch)
c
        xa1=ctriainfo(1,4,ipatch)
        ya1=ctriainfo(2,4,ipatch)
        za1=ctriainfo(3,4,ipatch)
c
        xa2=ctriainfo(1,5,ipatch)
        ya2=ctriainfo(2,5,ipatch)
        za2=ctriainfo(3,5,ipatch)
c
        xb1=ctriainfo(1,6,ipatch)
        yb1=ctriainfo(2,6,ipatch)
        zb1=ctriainfo(3,6,ipatch)
c
        xb2=ctriainfo(1,7,ipatch)
        yb2=ctriainfo(2,7,ipatch)
        zb2=ctriainfo(3,7,ipatch)
c
        xc1=ctriainfo(1,8,ipatch)
        yc1=ctriainfo(2,8,ipatch)
        zc1=ctriainfo(3,8,ipatch)
c
        xc2=ctriainfo(1,9,ipatch)
        yc2=ctriainfo(2,9,ipatch)
        zc2=ctriainfo(3,9,ipatch)
c
        xm=ctriainfo(1,10,ipatch)
        ym=ctriainfo(2,10,ipatch)
        zm=ctriainfo(3,10,ipatch)
c
c
ccc        write(*,*) x0,x1,x2,xa1,xa2,xb1,xb2,xc1,xc2,xm
ccc        write(*,*) y0,y1,y2,ya1,ya2,yb1,yb2,yc1,yc2,ym
ccc        write(*,*) z0,z1,z2,za1,za2,zb1,zb2,zc1,zc2,zm
c
c       ... process the geometry, return the point location and the
c       derivatives with respect to u and v
c
c
        x=x0+
     $     u*(-11/2d0*x0+x1+9*xa1-9/2d0*xa2)+
     $     v*(-11/2d0*x0+x2+9*xc2-9/2d0*xc1)+
     $     u**2*(18*xa2+9*x0-9/2d0*x1-45/2d0*xa1)+
     $     2*u*v*(9/4d0*(xa2+xc1-xb2-xb1))+
     $     2*u*v*(9*x0-45/4d0*(xc2+xa1)+27/2d0*xm)+
     $     v**2*(18*xc1+9*x0-9/2d0*x2-45/2d0*xc2)+
     $     u**3*(9/2d0*(x1-x0)+27/2d0*(xa1-xa2))+
     $     3*u**2*v*(9*(xa1-xm)+9/2d0*(xc2+xb1-xa2-x0))+
     $     3*u*v**2*(9*(xc2-xm)+9/2d0*(xa1+xb2-xc1-x0))+
     $     v**3*(9/2d0*(x2-x0)+27/2d0*(xc2-xc1)) 
c
        y=y0+
     $     u*(-11/2d0*y0+y1+9*ya1-9/2d0*ya2)+
     $     v*(-11/2d0*y0+y2+9*yc2-9/2d0*yc1)+
     $     u**2*(18*ya2+9*y0-9/2d0*y1-45/2d0*ya1)+
     $     2*u*v*(9/4d0*(ya2+yc1-yb2-yb1))+
     $     2*u*v*(9*y0-45/4d0*(yc2+ya1)+27/2d0*ym)+
     $     v**2*(18*yc1+9*y0-9/2d0*y2-45/2d0*yc2)+
     $     u**3*(9/2d0*(y1-y0)+27/2d0*(ya1-ya2))+
     $     3*u**2*v*(9*(ya1-ym)+9/2d0*(yc2+yb1-ya2-y0))+
     $     3*u*v**2*(9*(yc2-ym)+9/2d0*(ya1+yb2-yc1-y0))+
     $     v**3*(9/2d0*(y2-y0)+27/2d0*(yc2-yc1)) 
c
        z=z0+
     $     u*(-11/2d0*z0+z1+9*za1-9/2d0*za2)+
     $     v*(-11/2d0*z0+z2+9*zc2-9/2d0*zc1)+
     $     u**2*(18*za2+9*z0-9/2d0*z1-45/2d0*za1)+
     $     2*u*v*(9/4d0*(za2+zc1-zb2-zb1))+
     $     2*u*v*(9*z0-45/4d0*(zc2+za1)+27/2d0*zm)+
     $     v**2*(18*zc1+9*z0-9/2d0*z2-45/2d0*zc2)+
     $     u**3*(9/2d0*(z1-z0)+27/2d0*(za1-za2))+
     $     3*u**2*v*(9*(za1-zm)+9/2d0*(zc2+zb1-za2-z0))+
     $     3*u*v**2*(9*(zc2-zm)+9/2d0*(za1+zb2-zc1-z0))+
     $     v**3*(9/2d0*(z2-z0)+27/2d0*(zc2-zc1)) 
c
        xyz(1)=x
        xyz(2)=y
        xyz(3)=z
c
        dxyzduv(1,1)=
     $     (-11/2d0*x0+x1+9*xa1-9/2d0*xa2)+
     $     2*u*(18*xa2+9*x0-9/2d0*x1-45/2d0*xa1)+
     $     2*v*(9/4d0*(xa2+xc1-xb2-xb1))+
     $     2*v*(9*x0-45/4d0*(xc2+xa1)+27/2d0*xm)+
     $     3*u**2*(9/2d0*(x1-x0)+27/2d0*(xa1-xa2))+
     $     3*2*u*v*(9*(xa1-xm)+9/2d0*(xc2+xb1-xa2-x0))+
     $     3*v**2*(9*(xc2-xm)+9/2d0*(xa1+xb2-xc1-x0))
c
        dxyzduv(1,2)=
     $     (-11/2d0*x0+x2+9*xc2-9/2d0*xc1)+
     $     2*u*(9/4d0*(xa2+xc1-xb2-xb1))+
     $     2*u*(9*x0-45/4d0*(xc2+xa1)+27/2d0*xm)+
     $     2*v*(18*xc1+9*x0-9/2d0*x2-45/2d0*xc2)+
     $     3*u**2*(9*(xa1-xm)+9/2d0*(xc2+xb1-xa2-x0))+
     $     3*u*2*v*(9*(xc2-xm)+9/2d0*(xa1+xb2-xc1-x0))+
     $     3*v**2*(9/2d0*(x2-x0)+27/2d0*(xc2-xc1)) 
c
        dxyzduv(2,1)=
     $     (-11/2d0*y0+y1+9*ya1-9/2d0*ya2)+
     $     2*u*(18*ya2+9*y0-9/2d0*y1-45/2d0*ya1)+
     $     2*v*(9/4d0*(ya2+yc1-yb2-yb1))+
     $     2*v*(9*y0-45/4d0*(yc2+ya1)+27/2d0*ym)+
     $     3*u**2*(9/2d0*(y1-y0)+27/2d0*(ya1-ya2))+
     $     3*2*u*v*(9*(ya1-ym)+9/2d0*(yc2+yb1-ya2-y0))+
     $     3*v**2*(9*(yc2-ym)+9/2d0*(ya1+yb2-yc1-y0))
c
        dxyzduv(2,2)=
     $     (-11/2d0*y0+y2+9*yc2-9/2d0*yc1)+
     $     2*u*(9/4d0*(ya2+yc1-yb2-yb1))+
     $     2*u*(9*y0-45/4d0*(yc2+ya1)+27/2d0*ym)+
     $     2*v*(18*yc1+9*y0-9/2d0*y2-45/2d0*yc2)+
     $     3*u**2*(9*(ya1-ym)+9/2d0*(yc2+yb1-ya2-y0))+
     $     3*u*2*v*(9*(yc2-ym)+9/2d0*(ya1+yb2-yc1-y0))+
     $     3*v**2*(9/2d0*(y2-y0)+27/2d0*(yc2-yc1)) 
c
        dxyzduv(3,1)=
     $     (-11/2d0*z0+z1+9*za1-9/2d0*za2)+
     $     2*u*(18*za2+9*z0-9/2d0*z1-45/2d0*za1)+
     $     2*v*(9/4d0*(za2+zc1-zb2-zb1))+
     $     2*v*(9*z0-45/4d0*(zc2+za1)+27/2d0*zm)+
     $     3*u**2*(9/2d0*(z1-z0)+27/2d0*(za1-za2))+
     $     3*2*u*v*(9*(za1-zm)+9/2d0*(zc2+zb1-za2-z0))+
     $     3*v**2*(9*(zc2-zm)+9/2d0*(za1+zb2-zc1-z0))
c
        dxyzduv(3,2)=
     $     (-11/2d0*z0+z2+9*zc2-9/2d0*zc1)+
     $     2*u*(9/4d0*(za2+zc1-zb2-zb1))+
     $     2*u*(9*z0-45/4d0*(zc2+za1)+27/2d0*zm)+
     $     2*v*(18*zc1+9*z0-9/2d0*z2-45/2d0*zc2)+
     $     3*u**2*(9*(za1-zm)+9/2d0*(zc2+zb1-za2-z0))+
     $     3*u*2*v*(9*(zc2-zm)+9/2d0*(za1+zb2-zc1-z0))+
     $     3*v**2*(9/2d0*(z2-z0)+27/2d0*(zc2-zc1)) 
c
        return
        end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine rpatchpnt
     $     (ipatch,u,v,xyz,dxyzduv,par1,patchpnt,ipatchinfo,refineinfo)
        implicit real *8 (a-h,o-z)
        dimension xyz(3),dxyzduv(3,2)        
        dimension ipatchinfo(1),refineinfo(4,1)
        external patchpnt
c
c       This subroutine provides on-fly refinement routines for arbitrary
c       patches
c
c
c       ... refine the triangulation
c       
        jpatch=ipatchinfo(ipatch)
        uref=u*refineinfo(1,ipatch)+refineinfo(2,ipatch)
        vref=v*refineinfo(3,ipatch)+refineinfo(4,ipatch)
c
        call patchpnt(jpatch,uref,vref,xyz,dxyzduv,par1,
     $     xpar2,xpar3,xpar4)
c
c       ... adjust all derivatives
c       
        dxyzduv(1,1)=dxyzduv(1,1)*refineinfo(1,ipatch)
        dxyzduv(2,1)=dxyzduv(2,1)*refineinfo(1,ipatch)
        dxyzduv(3,1)=dxyzduv(3,1)*refineinfo(1,ipatch)
c
        dxyzduv(1,2)=dxyzduv(1,2)*refineinfo(3,ipatch)
        dxyzduv(2,2)=dxyzduv(2,2)*refineinfo(3,ipatch)
        dxyzduv(3,2)=dxyzduv(3,2)*refineinfo(3,ipatch)
c
        return
        end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine spatchpnt
     $     (ipatch,u,v,xyz,dxyzduv,triainfo,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension xyz(3),dxyzduv(3,2)        
        dimension triainfo(3,3,1)
c
c       This subroutine maps a standard simplex triangle (0,0), (1,0), (0,1)
c       into R^3
c
c       Contruct an analytical triangulation of a unit sphere.
c
c       ... setup a flat triangle in R^3
c
        x0=triainfo(1,1,ipatch)
        y0=triainfo(2,1,ipatch)
        z0=triainfo(3,1,ipatch)
c
        x1=triainfo(1,2,ipatch)
        y1=triainfo(2,2,ipatch)
        z1=triainfo(3,2,ipatch)
c
        x2=triainfo(1,3,ipatch)
        y2=triainfo(2,3,ipatch)
        z2=triainfo(3,3,ipatch)
c
c
c       ... process the geometry, return the point location and the
c       derivatives with respect to u and v
c
        x=x0+u*(x1-x0)+v*(x2-x0)
        y=y0+u*(y1-y0)+v*(y2-y0)
        z=z0+u*(z1-z0)+v*(z2-z0)
c
c
c       ... map the constructed flat patch onto the unit sphere
c
        r=sqrt(x**2+y**2+z**2)
c
        xyz(1)=x/r
        xyz(2)=y/r
        xyz(3)=z/r
c
        ru=x*(x1-x0)+y*(y1-y0)+z*(z1-z0)
        rv=x*(x2-x0)+y*(y2-y0)+z*(z2-z0)
c
        dxyzduv(1,1)=(x1-x0)/r-x/r**3*ru
        dxyzduv(1,2)=(x2-x0)/r-x/r**3*rv
c
        dxyzduv(2,1)=(y1-y0)/r-y/r**3*ru
        dxyzduv(2,2)=(y2-y0)/r-y/r**3*rv
c
        dxyzduv(3,1)=(z1-z0)/r-z/r**3*ru
        dxyzduv(3,2)=(z2-z0)/r-z/r**3*rv
c
        sx=1.0d0
        sy=1.0d0
        sz=1.0d0
c
        xyz(1)=xyz(1)*sx
        xyz(2)=xyz(2)*sy
        xyz(3)=xyz(3)*sz
c
        do i=1,2
        dxyzduv(1,i)=dxyzduv(1,i)*sx
        dxyzduv(2,i)=dxyzduv(2,i)*sy
        dxyzduv(3,i)=dxyzduv(3,i)*sz
        enddo        
c
        return
        end
c
c
c
c
c
        subroutine tpatchpnt
     $     (ipatch,u,v,xyz,dxyzduv,triainfo,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension xyz(3),dxyzduv(3,2)        
        dimension triainfo(3,3,1)
c
c       Contruct an analytical triangulation of a unit torus.
c
c       ... setup a flat triangle in R^3
c
        x0=triainfo(1,1,ipatch)
        y0=triainfo(2,1,ipatch)
        z0=triainfo(3,1,ipatch)
c
        x1=triainfo(1,2,ipatch)
        y1=triainfo(2,2,ipatch)
        z1=triainfo(3,2,ipatch)
c
        x2=triainfo(1,3,ipatch)
        y2=triainfo(2,3,ipatch)
        z2=triainfo(3,3,ipatch)
c
c
c       ... process the geometry, return the point location and the
c       derivatives with respect to u and v
c
        x=x0+u*(x1-x0)+v*(x2-x0)
        y=y0+u*(y1-y0)+v*(y2-y0)
        z=z0+u*(z1-z0)+v*(z2-z0)
c
c
c       ... make it a torus, we'll assume that triangle is in R^2
c       [0..2*pi] x [0..2*pi] -> torus
c
        a=3
        b=2
c
        s=x
        t=y
c
c       ... map the constructed flat patch onto the torus
c
        xyz(1)=cos(s)*(a+b*cos(t))    
        xyz(2)=sin(s)*(a+b*cos(t))
        xyz(3)=b*sin(t)
c
c       ... shift to the origin
c        
        xyz(1)=xyz(1)-a
c
        dxyzduv(1,1)=-sin(s)*(a+b*cos(t))
        dxyzduv(1,2)=cos(s)*(-b*sin(t))
c
        dxyzduv(2,1)=cos(s)*(a+b*cos(t))
        dxyzduv(2,2)=sin(s)*(-b*sin(t))
c
        dxyzduv(3,1)=0
        dxyzduv(3,2)=b*cos(t)
c
c       
c       ... the following assumes that the triangulation of
c       [0..2*pi] x [0..2*pi] contains exactly two patches
c
        done=1
        pi=4*atan(done)
c
        dxyzduv(1,1)=dxyzduv(1,1)*2*pi
        dxyzduv(1,2)=dxyzduv(1,2)*2*pi
c
        dxyzduv(2,1)=dxyzduv(2,1)*2*pi
        dxyzduv(2,2)=dxyzduv(2,2)*2*pi
c
        dxyzduv(3,1)=dxyzduv(3,1)*2*pi
        dxyzduv(3,2)=dxyzduv(3,2)*2*pi
c
        return
c
        sx=1.2d0
        sy=1.1d0
        sz=1.0d0
c
        xyz(1)=xyz(1)*sx
        xyz(2)=xyz(2)*sy
        xyz(3)=xyz(3)*sz
c
        do i=1,2
        dxyzduv(1,i)=dxyzduv(1,i)*sx
        dxyzduv(2,i)=dxyzduv(2,i)*sy
        dxyzduv(3,i)=dxyzduv(3,i)*sz
        enddo        
c
        return
        end
c
c
c
c
c
        subroutine wpatchpnt
     $     (ipatch,u,v,xyz,dxyzduv,par1,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension xyz(3),dxyzduv(3,2)        
c
c       This subroutine maps a standard simplex triangle (0,0), (1,0), (0,1)
c       into R^3. This is just an example, one triangle only.
c
c       ... setup a flat triangle in R^3
c
        x0=1
        y0=0
        z0=0
c
        x1=0
        y1=1
        z1=0
c
        x2=0
        y2=0
        z2=1
c
c       ... process the geometry, return the point location and the
c       derivatives with respect to u and v
c
        x=x0+u*(x1-x0)+v*(x2-x0)
        y=y0+u*(y1-y0)+v*(y2-y0)
        z=z0+u*(z1-z0)+v*(z2-z0)
c
        xyz(1)=x
        xyz(2)=y
        xyz(3)=z
c
        dxyzduv(1,1)=x1-x0
        dxyzduv(1,2)=x2-x0
c
        dxyzduv(2,1)=y1-y0
        dxyzduv(2,2)=y2-y0
c
        dxyzduv(3,1)=z1-z0
        dxyzduv(3,2)=z2-z0
c
        return
        end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       the extended patch descriptor routines 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine fpatchpnt_ext
     $     (ipatch,u,v,xyz,dxyzduv,triainfo,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension xyz(3),dxyzduv(3,5)        
        dimension triainfo(3,3,1)
c
c       This subroutine maps a standard simplex triangle (0,0), (1,0), (0,1)
c       into R^3
c
c       ... setup a flat triangle in R^3
c
c
c            2     
c          .   . 
c         .     .
c        .       .
c       0 .. . .. 1
c
c
        x0=triainfo(1,1,ipatch)
        y0=triainfo(2,1,ipatch)
        z0=triainfo(3,1,ipatch)
c
        x1=triainfo(1,2,ipatch)
        y1=triainfo(2,2,ipatch)
        z1=triainfo(3,2,ipatch)
c
        x2=triainfo(1,3,ipatch)
        y2=triainfo(2,3,ipatch)
        z2=triainfo(3,3,ipatch)
c
c
c       ... process the geometry, return the point location and the
c       derivatives with respect to u and v
c
        x=x0+u*(x1-x0)+v*(x2-x0)
        y=y0+u*(y1-y0)+v*(y2-y0)
        z=z0+u*(z1-z0)+v*(z2-z0)
c
        xyz(1)=x
        xyz(2)=y
        xyz(3)=z
c
        dxyzduv(1,1)=x1-x0
        dxyzduv(1,2)=x2-x0
c
        dxyzduv(2,1)=y1-y0
        dxyzduv(2,2)=y2-y0
c
        dxyzduv(3,1)=z1-z0
        dxyzduv(3,2)=z2-z0
c
c
c       ... second derivatives with respect to u and v
c
        dxyzduv(1,3)=0
        dxyzduv(1,4)=0
        dxyzduv(1,5)=0
c
        dxyzduv(2,3)=0
        dxyzduv(2,4)=0
        dxyzduv(2,5)=0
c
        dxyzduv(3,3)=0
        dxyzduv(3,4)=0
        dxyzduv(3,5)=0
c
        return
        end
c
c
c
c
c
        subroutine qpatchpnt_ext
     $     (ipatch,u,v,xyz,dxyzduv,qtriainfo,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension xyz(3),dxyzduv(3,5)        
        dimension qtriainfo(3,6,1)
c
c       This subroutine maps a standard simplex triangle (0,0), (1,0), (0,1)
c       into R^3
c
c
c       ... setup a quadratic triangle in R^3
c
c
c            2     
c          .   . 
c         C     B
c        .       .
c       0 .. A .. 1
c
c
        x0=qtriainfo(1,1,ipatch)
        y0=qtriainfo(2,1,ipatch)
        z0=qtriainfo(3,1,ipatch)
c
        x1=qtriainfo(1,2,ipatch)
        y1=qtriainfo(2,2,ipatch)
        z1=qtriainfo(3,2,ipatch)
c
        x2=qtriainfo(1,3,ipatch)
        y2=qtriainfo(2,3,ipatch)
        z2=qtriainfo(3,3,ipatch)
c
        xa=qtriainfo(1,4,ipatch)
        ya=qtriainfo(2,4,ipatch)
        za=qtriainfo(3,4,ipatch)
c
        xb=qtriainfo(1,5,ipatch)
        yb=qtriainfo(2,5,ipatch)
        zb=qtriainfo(3,5,ipatch)
c
        xc=qtriainfo(1,6,ipatch)
        yc=qtriainfo(2,6,ipatch)
        zc=qtriainfo(3,6,ipatch)
c
c
c       ... process the geometry, return the point location and the
c       derivatives with respect to u and v
c
c
        x=x0+u*(-3*x0+4*xa-x1)+v*(-3*x0+4*xc-x2)+
     $     2*(u**2*(x0-2*xa+x1)+u*v*2*(x0-xa+xb-xc)+v**2*(x0-2*xc+x2))
c
        y=y0+u*(-3*y0+4*ya-y1)+v*(-3*y0+4*yc-y2)+
     $     2*(u**2*(y0-2*ya+y1)+u*v*2*(y0-ya+yb-yc)+v**2*(y0-2*yc+y2))
c
        z=z0+u*(-3*z0+4*za-z1)+v*(-3*z0+4*zc-z2)+
     $     2*(u**2*(z0-2*za+z1)+u*v*2*(z0-za+zb-zc)+v**2*(z0-2*zc+z2))
c
        xyz(1)=x
        xyz(2)=y
        xyz(3)=z
c
        dxyzduv(1,1)=(-3*x0+4*xa-x1)+4*(u*(x0-2*xa+x1)+v*(x0-xa+xb-xc))
        dxyzduv(1,2)=(-3*x0+4*xc-x2)+4*(u*(x0-xa+xb-xc)+v*(x0-2*xc+x2))
c
        dxyzduv(2,1)=(-3*y0+4*ya-y1)+4*(u*(y0-2*ya+y1)+v*(y0-ya+yb-yc))
        dxyzduv(2,2)=(-3*y0+4*yc-y2)+4*(u*(y0-ya+yb-yc)+v*(y0-2*yc+y2))
c
        dxyzduv(3,1)=(-3*z0+4*za-z1)+4*(u*(z0-2*za+z1)+v*(z0-za+zb-zc))
        dxyzduv(3,2)=(-3*z0+4*zc-z2)+4*(u*(z0-za+zb-zc)+v*(z0-2*zc+z2))
c
c
c       ... second derivatives with respect to u and v
c
        dxyzduv(1,3)=4*(x0-2*xa+x1)
        dxyzduv(1,4)=4*(x0-xa+xb-xc)
        dxyzduv(1,5)=4*(x0-2*xc+x2)
c
        dxyzduv(2,3)=4*(y0-2*ya+y1)
        dxyzduv(2,4)=4*(y0-ya+yb-yc)
        dxyzduv(2,5)=4*(y0-2*yc+y2)
c
        dxyzduv(3,3)=4*(z0-2*za+z1)
        dxyzduv(3,4)=4*(z0-za+zb-zc)
        dxyzduv(3,5)=4*(z0-2*zc+z2)
c
        return
        end
c
c
c
c
c
        subroutine rpatchpnt_ext
     $     (ipatch,u,v,xyz,dxyzduv,par1,patchpnt,ipatchinfo,refineinfo)
        implicit real *8 (a-h,o-z)
        dimension xyz(3),dxyzduv(3,5)        
        dimension ipatchinfo(1),refineinfo(4,1)
        external patchpnt
c
c       This subroutine provides on-fly refinement routines for arbitrary
c       patches
c
c
c       ... refine the triangulation
c       
        jpatch=ipatchinfo(ipatch)
        uref=u*refineinfo(1,ipatch)+refineinfo(2,ipatch)
        vref=v*refineinfo(3,ipatch)+refineinfo(4,ipatch)
c
        call patchpnt(jpatch,uref,vref,xyz,dxyzduv,par1,
     $     xpar2,xpar3,xpar4)
c
c       ... adjust all derivatives
c       
        dxyzduv(1,1)=dxyzduv(1,1)*refineinfo(1,ipatch)
        dxyzduv(2,1)=dxyzduv(2,1)*refineinfo(1,ipatch)
        dxyzduv(3,1)=dxyzduv(3,1)*refineinfo(1,ipatch)
c
        dxyzduv(1,2)=dxyzduv(1,2)*refineinfo(3,ipatch)
        dxyzduv(2,2)=dxyzduv(2,2)*refineinfo(3,ipatch)
        dxyzduv(3,2)=dxyzduv(3,2)*refineinfo(3,ipatch)
c
        dxyzduv(1,3)=dxyzduv(1,3)*refineinfo(1,ipatch)
        dxyzduv(1,3)=dxyzduv(1,3)*refineinfo(1,ipatch)
        dxyzduv(1,4)=dxyzduv(1,4)*refineinfo(1,ipatch)
        dxyzduv(1,4)=dxyzduv(1,4)*refineinfo(3,ipatch)
        dxyzduv(1,5)=dxyzduv(1,5)*refineinfo(3,ipatch)
        dxyzduv(1,5)=dxyzduv(1,5)*refineinfo(3,ipatch)
c
        dxyzduv(2,3)=dxyzduv(2,3)*refineinfo(1,ipatch)
        dxyzduv(2,3)=dxyzduv(2,3)*refineinfo(1,ipatch)
        dxyzduv(2,4)=dxyzduv(2,4)*refineinfo(1,ipatch)
        dxyzduv(2,4)=dxyzduv(2,4)*refineinfo(3,ipatch)
        dxyzduv(2,5)=dxyzduv(2,5)*refineinfo(3,ipatch)
        dxyzduv(2,5)=dxyzduv(2,5)*refineinfo(3,ipatch)
c
        dxyzduv(3,3)=dxyzduv(3,3)*refineinfo(1,ipatch)
        dxyzduv(3,3)=dxyzduv(3,3)*refineinfo(1,ipatch)
        dxyzduv(3,4)=dxyzduv(3,4)*refineinfo(1,ipatch)
        dxyzduv(3,4)=dxyzduv(3,4)*refineinfo(3,ipatch)
        dxyzduv(3,5)=dxyzduv(3,5)*refineinfo(3,ipatch)
        dxyzduv(3,5)=dxyzduv(3,5)*refineinfo(3,ipatch)
c
        return
        end
c
c
c
        subroutine spatchpnt_ext
     $     (ipatch,u,v,xyz,dxyzduv,triainfo,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension xyz(3),dxyzduv(3,5)        
        dimension triainfo(3,3,1)
c
c       This subroutine maps a standard simplex triangle (0,0), (1,0), (0,1)
c       into R^3
c
c       Contruct an analytical triangulation of a unit sphere.
c
c       ... setup a flat triangle in R^3
c
        x0=triainfo(1,1,ipatch)
        y0=triainfo(2,1,ipatch)
        z0=triainfo(3,1,ipatch)
c
        x1=triainfo(1,2,ipatch)
        y1=triainfo(2,2,ipatch)
        z1=triainfo(3,2,ipatch)
c
        x2=triainfo(1,3,ipatch)
        y2=triainfo(2,3,ipatch)
        z2=triainfo(3,3,ipatch)
c
c
c       ... process the geometry, return the point location and the
c       derivatives with respect to u and v
c
        x=x0+u*(x1-x0)+v*(x2-x0)
        y=y0+u*(y1-y0)+v*(y2-y0)
        z=z0+u*(z1-z0)+v*(z2-z0)
c
c
c       ... map the constructed flat patch onto the unit sphere
c
        r=sqrt(x**2+y**2+z**2)
c
        xyz(1)=x/r
        xyz(2)=y/r
        xyz(3)=z/r
c
        ru=x*(x1-x0)+y*(y1-y0)+z*(z1-z0)
        rv=x*(x2-x0)+y*(y2-y0)+z*(z2-z0)
c
        dxyzduv(1,1)=(x1-x0)/r-x/r**3*ru
        dxyzduv(1,2)=(x2-x0)/r-x/r**3*rv
c
        dxyzduv(2,1)=(y1-y0)/r-y/r**3*ru
        dxyzduv(2,2)=(y2-y0)/r-y/r**3*rv
c
        dxyzduv(3,1)=(z1-z0)/r-z/r**3*ru
        dxyzduv(3,2)=(z2-z0)/r-z/r**3*rv
c
        ruu=(x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0)
        ruv=(x1-x0)*(x2-x0)+(y1-y0)*(y2-y0)+(z1-z0)*(z2-z0)
        rvv=(x2-x0)*(x2-x0)+(y2-y0)*(y2-y0)+(z2-z0)*(z2-z0)
c
        dxyzduv(1,3)=-(x1-x0)/r**3*ru-(x1-x0)/r**3*ru
     $     -x*3/r**2*(-ru/r**3)*ru-x/r**3*ruu
        dxyzduv(1,4)=-(x1-x0)/r**3*rv-(x2-x0)/r**3*ru
     $     -x*3/r**2*(-rv/r**3)*ru-x/r**3*ruv
        dxyzduv(1,5)=-(x2-x0)/r**3*rv-(x2-x0)/r**3*rv
     $     -x*3/r**2*(-rv/r**3)*rv-x/r**3*rvv
c
        dxyzduv(2,3)=-(y1-y0)/r**3*ru-(y1-y0)/r**3*ru
     $     -y*3/r**2*(-ru/r**3)*ru-y/r**3*ruu
        dxyzduv(2,4)=-(y1-y0)/r**3*rv-(y2-y0)/r**3*ru
     $     -y*3/r**2*(-rv/r**3)*ru-y/r**3*ruv
        dxyzduv(2,5)=-(y2-y0)/r**3*rv-(y2-y0)/r**3*rv
     $     -y*3/r**2*(-rv/r**3)*rv-y/r**3*rvv
c
        dxyzduv(3,3)=-(z1-z0)/r**3*ru-(z1-z0)/r**3*ru
     $     -z*3/r**2*(-ru/r**3)*ru-z/r**3*ruu
        dxyzduv(3,4)=-(z1-z0)/r**3*rv-(z2-z0)/r**3*ru
     $     -z*3/r**2*(-rv/r**3)*ru-z/r**3*ruv
        dxyzduv(3,5)=-(z2-z0)/r**3*rv-(z2-z0)/r**3*rv
     $     -z*3/r**2*(-rv/r**3)*rv-z/r**3*rvv
c
        sx=1.0d0
        sy=1.0d0
        sz=1.0d0
c
        xyz(1)=xyz(1)*sx
        xyz(2)=xyz(2)*sy
        xyz(3)=xyz(3)*sz
c
        do i=1,2
        dxyzduv(1,i)=dxyzduv(1,i)*sx
        dxyzduv(2,i)=dxyzduv(2,i)*sy
        dxyzduv(3,i)=dxyzduv(3,i)*sz
        enddo        
c
        do i=3,5
        dxyzduv(1,i)=dxyzduv(1,i)*sx
        dxyzduv(2,i)=dxyzduv(2,i)*sy
        dxyzduv(3,i)=dxyzduv(3,i)*sz
        enddo        
c
        return
        end
c
c
c
c
c
        subroutine tpatchpnt_ext
     $     (ipatch,u,v,xyz,dxyzduv,triainfo,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension xyz(3),dxyzduv(3,5)        
        dimension triainfo(3,3,1)
c
c       Contruct an analytical triangulation of a unit torus.
c
c       ... setup a flat triangle in R^3
c
        x0=triainfo(1,1,ipatch)
        y0=triainfo(2,1,ipatch)
        z0=triainfo(3,1,ipatch)
c
        x1=triainfo(1,2,ipatch)
        y1=triainfo(2,2,ipatch)
        z1=triainfo(3,2,ipatch)
c
        x2=triainfo(1,3,ipatch)
        y2=triainfo(2,3,ipatch)
        z2=triainfo(3,3,ipatch)
c
c
c       ... process the geometry, return the point location and the
c       derivatives with respect to u and v
c
        x=x0+u*(x1-x0)+v*(x2-x0)
        y=y0+u*(y1-y0)+v*(y2-y0)
        z=z0+u*(z1-z0)+v*(z2-z0)
c
c
c       ... make it a torus, we'll assume that triangle is in R^2
c       [0..2*pi] x [0..2*pi] -> torus
c
        a=3
        b=2
c
        s=x
        t=y
c
c       ... map the constructed flat patch onto the torus
c
        xyz(1)=cos(s)*(a+b*cos(t))    
        xyz(2)=sin(s)*(a+b*cos(t))
        xyz(3)=b*sin(t)
c
c       ... shift to the origin
c        
        xyz(1)=xyz(1)-a
c
        dxyzduv(1,1)=-sin(s)*(a+b*cos(t))
        dxyzduv(1,2)=cos(s)*(-b*sin(t))
c
        dxyzduv(2,1)=cos(s)*(a+b*cos(t))
        dxyzduv(2,2)=sin(s)*(-b*sin(t))
c
        dxyzduv(3,1)=0
        dxyzduv(3,2)=b*cos(t)
c
c
        dxyzduv(1,3)=-cos(s)*(a+b*cos(t))
        dxyzduv(1,4)=sin(s)*(-b*sin(t))
        dxyzduv(1,5)=cos(s)*(-b*cos(t))
c
        dxyzduv(2,3)=-sin(s)*(a+b*cos(t))
        dxyzduv(2,4)=cos(s)*(-b*sin(t))
        dxyzduv(2,5)=sin(s)*(-b*cos(t))
c
        dxyzduv(3,3)=0
        dxyzduv(3,4)=0
        dxyzduv(3,5)=-b*sin(t)
c
c       ... the following assumes that the triangulation of
c       [0..2*pi] x [0..2*pi] contains exactly two patches
c
        done=1
        pi=4*atan(done)
c
        do i=1,2
        dxyzduv(1,i)=dxyzduv(1,i)*2*pi
        dxyzduv(2,i)=dxyzduv(2,i)*2*pi
        dxyzduv(3,i)=dxyzduv(3,i)*2*pi
        enddo        
c
        do i=3,5
        dxyzduv(1,i)=dxyzduv(1,i)*(2*pi)**2
        dxyzduv(2,i)=dxyzduv(2,i)*(2*pi)**2
        dxyzduv(3,i)=dxyzduv(3,i)*(2*pi)**2
        enddo        
c
        return
c
        sx=1.2d0
        sy=1.1d0
        sz=1.0d0
c
        xyz(1)=xyz(1)*sx
        xyz(2)=xyz(2)*sy
        xyz(3)=xyz(3)*sz
c
        do i=1,2
        dxyzduv(1,i)=dxyzduv(1,i)*sx
        dxyzduv(2,i)=dxyzduv(2,i)*sy
        dxyzduv(3,i)=dxyzduv(3,i)*sz
        enddo        
c
        do i=3,5
        dxyzduv(1,i)=dxyzduv(1,i)*sx
        dxyzduv(2,i)=dxyzduv(2,i)*sy
        dxyzduv(3,i)=dxyzduv(3,i)*sz
        enddo        
c
        return
        end
c
c
c
c
c
        subroutine wpatchpnt_ext
     $     (ipatch,u,v,xyz,dxyzduv,par1,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension xyz(3),dxyzduv(3,5)        
c
c       This subroutine maps a standard simplex triangle (0,0), (1,0), (0,1)
c       into R^3. This is just an example, one triangle only.
c
c       ... setup a flat triangle in R^3
c
        x0=1
        y0=0
        z0=0
c
        x1=0
        y1=1
        z1=0
c
        x2=0
        y2=0
        z2=1
c
c       ... process the geometry, return the point location and the
c       derivatives with respect to u and v
c
        x=x0+u*(x1-x0)+v*(x2-x0)
        y=y0+u*(y1-y0)+v*(y2-y0)
        z=z0+u*(z1-z0)+v*(z2-z0)
c
        xyz(1)=x
        xyz(2)=y
        xyz(3)=z
c
        dxyzduv(1,1)=x1-x0
        dxyzduv(1,2)=x2-x0
c
        dxyzduv(2,1)=y1-y0
        dxyzduv(2,2)=y2-y0
c
        dxyzduv(3,1)=z1-z0
        dxyzduv(3,2)=z2-z0
c
c
c       ... second derivatives with respect to u and v
c
        dxyzduv(1,3)=0
        dxyzduv(1,4)=0
        dxyzduv(1,5)=0
c
        dxyzduv(2,3)=0
        dxyzduv(2,4)=0
        dxyzduv(2,5)=0
c
        dxyzduv(3,3)=0
        dxyzduv(3,4)=0
        dxyzduv(3,5)=0
c
        return
        end
c
c
c
c
c
        subroutine patchgeo_ext(
     $     patchpnt,ipatch,u,v,par1,par2,par3,par4,
     $     xyz,dxyzduv,ds,xyznorm,xyztang1,xyztang2,crvm,crvg)
        implicit real *8 (a-h,o-z)
        external patchpnt
        dimension xyz(3),dxyzduv(3,5)
        dimension xyznorm(3),xyztang1(3),xyztang2(3)
        dimension g(2,2),b(2,2),sw(2,2)
c
c       This subroutine maps the standard simplex patch ipatch 
c       at the location (u,v) and returns the location of the point in R^3, 
c       the derivatives with respect to parametrization, the area element,
c       the normal, and two standard orthonormal tangent vectors,
c       plus the mean and gaussian curvatures
c
c       Input parameters:
c
c       patchpnt: external: must be of the form
c              patchpnt(ipatch,u,v,xyz,dxyzduv,par1,par2,par3,par4)
c       ipatch: integer: the index of the patch
c       u,v: real *8: parametrization of the standard simplex triangle
c       par1,par2,par3,par4: extra parameters
c
c       Output parameters:
c       
c       xyz: real*8(3): the location of the point in R^3
c       dxyzduv: real*8(3,5): the derivatives of coordinate functions with
c             respect to (u,v) parameterization
c       ds: real *8: the area element
c       xyznorm: real*8(3): the normal vector
c       xyztang1: real*8(3): the first standard tangent vector
c       xyztang2: real*8(3): the second standard tangent vector
c       crvm: real *8: mean curvature at location xyz
c       crvg: real *8: gaussian curvature at location xyz
c
c       ... retrieve a point on the patch
c
        call patchpnt(ipatch,u,v,xyz,dxyzduv,par1,par2,par3,par4)
c
ccc        call prin2('xyz=*',xyz,3)
ccc        call prin2('dxyzduv=*',dxyzduv(1,1),3)
ccc        call prin2('dxyzduv=*',dxyzduv(1,2),3)
c
c       ... find the normal
c
        call cross_prod3d(dxyzduv(1,1),dxyzduv(1,2),xyznorm)
c
ccc        call prin2('after cross_prod3: xyznorm=*',xyznorm,3)
c
c       ... find the area element
c
        d=sqrt(xyznorm(1)**2+xyznorm(2)**2+xyznorm(3)**2)
c
        ds=d
c
ccc        call prin2('ds=*',ds,1)
c
        xyznorm(1)=xyznorm(1)/ds
        xyznorm(2)=xyznorm(2)/ds
        xyznorm(3)=xyznorm(3)/ds
c
c
c       ... construct the first tangent vector
c       
        xyztang1(1)=dxyzduv(1,1)
        xyztang1(2)=dxyzduv(2,1)
        xyztang1(3)=dxyzduv(3,1)
c
        d=sqrt(xyztang1(1)**2+xyztang1(2)**2+xyztang1(3)**2)
c
        xyztang1(1)=xyztang1(1)/d
        xyztang1(2)=xyztang1(2)/d
        xyztang1(3)=xyztang1(3)/d        
c
c
c       ... construct the second tangent vector
c       
        call cross_prod3d(xyznorm,xyztang1,xyztang2)
c
ccc        call prin2('xyznorm=*',xyznorm,3)
ccc        call prin2('xyztang1=*',xyztang1,3)
ccc        call prin2('xyztang2=*',xyztang2,3)
c
c
c       ... construct the coefficients of the first fundamental form
c
        call dot_prod3d(dxyzduv(1,1),dxyzduv(1,1),g11)
        call dot_prod3d(dxyzduv(1,1),dxyzduv(1,2),g12)
        call dot_prod3d(dxyzduv(1,2),dxyzduv(1,1),g21)
        call dot_prod3d(dxyzduv(1,2),dxyzduv(1,2),g22)
c
        g(1,1)=g11
        g(1,2)=g12
        g(2,1)=g21
        g(2,2)=g22
c
ccc        call prin2('g=*',g,4)
c
        detg=g11*g22-g12*g21
c
ccc        call prin2('sqrt(det|g|)=*',sqrt(detg),1)
ccc        call prin2('and ds=*',ds,1)
c
c
c       ... construct the coefficients of the second fundamental form
c
        call dot_prod3d(dxyzduv(1,3),xyznorm,b11)
        call dot_prod3d(dxyzduv(1,4),xyznorm,b12)
        call dot_prod3d(dxyzduv(1,4),xyznorm,b21)
        call dot_prod3d(dxyzduv(1,5),xyznorm,b22)
c
        b(1,1)=b11
        b(1,2)=b12
        b(2,1)=b21
        b(2,2)=b22
c
ccc        call prin2('b=*',b,4)
c
        detb=b11*b22-b12*b21
c
cc        call prin2('det|g|=*',detg,1)
cc        call prin2('det|b|=*',detb,1)
c
        crvg=detb/detg
c
cc        call prin2('gaussian curvature: detb/detg=*',crvg,1)
c
        crvm=(b22*g11-b12*g21-b21*g12+b11*g22)/2/detg
c
cc        call prin2('mean curvature=*',crvm,1)
c
        crv1=crvm+sqrt(abs((crvm)**2-crvg))
        crv2=crvm-sqrt(abs((crvm)**2-crvg))
c
cc        call prin2('principal curvature k1=*',crv1,1)
cc        call prin2('principal curvature k2=*',crv2,1)
c
        return
c
c       Weingarter equations for dndu and dndv (shape operator)
c
c       n_u = S( x_u )
c       n_v = S( x_v )
c
        sw(1,1)=(b12*g12-b11*g22)/detg
        sw(1,2)=(b11*g21-b12*g11)/detg
        sw(2,1)=(b22*g12-b21*g22)/detg
        sw(2,2)=(b21*g21-b22*g11)/detg
c
ccc        call prin2('sw=*',sw,4)
c
        return
        end
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       the discretization routines in R^3
c
c       Geometry 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       
c
c       
        subroutine patchgeo(patchpnt,ipatch,u,v,par1,par2,par3,par4,
     $     xyz,dxyzduv,ds,xyznorm,xyztang1,xyztang2)
        implicit real *8 (a-h,o-z)
        external patchpnt
        dimension xyz(3),dxyzduv(3,2)
        dimension xyznorm(3),xyztang1(3),xyztang2(3)
c
c       This subroutine maps the standard simplex patch ipatch 
c       at the location (u,v) and returns the location of the point in R^3, 
c       the derivatives with respect to parametrization, the area element,
c       the normal, and two standard orthonormal tangent vectors.
c
c       Input parameters:
c
c       patchpnt: external: must be of the form
c                 patchpnt(ipatch,u,v,xyz,dxyzduv,par1,par2,par3,par4)
c       ipatch: integer: the index of the patch
c       u,v: real *8: parametrization of the standard simplex triangle
c       par1,par2,par3,par4: extra parameters
c
c       Output parameters:
c       
c       xyz: real*8(3): the location of the point in R^3
c       dxyzduv: real*8(3,2): the derivatives of coordinate functions with
c             respect to (u,v) parameterization
c       ds: real *8: the area element
c       xyznorm: real*8(3): the normal vector
c       xyztang1: real*8(3): the first standard tangent vector
c       xyztang2: real*8(3): the second standard tangent vector
c
c
c       ... retrieve a point on the patch
c
        call patchpnt(ipatch,u,v,xyz,dxyzduv,par1,par2,par3,par4)
c
ccc        call prin2('xyz=*',xyz,3)
ccc        call prin2('dxyzduv=*',dxyzduv(1,1),3)
ccc        call prin2('dxyzduv=*',dxyzduv(1,2),3)
c
c       ... find the normal
c
        call cross_prod3d(dxyzduv(1,1),dxyzduv(1,2),xyznorm)
c
ccc        call prin2('after cross_prod3: xyznorm=*',xyznorm,3)
c
c       ... find the area element
c
        d=sqrt(xyznorm(1)**2+xyznorm(2)**2+xyznorm(3)**2)
c
        ds=d
c
ccc        call prin2('ds=*',ds,1)
c
        xyznorm(1)=xyznorm(1)/ds
        xyznorm(2)=xyznorm(2)/ds
        xyznorm(3)=xyznorm(3)/ds
c
c
c       ... construct the first tangent vector
c       
        xyztang1(1)=dxyzduv(1,1)
        xyztang1(2)=dxyzduv(2,1)
        xyztang1(3)=dxyzduv(3,1)
c
        d=sqrt(xyztang1(1)**2+xyztang1(2)**2+xyztang1(3)**2)
c
        xyztang1(1)=xyztang1(1)/d
        xyztang1(2)=xyztang1(2)/d
        xyztang1(3)=xyztang1(3)/d        
c
c
c       ... construct the second tangent vector
c       
        call cross_prod3d(xyznorm,xyztang1,xyztang2)
c
ccc        call prin2('xyznorm=*',xyznorm,3)
ccc        call prin2('xyztang1=*',xyztang1,3)
ccc        call prin2('xyztang2=*',xyztang2,3)
c
        return
        end
c
c
c
c
c
        subroutine patchinfo(xyz,xyznorm,xyztang1,xyztang2,xyzinfo)
        implicit real *8 (a-h,o-z)
        dimension xyz(3),xyznorm(3),xyztang1(3),xyztang2(3)
        dimension xyzinfo(12)
c
c       This subroutine compresses the geometry information into
c       the standard linear array.
c
c       Input parameters:
c       
c       xyz: real*8(3): the location of the point in R^3
c       xyznorm: real*8(3): the normal vector
c       xyztang1: real*8(3): the first standard tangent vector
c       xyztang2: real*8(3): the second standard tangent vector
c
c       Output parameters:
c
c       xyzinfo: real*8(12): the compressed geometry structure
c
c       ... compress geometry information
c
        xyzinfo(1)=xyz(1)
        xyzinfo(2)=xyz(2)
        xyzinfo(3)=xyz(3)
c
        xyzinfo(4)=xyznorm(1)
        xyzinfo(5)=xyznorm(2)
        xyzinfo(6)=xyznorm(3)
c
        xyzinfo(7)=xyztang1(1)
        xyzinfo(8)=xyztang1(2)
        xyzinfo(9)=xyztang1(3)
c
        xyzinfo(10)=xyztang2(1)
        xyzinfo(11)=xyztang2(2)
        xyzinfo(12)=xyztang2(3)        
c
        return
        end
c
c
c
c
c
        subroutine patchallpnts(npatches,patchpnt,par1,par2,par3,par4,
     $     npols,us,vs,ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts)
        implicit real *8 (a-h,o-z)
        external patchpnt
        dimension us(1),vs(1)
        dimension xyz(3),dxyzduv(3,2)
        dimension xyznorm(3),xyztang1(3),xyztang2(3)
c
        dimension ixyzs(2,1),xyzs(3,1),xyznorms(3,1)
        dimension xyztang1s(3,1),xyztang2s(3,1)
c
c       This subroutine return all points, normals and tangents from
c       geometry descriptor
c
c       Input parameters:
c
c       npatches: integer: the number of patches
c       patchpnt: external: must be of the form
c                 patchpnt(ipatch,u,v,xyz,dxyzduv,par1,par2,par3,par4)
c       par1,par2,par3,par4: extra parameters
c       npols: integer: the total number of polynomials for each patch
c       us: real *8(1): local u-discretization points for each patch
c       vs: real *8(1): local v-discretization points for each patch
c
c       Output parameters:
c       
c       ixyzs: integer(2,npts): the index array for each patch
c       xyzs: real*8(3,npts): the location of the points in R^3
c       xyznorms: real*8(3,npts): the normal vectors
c       xyztang1s: real*8(3,npts): the first standard tangent vectors
c       xyztang2s: real*8(3,npts): the second standard tangent vectors
c       npts: integer: the total number of points in discretization
c
c
        npts=0
        kk=1
c
        do 1400 ipatch=1,npatches
        ixyzs(1,ipatch)=kk
        ixyzs(2,ipatch)=npols
        do 1200 i=1,npols
c
        u=us(i)
        v=vs(i)
c
        call patchgeo(patchpnt,ipatch,u,v,
     $     par1,par2,par3,par4,
     $     xyz,dxyzduv,ds,xyznorm,xyztang1,xyztang2)
c
        xyzs(1,kk)=xyz(1)
        xyzs(2,kk)=xyz(2)
        xyzs(3,kk)=xyz(3)
c
        xyznorms(1,kk)=xyznorm(1)
        xyznorms(2,kk)=xyznorm(2)
        xyznorms(3,kk)=xyznorm(3)
c
        xyztang1s(1,kk)=xyztang1(1)
        xyztang1s(2,kk)=xyztang1(2)
        xyztang1s(3,kk)=xyztang1(3)
c
        xyztang2s(1,kk)=xyztang2(1)
        xyztang2s(2,kk)=xyztang2(2)
        xyztang2s(3,kk)=xyztang2(3)
c
        npts=npts+1
        kk=kk+1
c
 1200   continue
 1400   continue
c
c
        return
        end
c
c
c
c
c
        subroutine patchallwhts(npatches,patchpnt,par1,par2,par3,par4,
     $     npols,us,vs,ws,whts,npts)
        implicit real *8 (a-h,o-z)
        external patchpnt
        dimension us(1),vs(1),ws(1)
        dimension xyz(3),dxyzduv(3,2)
        dimension xyznorm(3),xyztang1(3),xyztang2(3)
c
c       This subroutine return the discretization weights
c
c       Input parameters:
c
c       npatches: integer: the number of patches
c       patchpnt: external: must be of the form
c                 patchpnt(ipatch,u,v,xyz,dxyzduv,par1,par2,par3,par4)
c       par1,par2,par3,par4: extra parameters
c       npols: integer: the total number of polynomials for each patch
c       us: real *8(npols): local u-discretization points for each patch
c       vs: real *8(npols): local v-discretization points for each patch
c       ws: real *8(npols): local integration weights for each patch
c
c       Output parameters:
c       
c       whts: real*8(npts): the discretization weights
c       npts: integer: the total number of points in discretization
c
c
        dimension whts(1)
c
        npts=0
        kk=1
c
        do 1400 ipatch=1,npatches
        do 1200 i=1,npols
c
        u=us(i)
        v=vs(i)
c
        call patchgeo(patchpnt,ipatch,u,v,
     $     par1,par2,par3,par4,
     $     xyz,dxyzduv,ds,xyznorm,xyztang1,xyztang2)
c
        whts(kk)=ds*ws(i)
c
        npts=npts+1
        kk=kk+1
c
 1200   continue
 1400   continue
c
c
        return
        end
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       the discretization routines in R^3
c
c       Matrix generation
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       
c
c       
        subroutine patchmatc(npatches,patchpnt,par1,par2,par3,par4,
     $     norder,npols,us,vs,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     interact,par5,par6,par7,par8,
     $     cmatr,w,lw,lused,ier)
        implicit real *8 (a-h,o-z)
c
c
c       This subroutine generates the (complex) matrix of interactions
c       on a user-defined geometry, described by patchpnt subroutine.
c
c       Input parameters:
c
c       npatches: integer: the number of patches
c       patchpnt: external: must be of the form
c                 patchpnt(ipatch,u,v,xyz,dxyzduv,par1,par2,par3,par4)
c       par1,par2,par3,par4: extra parameters
c       norder: integer: the order of interpolation
c       npols: integer: the total number of polynomials for each patch
c       us: real *8(1): local u-discretization points for each patch
c       vs: real *8(1): local v-discretization points for each patch
c       umatr: real *8(npols,npols): first interpolation matrix
c       vmatr: real *8(npols,npols): second interpolation matrix
c       interact: external: must be of the form
c            interact(srcinfo,targinfo,cout,par5,par6,par7,par8)
c
c       Output parameters:
c       
c       cmatr: complex*16(npts,npts): the interaction matrix
c
c       lused: integer: the total 
c       ier: integer: the error code
c
c       Work arrays:
c
c       w: real *8(1): work array, must be at least 2*10000 real *8 elements
c       lw: integer: the length of the work array
c
c
        external patchpnt,interact
        dimension us(1),vs(1)
        dimension xyz(3),dxyzduv(3,2)
        dimension xyznorm(3),xyztang1(3),xyztang2(3)
c
        dimension ixyzs(2,1),xyzs(3,1),xyznorms(3,1)
        dimension xyztang1s(3,1),xyztang2s(3,1)
c
        dimension w(1)
        complex *16 cmatr(npts,npts)
c
        ier=0
c
c       ... allocate work arrays
c        
c       ... max 100 points per patch
c
        nmax=100
c
        itmatr=1
        ltmatr=2*nmax*nmax
c
        lused7=ltmatr
c
        call patchmatc0(npatches,patchpnt,par1,par2,par3,par4,
     $     norder,npols,us,vs,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     interact,par5,par6,par7,par8,
     $     cmatr,w(itmatr),w(1+lused7),lw-lused7,lused8,ier)        
c
        lused=lused7+lused8
c
        return
        end
c
c
c
c
c
        subroutine patchmatc0(npatches,patchpnt,par1,par2,par3,par4,
     $     norder,npols,us,vs,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     interact,par5,par6,par7,par8,
     $     cmatr,tmatr,w,lw,lused,ier)
        implicit real *8 (a-h,o-z)
        external patchpnt,interact
        dimension us(1),vs(1)
        dimension xyz(3),dxyzduv(3,2)
        dimension xyznorm(3),xyztang1(3),xyztang2(3)
c
        dimension ixyzs(2,1),xyzs(3,1),xyznorms(3,1)
        dimension xyztang1s(3,1),xyztang2s(3,1)
c
        dimension w(1)
        real *8, allocatable :: w_omp(:)
        complex *16 tmatr(1)
        complex *16 cmatr(npts,npts)
        complex *16 tmatr_omp(10000)
c
c
c       ... construct the off-diagonal blocks
c
ccc        do 1400 k=1,1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ii,ipols,jj,jpols,tmatr_omp,lused,ier,w_omp) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 1400 j=1,npatches
c
        allocate(w_omp(2000000))
c
ccc        call prinf('od,j=*',j,1)
c
        do 1200 i=1,npatches
        if ( i .eq. j ) goto 1200
c        
c       ... (i,j), i index - target, j index - source
c
        ii=ixyzs(1,i)
        jj=ixyzs(1,j)
        ipols=ixyzs(2,i)
        jpols=ixyzs(2,j)
c
cc        call prinf('i=*',i,1)
cc        call prinf('j=*',j,1)
c
        call patchmatc_od(i,ipols,j,jpols,
     $     npatches,patchpnt,par1,par2,par3,par4,
     $     norder,npols,us,vs,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     interact,par5,par6,par7,par8,
     $     tmatr_omp,w_omp,lw,lused,ier)
c
        call patchsubcpy(npts,cmatr,tmatr_omp,ii,ipols,jj,jpols)
c        
 1200   continue
        deallocate(w_omp)
 1400   continue        
C$OMP END PARALLEL DO
c
c
ccc        return
c
c       
c       ... construct the diagonal (self interaction) blocks
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,ii,ipols,jj,jpols,tmatr_omp,lused,ier,w_omp) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 2200 i=1,npatches
c
        allocate(w_omp(2000000))
c
c       ... (i,i) 
c
        ii=ixyzs(1,i)
        ipols=ixyzs(2,i)
c
ccc        call prinf('dd,i=*',i,1)
c
        call patchmatc_dd(i,ipols,i,ipols,
     $     npatches,patchpnt,par1,par2,par3,par4,
     $     norder,npols,us,vs,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     interact,par5,par6,par7,par8,
     $     tmatr_omp,w_omp,lw,lused,ier)
c
        call patchsubcpy(npts,cmatr,tmatr_omp,ii,ipols,ii,ipols)
c              
        deallocate(w_omp)
 2200   continue
C$OMP END PARALLEL DO
c
        return
        end
c
c
c
c
c
        subroutine patchsubcpy(npts,cmatr,tmatr,ii,ipols,jj,jpols)
        implicit real *8 (a-h,o-z)
c
        complex *16 cmatr(npts,npts)
        complex *16 tmatr(ipols,jpols)
c
        do 1400 j=1,jpols
        do 1200 i=1,ipols
c
        cmatr(ii+i-1,jj+j-1)=tmatr(i,j)
 1200   continue
 1400   continue
c       
        return
        end
c
c
c
c
c
        subroutine patchmatc_od(ipatch,ipols,jpatch,jpols,
     $     npatches,patchpnt,par1,par2,par3,par4,
     $     norder,npols,us,vs,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     interact,par5,par6,par7,par8,
     $     tmatr,w,lw,lused,ier)
        implicit real *8 (a-h,o-z)
c
c       ... generate the off-diagonal block of interaction matrix
c
        external patchpnt,interact
        dimension us(1),vs(1)
        dimension xyz(3),dxyzduv(3,2)
        dimension xyznorm(3),xyztang1(3),xyztang2(3)
c
        dimension ixyzs(2,1),xyzs(3,1),xyznorms(3,1)
        dimension xyztang1s(3,1),xyztang2s(3,1)
c
        dimension targinfo(12),info(20)
        dimension xpar1(20),xpar2(20)
c
        dimension w(1)
        complex *16 tmatr(ipols,jpols)
        complex *16 cvals(1000),coefs(1000)
c       
        dimension vert1(2),vert2(2),vert3(2)
        external patchfun3
c
ccc        write(*,*) '.', ipatch, ipols, jpatch, jpols
c
        ii=ixyzs(1,ipatch)
        jj=ixyzs(1,jpatch)
c
ccc        write(*,*) '.', ii, jj
c
c       ... construct one off-diagonal block via collocation
c       
c       ... (i,j), i index - target, j index - source
c
        if( ipols .ne. npols ) then
        write(*,*) 'ipols .ne. npols'
        endif
        if( jpols .ne. npols ) then
        write(*,*) 'jpols .ne. npols'
        endif
c        
        do 1400 i=1,npols
c
c       ... on j-th triangle integrate all basis functions multiplied 
c       by interaction function at the target point us(i),vs(i)
c
c       ... first, initialize function to be integrated
c
        xpar1(1)=norder
        xpar1(2)=npols
        xpar1(3)=jpatch
c
        call patchinfo(xyzs(1,ii+i-1),xyznorms(1,ii+i-1),
     $     xyztang1s(1,ii+i-1),xyztang2s(1,ii+i-1),targinfo)
c
cccc        call prin2('targinfo=*',targinfo,12)
        do j=1,12
        xpar2(j)=targinfo(j)
        enddo
c
c       ... then, call adaptive gaussian integration routine 
c       
        m=6
ccc        eps=1d-12
ccc        eps=1d-10
ccc        eps=1d-8
        eps=1d-6
c        eps=1d-4
c
c        m=3
c        eps=1d-3
c
        m=12
        eps=1d-12
c
        m=16
        eps=1d-16
c
        iquadtype=1
        nrec=20
c
c        iquadtype=2
c
        vert1(1)=0
        vert1(2)=0
        vert2(1)=1
        vert2(2)=0
        vert3(1)=0
        vert3(2)=1
c
        if( iquadtype .eq. 1 )        
     $     call c28triaadam(ier,vert1,vert2,vert3,patchfun3,npols,
     $     patchpnt,par1,par2,par3,par4,
     $     interact,par5,par6,par7,par8,xpar1,xpar2,
     $     m,eps,coefs,maxrec,numfunev,w,nrec,jpatch,targinfo,info)
c
        if( iquadtype .eq. 2 )
     $     call c29triaadam(ier,vert1,vert2,vert3,patchfun3,npols,
     $     patchpnt,par1,par2,par3,par4,
     $     interact,par5,par6,par7,par8,xpar1,xpar2,
     $     m,eps,coefs,maxrec,numfunev,w)
c       
        if( ier .eq. 8 ) then
        write(*,*) 'maximum recursion depth of 200 has been reached'
        write(*,*) 'abort'
        stop
        endif
c
ccc        call prinf('numfunev=*',numfunev,1)
ccc        call prinf('info=*',info,2)
ccc        call prin2('coefs=*',coefs,2*npols)
c
c
c       ... finally, convert the linear form of integral values to the
c       pointwise interation matrix, we will need umatr and vmatr for
c       this operation
c
        call patchcoefs2cvals(npols,umatr,vmatr,coefs,cvals)
c
ccc        call prin2('cvals=*',cvals,2*npols)
c
        do 1200 j=1,npols
        tmatr(i,j)=cvals(j)
 1200   continue
c
 1400   continue        
c
c
        return
        end
c
c
c
c
c
        subroutine patchmatc_dd(ipatch,ipols,jpatch,jpols,
     $     npatches,patchpnt,par1,par2,par3,par4,
     $     norder,npols,us,vs,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     interact,par5,par6,par7,par8,
     $     tmatr,w,lw,lused,ier)
        implicit real *8 (a-h,o-z)
c
c       ... generate the diagonal block of interaction matrix,
c       self-interaction
c
        external patchpnt,interact
        dimension us(1),vs(1)
        dimension xyz(3),dxyzduv(3,2)
        dimension xyznorm(3),xyztang1(3),xyztang2(3)
c
        dimension ixyzs(2,1),xyzs(3,1),xyznorms(3,1)
        dimension xyztang1s(3,1),xyztang2s(3,1)
c
        dimension targinfo(12),info(20)
        dimension xpar1(20),xpar2(20)
c
        dimension w(1)
        complex *16 tmatr(ipols,jpols)
        complex *16 cvals(1000),coefs(1000)
        complex *16 coefs1(1000),coefs2(1000),coefs3(1000)
c
        dimension xs(10000),ys(10000),ws(10000)
c       
        dimension vert1(2),vert2(2),vert3(2)
        external patchfun3
c
c
ccc        write(*,*) '.', ipatch, ipols, jpatch, jpols
c
        ii=ixyzs(1,ipatch)
        jj=ixyzs(1,jpatch)
c
ccc        write(*,*) '.', ii, jj
c
c       ... construct one off-diagonal block via collocation
c       
c       ... (i,j), i index - target, j index - source
c
        if( ipols .ne. npols ) then
        write(*,*) 'ipols .ne. npols'
        endif
        if( jpols .ne. npols ) then
        write(*,*) 'jpols .ne. npols'
        endif
c        
        do 1400 i=1,npols
c
c       ... on j-th triangle integrate all basis functions multiplied 
c       by interaction function at the target point us(i),vs(i)
c
c       ... first, initialize function to be integrated
c
        xpar1(1)=norder
        xpar1(2)=npols
        xpar1(3)=jpatch
c
        call patchinfo(xyzs(1,ii+i-1),xyznorms(1,ii+i-1),
     $     xyztang1s(1,ii+i-1),xyztang2s(1,ii+i-1),targinfo)
c
cccc        call prin2('targinfo=*',targinfo,12)
        do j=1,12
        xpar2(j)=targinfo(j)
        enddo
c
cccc        write(*,*) i,us(i),vs(i)
c
c       ... then, call adaptive gaussian integration routine 
c
        m=6
ccc        eps=1d-12
ccc        eps=1d-10
ccc        eps=1d-8
        eps=1d-6
c        eps=1d-4
c
c        m=3
c        eps=1d-3
c
        m=12
        eps=1d-12
c
        iquadtype=1
        nrec=20
c
        iquadtype=2
c
        if( norder .eq. 0 ) iquadtype=3
        if( norder .eq. 1 ) iquadtype=3
        if( norder .eq. 2 ) iquadtype=3
        if( norder .eq. 3 ) iquadtype=3
ccc        if( norder .eq. 6 ) iquadtype=3
c
        if( iquadtype .eq. 3 ) then 
c
        vert1(1)=0
        vert1(2)=0
        vert2(1)=1
        vert2(2)=0
        vert3(1)=0
        vert3(2)=1
c
        itype=2
        inode=i
        call triaselfquad
     $     (ier,itype,norder,inode,xs,ys,ws,ns,x0,y0)
c
        do j=1,npols
        coefs(j)=0
        enddo
c
        do k=1,ns
        call patchfun3(xs(k),ys(k),patchpnt,par1,par2,par3,par4,
     $     interact,par5,par6,par7,par8,xpar1,xpar2,cvals)
        do j=1,npols
        coefs(j)=coefs(j)+ws(k)*cvals(j)
        enddo
        enddo

        endif
c
c
        if( iquadtype .eq. 1 .or. iquadtype .eq. 2 ) then 
c
        vert1(1)=us(i)
        vert1(2)=vs(i)
        vert2(1)=0
        vert2(2)=0
        vert3(1)=0
        vert3(2)=1
c
        if( iquadtype .eq. 1 )
     $     call c28triaadam(ier,vert1,vert2,vert3,patchfun3,npols,
     $     patchpnt,par1,par2,par3,par4,
     $     interact,par5,par6,par7,par8,xpar1,xpar2,
     $     m,eps,coefs1,maxrec,numfunev,w,nrec,jpatch,targinfo,info)
c
        if( iquadtype .eq. 2 )
     $     call c29triaadam(ier,vert1,vert2,vert3,patchfun3,npols,
     $     patchpnt,par1,par2,par3,par4,
     $     interact,par5,par6,par7,par8,xpar1,xpar2,
     $     m,eps,coefs1,maxrec,numfunev,w)
c
ccc        call prinf('numfunev=*',numfunev,1)
ccc        call prinf('info=*',info,2)
c       
        if( ier .eq. 8 ) then
        write(*,*) 'maximum recursion depth of 200 has been reached'
        write(*,*) 'abort'
        stop
        endif
c
        vert1(1)=us(i)
        vert1(2)=vs(i)
        vert2(1)=0
        vert2(2)=1
        vert3(1)=1
        vert3(2)=0
c
c
        if( iquadtype .eq. 1 )
     $     call c28triaadam(ier,vert1,vert2,vert3,patchfun3,npols,
     $     patchpnt,par1,par2,par3,par4,
     $     interact,par5,par6,par7,par8,xpar1,xpar2,
     $     m,eps,coefs2,maxrec,numfunev,w,nrec,jpatch,targinfo,info)
c
        if( iquadtype .eq. 2 )
     $     call c29triaadam(ier,vert1,vert2,vert3,patchfun3,npols,
     $     patchpnt,par1,par2,par3,par4,
     $     interact,par5,par6,par7,par8,xpar1,xpar2,
     $     m,eps,coefs2,maxrec,numfunev,w)
c
ccc        call prinf('numfunev=*',numfunev,1)
ccc        call prinf('info=*',info,2)
c
        if( ier .eq. 8 ) then
        write(*,*) 'maximum recursion depth of 200 has been reached'
        write(*,*) 'abort'
        stop
        endif
c
        vert1(1)=us(i)
        vert1(2)=vs(i)
        vert2(1)=1
        vert2(2)=0
        vert3(1)=0
        vert3(2)=0
c
        if( iquadtype .eq. 1 )
     $     call c28triaadam(ier,vert1,vert2,vert3,patchfun3,npols,
     $     patchpnt,par1,par2,par3,par4,
     $     interact,par5,par6,par7,par8,xpar1,xpar2,
     $     m,eps,coefs3,maxrec,numfunev,w,nrec,jpatch,targinfo,info)
c
        if( iquadtype .eq. 2 )
     $     call c29triaadam(ier,vert1,vert2,vert3,patchfun3,npols,
     $     patchpnt,par1,par2,par3,par4,
     $     interact,par5,par6,par7,par8,xpar1,xpar2,
     $     m,eps,coefs3,maxrec,numfunev,w)
c
ccc        call prinf('numfunev=*',numfunev,1)
ccc        call prinf('info=*',info,2)
c
        if( ier .eq. 8 ) then
        write(*,*) 'maximum recursion depth of 200 has been reached'
        write(*,*) 'abort'
        stop
        endif
c
        do j=1,npols
        coefs(j)=coefs1(j)+coefs2(j)+coefs3(j)
        enddo
c
        endif
c

c
ccc        call prin2('coefs=*',coefs,2*npols)
c
c       ... finally, convert the linear form of integral values to the
c       pointwise interation matrix, we will need umatr and vmatr for
c       this operation
c
        call patchcoefs2cvals(npols,umatr,vmatr,coefs,cvals)
c
ccc        call prin2('cvals=*',cvals,2*npols)
c
        do 1200 j=1,npols
        tmatr(i,j)=cvals(j)
 1200   continue
c
 1400   continue        
c
c
        return
        end
c
c
c
c
c
        subroutine patchcoefs2cvals(npols,umatr,vmatr,coefs,cvals)
        implicit real *8 (a-h,o-z)
        dimension umatr(npols,npols),vmatr(npols,npols)
        complex *16 cvals(1000),coefs(1000),cd
c
        do 1400 i=1,npols
        cd=0
        do 1200 j=1,npols
        cd=cd+umatr(j,i)*coefs(j)
 1200   continue
        cvals(i)=cd
 1400   continue
c
        return
        end
c
c
c
c
c
        subroutine patchfun3(u,v,patchpnt,par1,par2,par3,par4,
     $     interact,par5,par6,par7,par8,xpar1,xpar2,cvals)
        implicit real *8 (a-h,o-z)
c
c       ... provides functions to be integrated 
c      must be initialized via arrays xpar1 and xpar2
c
c
        external patchpnt,interact
        complex *16 cvals(1),cout
        dimension pols(1000)
        dimension xyz(3),dxyzduv(3,2)
        dimension xyznorm(3),xyztang1(3),xyztang2(3)
        dimension srcinfo(12),targinfo(12)
        dimension xpar1(20),xpar2(20)
c
c       ... initialize the function evaluator
c
        norder=xpar1(1)
        npols=xpar1(2)
        jpatch=xpar1(3)
c
        do i=1,12
        targinfo(i)=xpar2(i)
        enddo
c
c
ccc        ntimes=ntimes+1
ccc        if( mod(ntimes,1000) .eq. 0 ) call prinf('ntimes=*',ntimes,1)
c
        call ortho2sipols(u,v,norder,pols)
c
ccc        call prin2('pols=*',pols,npols)
c
        call patchgeo(patchpnt,jpatch,u,v,
     $     par1,par2,par3,par4,
     $     xyz,dxyzduv,ds,xyznorm,xyztang1,xyztang2)
c
ccc        call prin2('xyz=*',xyz,3)
c
        call patchinfo(xyz,xyznorm,xyztang1,xyztang2,srcinfo)
c
        call interact(srcinfo,targinfo,cout,par5,par6,par7,par8)
c
        do 1200 i=1,npols
        cvals(i)=cout*pols(i)*ds
ccc        cvals(i)=pols(i)*ds
 1200   continue
c
ccc        call prin2('cvals=*',cvals,2*npols)
c
        return
        end
c
c
c
c
c
