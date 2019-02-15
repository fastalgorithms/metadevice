c
c
c
c
c
        subroutine em3multa(a,par1,par2,x,y,n)
        implicit real *8 (a-h,o-z)
        complex *16 a(n,n),x(n),y(n)
c
ccc$OMP PARALLEL DO DEFAULT(SHARED)
cc        do i=1,n
cc          y(i)=0
cc          do j=1,n
cc            y(i)=y(i)+a(i,j)*x(j)
cc          enddo
cc        enddo
ccc$OMP END PARALLEL DO
c
c       call the BLAS2 routine
c
        call zmatvec(n, n, a, x, y)        
c
        return
        end
c
c
c
c
c
        subroutine em3multb(a,par1,par2,x,y,n)
        implicit real *8 (a-h,o-z)
        complex *16 a(n,n),x(n),y(n)
c
        do 1200 i=1,n
        y(i)=0
        do 1100 j=1,n
        y(i)=y(i)+conjg(a(j,i))*x(j)
 1100   continue
 1200   continue
c
        return
        end
c
c
c
c
c

        subroutine em3debug(npts,cmatr,whts)
        implicit real *8 (a-h,o-z)
        complex *16 cmatr(npts,npts)
        dimension whts(1)
c
        do i=1,3
        call prinf('i=*',i,1)
        call prin2('cmatr=*',cmatr(i,i)*whts(i),2)
        enddo
c
        stop
        return
        end
c
c
c
c
c
        subroutine em3subcpy(npts,cmatr,ipts2,jpts2,cmatr2,ii,jj)
        implicit real *8 (a-h,o-z)
c
        complex *16 cmatr(npts,npts)
        complex *16 cmatr2(ipts2,jpts2)
c
        do 1400 j=1,jpts2
        do 1200 i=1,ipts2
c
        cmatr(ii+i-1,jj+j-1)=cmatr2(i,j)
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
        subroutine em3submulcpy(npts,cmatr,ipts2,jpts2,cmatr2,ii,jj,cd)
        implicit real *8 (a-h,o-z)
c
        complex *16 cmatr(npts,npts)
        complex *16 cmatr2(ipts2,jpts2),cd
c
        do 1400 j=1,jpts2
        do 1200 i=1,ipts2
c
        cmatr(ii+i-1,jj+j-1)=cmatr2(i,j)*cd
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
        subroutine em3submuladd(npts,cmatr,ipts2,jpts2,cmatr2,ii,jj,cd)
        implicit real *8 (a-h,o-z)
c
        complex *16 cmatr(npts,npts)
        complex *16 cmatr2(ipts2,jpts2),cd
c
        do 1400 j=1,jpts2
        do 1200 i=1,ipts2
c
        cmatr(ii+i-1,jj+j-1)=cmatr(ii+i-1,jj+j-1)+cmatr2(i,j)*cd
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
        subroutine em3submul(n,m,cmatr,cd)
        implicit real *8 (a-h,o-z)
c
        complex *16 cmatr(n,m),cd
c
        do 1400 j=1,m
        do 1200 i=1,n
c
        cmatr(i,j)=cmatr(i,j)*cd
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
        subroutine em3subadd(npts,cmatr,ipts2,jpts2,cmatr2,ii,jj)
        implicit real *8 (a-h,o-z)
c
        complex *16 cmatr(npts,npts)
        complex *16 cmatr2(ipts2,jpts2)
c
        do 1400 j=1,jpts2
        do 1200 i=1,ipts2
c
        cmatr(ii+i-1,jj+j-1)=cmatr(ii+i-1,jj+j-1)+cmatr2(i,j)
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
        subroutine patchdiag(cmatr,npts,d)
        implicit real *8 (a-h,o-z)
        complex *16 cmatr(npts,npts)
c
        do 1200 i=1,npts
        cmatr(i,i)=cmatr(i,i)+d
ccc        cmatr(i,i)=cmatr(i,i)-d
 1200   continue
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
        subroutine em3getrhs(rk,source,xyzinfo,rhs,npts)
        implicit real *8 (a-h,o-z)
        dimension source(3),xyzinfo(12,1)
        complex *16 rhs(npts,4),rk,ima
        complex *16 evec(3),hvec(3),cjvec(3)
        complex *16 cvec(3)
        dimension xyz(3),xyznorm(3),xyztang1(3),xyztang2(3)
        data ima/(0.0d0,1.0d0)/
c
        cjvec(1)=1
        cjvec(2)=1
        cjvec(3)=1
c       
        do 1200 i=1,npts
c
        xyz(1)=xyzinfo(1,i)-source(1)
        xyz(2)=xyzinfo(2,i)-source(2)
        xyz(3)=xyzinfo(3,i)-source(3)
c
        call dipole3e(rk,xyz,cjvec,evec,hvec)
c
        xyznorm(1)=xyzinfo(4,i)
        xyznorm(2)=xyzinfo(5,i)
        xyznorm(3)=xyzinfo(6,i)
c
ccc        call cross_prod3d_dcc(xyznorm,hvec,cvec)
c
        cvec(1)=hvec(1)
        cvec(2)=hvec(2)
        cvec(3)=hvec(3)
c
        xyztang1(1)=xyzinfo(7,i)
        xyztang1(2)=xyzinfo(8,i)
        xyztang1(3)=xyzinfo(9,i)
        xyztang2(1)=xyzinfo(10,i)
        xyztang2(2)=xyzinfo(11,i)
        xyztang2(3)=xyzinfo(12,i)
        call dot_prod3d_dcc(xyztang1,cvec,rhs(i,1))
        call dot_prod3d_dcc(xyztang2,cvec,rhs(i,2))
c
ccc        call cross_prod3d_dcc(xyznorm,evec,cvec)
c
        cvec(1)=evec(1)
        cvec(2)=evec(2)
        cvec(3)=evec(3)
c
        xyztang1(1)=xyzinfo(7,i)
        xyztang1(2)=xyzinfo(8,i)
        xyztang1(3)=xyzinfo(9,i)
        xyztang2(1)=xyzinfo(10,i)
        xyztang2(2)=xyzinfo(11,i)
        xyztang2(3)=xyzinfo(12,i)
        call dot_prod3d_dcc(xyztang1,cvec,rhs(i,3))
        call dot_prod3d_dcc(xyztang2,cvec,rhs(i,4))
 1200   continue
c        
        return
        end
c
c
c
c
c
        subroutine em3getrhs2(rk,ceps,cmu,source,xyzinfo,rhs,npts)
        implicit real *8 (a-h,o-z)
        dimension source(3),xyzinfo(12,1)
        complex *16 rhs(npts,4),rk,ima,ceps,cmu
        complex *16 evec(3),hvec(3),cjvec(3)
        complex *16 cvec(3)
        dimension xyz(3),xyznorm(3),xyztang1(3),xyztang2(3)
        data ima/(0.0d0,1.0d0)/
c
        cjvec(1)=1
        cjvec(2)=1
        cjvec(3)=1
c       
        do 1200 i=1,npts
c
        xyz(1)=xyzinfo(1,i)-source(1)
        xyz(2)=xyzinfo(2,i)-source(2)
        xyz(3)=xyzinfo(3,i)-source(3)
c
        call dipole3eimp(rk,ceps,cmu,xyz,cjvec,evec,hvec)
c
        xyznorm(1)=xyzinfo(4,i)
        xyznorm(2)=xyzinfo(5,i)
        xyznorm(3)=xyzinfo(6,i)
c
ccc        call cross_prod3d_dcc(xyznorm,hvec,cvec)
c
        cvec(1)=hvec(1)
        cvec(2)=hvec(2)
        cvec(3)=hvec(3)
c
        xyztang1(1)=xyzinfo(7,i)
        xyztang1(2)=xyzinfo(8,i)
        xyztang1(3)=xyzinfo(9,i)
        xyztang2(1)=xyzinfo(10,i)
        xyztang2(2)=xyzinfo(11,i)
        xyztang2(3)=xyzinfo(12,i)
        call dot_prod3d_dcc(xyztang1,cvec,rhs(i,1))
        call dot_prod3d_dcc(xyztang2,cvec,rhs(i,2))
c
        cvec(1)=evec(1)
        cvec(2)=evec(2)
        cvec(3)=evec(3)
c
        xyztang1(1)=xyzinfo(7,i)
        xyztang1(2)=xyzinfo(8,i)
        xyztang1(3)=xyzinfo(9,i)
        xyztang2(1)=xyzinfo(10,i)
        xyztang2(2)=xyzinfo(11,i)
        xyztang2(3)=xyzinfo(12,i)
        call dot_prod3d_dcc(xyztang1,cvec,rhs(i,3))
        call dot_prod3d_dcc(xyztang2,cvec,rhs(i,4))
 1200   continue
c        
        return
        end
c
c
c
c
c
        subroutine em3getrhs3(rk,ceps,cmu,
     $     cjvec,cmvec,source,xyzinfo,rhs,npts)
        implicit real *8 (a-h,o-z)
        dimension source(3),xyzinfo(12,1)
        complex *16 rhs(npts,4),rk,ima,ceps,cmu
        complex *16 evec(3),hvec(3),cjvec(3),cmvec(3)
        complex *16 evec1(3),hvec1(3)
        complex *16 evec2(3),hvec2(3)
        complex *16 cvec(3)
        dimension xyz(3),xyznorm(3),xyztang1(3),xyztang2(3)
        data ima/(0.0d0,1.0d0)/
c
c
        do 1200 i=1,npts
c
        xyz(1)=xyzinfo(1,i)-source(1)
        xyz(2)=xyzinfo(2,i)-source(2)
        xyz(3)=xyzinfo(3,i)-source(3)
        call dipole3eimp(rk,ceps,cmu,xyz,cjvec,evec1,hvec1)
        call dipole3mimp(rk,ceps,cmu,xyz,cmvec,evec2,hvec2)
        evec(1)=evec1(1)+evec2(1)
        evec(2)=evec1(2)+evec2(2)
        evec(3)=evec1(3)+evec2(3)
        hvec(1)=hvec1(1)+hvec2(1)
        hvec(2)=hvec1(2)+hvec2(2)
        hvec(3)=hvec1(3)+hvec2(3)
c
        if( 1 .eq. 2 ) then
c
        call emplanew(rk,target,cjvec,cmvec,evec,hvec)
c
c       ... plane wave is coming from exterior domain, flip the rhs
c       
        evec(1)=-evec(1)
        evec(2)=-evec(2)
        evec(3)=-evec(3)
        hvec(1)=-hvec(1)
        hvec(2)=-hvec(2)
        hvec(3)=-hvec(3)
c
        endif

c
        xyznorm(1)=xyzinfo(4,i)
        xyznorm(2)=xyzinfo(5,i)
        xyznorm(3)=xyzinfo(6,i)
c
        call cross_prod3d_dcc(xyznorm,hvec,cvec)
c
c        cvec(1)=hvec(1)
c        cvec(2)=hvec(2)
c        cvec(3)=hvec(3)
c
        xyztang1(1)=xyzinfo(7,i)
        xyztang1(2)=xyzinfo(8,i)
        xyztang1(3)=xyzinfo(9,i)
        xyztang2(1)=xyzinfo(10,i)
        xyztang2(2)=xyzinfo(11,i)
        xyztang2(3)=xyzinfo(12,i)
        call dot_prod3d_dcc(xyztang1,cvec,rhs(i,1))
        call dot_prod3d_dcc(xyztang2,cvec,rhs(i,2))
c
        call cross_prod3d_dcc(xyznorm,evec,cvec)
c
c        cvec(1)=evec(1)
c        cvec(2)=evec(2)
c        cvec(3)=evec(3)
c
        xyztang1(1)=xyzinfo(7,i)
        xyztang1(2)=xyzinfo(8,i)
        xyztang1(3)=xyzinfo(9,i)
        xyztang2(1)=xyzinfo(10,i)
        xyztang2(2)=xyzinfo(11,i)
        xyztang2(3)=xyzinfo(12,i)
        call dot_prod3d_dcc(xyztang1,cvec,rhs(i,3))
        call dot_prod3d_dcc(xyztang2,cvec,rhs(i,4))
 1200   continue
c  
        return
        end
c
c
c
c
c
        subroutine em3soleva(rk,target,xyzinfo,sol,whts,npts,evec,hvec)
        implicit real *8 (a-h,o-z)
        dimension target(3),xyzinfo(12,1),whts(1)
        complex *16 sol(npts,4),rk,ima
        complex *16 evec(3),hvec(3),cjvec(3),cmvec(3)
        complex *16 evec0(3),hvec0(3)
        dimension xyz(3),xyznorm(3),xyztang1(3),xyztang2(3)
        data ima/(0.0d0,1.0d0)/
c
        evec(1)=0
        evec(2)=0
        evec(3)=0
c
        hvec(1)=0
        hvec(2)=0
        hvec(3)=0
c
        do i=1,npts
c
          xyz(1)=target(1)-xyzinfo(1,i)
          xyz(2)=target(2)-xyzinfo(2,i)
          xyz(3)=target(3)-xyzinfo(3,i)
c
          xyznorm(1)=xyzinfo(4,i)
          xyznorm(2)=xyzinfo(5,i)
          xyznorm(3)=xyzinfo(6,i)
c
          xyztang1(1)=xyzinfo(7,i)
          xyztang1(2)=xyzinfo(8,i)
          xyztang1(3)=xyzinfo(9,i)
c
          xyztang2(1)=xyzinfo(10,i)
          xyztang2(2)=xyzinfo(11,i)
          xyztang2(3)=xyzinfo(12,i)
c
          cjvec(1)=sol(i,1)*xyztang1(1)+sol(i,2)*xyztang2(1)
          cjvec(2)=sol(i,1)*xyztang1(2)+sol(i,2)*xyztang2(2)
          cjvec(3)=sol(i,1)*xyztang1(3)+sol(i,2)*xyztang2(3)
c
          call dipole3e(rk,xyz,cjvec,evec0,hvec0)
c
          evec(1)=evec(1)+whts(i)*evec0(1)
          evec(2)=evec(2)+whts(i)*evec0(2)
          evec(3)=evec(3)+whts(i)*evec0(3)
c
          hvec(1)=hvec(1)+whts(i)*hvec0(1)
          hvec(2)=hvec(2)+whts(i)*hvec0(2)
          hvec(3)=hvec(3)+whts(i)*hvec0(3)
c
          cmvec(1)=sol(i,3)*xyztang1(1)+sol(i,4)*xyztang2(1)
          cmvec(2)=sol(i,3)*xyztang1(2)+sol(i,4)*xyztang2(2)
          cmvec(3)=sol(i,3)*xyztang1(3)+sol(i,4)*xyztang2(3)
c
          call dipole3m(rk,xyz,cmvec,evec0,hvec0)
c
          evec(1)=evec(1)+whts(i)*evec0(1)
          evec(2)=evec(2)+whts(i)*evec0(2)
          evec(3)=evec(3)+whts(i)*evec0(3)
c
          hvec(1)=hvec(1)+whts(i)*hvec0(1)
          hvec(2)=hvec(2)+whts(i)*hvec0(2)
          hvec(3)=hvec(3)+whts(i)*hvec0(3)
c
        enddo
c        
        return
        end
c
c
c
c
c
        subroutine em3soleva2(rk,ceps,cmu,
     $     target,xyzinfo,sol,whts,npts,evec,hvec)
        implicit real *8 (a-h,o-z)
        dimension target(3),xyzinfo(12,1),whts(1)
        complex *16 sol(npts,4),rk,ima,ceps,cmu
        complex *16 evec(3),hvec(3),cjvec(3),cmvec(3)
        complex *16 evec0(3),hvec0(3)
        dimension xyz(3),xyznorm(3),xyztang1(3),xyztang2(3)
        data ima/(0.0d0,1.0d0)/
c
        evec(1)=0
        evec(2)=0
        evec(3)=0
c
        hvec(1)=0
        hvec(2)=0
        hvec(3)=0
c
        do 1200 i=1,npts
c
        xyz(1)=target(1)-xyzinfo(1,i)
        xyz(2)=target(2)-xyzinfo(2,i)
        xyz(3)=target(3)-xyzinfo(3,i)
c
        xyznorm(1)=xyzinfo(4,i)
        xyznorm(2)=xyzinfo(5,i)
        xyznorm(3)=xyzinfo(6,i)
c
        xyztang1(1)=xyzinfo(7,i)
        xyztang1(2)=xyzinfo(8,i)
        xyztang1(3)=xyzinfo(9,i)
c
        xyztang2(1)=xyzinfo(10,i)
        xyztang2(2)=xyzinfo(11,i)
        xyztang2(3)=xyzinfo(12,i)
c
        cjvec(1)=sol(i,1)*xyztang1(1)+sol(i,2)*xyztang2(1)
        cjvec(2)=sol(i,1)*xyztang1(2)+sol(i,2)*xyztang2(2)
        cjvec(3)=sol(i,1)*xyztang1(3)+sol(i,2)*xyztang2(3)
c       
        cjvec(1)=cjvec(1)*ceps
        cjvec(2)=cjvec(2)*ceps
        cjvec(3)=cjvec(3)*ceps
c       
        call dipole3eimp(rk,ceps,cmu,xyz,cjvec,evec0,hvec0)
c
        evec(1)=evec(1)+whts(i)*evec0(1)
        evec(2)=evec(2)+whts(i)*evec0(2)
        evec(3)=evec(3)+whts(i)*evec0(3)
c
        hvec(1)=hvec(1)+whts(i)*hvec0(1)
        hvec(2)=hvec(2)+whts(i)*hvec0(2)
        hvec(3)=hvec(3)+whts(i)*hvec0(3)
c
        cmvec(1)=sol(i,3)*xyztang1(1)+sol(i,4)*xyztang2(1)
        cmvec(2)=sol(i,3)*xyztang1(2)+sol(i,4)*xyztang2(2)
        cmvec(3)=sol(i,3)*xyztang1(3)+sol(i,4)*xyztang2(3)
c       
        cmvec(1)=cmvec(1)*cmu
        cmvec(2)=cmvec(2)*cmu
        cmvec(3)=cmvec(3)*cmu
c
        call dipole3mimp(rk,ceps,cmu,xyz,cmvec,evec0,hvec0)
c
        evec(1)=evec(1)+whts(i)*evec0(1)
        evec(2)=evec(2)+whts(i)*evec0(2)
        evec(3)=evec(3)+whts(i)*evec0(3)
c
        hvec(1)=hvec(1)+whts(i)*hvec0(1)
        hvec(2)=hvec(2)+whts(i)*hvec0(2)
        hvec(3)=hvec(3)+whts(i)*hvec0(3)
c
 1200   continue
c        
        return
        end
c
c
c
c
c
        subroutine em3direva(rk,source,target,evec,hvec)
        implicit real *8 (a-h,o-z)
        dimension source(3),target(3)
        complex *16 rk,ima
        complex *16 evec(3),hvec(3),cjvec(3)
        dimension xyz(3),xyznorm(3),xyztang1(3),xyztang2(3)
        data ima/(0.0d0,1.0d0)/
c
        cjvec(1)=1
        cjvec(2)=1
        cjvec(3)=1
c       
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        call dipole3e(rk,xyz,cjvec,evec,hvec)
c
        return
        end
c
c
c
c
c
        subroutine em3direva2(rk,ceps,cmu,source,target,evec,hvec)
        implicit real *8 (a-h,o-z)
        dimension source(3),target(3)
        complex *16 rk,ima,ceps,cmu
        complex *16 evec(3),hvec(3),cjvec(3)
        dimension xyz(3),xyznorm(3),xyztang1(3),xyztang2(3)
        data ima/(0.0d0,1.0d0)/
c
        cjvec(1)=1
        cjvec(2)=1
        cjvec(3)=1
c       
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        call dipole3eimp(rk,ceps,cmu,xyz,cjvec,evec,hvec)
c
        return
        end
c
c
c
c
c
        subroutine em3direva3(rk,ceps,cmu,
     $     cjvec,cmvec,source,target,evec,hvec)
        implicit real *8 (a-h,o-z)
        dimension source(3),target(3)
        complex *16 rk,ima,ceps,cmu
        complex *16 evec(3),hvec(3),cjvec(3),cmvec(3)
        complex *16 evec1(3),hvec1(3)
        complex *16 evec2(3),hvec2(3)
        dimension xyz(3),xyznorm(3),xyztang1(3),xyztang2(3)
        data ima/(0.0d0,1.0d0)/
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
        call dipole3eimp(rk,ceps,cmu,xyz,cjvec,evec1,hvec1)
        call dipole3mimp(rk,ceps,cmu,xyz,cmvec,evec2,hvec2)
        evec(1)=evec1(1)+evec2(1)
        evec(2)=evec1(2)+evec2(2)
        evec(3)=evec1(3)+evec2(3)
        hvec(1)=hvec1(1)+hvec2(1)
        hvec(2)=hvec1(2)+hvec2(2)
        hvec(3)=hvec1(3)+hvec2(3)
c
        if( 1 .eq. 2 ) then
        call emplanew(rk,target,cjvec,cmvec,evec,hvec)
        endif
c
        return
        end
c
c
c
c
c
        subroutine dipole3eimp(rk,ceps,cmu,xyz,cjvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic electric dipole cjvec located at the origin
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        dimension xyz(3)
        complex *16 cjvec(3),evec(3),hvec(3),rk,ima
        complex *16 ceps,cmu,cimp,cjvec0(3)
c
        data ima/(0.0d0,1.0d0)/
        save
c
        cjvec0(1)=cjvec(1)*sqrt(cmu)
        cjvec0(2)=cjvec(2)*sqrt(cmu)
        cjvec0(3)=cjvec(3)*sqrt(cmu)
c
        call green3e(rk,xyz,cjvec0,evec)
        evec(1)=evec(1)*(ima*rk) /sqrt(ceps)
        evec(2)=evec(2)*(ima*rk) /sqrt(ceps)
        evec(3)=evec(3)*(ima*rk) /sqrt(ceps)
c       
        call green3m(rk,xyz,cjvec0,hvec)
        hvec(1)=hvec(1) /sqrt(cmu)
        hvec(2)=hvec(2) /sqrt(cmu)
        hvec(3)=hvec(3) /sqrt(cmu)
c
        return
        end
c
c
c
c
c
        subroutine dipole3mimp(rk,ceps,cmu,xyz,cmvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic magnetic dipole cmvec located at the origin
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        dimension xyz(3)
        complex *16 cmvec(3),evec(3),hvec(3),rk,ima
        complex *16 ceps,cmu,cimp,cmvec0(3)
c
        data ima/(0.0d0,1.0d0)/
        save
c
        cmvec0(1)=cmvec(1)*sqrt(ceps)
        cmvec0(2)=cmvec(2)*sqrt(ceps)
        cmvec0(3)=cmvec(3)*sqrt(ceps)
c
        call green3m(rk,xyz,cmvec0,evec)
        evec(1)=evec(1) /sqrt(ceps)
        evec(2)=evec(2) /sqrt(ceps)
        evec(3)=evec(3) /sqrt(ceps)
c
        call green3e(rk,xyz,cmvec0,hvec)
        hvec(1)=hvec(1)*(-ima*rk) /sqrt(cmu)
        hvec(2)=hvec(2)*(-ima*rk) /sqrt(cmu)
        hvec(3)=hvec(3)*(-ima*rk) /sqrt(cmu)
c       
        return
        end
c
c
c
c
c
        subroutine dipole3imp(ceps,cmu,evec,hvec)
        implicit real *8 (a-h,o-z)
        complex *16 ceps,cmu,evec(3),hvec(3)
        evec(1)=evec(1)/sqrt(ceps)
        evec(2)=evec(2)/sqrt(ceps)
        evec(3)=evec(3)/sqrt(ceps)
        hvec(1)=hvec(1)/sqrt(cmu)
        hvec(2)=hvec(2)/sqrt(cmu)
        hvec(3)=hvec(3)/sqrt(cmu)
        return
        end
c
c
c
        subroutine dot_prod3d_dcc(x,y,d)
        implicit real *8 (a-h,o-z)
        dimension x(3)
        complex *16 y(3),d
c
c       d = x \dot y
c
        d=x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
c
        return
        end
c
c
c
c
c
        subroutine cross_prod3d_dcc(x,y,z)
        implicit real *8 (a-h,o-z)
        dimension x(3)
        complex *16 y(3),z(3)
c
c       z = x \cross y
c
        z(1)=x(2)*y(3)-x(3)*y(2)
        z(2)=x(3)*y(1)-x(1)*y(3)
        z(3)=x(1)*y(2)-x(2)*y(1)
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
        subroutine h3getrhs(rk,source,xyzs,rhs,npts)
        implicit real *8 (a-h,o-z)
        dimension source(3),xyzs(3,1)
        complex *16 rhs(npts),rk,ima
        data ima/(0.0d0,1.0d0)/
c
        do 1200 i=1,npts
        dx=xyzs(1,i)-source(1)
        dy=xyzs(2,i)-source(2)
        dz=xyzs(3,i)-source(3)
        r=sqrt(dx**2+dy**2+dz**2)
        rhs(i)=exp(ima*rk*r)/r
 1200   continue
c        
        return
        end
c
c
c
c
c
        subroutine h3getrhsneu(rk,source,xyzs,xyznorms,rhs,npts)
        implicit real *8 (a-h,o-z)
        dimension source(3),xyzs(3,1),xyznorms(3,1)
        complex *16 rhs(npts),rk,ima,cd
        data ima/(0.0d0,1.0d0)/
c
        do 1200 i=1,npts
        dx=xyzs(1,i)-source(1)
        dy=xyzs(2,i)-source(2)
        dz=xyzs(3,i)-source(3)
        r=sqrt(dx**2+dy**2+dz**2)
        cd=dx*xyznorms(1,i)+dy*xyznorms(2,i)+dz*xyznorms(3,i)
        cd=cd*(1-ima*rk*r)
        rhs(i)=cd*exp(ima*rk*r)/r**3
 1200   continue
c        
        return
        end
c
c
c
c
c
        subroutine h3soleva(rk,target,xyzs,sol,whts,npts,cpot)
        implicit real *8 (a-h,o-z)
        dimension target(3),xyzs(3,1),whts(1)
        complex *16 sol(npts),cpot,rk,ima
        data ima/(0.0d0,1.0d0)/
c
        cpot=0
        do 1200 i=1,npts
        dx=target(1)-xyzs(1,i)
        dy=target(2)-xyzs(2,i)
        dz=target(3)-xyzs(3,i)
        r=sqrt(dx**2+dy**2+dz**2)
        cpot=cpot+exp(ima*rk*r)/r*sol(i)*whts(i)
 1200   continue
c        
        return
        end
c
c
c
c
c
        subroutine h3direva(rk,source,target,cpot)
        implicit real *8 (a-h,o-z)
        dimension source(3),target(3)
        complex *16 cpot,rk,ima
        data ima/(0.0d0,1.0d0)/
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        dz=target(3)-source(3)
        r=sqrt(dx**2+dy**2+dz**2)
        cpot=exp(ima*rk*r)/r
c        
        return
        end
c
c
c
c
c
        subroutine em3mpzero(ampole,nterms)
        implicit real *8 (a-h,o-z)
        complex *16 ampole(0:nterms,-nterms:nterms)
c
        do n=0,nterms
        do m=-n,n
        ampole(n,m)=0
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine em3mpclear(ampole,nterms)
        implicit real *8 (a-h,o-z)
        complex *16 ampole(0:nterms,-nterms:nterms)
c
        do n=0,nterms
        do m=-nterms,nterms
        ampole(n,m)=0
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine em3mpset(ampole,nterms,n,m,cd)
        implicit real *8 (a-h,o-z)
        complex *16 ampole(0:nterms,-nterms:nterms),cd
c
        ampole(n,m)=cd
c
        return
        end
c
c
c
c
c

        subroutine em3sphlin(ampole,nterms,cvec)
        implicit real *8 (a-h,o-z)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 cvec(1)
c
c       ... compress multipole expansion storage into linear array
c
        kk=0
        do n=0,nterms
        do m=-n,n
        kk=kk+1
        cvec(kk)=ampole(n,m)
        enddo
        enddo
c
ccc        call prinf('kk=*',kk,1)
c
        return
        end
c
c
c
c
c
        subroutine em3linsph(ampole,nterms,cvec)
        implicit real *8 (a-h,o-z)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 cvec(1)
c
c       ... unroll linear storage array into multipole expansion 
c
        kk=0
        do n=0,nterms
        do m=-n,n
        kk=kk+1
        ampole(n,m)=cvec(kk)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine em3_get_inc_field_loc
     $     (rk,ampole,bmpole,nterms,xyz,evec,hvec)
        implicit real *8 (a-h,o-z)
        complex *16 rk
        dimension center(3),xyz(3)
c
        complex *16 ceps,cmu
        complex *16 evec(3),hvec(3)
        complex *16 cjvec(3),cmvec(3)
        complex *16 source_cjvec(3),source_cmvec(3)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
        complex *16 cd
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
        save
c
        call em3taeval
     $     (rk,center,ampole,bmpole,nterms,xyz,evec,hvec)
c
c       ... spherical wave is coming from exterior domain, flip the rhs
c       
        evec(1)=-evec(1)
        evec(2)=-evec(2)
        evec(3)=-evec(3)
        hvec(1)=-hvec(1)
        hvec(2)=-hvec(2)
        hvec(3)=-hvec(3)
c        
        return
        end
c
c
c
c
c
        subroutine em3getrhs4(rk,ceps,cmu,
     $     itype_mp,n,m,nterms,xyzinfo,rhs,npts)
        implicit real *8 (a-h,o-z)
        dimension source(3),xyzinfo(12,1)
        complex *16 rhs(npts,4),rk,ima,ceps,cmu,cjvec(3),cmvec(3)
        complex *16 evec(3),hvec(3),source_cjvec(3),source_cmvec(3)
        complex *16 cvec(3)
        dimension xyz(3),xyznorm(3),xyztang1(3),xyztang2(3)
 
        complex *16, allocatable :: ampole(:,:)
        complex *16, allocatable :: bmpole(:,:)

        complex *16 cd
        dimension center(3)
        data ima/(0.0d0,1.0d0)/
c
c
        allocate( ampole(0:nterms,-nterms:nterms) ) 
        allocate( bmpole(0:nterms,-nterms:nterms) ) 
c
c
        center(1)=0
        center(2)=0
        center(3)=0
        call em3mpzero(ampole,nterms)
        call em3mpzero(bmpole,nterms)
c
c        call prinf('inside getrhs4, itype_mp=*',itype_mp,1)
c        call prinf('inside getrhs4, nterms=*',nterms,1)
c        call prinf('inside getrhs4, n=*',n,1)
c        call prinf('inside getrhs4, m=*',m,1)
        
        cd=1
        if( itype_mp .eq. 1 ) call em3mpset(ampole,nterms,n,m,cd)
        if( itype_mp .eq. 2 ) call em3mpset(bmpole,nterms,n,m,cd)
c
c
c
        do 1200 i=1,npts
c
        xyz(1)=xyzinfo(1,i)-source(1)
        xyz(2)=xyzinfo(2,i)-source(2)
        xyz(3)=xyzinfo(3,i)-source(3)
c
        call em3_get_inc_field_loc
     $     (rk,ampole,bmpole,nterms,xyz,evec,hvec)
c
        xyznorm(1)=xyzinfo(4,i)
        xyznorm(2)=xyzinfo(5,i)
        xyznorm(3)=xyzinfo(6,i)
c
        call cross_prod3d_dcc(xyznorm,hvec,cvec)
c
c        cvec(1)=hvec(1)
c        cvec(2)=hvec(2)
c        cvec(3)=hvec(3)
c
        xyztang1(1)=xyzinfo(7,i)
        xyztang1(2)=xyzinfo(8,i)
        xyztang1(3)=xyzinfo(9,i)
        xyztang2(1)=xyzinfo(10,i)
        xyztang2(2)=xyzinfo(11,i)
        xyztang2(3)=xyzinfo(12,i)
        call dot_prod3d_dcc(xyztang1,cvec,rhs(i,1))
        call dot_prod3d_dcc(xyztang2,cvec,rhs(i,2))
c
        call cross_prod3d_dcc(xyznorm,evec,cvec)
c
c        cvec(1)=evec(1)
c        cvec(2)=evec(2)
c        cvec(3)=evec(3)
c
        xyztang1(1)=xyzinfo(7,i)
        xyztang1(2)=xyzinfo(8,i)
        xyztang1(3)=xyzinfo(9,i)
        xyztang2(1)=xyzinfo(10,i)
        xyztang2(2)=xyzinfo(11,i)
        xyztang2(3)=xyzinfo(12,i)
        call dot_prod3d_dcc(xyztang1,cvec,rhs(i,3))
        call dot_prod3d_dcc(xyztang2,cvec,rhs(i,4))
 1200   continue
c        
        return
        end
c
c
c
c
c
        subroutine em3soleva2a
     $     (xyzinfo,sol,whts,npts,sol_cjvecs,sol_cmvecs)
        implicit real *8 (a-h,o-z)
        dimension target(3),xyzinfo(12,1),whts(1)
        complex *16 sol(npts,4),rk,ima,ceps,cmu
        complex *16 evec(3),hvec(3),cjvec(3),cmvec(3)
        complex *16 evec0(3),hvec0(3)
        dimension xyz(3),xyznorm(3),xyztang1(3),xyztang2(3)
        complex *16 sol_cjvecs(3,1)
        complex *16 sol_cmvecs(3,1)
        data ima/(0.0d0,1.0d0)/
c
c
        do 1200 i=1,npts
c
        xyz(1)=target(1)-xyzinfo(1,i)
        xyz(2)=target(2)-xyzinfo(2,i)
        xyz(3)=target(3)-xyzinfo(3,i)
c
        xyznorm(1)=xyzinfo(4,i)
        xyznorm(2)=xyzinfo(5,i)
        xyznorm(3)=xyzinfo(6,i)
c
        xyztang1(1)=xyzinfo(7,i)
        xyztang1(2)=xyzinfo(8,i)
        xyztang1(3)=xyzinfo(9,i)
c
        xyztang2(1)=xyzinfo(10,i)
        xyztang2(2)=xyzinfo(11,i)
        xyztang2(3)=xyzinfo(12,i)
c
        cjvec(1)=sol(i,1)*xyztang1(1)+sol(i,2)*xyztang2(1)
        cjvec(2)=sol(i,1)*xyztang1(2)+sol(i,2)*xyztang2(2)
        cjvec(3)=sol(i,1)*xyztang1(3)+sol(i,2)*xyztang2(3)
c       
        cjvec(1)=cjvec(1)
        cjvec(2)=cjvec(2)
        cjvec(3)=cjvec(3)
c
        sol_cjvecs(1,i)=cjvec(1)
        sol_cjvecs(2,i)=cjvec(2)
        sol_cjvecs(3,i)=cjvec(3)
c       
c
        cmvec(1)=sol(i,3)*xyztang1(1)+sol(i,4)*xyztang2(1)
        cmvec(2)=sol(i,3)*xyztang1(2)+sol(i,4)*xyztang2(2)
        cmvec(3)=sol(i,3)*xyztang1(3)+sol(i,4)*xyztang2(3)
c       
        cmvec(1)=cmvec(1)
        cmvec(2)=cmvec(2)
        cmvec(3)=cmvec(3)
c
        sol_cmvecs(1,i)=cmvec(1)
        sol_cmvecs(2,i)=cmvec(2)
        sol_cmvecs(3,i)=cmvec(3)
c       
c
 1200   continue
c        
        return
        end
c
c
c
c
c
