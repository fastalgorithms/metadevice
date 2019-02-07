c
c
c
        implicit real *8 (a-h,o-z)
        dimension xyz(3),dxyzduv(3,2)
        dimension xyznorm(3),xyztang1(3),xyztang2(3)
        dimension g(2,2)
c    
        dimension usout(100),vsout(100),wsout(100)
        dimension umatr(10000),vmatr(10000)
c
        allocatable :: verts(:,:),ifaces(:,:),iqfaces(:,:)
        allocatable :: triainfo(:,:,:)
        allocatable :: qtriainfo(:,:,:)
        allocatable :: ipatchinfo(:),refineinfo(:,:)
c
        dimension scale_geo(3),shift_geo(3)
c
        allocatable :: ixyzs(:,:),xyzs(:,:)
        allocatable :: xyznorms(:,:)
        allocatable :: xyztang1s(:,:),xyztang2s(:,:)
        allocatable :: xyzinfo(:,:)
        allocatable :: whts(:)
c
        external fpatchpnt,qpatchpnt,cpatchpnt
        external spatchpnt,wpatchpnt,rpatchpnt
c
        external lfinter1,lfinter2,lfinter3,lfinter4
        external hfinter1,hfinter2,hfinter3,hfinter4
c
        external eminter1,eminter3,eminter4
        external eminter1h,eminter3h
        external eminter1n,eminter3n
c
        external em3multa
c
        complex *16 rk
        complex *16 ima
c
        complex *16, allocatable :: cmatr0(:)
c
        complex *16, allocatable :: cmatr(:)
        complex *16, allocatable :: rhs(:)
        complex *16, allocatable :: sol(:)
c
        complex *16 cd
c
        allocatable :: w(:)
c
        dimension source(3),target(3)
        complex *16 source_cjvec(3),source_cmvec(3)
        complex *16 cpot,cpot0
c
        dimension info(2)
        complex *16 evec(3),hvec(3)
        complex *16 evec0(3),hvec0(3)
        complex *16 evec1(3),hvec1(3)
c       
        complex *16 ceps(10),cmus(10)
        complex *16 rk_id,ceps_id,cmus_id,cems_id
c
        dimension errs(10000)
c
        character*256 config
        character*256 filename_geo
        character*256 filename_out
c
        complex *16, allocatable :: ampole(:,:)
        complex *16, allocatable :: bmpole(:,:)
        dimension center(3)
c
        complex *16, allocatable :: sol_cjvecs(:,:)
        complex *16, allocatable :: sol_cmvecs(:,:)
c


        data ima/(0.0d0,1.0d0)/
c
c
c
c       SET ALL PARAMETERS
c        
        call prini(6,13)
c
c       ... get configuration filename_geo
c
        call getarg(1,config)
        call prina('config=*',config,78)
c
        ir = 10
        open(unit=ir,file=config)
c
        read(ir,*) igeom
        read(ir,*) filename_geo
        filename_geo=trim(adjustl(filename_geo))
        call prina('filename_geo=*',filename_geo,78)
c
        read(ir,*) scale_geo(1), scale_geo(2), scale_geo(3), 
     $     shift_geo(1), shift_geo(2), shift_geo(3)
c
        read(ir,*) radius
        read(ir,*) wavelength, dreal_n, dimag_n
        read(ir,*) nterms
        read(ir,*) itype_solve, eps, numit
        
        read(ir,*) filename_out        
        filename_out=trim(adjustl(filename_out))
        call prina('filename_out=*',filename_out,78)
c
        open(unit=72,file=filename_out)
c
        done=1
        pi=4*atan(done)
c
cccc        rk=1.0d0
cccc        rk=1.0d0*pi
cccc        rk=1.0d0*pi*2
c
c        scale_geo(1)=25d0/2
c        scale_geo(2)=25d0/2
c        scale_geo(3)=75d0/2
c        shift_geo(1)=0
c        shift_geo(2)=0
c        shift_geo(3)=0
c
c        radius=60d0
c
c        wavelength = 1000d0
c
c        omega=(2*pi)/wavelength
c        ceps(1)=1
c        cmus(1)=1
c        ceps(2)=(0.2283+6.4700*ima)**2
c        cmus(2)=1
c
c        nterms=3
c
c
        omega=(2*pi)/wavelength
        ceps(1)=1
        cmus(1)=1
c
        call prin2('omega=*',omega,1)
        call prin2('ceps=*',ceps,2*1)
        call prin2('cmus=*',cmus,2*11)
c
        call prin2('wavelength=*',wavelength,1)
        call prin2('scale_geo=*',scale_geo,3)
        call prin2('shift_geo=*',shift_geo,3)
        call prinf('nterms=*',nterms,1)
        call prinf('itype_solve=*',itype_solve,1)
        call prin2('eps=*',eps,1)
        call prinf('numit=*',numit,1)
c
        call prinf('============================*',i,0)
c
        rk=omega*sqrt(cmus(1))*sqrt(ceps(1))
        call prin2('rk=*',rk,2)
c       
c
c       ... retrieve the interpolation nodes
c
        itype=1
        norder=1
        call ortho2siexps(itype,norder,npols,usout,vsout,
     1     umatr,vmatr,wsout)
c
        call prinf('norder=*',norder,1)
        call prinf('npols=*',npols,1)
        call prin2('usout=*',usout,npols)
        call prin2('vsout=*',vsout,npols)
        call prin2('wsout=*',wsout,npols)
c
        d=0
        do 1100 i=1,npols
        d=d+wsout(i)
 1100   continue
c
        call prin2('sum of weights=*',d,1)
c
c
        itest=2
        call prinf('itest=*',itest,1)
c
c       
c       ... allocate work arrays for geometry descriptor
c
        max_verts = 100000
        max_faces = 100000
        max_tri = 100000
c
        allocate( verts(3,max_verts) ) 
        allocate( ifaces(3,max_faces) )
        allocate( triainfo(3,3,max_tri) )
c
        allocate( iqfaces(6,max_faces) )
        allocate( qtriainfo(3,6,max_tri) )
c
        max_ref_tri = 100000*25
c
        allocate( ipatchinfo(max_ref_tri), refineinfo(4,max_ref_tri) )
c
c
cccc        igeom = 3
        call prinf('igeom=*',igeom,1)
c
        if( igeom .eq. 1 ) then
c
c       ... retrieve a triangulation from a library
c
        itype=4
        call rsolid(itype,verts,nverts,ifaces,nfaces)
c
        noversamp=4
c
        endif
c
        if( igeom .eq. 2 ) then
c
c       ... retrieve a triangulation from a file
c
        ir = 17
        open (unit = ir,file=filename_geo)
ccc        open (unit = ir,file='sphere320.a.tri')
ccc        open (unit = ir,file='../../Data/sphere720.a.tri')
        call atri3init(iergeom,ir,nverts,nfaces)
        close(ir)
        open (unit = ir,file=filename_geo)
ccc        open (unit = ir,file='sphere320.a.tri')
ccc        open (unit = ir,file='../../Data/sphere720.a.tri')
        call atriread3(iergeom,ir,verts,nverts,ifaces,nfaces)
c
        noversamp=1
c
        close(ir)
c
        endif
c
        if( igeom .eq. 3 ) then
c
c       ... retrieve a triangulation from a file (quadratic triangles)
c
        ir = 17
        open (unit = ir,file=filename_geo)
ccc        open (unit = ir,file='sphere320.q.tri')
ccc        open (unit = ir,file='../../Data/sphere720.q.tri')
        call qtri3init(iergeom,ir,nverts,nfaces)
        close(ir)
        open (unit = ir,file=filename_geo)
ccc        open (unit = ir,file='sphere320.q.tri')
ccc        open (unit = ir,file='../../Data/sphere720.q.tri')
        call qtriread3(iergeom,ir,verts,nverts,iqfaces,nfaces)
c
        noversamp=1
c
        close(ir)
c
        endif
c
c
c
c
c       ... optional, resize and shift
c
        do i=1,nverts
        verts(1,i)=verts(1,i)*scale_geo(1) + shift_geo(1)
        verts(2,i)=verts(2,i)*scale_geo(2) + shift_geo(2)
        verts(3,i)=verts(3,i)*scale_geo(3) + shift_geo(3)
        enddo
c
c
        if( igeom .eq. 1 .or. igeom .eq. 2 ) then
c
        call genqtriainfo_flat(verts,nverts,ifaces,nfaces,qtriainfo)
        npatches=nfaces
c        
        endif
c
        if( igeom .eq. 3 ) then
c
        call genqtriainfo(verts,nverts,iqfaces,nfaces,qtriainfo)
        npatches=nfaces
c        
        endif
c
c
ccc        call prin2('verts=*',verts,3*nverts)
ccc        call prinf('iqfaces=*',iqfaces,6*nfaces)
ccc        call prin2('qtriainfo=*',qtriainfo,3*6*npatches)
c
c
c       ... refine the triangulation
c
cccc        noversamp=1
        call genrefineinfo(noversamp,npatches,
     $     npatchesout,ipatchinfo,refineinfo)
c
        npatches=npatchesout
c
        call prinf('after oversampling, npatches=*',npatches,1)
ccc        call prinf('ipatchinfo=*',ipatchinfo,npatches)
ccc        call prin2('refineinfo=*',refineinfo,4*npatches)
c
c
c       ... map the interpolation points on the patch into R^3 
c        
cccc        ipatch=1
c
        npts=npols*npatches
        call prinf('npts=*',npts,1)
        call prinf('npatches=*',npatches,1)
        call prinf('npols=*',npols,1)
c
c
c
c       ... allocate work arrays for discretization 
c
        allocate( ixyzs(2,npts), xyzs(3,npts) )
        allocate( xyznorms(3,npts) )
        allocate( xyztang1s(3,npts), xyztang2s(3,npts) )
        allocate( xyzinfo(12,npts) )
        allocate( whts(npts) )
c
c       ... allocate temporary matrix for patchmatc discretizer
c       
        allocate( cmatr0(npts*npts) )
c
c       ... allocate work arrays for the solver
c
        allocate( cmatr(2*npts*2*npts) )
        allocate( rhs(2*npts) ) 
        allocate( sol(2*npts) )
c
c
c
c       ... plot all discretization nodes and normals
c
c       
ccc        call prinf('============================*',i,0)
c
c
        do 1400 j=1,npatches
        do 1200 i=1,npols
c
        ipatch=j
ccc        call prinf('ipatch=*',ipatch,1)
c
        u=usout(i)
        v=vsout(i)
c
        call patchgeo(rpatchpnt,ipatch,u,v,
     $     qtriainfo,qpatchpnt,ipatchinfo,refineinfo,
     $     xyz,dxyzduv,ds,xyznorm,xyztang1,xyztang2)
c
        write(16,*) xyz(1),xyz(2),xyz(3)
        write(17,*) xyz(1),xyz(2),xyz(3)
        write(17,*) 
     $     xyz(1)+xyznorm(1)/3,
     $     xyz(2)+xyznorm(2)/3,
     $     xyz(3)+xyznorm(3)/3
        write(17,*) 
        write(17,*) 
c
        write(18,*) xyz(1),xyz(2),xyz(3)
        write(18,*) 
     $     xyz(1)+xyztang1(1)/3,
     $     xyz(2)+xyztang1(2)/3,
     $     xyz(3)+xyztang1(3)/3
        write(18,*) 
        write(18,*) 
        write(19,*) xyz(1),xyz(2),xyz(3)
        write(19,*) 
     $     xyz(1)+xyztang2(1)/3,
     $     xyz(2)+xyztang2(2)/3,
     $     xyz(3)+xyztang2(3)/3
        write(19,*) 
        write(19,*) 
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
cc        call prin2('g=*',g,4)
c
        detg=g11*g22-g12*g21
c
cc        call prin2('sqrt(det|g|)=*',sqrt(detg),1)
cc        call prin2('and ds=*',ds,1)
c
 1200   continue
 1400   continue
c
c
c       ... generate all discretization nodes and normals
c
c       
        call prinf('============================*',i,0)
c
        call patchallpnts(npatches,rpatchpnt,
     $     qtriainfo,qpatchpnt,ipatchinfo,refineinfo,
     $     npols,usout,vsout,ixyzs,xyzs,xyznorms,
     $     xyztang1s,xyztang2s,npts)
c
        call prinf('npts=*',npts,1)
        call prinf('npatches=*',npatches,1)
        call prinf('npols=*',npols,1)
ccc        call prinf('ixyzs=*',ixyzs,2*npatches)
ccc        call prin2('xyzs=*',xyzs,3*npatches)
ccc        call prin2('xyznorms=*',xyznorms,3*npatches)
c
        call patchallwhts(npatches,rpatchpnt,
     $     qtriainfo,qpatchpnt,ipatchinfo,refineinfo,
     $     npols,usout,vsout,wsout,whts,npts)
c
c
c
c       ... call patchmatc discretizer
c
        lw=2 000 000
        allocate( w(lw) )
c       
        call prinf('============================*',i,0)
c
        t1=second()
C$        t1=omp_get_wtime()
c
        id=1
        rk_id=omega*sqrt(cmus(id))*sqrt(ceps(id))
        ceps_id=ceps(id)
        cmus_id=cmus(id)
        cems_id=sqrt(cmus(id))*sqrt(ceps(id))
c
        info(1)=1
        info(2)=1
        call prinf('tangential component: H(11)=*',id,1)
        call patchmatc(npatches,rpatchpnt,
     $     qtriainfo,qpatchpnt,ipatchinfo,refineinfo,
     $     norder,npols,usout,vsout,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     eminter3n,rk_id,info,par7,par8,
     $     cmatr0,w,lw,lused,ier)
c
        call patchdiag(cmatr0,npts,-2*pi)
        call em3submulcpy
     $     (2*npts,cmatr,npts,npts,cmatr0,1,1,ceps_id)
c
        info(1)=1
        info(2)=2
        call prinf('tangential component: H(12)=*',id,1)
        call patchmatc(npatches,rpatchpnt,
     $     qtriainfo,qpatchpnt,ipatchinfo,refineinfo,
     $     norder,npols,usout,vsout,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     eminter3n,rk_id,info,par7,par8,
     $     cmatr0,w,lw,lused,ier)
c
        call em3submulcpy
     $     (2*npts,cmatr,npts,npts,cmatr0,npts+1,1,ceps_id)
c
        info(1)=2
        info(2)=1
        call prinf('tangential component: H(21)=*',id,1)
        call patchmatc(npatches,rpatchpnt,
     $     qtriainfo,qpatchpnt,ipatchinfo,refineinfo,
     $     norder,npols,usout,vsout,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     eminter3n,rk_id,info,par7,par8,
     $     cmatr0,w,lw,lused,ier)
c
        call em3submulcpy
     $     (2*npts,cmatr,npts,npts,cmatr0,1,npts+1,ceps_id)
c
        info(1)=2
        info(2)=2
        call prinf('tangential component: H(22)=*',id,1)
        call patchmatc(npatches,rpatchpnt,
     $     qtriainfo,qpatchpnt,ipatchinfo,refineinfo,
     $     norder,npols,usout,vsout,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     eminter3n,rk_id,info,par7,par8,
     $     cmatr0,w,lw,lused,ier)
c
        call patchdiag(cmatr0,npts,-2*pi)
        call em3submulcpy
     $     (2*npts,cmatr,npts,npts,cmatr0,npts+1,npts+1,ceps_id)
c
        t2=second()
C$        t2=omp_get_wtime()
c
        call prinf('after patchmatc, ier=*',ier,1)
        call prin2('in patchmatc, time=*',t2-t1,1)
c
cccc        call em3debug(4*npts,cmatr,whts)
c
        do i=1,npts
        call patchinfo(xyzs(1,i),xyznorms(1,i),
     $     xyztang1s(1,i),xyztang2s(1,i),xyzinfo(1,i))
        enddo
c
ccc        call prin2('whts=*',whts,npts)
c
        cd=0
        do i=1,npts
        cd=cd+whts(i)
        enddo
c
        call prin2('sum, whts=*',cd,2)
        call prin2('whts-4 pi=*',cd-4*pi,2)
c
c
        ifexterior=1
c
c
        goto 6000
c
c
c       .. simple test for MFIE solver
c
c       .. construct the right hand side
c
c
        source_cjvec(1)=1
        source_cjvec(2)=1
        source_cjvec(3)=1
        source_cmvec(1)=0
        source_cmvec(2)=0
        source_cmvec(3)=0
c
c
        if( ifexterior .eq. 0 ) then
c
        source(1)=0.2
        source(2)=-0.1
        source(3)=0.3
c
        target(1)=10
        target(2)=20
        target(3)=-30

        endif

        if( ifexterior .eq. 1 ) then
c
        target(1)=0.2
        target(2)=-0.1
        target(3)=0.3
c
        source(1)=10
        source(2)=20
        source(3)=-30
c
        endif
c
c
c
        id = 1
        rk_id=omega*sqrt(cmus(id))*sqrt(ceps(id))
        ceps_id=ceps(id)
        cmus_id=cmus(id)
        cems_id=sqrt(cmus(id))*sqrt(ceps(id))
        call em3getrhs3(rk_id,ceps_id,cmus_id,
     $     source_cjvec,source_cmvec,source,xyzinfo,rhs,npts)
        if( ifexterior .eq. 0 ) then
        do i=1,npts*2
        rhs(i)=-rhs(i)
        enddo
        endif
        
c
c
c
c       call the solver
c
c
        call prinf('npts=*',npts,1)
c
ccc        itype_solve = 2
c
c
        if( itype_solve .eq. 1 ) then
c
        if( allocated(w) ) deallocate(w)
        lw = 2*(2*npts)**2+(2*npts)+10
        lw = 2*lw
        allocate( w(lw), stat=ier )
ccc        write(*,*) 'ier=',ier
c
        call prinf('entering cqrdecom, 2 x npts=*',2*npts,1)
c
        t1=second()
C$        t1=omp_get_wtime()
        call cqrdecom(cmatr,2*npts,w,rcond)
        t2=second()
C$        t2=omp_get_wtime()
c
        call prin2('after cqrdecom, rcond=*',rcond,1)
        call prin2('in cqrdecom, time=*',t2-t1,1)
c
        t1=second()
C$        t1=omp_get_wtime()
        call cqrsolve(2*npts,w,rhs,sol)
        t2=second()
C$        t2=omp_get_wtime()
c
ccc        call prin2('after cqrsolve, rhs=*',rhs,2*npts)
ccc        call prin2('after cqrsolve, sol=*',sol,2*npts)
        call prin2('in cqrsolve, time=*',t2-t1,1)
c
        endif
c
c
        if( itype_solve .eq. 2 ) then
c
        call prinf('entering cgmres, 2 x npts=*',2*npts,1)
c
        t1=second()
C$        t1=omp_get_wtime()
ccc        eps=1e-5
ccc        numit=40
        ngmrec=numit
c
        if( allocated(w) ) deallocate(w)
        allocate( w( 2*(ngmrec*2+4)*(2*npts) ) )
c 
        call cgmres(ier,2*npts,cmatr,
     $     em3multa,par1,par2,rhs,eps,numit,
     1     sol,niter,errs,ngmrec,w)
        t2=second()
C$        t2=omp_get_wtime()
c
        call prinf('after cgmres, ier=*',ier,1)
        call prinf('after cgmres, niter=*',niter,1)
        call prin2('after cgmres, errs=*',errs,niter)
        call prin2('in cgmres, time=*',t2-t1,1)
c
        endif
c
c
c
        if( itype_solve .eq. 3 ) then
c
        call prinf('entering cbicgstab, 2 x npts=*',2*npts,1)
c
        t1=second()
C$        t1=omp_get_wtime()
ccc        eps=1e-5
ccc        numit=40
c
        if( allocated(w) ) deallocate(w)
        allocate( w( 2*11*(2*npts) ) )
c 
        call cbicgstab(ier,2*npts,cmatr,
     $     em3multa,par1,par2,rhs,eps,numit,
     1     sol,niter,errs,w)
        t2=second()
C$        t2=omp_get_wtime()
c
        call prinf('after cbicgstab, ier=*',ier,1)
        call prinf('after cbicgstab, niter=*',niter,1)
        call prin2('after cbicgstab, errs=*',errs,niter)
        call prin2('in cbicgstab, time=*',t2-t1,1)
c
        endif
c
c
c
c
c       ... check the solution 
c
        id = 1
        rk_id=omega*sqrt(cmus(id))*sqrt(ceps(id))
        ceps_id=ceps(id)
        cmus_id=cmus(id)
        cems_id=sqrt(cmus(id))*sqrt(ceps(id))
        call em3soleva2(rk_id,ceps_id,cmus_id,
     $     target,xyzinfo,sol,whts,npts,evec,hvec)
c
        call prin2('target=*',target,3)
        call prin2('evec=*',evec,6)
        call prin2('hvec=*',hvec,6)
c
ccc        call em3direva(rk,source,target,evec0,hvec0)
c
        id = 1
        rk_id=omega*sqrt(cmus(id))*sqrt(ceps(id))
        ceps_id=ceps(id)
        cmus_id=cmus(id)
        cems_id=sqrt(cmus(id))*sqrt(ceps(id))
        call em3direva3(rk_id,ceps_id,cmus_id,
     $     source_cjvec,source_cmvec,source,target,evec0,hvec0)
c
        call prin2('directly, evec0=*',evec0,6)
        call prin2('directly, hvec0=*',hvec0,6)
c
        do i=1,3
        evec1(1)=(evec(1)-evec0(1))/evec0(1)
        evec1(2)=(evec(2)-evec0(2))/evec0(2)
        evec1(3)=(evec(3)-evec0(3))/evec0(3)
        hvec1(1)=(hvec(1)-hvec0(1))/hvec0(1)
        hvec1(2)=(hvec(2)-hvec0(2))/hvec0(2)
        hvec1(3)=(hvec(3)-hvec0(3))/hvec0(3)
        enddo

        call prin2('relative error, E=*',evec1,6)
        call prin2('relative error, H=*',hvec1,6)
c
ccc        stop
c
        do i=1,3
        evec1(1)=(evec(1)-evec0(1))
        evec1(2)=(evec(2)-evec0(2))
        evec1(3)=(evec(3)-evec0(3))
        hvec1(1)=(hvec(1)-hvec0(1))
        hvec1(2)=(hvec(2)-hvec0(2))
        hvec1(3)=(hvec(3)-hvec0(3))
        enddo

        call prin2('absolute error, E=*',evec1,6)
        call prin2('absolute error, H=*',hvec1,6)
c
ccc        stop

 6000   continue

c
c
c
c       ... call MFIE solver for multiple incoming fields
c
c

        
        allocate( ampole(0:nterms,-nterms:nterms) ) 
        allocate( bmpole(0:nterms,-nterms:nterms) ) 

        allocate( sol_cjvecs(3,2*npts) )
        allocate( sol_cmvecs(3,2*npts) )
c

        id = 1
        rk_id=omega*sqrt(cmus(id))*sqrt(ceps(id))
        ceps_id=ceps(id)
        cmus_id=cmus(id)
        cems_id=sqrt(cmus(id))*sqrt(ceps(id))
c

        call prini(0,13)
c
c
c       ... call Muller solver for multiple incoming fields
c
c
        call xprini(0,72)
c
ccc        call xprin2('# gold ellipsoid*',i,0)
ccc        call xprin2('omega=(2 pi)/wavelength=*',omega,1)
        call xprin2('omega=*',omega,1)
        call xprin2('rk=*',rk,2)
        call xprin2('radius=*',radius,1)


        call xprinf('nterms=*',nterms,1)


        do 3200 itype_mp=1,2
        do 3100 n=1,nterms
        do 3000 m=-n,n

        call xprinf('itype_mp=*',itype_mp,1)
        call xprinf('n=*',n,1)
        call xprinf('m=*',m,1)

        call em3getrhs4(rk_id,ceps_id,cmus_id,
     $     itype_mp,n,m,nterms,xyzinfo,rhs,npts)
        if( ifexterior .eq. 0 ) then
        do i=1,npts*2
        rhs(i)=-rhs(i)
        enddo
        endif
        
c
        call prinf('npts=*',npts,1)
c
c       ... add the diagonal term
c
ccc        itype_solve = 2
c
        if( itype_solve .eq. 1 ) then
c
        call prinf('entering cqrdecom, 2 x npts=*',2*npts,1)
c
        if( allocated(w) ) deallocate(w)
        lw = 2*(2*npts)**2+(2*npts)+10
        lw = 2*lw
        allocate( w(lw), stat=ier )
ccc        write(*,*) 'ier=',ier
c
        t1=second()
C$        t1=omp_get_wtime()
        call cqrdecom(cmatr,4*npts,w,rcond)
        t2=second()
C$        t2=omp_get_wtime()
c
        call prin2('after cqrdecom, rcond=*',rcond,1)
        call prin2('in cqrdecom, time=*',t2-t1,1)
c
        t1=second()
C$        t1=omp_get_wtime()
        call cqrsolve(4*npts,w,rhs,sol)
        t2=second()
C$        t2=omp_get_wtime()
c
ccc        call prin2('after cqrsolve, rhs=*',rhs,2*2*npts)
ccc        call prin2('after cqrsolve, sol=*',sol,2*2*npts)
        call prin2('in cqrsolve, time=*',t2-t1,1)
c
        endif
c
c
        if( itype_solve .eq. 2 ) then
c
        call prinf('entering cgmres, 2 x npts=*',2*npts,1)
c
        t1=second()
C$        t1=omp_get_wtime()
ccc        eps=1e-5
ccc        numit=40
        ngmrec=numit
c
        if( allocated(w) ) deallocate(w)
        allocate( w( 2*(ngmrec*2+4)*(2*npts) ) )
c 
        call cgmres(ier,2*npts,cmatr,
     $     em3multa,par1,par2,rhs,eps,numit,
     1     sol,niter,errs,ngmrec,w)
        t2=second()
C$        t2=omp_get_wtime()
c
        call prinf('after cgmres, ier=*',ier,1)
        call prinf('after cgmres, niter=*',niter,1)
        call prin2('after cgmres, errs=*',errs,niter)
        call prin2('in cgmres, time=*',t2-t1,1)
c
        endif
c
c
        if( itype_solve .eq. 3 ) then
c
        call prinf('entering cbicgstab, 2 x npts=*',2*npts,1)
c
        t1=second()
C$        t1=omp_get_wtime()
ccc        eps=1e-5
ccc        numit=60
c
        if( allocated(w) ) deallocate(w)
        allocate( w( 2*11*(2*npts) ) )
c 
        call cbicgstab(ier,2*npts,cmatr,
     $     em3multa,par1,par2,rhs,eps,numit,
     1     sol,niter,errs,w)
        t2=second()
C$        t2=omp_get_wtime()
c
        call prinf('after cbicgstab, ier=*',ier,1)
        call prinf('after cbicgstab, niter=*',niter,1)
        call prin2('after cbicgstab, errs=*',errs,niter)
        call prin2('in cbicgstab, time=*',t2-t1,1)
c
        endif
c
c
c
c
ccc        call em3soleva(rk,target,xyzinfo,sol,whts,npts,evec,hvec)
c
        if( 1.eq.2 ) then
        id = 1
        rk_id=omega*sqrt(cmus(id))*sqrt(ceps(id))
        ceps_id=ceps(id)
        cmus_id=cmus(id)
        cems_id=sqrt(cmus(id))*sqrt(ceps(id))
        call em3soleva2(rk_id,ceps_id,cmus_id,
     $     target,xyzinfo,sol,whts,npts,evec,hvec)
c
        call prin2('target=*',target,3)
        call prin2('evec=*',evec,6)
        call prin2('hvec=*',hvec,6)
        endif
c
c
c       ... evaluate the outgoing multipole expansion
c
        call em3soleva2a(
     $     xyzinfo,sol,whts,npts,sol_cjvecs,sol_cmvecs)
c
ccc        call prin2('sol_cjvecs=*',sol_cjvecs,3*npts)
ccc        call prin2('sol_cmvecs=*',sol_cmvecs,3*npts)
c
        do i=1,npts
        sol_cjvecs(1,i)=sol_cjvecs(1,i)*whts(i)
        sol_cjvecs(2,i)=sol_cjvecs(2,i)*whts(i)
        sol_cjvecs(3,i)=sol_cjvecs(3,i)*whts(i)
        sol_cmvecs(1,i)=sol_cmvecs(1,i)*whts(i)
        sol_cmvecs(2,i)=sol_cmvecs(2,i)*whts(i)
        sol_cmvecs(3,i)=sol_cmvecs(3,i)*whts(i)
        enddo
c
c
        center(1)=0
        center(2)=0
        center(3)=0
c
        call em3formmp
     $     (rk,xyzs,sol_cjvecs,sol_cmvecs,
     $     npts,center,ampole,bmpole,nterms)
c
        call em3sphlin(ampole,nterms,w)
        call xprin2('outgoing ampole=*',w,2*(nterms+1)**2)
        call em3sphlin(bmpole,nterms,w)
        call xprin2('outgoing bmpole=*',w,2*(nterms+1)**2)
c
 3000   continue
 3100   continue
 3200   continue





        stop
        end
c
c
c
c
c
        subroutine em3multa(a,par1,par2,x,y,n)
        implicit real *8 (a-h,o-z)
        complex *16 a(n,n),x(n),y(n)
c
c$OMP PARALLEL DO DEFAULT(SHARED)
        do 1200 i=1,n
        y(i)=0
        do 1100 j=1,n
        y(i)=y(i)+a(i,j)*x(j)
 1100   continue
 1200   continue
c$OMP END PARALLEL DO
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
        complex *16 rhs(npts,2),rk,ima
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
        complex *16 rhs(npts,2),rk,ima,ceps,cmu
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
        complex *16 rhs(npts,2),rk,ima,ceps,cmu
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
        complex *16 sol(npts,2),rk,ima
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
 1200   continue
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
        complex *16 sol(npts,2),rk,ima,ceps,cmu
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
        complex *16 rhs(npts,2),rk,ima,ceps,cmu,cjvec(3),cmvec(3)
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
        complex *16 sol(npts,2),rk,ima,ceps,cmu
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
        sol_cmvecs(1,i)=0
        sol_cmvecs(2,i)=0
        sol_cmvecs(3,i)=0
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
