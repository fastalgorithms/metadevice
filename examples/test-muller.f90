program test_muller

  implicit double precision (a-h,o-z)

  dimension xyz(3),dxyzduv(3,2)
  dimension xyznorm(3),xyztang1(3),xyztang2(3)
  dimension g(2,2)

  dimension usout(100),vsout(100),wsout(100), uvs(2,100)
  dimension umatr(10000),vmatr(10000)

  allocatable :: verts(:,:),ifaces(:,:),iqfaces(:,:)
  allocatable :: triainfo(:,:,:)
  allocatable :: qtriainfo(:,:,:)
  allocatable :: ipatchinfo(:),refineinfo(:,:)

  dimension scale_geo(3),shift_geo(3)

  allocatable :: ixyzs(:,:),xyzs(:,:)
  allocatable :: xyznorms(:,:)
  allocatable :: xyztang1s(:,:),xyztang2s(:,:)
  allocatable :: xyzinfo(:,:)
  allocatable :: whts(:)

  external fpatchpnt,qpatchpnt,cpatchpnt
  external spatchpnt,wpatchpnt,rpatchpnt, epatchpnt

  external lfinter1,lfinter2,lfinter3,lfinter4
  external hfinter1,hfinter2,hfinter3,hfinter4

  external eminter1,eminter3,eminter4
  external eminter1h,eminter3h
  external eminter1n,eminter3n

  external em3multa

  double complex :: rk, ima

  double complex, allocatable :: cmatr0(:)
  double complex, allocatable :: cmatr(:)
  double complex, allocatable :: rhs(:)
  double complex, allocatable :: sol(:)

  double complex cd

  allocatable :: w(:)

  dimension source(3),targ(3)
  double complex source_cjvec(3),source_cmvec(3)
  double complex cpot,cpot0

  dimension info(2)
  double complex evec(3),hvec(3)
  double complex evec0(3),hvec0(3)
  double complex evec1(3),hvec1(3)
  double complex :: evec2(3),hvec2(3)

  double complex ceps(10),cmus(10)
  double complex rk_id,ceps_id,cmus_id,cems_id

  dimension errs(10000)

  character*256 config
  character*256 filename_geo
  character*256 filename_out

  double complex, allocatable :: ampole(:,:)
  double complex, allocatable :: bmpole(:,:)
  dimension center(3)

  double complex, allocatable :: sol_cjvecs(:,:)
  double complex, allocatable :: sol_cmvecs(:,:)


  ima = (0,1)
  
  !
  ! SET ALL PARAMETERS
  !        
  call prini(6,13)

  !
  ! ... get configuration filename_geo
  !
  ! call getarg(1,config)
  ! call prina('config file =*',config,78)

  ir = 10
  !open(unit=ir,file=config)
  open(unit=ir, file='config_muller.txt')

  read(ir,*) igeom
  igeom = 1
  
  read(ir,'(a)') filename_geo
  filename_geo=trim(adjustl(filename_geo))
  call prina('filename_geo=*',filename_geo,78)

  read(ir,*) scale_geo(1), scale_geo(2), scale_geo(3),  &
      shift_geo(1), shift_geo(2), shift_geo(3)

  read(ir,*) radius
  read(ir,*) wavelength, dreal_n, dimag_n
  read(ir,*) nterms
  read(ir,*) itype_solve, eps, numit
        
  read(ir,*) filename_out        
  filename_out=trim(adjustl(filename_out))
  call prina('filename_out=*',filename_out,78)
  

  done=1
  pi=4*atan(done)

  !
  ! . . . override the parameters if need be
  !
  ! rk=1.0d0
  ! rk=1.0d0*pi
  ! rk=1.0d0*pi*2

  ! scale_geo(1)=25d0/2
  ! scale_geo(2)=25d0/2
  ! scale_geo(3)=75d0/2
  ! shift_geo(1)=0
  ! shift_geo(2)=0
  ! shift_geo(3)=0

  ! radius=60d0

  ! wavelength = 1000d0
  ! nterms=3

  omega=(2*pi)/wavelength
  ceps(1)=1
  cmus(1)=1
  ceps(2)=(dreal_n+ima*dimag_n)**2
  cmus(2)=1

  call prin2('omega=*',omega,1)
  call prin2('ceps=*',ceps,2*2)
  call prin2('cmus=*',cmus,2*2)

  call prin2('wavelength=*',wavelength,1)
  call prin2('Re(n)=*',dreal_n,1)
  call prin2('Im(n)=*',dimag_n,1)
  call prin2('scale_geo=*',scale_geo,3)
  call prin2('shift_geo=*',shift_geo,3)
  call prinf('nterms=*',nterms,1)
  call prinf('itype_solve=*',itype_solve,1)
  call prin2('eps=*',eps,1)
  call prinf('numit=*',numit,1)

  call prinf('============================*',i,0)

  rk=omega*sqrt(cmus(1))*sqrt(ceps(1))
  call prin2('rk=*',rk,2)

  !
  ! ... retrieve the interpolation nodes
  ! norder sets the order of discretization (i.e. the number of
  ! unknowns per triangle)
  !
  itype=1
  norder = 1
  call ortho2siexps(itype,norder,npols,usout,vsout, &
      umatr,vmatr,wsout)

  call prinf('norder=*',norder,1)
  call prinf('npols=*',npols,1)
  call prin2('usout=*',usout,npols)
  call prin2('vsout=*',vsout,npols)
  call prin2('wsout=*',wsout,npols)

  d=0
  do  i=1,npols
    d=d+wsout(i)
  end do

  call prin2('sum of weights=*',d,1)

  !       
  ! ... allocate work arrays for geometry descriptor
  !
  max_verts = 100000
  max_faces = 100000
  max_tri = 100000

  allocate( verts(3,max_verts) ) 
  allocate( ifaces(3,max_faces) )
  allocate( triainfo(3,3,max_tri) )

  allocate( iqfaces(6,max_faces) )
  allocate( qtriainfo(3,6,max_tri) )

  max_ref_tri = 100000*25
  allocate( ipatchinfo(max_ref_tri), refineinfo(4,max_ref_tri) )


  !cccc        igeom = 3
  call prinf('igeom=*',igeom,1)
  !c
  !c
  !c       ... compute a triangulation exactly
  !c
  if( igeom .eq. 1 ) then
    itype = 4
    call rsolid(itype,verts,nverts,ifaces,nfaces)
    noversamp = 5
  endif

  !c
  !c       ... retrieve a triangulation from a file
  !c
  if( igeom .eq. 2 ) then
    ir = 17
    open (unit = ir,file=filename_geo)
    call atri3init(iergeom,ir,nverts,nfaces)
    close(ir)

    open (unit = ir,file=filename_geo)
    call atriread3(iergeom,ir,verts,nverts,ifaces,nfaces)
    noversamp=1
    close(ir)
  endif

  !c
  !c
  !c       ... retrieve a triangulation from a file (quadratic triangles)
  !c
  if( igeom .eq. 3 ) then
    ir = 17
    open (unit = ir,file=filename_geo)
    call qtri3init(iergeom,ir,nverts,nfaces)
    close(ir)

    open (unit = ir,file=filename_geo)
    call qtriread3(iergeom,ir,verts,nverts,iqfaces,nfaces)
    noversamp=1
    close(ir)
  endif

  !
  !
  ! package up info into qtriainfo
  !
  
  if (igeom .eq. 1) then
    call gentriainfo(verts,nverts,ifaces,nfaces,qtriainfo)
    npatches=nfaces
  endif

  if (igeom .eq. 2) then
    do i=1,nverts
      verts(1,i)=verts(1,i)*scale_geo(1) + shift_geo(1)
      verts(2,i)=verts(2,i)*scale_geo(2) + shift_geo(2)
      verts(3,i)=verts(3,i)*scale_geo(3) + shift_geo(3)
    enddo
    call genqtriainfo_flat(verts,nverts,ifaces,nfaces,qtriainfo)
    npatches=nfaces
  endif

  if( igeom .eq. 3 ) then
    do i=1,nverts
      verts(1,i)=verts(1,i)*scale_geo(1) + shift_geo(1)
      verts(2,i)=verts(2,i)*scale_geo(2) + shift_geo(2)
      verts(3,i)=verts(3,i)*scale_geo(3) + shift_geo(3)
    enddo
    call genqtriainfo(verts,nverts,iqfaces,nfaces,qtriainfo)
    npatches=nfaces
  endif


  !
  ! ... refine the triangulation if so desired (by generating
  ! refinment info)
  !
  call genrefineinfo(noversamp,npatches, &
      npatchesout,ipatchinfo,refineinfo)

  npatches=npatchesout

  call prinf('after oversampling, npatches=*',npatches,1)


  !
  ! plot these triangles
  !
  iw = 100
  kover = 0
  call xtri_vtk_surf(iw, npatches, rpatchpnt, &
      qtriainfo,epatchpnt,ipatchinfo,refineinfo, kover, &
      'the nanofin discretization')


  !c
  !c       ... map the interpolation points on the patch into R^3 
  !c        
  npts=npols*npatches
  call prinf('npts=*',npts,1)
  call prinf('npatches=*',npatches,1)
  call prinf('npols=*',npols,1)


  !
  !  ... allocate work arrays for discretization 
  !
  allocate( ixyzs(2,npts), xyzs(3,npts) )
  allocate( xyznorms(3,npts) )
  allocate( xyztang1s(3,npts), xyztang2s(3,npts) )
  allocate( xyzinfo(12,npts) )
  allocate( whts(npts) )

  !
  ! ... allocate temporary matrix for patchmatc discretizer
  !       
  allocate( cmatr0(npts*npts) )

  !
  ! ... allocate work arrays for the solver
  !
  gb = 4*npts*4.0d0*npts*2*8/((2**10)*(2**10)*(2**10))
  call prin2('trying to allocate mem, gb = *', gb, 1)
  
  allocate( cmatr(4*npts*4*npts) )
  allocate( rhs(4*npts) ) 
  allocate( sol(4*npts) )


  
  !
  ! ... generate all discretization nodes and normals
  !
  call prinf('============================*',i,0)

  call patchallpnts(npatches,rpatchpnt, &
      qtriainfo,epatchpnt,ipatchinfo,refineinfo, &
      npols,usout,vsout,ixyzs,xyzs,xyznorms, &
      xyztang1s,xyztang2s,npts)
  
  call prinf('npts=*',npts,1)
  call prinf('npatches=*',npatches,1)
  call prinf('npols=*',npols,1)
  ! call prinf('ixyzs=*',ixyzs,2*npatches)
  ! call prin2('xyzs=*',xyzs,3*npatches)
  ! call prin2('xyznorms=*',xyznorms,3*npatches)

  call patchallwhts(npatches,rpatchpnt, &
      qtriainfo,epatchpnt,ipatchinfo,refineinfo, &
      npols,usout,vsout,wsout,whts,npts)


  !
  ! ... call patchmatc discretizer
  !
  !
  ! note that the order of the self-quadrature is set by the paraemter
  ! norder0 in the subroutine patchmatc in the file patchmatc4.f, and
  ! is not accessible outside of that routine
  !
  ! furthermore, the order and precision of the adpative discretizer
  ! for the off-diagonal blocks is set by m and eps in the subroutine
  ! patchmatc_od in the file patchmatc4.f
  !
  call prinf('============================*',i,0)
  call prinf('... assembling system matrix*',i,0)



  !call build_muller_matrix(info, npatches, rpatchpnt, qtriainfo, &
  !    epatchpnt, ipatchinfo, refineinfo, &
  !    norder, npols, usout, vsout, umatr, vmatr, &
  !    ixyzs, xyzs, xyznorms, xyztang1s, xyztang2s, npts, &
  !    cmatr, ier)


  
  lw=2000000
  allocate( w(lw) )


  
  id=1
  rk_id=omega*sqrt(cmus(id))*sqrt(ceps(id))
  ceps_id=ceps(id)
  cmus_id=cmus(id)
  cems_id=sqrt(cmus(id))*sqrt(ceps(id))

  info(1)=1
  info(2)=1
  call prinf('tangential component: H(11)=*',id,1)

  !do i =1,npols
  !  uvs(1,i) = usout(i)
  !  uvs(2,i) = vsout(i)
  !end do
  
  call patchmatc(npatches,rpatchpnt, &
      qtriainfo,epatchpnt,ipatchinfo,refineinfo, &
      norder,npols,usout,vsout,umatr,vmatr, &
      ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts, &
      eminter3n,rk_id,info,par7,par8, &
      cmatr0,w,lw,lused,ier)

  call patchdiag(cmatr0,npts,+2*pi)
  call em3submulcpy &
      (4*npts,cmatr,npts,npts,cmatr0,1,1,ceps_id)
  call em3submulcpy &
      (4*npts,cmatr,npts,npts,cmatr0,2*npts+1,2*npts+1,cmus_id)

  call prinf('tangential component: E(11)=*',id,1)
  call patchmatc(npatches,rpatchpnt, &
      qtriainfo,epatchpnt,ipatchinfo,refineinfo, &
      norder,npols,usout,vsout,umatr,vmatr, &
      ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts, &
      eminter1n,rk_id,info,par7,par8, &
      cmatr0,w,lw,lused,ier)

  call em3submulcpy &
      (4*npts,cmatr,npts,npts,cmatr0,2*npts+1,1,+cems_id)
  call em3submulcpy &
      (4*npts,cmatr,npts,npts,cmatr0,1,2*npts+1,-cems_id)

  info(1)=1
  info(2)=2
  call prinf('tangential component: H(12)=*',id,1)
  call patchmatc(npatches,rpatchpnt, &
      qtriainfo,epatchpnt,ipatchinfo,refineinfo, &
      norder,npols,usout,vsout,umatr,vmatr, &
      ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts, &
      eminter3n,rk_id,info,par7,par8, &
      cmatr0,w,lw,lused,ier)

  call em3submulcpy &
      (4*npts,cmatr,npts,npts,cmatr0,npts+1,1,ceps_id)
  call em3submulcpy &
      (4*npts,cmatr,npts,npts,cmatr0,3*npts+1,2*npts+1,cmus_id)

  call prinf('tangential component: E(12)=*',id,1)
  call patchmatc(npatches,rpatchpnt, &
      qtriainfo,epatchpnt,ipatchinfo,refineinfo, &
      norder,npols,usout,vsout,umatr,vmatr, &
      ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts, &
      eminter1n,rk_id,info,par7,par8, &
      cmatr0,w,lw,lused,ier)

  call em3submulcpy &
      (4*npts,cmatr,npts,npts,cmatr0,3*npts+1,1,+cems_id)
  call em3submulcpy &
      (4*npts,cmatr,npts,npts,cmatr0,npts+1,2*npts+1,-cems_id)

  info(1)=2
  info(2)=1
  call prinf('tangential component: H(21)=*',id,1)
  call patchmatc(npatches,rpatchpnt, &
      qtriainfo,epatchpnt,ipatchinfo,refineinfo, &
      norder,npols,usout,vsout,umatr,vmatr, &
      ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts, &
      eminter3n,rk_id,info,par7,par8, &
      cmatr0,w,lw,lused,ier)

  call em3submulcpy &
      (4*npts,cmatr,npts,npts,cmatr0,1,npts+1,ceps_id)
  call em3submulcpy &
      (4*npts,cmatr,npts,npts,cmatr0,2*npts+1,3*npts+1,cmus_id)

  call prinf('tangential component: E(21)=*',id,1)
  call patchmatc(npatches,rpatchpnt, &
      qtriainfo,epatchpnt,ipatchinfo,refineinfo, &
      norder,npols,usout,vsout,umatr,vmatr, &
      ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts, &
      eminter1n,rk_id,info,par7,par8, &
      cmatr0,w,lw,lused,ier)

  call em3submulcpy &
      (4*npts,cmatr,npts,npts,cmatr0,2*npts+1,npts+1,+cems_id)
  call em3submulcpy &
      (4*npts,cmatr,npts,npts,cmatr0,1,3*npts+1,-cems_id)

  info(1)=2
  info(2)=2
  call prinf('tangential component: H(22)=*',id,1)
  call patchmatc(npatches,rpatchpnt, &
      qtriainfo,epatchpnt,ipatchinfo,refineinfo, &
      norder,npols,usout,vsout,umatr,vmatr, &
      ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts, &
      eminter3n,rk_id,info,par7,par8, &
      cmatr0,w,lw,lused,ier)

  call patchdiag(cmatr0,npts,+2*pi)
  call em3submulcpy &
      (4*npts,cmatr,npts,npts,cmatr0,npts+1,npts+1,ceps_id)
  call em3submulcpy &
      (4*npts,cmatr,npts,npts,cmatr0,3*npts+1,3*npts+1,cmus_id)

  call prinf('tangential component: E(22)=*',id,1)
  call patchmatc(npatches,rpatchpnt, &
      qtriainfo,epatchpnt,ipatchinfo,refineinfo, &
      norder,npols,usout,vsout,umatr,vmatr, &
      ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts, &
      eminter1n,rk_id,info,par7,par8, &
      cmatr0,w,lw,lused,ier)

  call em3submulcpy &
      (4*npts,cmatr,npts,npts,cmatr0,3*npts+1,npts+1,+cems_id)
  call em3submulcpy &
      (4*npts,cmatr,npts,npts,cmatr0,npts+1,3*npts+1,-cems_id)


  id=2
  rk_id=omega*sqrt(cmus(id))*sqrt(ceps(id))
  ceps_id=ceps(id)
  cmus_id=cmus(id)
  cems_id=sqrt(cmus(id))*sqrt(ceps(id))

  info(1)=1
  info(2)=1
  call prinf('tangential component: H(11)=*',id,1)
  call patchmatc(npatches,rpatchpnt, &
      qtriainfo,epatchpnt,ipatchinfo,refineinfo, &
      norder,npols,usout,vsout,umatr,vmatr, &
      ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts, &
      eminter3n,rk_id,info,par7,par8, &
      cmatr0,w,lw,lused,ier)

  call patchdiag(cmatr0,npts,-2*pi)
  call em3submuladd &
      (4*npts,cmatr,npts,npts,cmatr0,1,1,-ceps_id)
  call em3submuladd &
      (4*npts,cmatr,npts,npts,cmatr0,2*npts+1,2*npts+1,-cmus_id)

  call prinf('tangential component: E(11)=*',id,1)
  call patchmatc(npatches,rpatchpnt, &
      qtriainfo,epatchpnt,ipatchinfo,refineinfo, &
      norder,npols,usout,vsout,umatr,vmatr, &
      ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts, &
      eminter1n,rk_id,info,par7,par8, &
      cmatr0,w,lw,lused,ier)

  call em3submuladd &
      (4*npts,cmatr,npts,npts,cmatr0,2*npts+1,1,-cems_id)
  call em3submuladd &
      (4*npts,cmatr,npts,npts,cmatr0,1,2*npts+1,+cems_id)

  info(1)=1
  info(2)=2
  call prinf('tangential component: H(12)=*',id,1)
  call patchmatc(npatches,rpatchpnt, &
      qtriainfo,epatchpnt,ipatchinfo,refineinfo, &
      norder,npols,usout,vsout,umatr,vmatr, &
      ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts, &
      eminter3n,rk_id,info,par7,par8, &
      cmatr0,w,lw,lused,ier)

  call em3submuladd &
      (4*npts,cmatr,npts,npts,cmatr0,npts+1,1,-ceps_id)
  call em3submuladd  &
      (4*npts,cmatr,npts,npts,cmatr0,3*npts+1,2*npts+1,-cmus_id)

  call prinf('tangential component: E(12)=*',id,1)
  call patchmatc(npatches,rpatchpnt, &
      qtriainfo,epatchpnt,ipatchinfo,refineinfo, &
      norder,npols,usout,vsout,umatr,vmatr, &
      ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts, &
      eminter1n,rk_id,info,par7,par8, &
      cmatr0,w,lw,lused,ier)

  call em3submuladd &
      (4*npts,cmatr,npts,npts,cmatr0,3*npts+1,1,-cems_id)
  call em3submuladd &
      (4*npts,cmatr,npts,npts,cmatr0,npts+1,2*npts+1,+cems_id)

  info(1)=2
  info(2)=1
  call prinf('tangential component: H(21)=*',id,1)
  call patchmatc(npatches,rpatchpnt, &
      qtriainfo,epatchpnt,ipatchinfo,refineinfo, &
      norder,npols,usout,vsout,umatr,vmatr, &
      ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts, &
      eminter3n,rk_id,info,par7,par8, &
      cmatr0,w,lw,lused,ier)

  call em3submuladd &
      (4*npts,cmatr,npts,npts,cmatr0,1,npts+1,-ceps_id)
  call em3submuladd &
      (4*npts,cmatr,npts,npts,cmatr0,2*npts+1,3*npts+1,-cmus_id)

  call prinf('tangential component: E(21)=*',id,1)
  call patchmatc(npatches,rpatchpnt, &
      qtriainfo,epatchpnt,ipatchinfo,refineinfo, &
      norder,npols,usout,vsout,umatr,vmatr, &
      ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts, &
      eminter1n,rk_id,info,par7,par8, &
      cmatr0,w,lw,lused,ier)

  call em3submuladd &
      (4*npts,cmatr,npts,npts,cmatr0,2*npts+1,npts+1,-cems_id)
  call em3submuladd &
      (4*npts,cmatr,npts,npts,cmatr0,1,3*npts+1,+cems_id)

  info(1)=2
  info(2)=2
  call prinf('tangential component: H(22)=*',id,1)
  call patchmatc(npatches,rpatchpnt, &
      qtriainfo,epatchpnt,ipatchinfo,refineinfo, &
      norder,npols,usout,vsout,umatr,vmatr, &
      ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts, &
      eminter3n,rk_id,info,par7,par8, &
      cmatr0,w,lw,lused,ier)

  call patchdiag(cmatr0,npts,-2*pi)
  call em3submuladd &
      (4*npts,cmatr,npts,npts,cmatr0,npts+1,npts+1,-ceps_id)
  call em3submuladd &
      (4*npts,cmatr,npts,npts,cmatr0,3*npts+1,3*npts+1,-cmus_id)

  call prinf('tangential component: E(22)=*',id,1)
  call patchmatc(npatches,rpatchpnt, &
      qtriainfo,epatchpnt,ipatchinfo,refineinfo, &
      norder,npols,usout,vsout,umatr,vmatr, &
      ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts, &
      eminter1n,rk_id,info,par7,par8, &
      cmatr0,w,lw,lused,ier)

  call em3submuladd &
      (4*npts,cmatr,npts,npts,cmatr0,3*npts+1,npts+1,-cems_id)
  call em3submuladd &
      (4*npts,cmatr,npts,npts,cmatr0,npts+1,3*npts+1,+cems_id)


  call prinf('after patchmatc, ier=*',ier,1)

  do i=1,npts
    call patchinfo(xyzs(1,i),xyznorms(1,i), &
        xyztang1s(1,i),xyztang2s(1,i),xyzinfo(1,i))
  enddo


  !
  ! . . . at this point the muller system matrix has been built
  !

  

  ifexterior = 1
  iftest = 1

  !
  ! . . . simple test for Muller solver
  !
  if (iftest .eq. 1) then

    !
    !.. construct the right hand side
    !
    source_cjvec(1)=1
    source_cjvec(2)=1
    source_cjvec(3)=1
    source_cmvec(1)=0
    source_cmvec(2)=0
    source_cmvec(3)=0


    ! place source inside, targ outside
    if( ifexterior .eq. 1 ) then
      source(1)=0.2
      source(2)=-0.1
      source(3)=0.3
      targ(1)=10
      targ(2)=30
      targ(3)=-30
    endif

    ! or place targ inside, source outside
    if( ifexterior .eq. 0 ) then
      targ(1)=0.2
      targ(2)=-0.1
      targ(3)=0.3
      source(1)=10
      source(2)=20
      source(3)=-30
    endif

    ! set the material parameters as need be for the exact solution
    if( ifexterior .eq. 1 ) id=1
    if( ifexterior .eq. 0 ) id=2    
    rk_id = omega*sqrt(cmus(id))*sqrt(ceps(id))
    ceps_id = ceps(id)
    cmus_id = cmus(id)
    cems_id = sqrt(cmus(id))*sqrt(ceps(id))

    call em3getrhs3(rk_id,ceps_id,cmus_id, &
        source_cjvec,source_cmvec,source,xyzinfo,rhs,npts)

    ! flip the sign on the rhs if need be
    if( ifexterior .eq. 0 ) then
      do i=1,npts*4
        rhs(i)=-rhs(i)
      enddo
    endif

    !
    ! call the solver
    !
    call prinf('npts=*',npts,1)

    !
    ! hardcode the gmres solver...
    !
    call prinf('entering cgmres, 4 x npts=*',4*npts,1)

    ngmrec=numit
    if( allocated(w) ) deallocate(w)
    allocate( w( 2*(ngmrec*2+4)*(4*npts) ) )

    call cgmres(ier,4*npts,cmatr, &
        em3multa,par1,par2,rhs,eps,numit, &
        sol,niter,errs,ngmrec,w)
    
    call prinf('after cgmres, ier=*',ier,1)
    call prinf('after cgmres, niter=*',niter,1)
    call prin2('after cgmres, errs=*',errs,niter)


    !
    ! ... check the solution
    !
    print *
    call prin2('targ=*',targ,3)
    call em3soleva2(rk_id,ceps_id,cmus_id, &
        targ,xyzinfo,sol,whts,npts,evec,hvec)

    call em3direva3(rk_id,ceps_id,cmus_id, &
        source_cjvec,source_cmvec,source,targ,evec0,hvec0)

    do i=1,3
      evec1(i)=(evec(i)-evec0(i))/evec0(i)
      hvec1(i)=(hvec(i)-hvec0(i))/hvec0(i)
      evec2(i) = abs(evec1(i))
      hvec2(i) = abs(hvec1(i))
    enddo

    etop = 0
    ebottom = 0
    htop = 0
    hbottom = 0
    do i = 1,3
      etop = etop + abs(evec(i)-evec0(i))**2
      ebottom = ebottom + abs(evec0(i))**2
      htop = htop + abs(hvec(i)-hvec0(i))**2
      hbottom = hbottom + abs(hvec0(i))**2
    enddo

    eerr = sqrt(etop/ebottom)
    herr = sqrt(htop/hbottom)


    print *
    call prin2('from solver, evec=*',evec,6)
    call prin2('directly, evec0=*',evec0,6)
    call prin2('relative error, E=*',evec1,6)
    call prin2('abs relative error, E=*',evec2,6)
    call prin2('norm of abs relative err, E = *', eerr, 1)

    print *
    call prin2('from solver, hvec=*',hvec,6)
    call prin2('directly, hvec0=*',hvec0,6)
    call prin2('abs relative error, H=*',hvec2,6)
    call prin2('norm of abs relative err, H = *', herr, 1)

    stop



  endif
      

        
        


  !
  ! . . . call Muller solver for multiple incoming fields and generate
  !    the corresponding scattering matrix acting on coefficients
  !

  open(unit=72,file=filename_out)

  allocate( ampole(0:nterms,-nterms:nterms) ) 
  allocate( bmpole(0:nterms,-nterms:nterms) ) 

  allocate( sol_cjvecs(3,4*npts) )
  allocate( sol_cmvecs(3,4*npts) )

  if( ifexterior .eq. 1 ) id=1
  if( ifexterior .eq. 0 ) id=2
  rk_id=omega*sqrt(cmus(id))*sqrt(ceps(id))
  ceps_id=ceps(id)
  cmus_id=cmus(id)
  cems_id=sqrt(cmus(id))*sqrt(ceps(id))

  ! turn off usual printing to the screen, setup printing to file
  call prini(0,13)
  call xprini(0,72)
  
  !
  ! ... call Muller solver for multiple incoming fields
  !
  call xprin2('omega=*',omega,1)
  call xprin2('rk=*',rk,2)
  call xprin2('radius=*',radius,1)
  call xprinf('nterms=*',nterms,1)

  !
  ! loop over mp type (a's or b's in the Mie expansions)
  ! loop over the orders
  ! loop over the degrees
  !
  do itype_mp=1,2
    do n=1,nterms
      do m=-n,n

        call xprinf('itype_mp=*',itype_mp,1)
        call xprinf('n=*',n,1)
        call xprinf('m=*',m,1)

        call em3getrhs4(rk_id,ceps_id,cmus_id, &
            itype_mp,n,m,nterms,xyzinfo,rhs,npts)
        if( ifexterior .eq. 0 ) then
          do i=1,npts*4
            rhs(i)=-rhs(i)
          enddo
        endif

        ! use GMRES to solve the system

        call prinf('entering cgmres, 4 x npts=*',4*npts,1)

        !eps=1e-5
        !numit=40
        ngmrec=numit

        if( allocated(w) ) deallocate(w)
        allocate( w( 2*(ngmrec*2+4)*(4*npts) ) )

        print *
        write(*,'(a,i2,a,i2,a,i2)') 'id = ', id, ', n = ', n, ', m = ', m
        call cgmres(ier,4*npts,cmatr, em3multa,par1,par2,rhs, &
            eps,numit, sol,niter,errs,ngmrec,w)

        call prinf('after cgmres, ier=*',ier,1)
        call prinf('after cgmres, niter=*',niter,1)
        call prin2('after cgmres, errs=*',errs,niter)



        !        ccc        call em3soleva(rk,targ,xyzinfo,sol,whts,npts,evec,hvec)

        if( 1.eq.2 ) then
          if( ifexterior .eq. 1 ) id=1
          if( ifexterior .eq. 0 ) id=2
          rk_id=omega*sqrt(cmus(id))*sqrt(ceps(id))
          ceps_id=ceps(id)
          cmus_id=cmus(id)
          cems_id=sqrt(cmus(id))*sqrt(ceps(id))
          call em3soleva2(rk_id,ceps_id,cmus_id, &
              targ,xyzinfo,sol,whts,npts,evec,hvec)

          call prin2('targ=*',targ,3)
          call prin2('evec=*',evec,6)
          call prin2('hvec=*',hvec,6)
        endif

        !
        ! ... evaluate the outgoing multipole expansion
        !
        call em3soleva2a(xyzinfo,sol,whts,npts,sol_cjvecs,sol_cmvecs)
        !call prin2('sol_cjvecs=*',sol_cjvecs,3*npts)
        !call prin2('sol_cmvecs=*',sol_cmvecs,3*npts)

        do i=1,npts
          sol_cjvecs(1,i)=sol_cjvecs(1,i)*whts(i)
          sol_cjvecs(2,i)=sol_cjvecs(2,i)*whts(i)
          sol_cjvecs(3,i)=sol_cjvecs(3,i)*whts(i)
          sol_cmvecs(1,i)=sol_cmvecs(1,i)*whts(i)
          sol_cmvecs(2,i)=sol_cmvecs(2,i)*whts(i)
          sol_cmvecs(3,i)=sol_cmvecs(3,i)*whts(i)
        enddo

        center(1)=0
        center(2)=0
        center(3)=0

        call em3formmp(rk,xyzs,sol_cjvecs,sol_cmvecs, &
            npts,center,ampole,bmpole,nterms)

        ! print the columns of the scattering matrix to file
        call em3sphlin(ampole,nterms,w)
        call xprin2('outgoing ampole=*',w,2*(nterms+1)**2)
        call em3sphlin(bmpole,nterms,w)
        call xprin2('outgoing bmpole=*',w,2*(nterms+1)**2)

      enddo
    enddo
  enddo





end program test_muller






