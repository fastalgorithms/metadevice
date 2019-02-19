subroutine patchmatc0(npatches,patchpnt,par1,par2,par3,par4, &
    norder,npols,us,vs,umatr,vmatr, &
    ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts, &
    interact,par5,par6,par7,par8, &
    cmatr, ier)
  implicit double precision (a-h,o-z)
  external patchpnt,interact
  double precision :: us(npols),vs(npols)
  dimension ixyzs(2,1),xyzs(3,1),xyznorms(3,1)
  dimension xyztang1s(3,1),xyztang2s(3,1)
  double complex :: cmatr(npts,npts)

  integer :: omp_get_thread_num
  double precision :: rad(1000000)
  double complex :: tmatr(100000), work(10000000)

  !
  ! ... construct the diagonal (self interaction) blocks
  !

  !
  ! get the quadrature
  !
  lrad = 1000000
  nquad = 12
  call radial_init(jer0, nquad, rad, lrad, lkeep)

  !$omp parallel do default(shared) &
  !$omp private(i, ii, ipols, tmatr, ier)
  do i=1,npatches
    ii = ixyzs(1,i)
    ipols = ixyzs(2,i)
    ier = 0
    call patchmatc_dd(i, ipols, i, ipols, &
        npatches, patchpnt, par1, par2, par3, par4, &
        norder, npols, us, vs, umatr, vmatr, rad, &
        ixyzs, xyzs, xyznorms, xyztang1s, xyztang2s, npts, &
        interact, par5, par6, par7, par8, &
        tmatr, ier)
    call patchsubcpy(npts, cmatr, tmatr, ii, ipols, ii, ipols)
  enddo
  !$omp end parallel do

  
  lw = 10000000
  lused = 0

  !
  ! ... construct the off-diagonal blocks
  !
  !$omp parallel do default(shared) &
  !$omp private(i, ii, ipols, j, jj, jpols, work, lused, tmatr, ier)
  do j=1,npatches
    do i=1,npatches

      ! ... (i,j), i index - target, j index - source
      if (i .ne. j) then
        ii=ixyzs(1,i)
        jj=ixyzs(1,j)
        ipols=ixyzs(2,i)
        jpols=ixyzs(2,j)
        call patchmatc_od(i,ipols,j,jpols, &
            npatches,patchpnt,par1,par2,par3,par4, &
            norder,npols,us,vs,umatr,vmatr, &
            ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts, &
            interact,par5,par6,par7,par8,&
            tmatr, work, lw, lused, ier)
        call patchsubcpy(npts, cmatr, tmatr, ii, ipols, jj, jpols)
      endif

    enddo

  enddo
  !$omp end parallel do
  



  return
end subroutine patchmatc0




        

subroutine patchmatc_dd(ipatch,ipols,jpatch,jpols,&
    npatches,patchpnt,par1,par2,par3,par4,&
    norder, npols,us,vs,umatr,vmatr, rad, &
    ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts, &
    interact,par5,par6,par7,par8, &
    tmatr, ier)
  implicit double precision (a-h,o-z)
  !
  !   ... generate the diagonal block of interaction matrix,
  !   self-interaction
  !
  integer :: ixyzs(2,1)
  double precision :: us(npols), vs(npols), rad(*)
  double precision :: umatr(npols,npols), vmatr(npols,npols)
  double precision :: xyzs(3,1), xyznorms(3,1)
  double precision :: xyztang1s(3,1),xyztang2s(3,1)
  double complex :: tmatr(ipols,jpols)

  external patchpnt, interact
  
  !double precision :: rad(1000000)
  double precision :: xyz(3), dxyzduv(3,10)
  double precision :: xyznorm(30), xyztang1(30), xyztang2(30)
  double precision :: targinfo(120), info(200)
  double precision :: xpar1(200), xpar2(200)
  double precision :: xs(10000), ys(10000), ws(10000)
  double precision :: vert1a(2,10)

  double complex :: cvals(1000), coefs(1000)


  ii = ixyzs(1,ipatch)



  ! print *, 'ipatch = ', ipatch
  ! print *, 'ipols = ', ipols 
  ! print *, 'jpatch = ', jpatch
  ! print *, 'npatches = ', npatches
  ! print *, 'norder = ', norder 
  ! print *, 'npols = ', npols

  ! print *, 'us = '
  ! do i = 1,npols
  !   print *, us(i)
  ! end do

  ! print *, 'vs = ' 
  ! do i = 1,npols
  !   print *, vs(i)
  ! end do

  ! print *, 'umatr = ' 
  ! do i = 1,npols
  !   do j = 1,npols
  !     !print *, umatr(i,j)
  !   end do
  ! end do

  ! !print *, 'vmatr = ' 
  ! do i = 1,npols
  !   do j = 1,npols
  !     !print *, vmatr(i,j)
  !   end do
  ! end do

  ! print *, 'npts = ', npts 
  

  if( ipols .ne. npols ) then
    write(*,*) 'ipols .ne. npols'
  endif
  if( jpols .ne. npols ) then
    write(*,*) 'jpols .ne. npols'
  endif

  
  
  
  ! get the quadrature
  
  ! lrad = 1000000
  ! allocate(rad(lrad))
  ! norder12 = 12
  ! print *, 'calling radial_init . . .'
  ! call radial_init(jer0, norder12, rad, lrad, lkeep)
  ! print *, 'loaded rad from radial_init'
  
  !
  !  ... on i-th triangle integrate all basis functions multiplied 
  !  by interaction function at the target point us(i),vs(i)
  !
  xpar1(1) = norder
  xpar1(2) = npols
  xpar1(3) = ipatch

  do i=1,npols
    !print *, 'scanning over npols, i = ', i

    ! ... first, initialize function to be integrated
    call patchinfo(xyzs(1,ii+i-1),xyznorms(1,ii+i-1), &
        xyztang1s(1,ii+i-1),xyztang2s(1,ii+i-1),targinfo)

    do j=1,12
      xpar2(j) = targinfo(j)
    enddo

    !
    ! Jim Bremer's quadratures
    !
    vert1a(1,1) = 0
    vert1a(2,1) = 0
    vert1a(1,2) = 1
    vert1a(2,2) = 0
    vert1a(1,3) = 0
    vert1a(2,3) = 1

    u=us(i)
    v=vs(i)
    call patchgeo(patchpnt,ipatch,u,v, &
        par1,par2,par3,par4, &
        xyz,dxyzduv,ds,xyznorm,xyztang1,xyztang2)

    x0=us(i)
    y0=vs(i)
    call self_quadrature_new(ier, rad, vert1a, x0, y0, dxyzduv, ns, &
        xs,ys,ws)
    !print *, 'after self_quad_new, ier = ', ier

    do j=1,npols
      coefs(j)=0
    enddo

    do k=1,ns
      call patchfun3(xs(k),ys(k),patchpnt,par1,par2,par3,par4, &
          interact,par5,par6,par7,par8,xpar1,xpar2,cvals)
      do j=1,npols
        coefs(j)=coefs(j)+ws(k)*cvals(j)
      enddo
    enddo

    !
    ! ...finally, convert the linear form of integral values to the
    ! pointwise interation matrix, we will need umatr and vmatr for
    ! this operation
    !
    call patchcoefs2cvals(npols, umatr, vmatr, coefs, cvals)

    
    do j=1,npols
      tmatr(i,j)=cvals(j)
    enddo

  enddo


  return
end subroutine patchmatc_dd
  




subroutine patchmatc_od(ipatch,ipols,jpatch,jpols, &
    npatches,patchpnt,par1,par2,par3,par4, &
    norder,npols,us,vs,umatr,vmatr, &
    ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts, &
    interact,par5,par6,par7,par8, &
    tmatr,w,lw,lused,ier)
  implicit real *8 (a-h,o-z)
  !
  ! ... generate the off-diagonal block of interaction matrix
  !
  external patchpnt,interact
  dimension us(1),vs(1)
  dimension ixyzs(2,1),xyzs(3,1),xyznorms(3,1)
  dimension xyztang1s(3,1),xyztang2s(3,1)

  dimension xyz(3),dxyzduv(3,20)
  dimension xyznorm(30),xyztang1(30),xyztang2(30)
  dimension targinfo(120),info(200)
  dimension xpar1(200),xpar2(200)

  dimension w(1)
  complex *16 tmatr(ipols,jpols)
  complex *16 cvals(1000),coefs(1000)

  dimension vert1(2),vert2(2),vert3(2)
  external patchfun3

  ii=ixyzs(1,ipatch)
  jj=ixyzs(1,jpatch)

  !
  ! ... construct one off-diagonal block via collocation       
  ! ... (i,j), i index - target, j index - source
  !
  if( ipols .ne. npols ) then
    write(*,*) 'ipols .ne. npols'
  endif
  if( jpols .ne. npols ) then
    write(*,*) 'jpols .ne. npols'
  endif

  do  i=1,npols
    !c
    !c       ... on j-th triangle integrate all basis functions multiplied 
    !c       by interaction function at the target point us(i),vs(i)
    !c
    !c       ... first, initialize function to be integrated
    !c
    xpar1(1)=norder
    xpar1(2)=npols
    xpar1(3)=jpatch

    call patchinfo(xyzs(1,ii+i-1),xyznorms(1,ii+i-1), &
        xyztang1s(1,ii+i-1),xyztang2s(1,ii+i-1),targinfo)

    do j=1,12
      xpar2(j)=targinfo(j)
    enddo
    !
    ! ... then, call adaptive gaussian integration routine 
    !       
    m=12
    eps=1d-12

    !iquadtype=1
    !iquadtype=2

    vert1(1)=0
    vert1(2)=0
    vert2(1)=1
    vert2(2)=0
    vert3(1)=0
    vert3(2)=1

    !if( iquadtype .eq. 1 )        
    nrec=20
    call c28triaadam(ier,vert1,vert2,vert3,patchfun3,npols, &
        patchpnt,par1,par2,par3,par4, &
        interact,par5,par6,par7,par8,xpar1,xpar2, &
        m,eps,coefs,maxrec,numfunev,w,nrec,jpatch,targinfo,info)

    !if( iquadtype .eq. 2 )
    !     call c29triaadam(ier,vert1,vert2,vert3,patchfun3,npols,
    !     patchpnt,par1,par2,par3,par4,
    !     interact,par5,par6,par7,par8,xpar1,xpar2,
    !     m,eps,coefs,maxrec,numfunev,w)

    if( ier .eq. 8 ) then
      write(*,*) 'maximum recursion depth of 200 has been reached'
      write(*,*) 'abort'
      stop
    endif

    ! 
    ! ... finally, convert the linear form of integral values to the
    ! pointwise interation matrix, we will need umatr and vmatr for
    ! this operation
    ! 
    call patchcoefs2cvals(npols,umatr,vmatr,coefs,cvals)
    do  j=1,npols
      tmatr(i,j)=cvals(j)
    end do

  end do


  return
end subroutine patchmatc_od
