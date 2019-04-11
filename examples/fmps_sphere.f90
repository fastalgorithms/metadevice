program fmps_sphere

  implicit double precision (a-h,o-z)

  !, outgt(:,:)

  ! double complex, allocatable :: ceps(:), cmu(:), zk(:)
  ! double complex, allocatable :: eveci(:,:), hveci(:,:), &
  !     evect(:,:), hvect(:,:)
  ! double complex, allocatable :: pynt(:,:), pyns(:,:)

  double precision, allocatable :: centers(:,:), radius(:), angles(:)
  double precision, allocatable :: rnodes(:,:), whts(:)
  double precision, allocatable :: targs(:,:)

  double complex :: ceps0, cmu0, zk0
  double complex :: rkvec(3), epol(3)
  double complex, parameter :: ima = (0,1)
  double complex, allocatable :: rmatr(:,:)
  double complex, allocatable :: evec1(:,:), hvec1(:,:)
  double complex, allocatable :: aimpole(:,:,:), bimpole(:,:,:)
  double complex, allocatable :: aompole(:,:,:), bompole(:,:,:)
  double complex, allocatable :: ampvec(:,:), bmpvec(:,:)
  double complex, allocatable :: rhs(:), sol(:), abvec(:)

  character *256 :: scatfile       

  ! !       double complex, allocatable :: raa(:,:), rab(:,:), & 
  ! !                           rba(:,:), rbb(:,:)
  ! double complex, allocatable :: ampout(:,:), bmpout(:,:), & 
  !     x(:),y(:)

  ! double precision, allocatable :: errs(:)



  call prini(6,13)

  done = 1
  pi = 4*atan(done)

  !
  ! load pre-computed scattering matrix
  !
  !       scatfile = '../fmps-3.1/src/muller/scatmatr_muller.txt'
  !       scatfile = './scatmatr_high/scatmatr_5_3.txt'
  scatfile = './scatmatrs/scatmatr_6_5.txt'
  print *, 'loading file: ', trim(scatfile)

  nterms = 5
  ncoefs = (nterms+1)**2
  allocate(rmatr(2*ncoefs,2*ncoefs))
  call load_scatmatr(scatfile, nterms, omega, radius0, rmatr, ier)

  call prinf('after loading scatmatr, nterms = *', nterms, 1)
  call prin2('after loading scatmatr, omega = *', omega, 1)
  call prin2('after loading scatmatr, radius = *', radius0, 1)
  !call prin2('after loading scatmatr, rmatr = *', rmatr, 4*ncoefs**2)
  call prinf('after loading scatmatr, ier = *', ier, 1)
  
  if (ier .ne. 0) then
    call prinf('error loading scattering matrix, ier = *', ier, 1)
    stop
  end if
  
  wavelength = 2*pi/omega
  call prin2('wavelength = *', wavelength, 1)

  ceps0 = 1
  cmu0 = 1
  zk0 = omega*sqrt(ceps0)*sqrt(cmu0)
  call prin2('eps0 = *', ceps0, 2)
  call prin2('mu0 = *', cmu0, 2)
  call prin2('zk0 = *', zk0, 2)


  !
  ! setup geometry confiuration
  !
  
  ns = 5
  nspheres = (2*ns+1)**2


  allocate(centers(3,nspheres))
  allocate(radius(nspheres))
  allocate(angles(nspheres))
  
  do i=1,nspheres
    radius(i) = radius0
  end do

  call prin2('radii = ', radius, nspheres)

  ! compute the centers on a grid, making sure the separation is
  ! commensurate with radius0
  kk=1
  do i=-ns,ns
    do j=-ns,ns
      centers(1,kk) = i*2*radius0
      centers(2,kk) = j*2*radius0
      centers(3,kk) = 0
      kk=kk+1
    end do
  end do

  call prin2('centers = *', centers, 3*nspheres)
  
  ! set some random angles of rotation
  do i=1,nspheres
    angles(i) = acos(cos(1.34d0*done*i**4))
  end do

  call prin2('angles = *', angles, nspheres)
  
  
  !
  ! generate quadrature nodes on the sphere for computing
  ! coefficients in spherical harmonic expansions
  !
  itype = 1
  nquad = nterms
  nphi = 2*nquad+1
  ntheta = nquad+1

  nnodes = nphi*ntheta
  allocate(rnodes(3,nnodes))
  allocate(whts(nnodes))
  call e3fgrid(itype,nterms,nphi,ntheta,rnodes,whts,nnodes)


  call prin2('quad nodes = *', rnodes, 3*nnodes)
  call prin2('quad whts = *', whts, nnodes)

  
  !
  ! compute the local expansions for an incoming planewave, one
  ! inclusion at a time
  !
  !
  allocate(evec1(3,nnodes))
  allocate(hvec1(3,nnodes))
  allocate(targs(3,nnodes))

  allocate(aimpole(0:nterms,-nterms:nterms,nspheres))
  allocate(bimpole(0:nterms,-nterms:nterms,nspheres))
  allocate(aompole(0:nterms,-nterms:nterms,nspheres))
  allocate(bompole(0:nterms,-nterms:nterms,nspheres))

  allocate( ampvec(0:nterms,-nterms:nterms) )
  allocate( bmpvec(0:nterms,-nterms:nterms) )

  !
  ! setup parameters for the planewave
  !
  do i=1,3
    rkvec(i) = 0
    epol(i) = 0
  end do

  rkvec(3) = -zk0
  epol(1) = 1.0d0

  !
  ! evaluate the planewaves and compute local expansions
  !
  do i=1,nspheres

    do j=1,nnodes
      do k = 1,3
        targs(k,j) = rnodes(k,j)*radius(i) + centers(k,i)
      end do
    end do

    do j = 1,nnodes
      call emplanearb(rkvec, epol, targs(1,j), evec1, hvec1)
    end do

    ! form local expansion coefs for incoming field at each center
    call em3ehformta(zk0, centers(1,i), radius(i), &
        evec1, hvec1, rnodes, whts, nphi, ntheta, &
        centers(1,i), ampvec, bmpvec, nterms) 

    do n=0,nterms
      do m=-n,n
        aimpole(n,m,i) = ampvec(n,m)
        bimpole(n,m,i) = bmpvec(n,m)
      end do
    end do


  end do


  !
  ! form the rhs 
  !
  
  n = 2*ncoefs*nspheres
  allocate(rhs(n), sol(n), abvec(n))

  call multi_sphlin(aimpole, bimpole, nterms, nspheres, abvec)

  ! apply 1scattering matrix to incoming field: 
  call scatmatr_multa_rot(rmatr, angles, nterms, nspheres, &
      abvec, rhs)

  do i=1,n
    rhs(i)=-rhs(i)
  end do

!   numit = 200
!   eps = 1e-10
!   allocate(errs(numit))



! !!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!! gmres solver
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   t1=second()


!   call cgmres_fmps(zk0,radius,centers,angles,nspheres, &
!       rnodes,whts,nphi,ntheta,nterms,rmatr, &
!       rhs,eps,numit,sol,niter,errs)


!   call multi_linsph(sol,nterms,nspheres,aompole,bompole)


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!! evaluate outgoing fields:
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   nx=200
!   ntgt = nx**2
!   allocate(outgt(3,ntgt))

!   !! x = [lx,rx]
!   lx = -1500
!   rx = 1500
!   hx = (rx-lx)/nx
!   f = 300

!   kk=1
!   do i=1,nx
!     do j=1,nx
!       outgt(1:3,kk)=(/lx+i*hx,lx+j*hx,f/)
!       kk=kk+1
!     end do
!   end do


! !!!!!    single point test case:       
!   !       ntgt = 1
!   !       allocate(outgt(3,ntgt))
!   !       outgt(1:3,1)=(/0.0d0,0.0d0,100d0/)

!   allocate(eveci(3,ntgt))
!   allocate(hveci(3,ntgt))
!   allocate(evect(3,ntgt))
!   allocate(hvect(3,ntgt))

!   !! evaluate incoming fields on target grid points:
!   call emplanearbtargeval(rkvec,epol, &
!       ntgt,outgt,eveci,hveci)

!   !! evaluate scattered fields from each sphere:
!   call em3dmpoletargeval(nspheres,nterms,ncoefs,omega, &
!       ceps0,cmu0,centers,centers,aompole,bompole, &
!       ntgt,outgt,evect,hvect)




! !!!! poynting vector of scattered fields:
!   allocate(pyns(3,ntgt))
!   allocate(pynt(3,ntgt))

!   do i=1,ntgt
!     call crossprod_cmplx(evect(1,i),hvect(1,i),pyns(1,i))
!   end do


! !!! total field = incoming + scattered 
!   do j=1,ntgt
!     do l=1,3
!       evect(l,j)=evect(l,j)+eveci(l,j)
!       hvect(l,j)=hvect(l,j)+hveci(l,j)         
!     end do

!     call crossprod_cmplx(evect(1,j),hvect(1,j),pynt(1,j))

!   end do



!   t2=second()
!   print *, 'FMPS timing=', t2-t1, ' seconds'




! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
! !!!!!!   write to file:
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

!   ir = 7
!   open(unit=ir, file='xgrid_121_5.txt')
!   do i=1,nx
!     write(ir,*) lx+i*hx 
!   end do
!   close(ir)


!   ir = 8 
!   open(unit=ir,file='total_E_121_5.txt')       
!   do j=1,ntgt
!     do l=1,3
!       write(ir,*) dreal(evect(l,j)), dimag(evect(l,j))
!     end do
!   end do
!   close(ir)     


!   ir = 9 
!   open(unit=ir,file='total_H_121_5.txt')

!   do j=1,ntgt
!     do l=1,3
!       write(ir,*) dreal(hvect(l,j)), dimag(hvect(l,j))
!     end do
!   end do
!   close(ir)     


!   ir = 10 
!   open(unit=ir,file='pynt_121_5.txt')

!   do j=1,ntgt
!     do l=1,3
!       write(ir,*) dreal(pynt(l,j)), dimag(pynt(l,j))
!     end do
!   end do
!   close(ir)     



!   ir = 11
!   open(unit=ir,file='pyns_121_5.txt')

!   do j=1,ntgt
!     do l=1,3
!       write(ir,*) dreal(pyns(l,j)), dimag(pyns(l,j))
!     end do
!   end do
!   close(ir)     




end program fmps_sphere






subroutine load_scatmatr(filename, nterms, omega, radius, rmatr, ier)
  implicit double precision (a-h,o-z)

  !
  ! cleaner version for loading the file
  ! each itype,n,m: one column in the scatrmatr
  !
  ! Input
  !   filename -
  ! Output
  !
  !
  
  character(*) :: filename
  character (len=72) :: str
  double complex :: rmatr(2*(nterms+1)**2,2*(nterms+1)**2)
  double complex :: zk

  ir = 777
  open (unit = ir, file = trim(filename))

  read (ir,*) str
  read (ir,*) omega

  read (ir,*) str
  read (ir,*) a,b
  zk = dcmplx(a,b)

  read (ir,*) str
  read (ir,*) radius

  read (ir,*) str
  read (ir,*) nn

  if (nterms.ne.nn) then
    print *, 'nterms not matched!!'
    ier = 1
    stop
  end if
  
  nmp = (nterms+1)**2


  kk=2
  do i=1,2
    do n=1,nterms
      do m=-n,n
        read (ir,*) str
        read (ir,*) itype

        read (ir,*) str
        read (ir,*) ni

        read (ir,*) str
        read (ir,*) mi

        if ( (itype.ne.i).or.(ni.ne.n).or.(mi.ne.m) ) then
          ier = 1
          print *, 'itype=', itype, 'i=', i
          print *, 'ni=', ni, 'n=', n
          print *, 'mi=', mi, 'm=', m

          return
        end if

        !! one column in the scattering matrix
        !! ampole:
        read (ir,*) str
        do j = 1,nmp
          read (ir,*) x,y
          rmatr(j,kk)=dcmplx(x,y)
        end do

        !! bmpole:
        read (ir,*) str
        do j = 1,nmp
          read (ir,*) x,y
          rmatr(j+nmp,kk)=dcmplx(x,y)
        end do

        kk = kk+1

      end do
    end do
  end do


  close(ir)
  ier = 0

  return
end subroutine load_scatmatr









subroutine crossprod_cmplx(x,y,z)
  implicit double precision (a-h,o-z)
  double complex :: x(3),y(3),z(3)

  ! z = x \cross y

  z(1)=x(2)*y(3)-x(3)*y(2)
  z(2)=x(3)*y(1)-x(1)*y(3)
  z(3)=x(1)*y(2)-x(2)*y(1)

  return
end subroutine crossprod_cmplx

