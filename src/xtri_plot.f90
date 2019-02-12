
!
!
! xtri plotting library (c) 2018 Mike O'Neil
! Contact: oneil@cims.nyu.edu
!
! This file contains plotting routines to accompany the xtri
! surface triangulation routines
!
!


subroutine xtri_vtk_surf(iw, ntri, xeval, &
    par1, par2, par3, par4, kover, title)
  implicit real *8 (a-h,o-z)
  real *8 :: xtriainfo(3,1)
  character (len=*) :: title
  external xeval

  real *8 :: uvs(2,1000), xyz(10), dxyzduv(3,20)
  real *8 :: usout(1000), vsout(1000), umatr(10000), vmatr(10000)
  real *8 :: wsout(1000)
  real *8, allocatable :: sigma(:,:)

  !
  ! Dump out a .vtk ASCII file for reading in paraview or other
  !
  ! Merely plot the triangulated surface using the triangle evalation
  ! routine xeval, which should have the calling sequence:
  !
  !      xeval(itri, u, v, xyz, dxyzduv, par1, par2, par3, par4)
  !
  ! This calling sequence is consistent with legacy "patchmatc"
  ! codes. When using the xtri data structure, usually
  ! xeval=xtri_eval, and the calling sequence and parameters will be:
  !
  !      xtri_eval(itri, u, v, xyz, dxyzduv, korder, xtriainfo, &
  !          par3, par4)
  !
  ! Each triangle is oversampled kover times into flat triangles.
  ! Setting kover=1 results in 4 times as many triangles, kover=2
  ! results in 16 times as many triangles, etc.
  !
  ! input:
  !   iw - plot number
  !   ntri - number of triangles
  !   xeval - subroutine evaluating points on the surface
  !   par1, par2, par3, par4 - parameters for xeval
  !   kover - number of times to oversample each triangle, keep in
  !      mind this results in a factor of 4 triangles each time
  !   title - title of the plot, example: 'plot title'
  !
  !

  vmax = -1000
  vmin = 1000

  itype=1
  norder = 1
  call ortho2siexps(itype, norder, npols, usout, vsout, &
      umatr, vmatr, wsout)


  allocate(sigma(npols,ntri))

  !
  ! evaluate z height
  !
  do itri = 1,ntri
    do j = 1,npols
      call xeval(itri, usout(j), vsout(j), xyz, dxyzduv, &
          par1, par2, par3, par4)
      sigma(j,itri) = xyz(3)
      if (xyz(3) .gt. vmax) vmax = xyz(3)
      if (xyz(3) .lt. vmin) vmin = xyz(3)
    end do
  end do

  !
  ! now call the regular plotting routine...
  !
  call xtri_vtk_plot(iw, ntri, xeval, par1, par2, &
      par3, par4, norder, sigma, kover, title)

  return
end subroutine xtri_vtk_surf





subroutine xtri_vtk_surf2(iw, ntri, xeval, &
    par1, par2, par3, par4, kover, title)
  implicit real *8 (a-h,o-z)
  real *8 :: xtriainfo(3,1)
  character (len=*) :: title
  external xeval

  real *8 :: uvs(2,1000), xyz(10), dxyzduv(3,20)
  real *8 :: usout(1000), vsout(1000), umatr(10000), vmatr(10000)
  real *8 :: wsout(1000)
  real *8, allocatable :: sigma(:,:), qtria(:,:,:)

  !
  ! Dump out a .vtk ASCII file for reading in paraview or other
  !
  ! Merely plot the triangulated surface using the triangle evalation
  ! routine xeval, which should have the calling sequence:
  !
  !      xeval(itri, u, v, xyz, dxyzduv, par1, par2, par3, par4)
  !
  ! This calling sequence is consistent with legacy "patchmatc"
  ! codes. When using the xtri data structure, usually
  ! xeval=xtri_eval, and the calling sequence and parameters will be:
  !
  !      xtri_eval(itri, u, v, xyz, dxyzduv, korder, xtriainfo, &
  !          par3, par4)
  !
  ! Each triangle is oversampled kover times into flat triangles.
  ! Setting kover=1 results in 4 times as many triangles, kover=2
  ! results in 16 times as many triangles, etc.
  !
  ! input:
  !   iw - plot number
  !   ntri - number of triangles
  !   xeval - subroutine evaluating points on the surface
  !   par1, par2, par3, par4 - parameters for xeval
  !   kover - number of times to oversample each triangle, keep in
  !      mind this results in a factor of 4 triangles each time
  !   title - title of the plot, example: 'plot title'
  !
  !

  vmax = -1000
  vmin = 1000

  itype=1
  norder = 1
  call ortho2siexps(itype, norder, npols, usout, vsout, &
      umatr, vmatr, wsout)


  allocate(sigma(npols,ntri))

  !
  ! evaluate z height
  !
  do itri = 1,ntri
    do j = 1,npols
      call xeval(itri, usout(j), vsout(j), xyz, dxyzduv, &
          par1, par2, par3, par4)
      sigma(j,itri) = xyz(3)
      if (xyz(3) .gt. vmax) vmax = xyz(3)
      if (xyz(3) .lt. vmin) vmin = xyz(3)
    end do
  end do

  !
  ! now call the regular plotting routine...
  !
  call xtri_vtk_plot2(iw, ntri, xeval, par1, par2, &
      par3, par4, norder, sigma, kover, title)

  return
end subroutine xtri_vtk_surf2





subroutine xtri_vtk_plot(iw, ntri, xeval, par1, par2, &
    par3, par4, norder, sigma, kover, title)
  implicit real *8 (a-h,o-z)
  real *8 :: sigma(*)
  character (len=*) :: title
  external xeval

  real *8 :: umatr(10000), vmatr(10000)
  real *8 :: usout(10000), vsout(10000), wsout(10000)
  real *8, allocatable :: sigmaout(:,:), xtriout(:,:,:)

  !
  ! Plot the real-valued function sigma on eaah of the triangles using
  ! the triangle evalation routine xeval, which should have the
  ! calling sequence:
  !
  !      xeval(itri, u, v, xyz, dxyzduv, par1, par2, par3, par4)
  !
  ! This calling sequence is consistent with legacy "patchmatc"
  ! codes. When using the xtri data structure, usually
  ! xeval=xtri_eval, and the calling sequence and parameters will be:
  !
  !      xtri_eval(itri, u, v, xyz, dxyzduv, korder, xtriainfo, &
  !          par3, par4)
  !
  ! Each triangle (and function sigma) is oversampled kover times into
  ! flat triangles.  Setting kover=1 results in 4 times as many
  ! triangles, kover=2 results in 16 times as many triangles, etc.
  ! contruct a script to plot the function sigma on the triangles
  !

  len = ntri*(4**kover)
  allocate(sigmaout(3,len))
  allocate(xtriout(3,3,len))

  itype=1
  call ortho2siexps(itype, norder, npols, usout, vsout, &
      umatr, vmatr, wsout)

  call xtri_flatten(ntri, xeval, par1, par2, &
    par3, par4, norder, sigma, npols, umatr, &
    kover, ntriout, xtriout, sigmaout)

  m = 1
  call xtri_vtk_flat_scalars(iw, ntriout, xtriout, m, sigmaout, &
      title)

  return
end subroutine xtri_vtk_plot





subroutine xtri_vtk_plot2(iw, ntri, xeval, par1, par2, &
    par3, par4, norder, sigma, kover, title)
  implicit real *8 (a-h,o-z)
  real *8 :: sigma(*)
  character (len=*) :: title
  external xeval

  real *8 :: umatr(10000), vmatr(10000)
  real *8 :: usout(10000), vsout(10000), wsout(10000)
  real *8, allocatable :: sigmaout(:,:), xtriout(:,:,:)

  !
  ! Plot the real-valued function sigma on eah of the triangles using
  ! the triangle evalation routine xeval, which should have the
  ! calling sequence:
  !
  !      xeval(itri, u, v, xyz, dxyzduv, par1, par2, par3, par4)
  !
  ! THIS ROUTINE PLOTS VTK QUADRATIC TRIANGLES
  !
  ! This calling sequence is consistent with legacy "patchmatc"
  ! codes. When using the xtri data structure, usually
  ! xeval=xtri_eval, and the calling sequence and parameters will be:
  !
  !      xtri_eval(itri, u, v, xyz, dxyzduv, korder, xtriainfo, &
  !          par3, par4)
  !
  ! Each triangle (and function sigma) is oversampled kover times into
  ! flat triangles.  Setting kover=1 results in 4 times as many
  ! triangles, kover=2 results in 16 times as many triangles, etc.
  ! contruct a script to plot the function sigma on the triangles
  !

  len = ntri*(4**kover)
  allocate(sigmaout(6,len))
  allocate(xtriout(3,6,len))

  itype=1
  call ortho2siexps(itype, norder, npols, usout, vsout, &
      umatr, vmatr, wsout)

  call xtri_quadratic(ntri, xeval, par1, par2, &
    par3, par4, norder, sigma, npols, umatr, &
    kover, ntriout, xtriout, sigmaout)

  m = 1
  call xtri_vtk_quadratic_scalars(iw, ntriout, xtriout, m, sigmaout, &
      title)

  return
end subroutine xtri_vtk_plot2





subroutine xtri_quadratic(ntri, xeval, par1, par2, &
    par3, par4, norder, sigma, npols, umatr, &
    kover, ntriout, xtriout, sigmaout)
  implicit real *8 (a-h,o-z)
  real *8 :: sigma(*), xtriout(3,6,1)
  real *8 :: sigmaout(6,*), umatr(npols,npols)
  external xeval

  real *8 :: uv0(2,6)
  real *8 :: xyz(3), dxyzduv(3,20)
  real *8 :: coefs(10000), pols(10000)
  real *8, allocatable :: uvs(:,:,:), centers(:,:)

  !
  ! Oversample each curvilinear triangle to QUADRATIC ones, korder times,
  ! and evaluate sigma at the nodes of each triangle - this routine
  ! is mainly for plotting purposes.
  !
  ! Input:
  !   ntri - number of original triangles
  !   xeval - subroutine evaluating the triangles, must have the
  !       calling sequence:
  !
  !          xeval(itri, u, v, xyz, dxyzduv, par1, par2, par3, par4)
  !   norder - order of sigma on each triangle
  !   npols - number of point at which sigma is sampled on each
  !       triangle
  !   umatr - matrix converting values of sigma to expansion coefs
  !   kover - number of times to oversample each triangle
  !
  ! Output:
  !   xtriout - output quadratic triangles (for VTK):
  !
  !           2
  !           | \
  !           5  4
  !           |    \
  !           0--3--1
  !
  !   sigmaout - values of sigma at the corners of each triangle in
  !       xtriout
  !
  !

  neach = 4**kover
  allocate(uvs(2,6,neach))
  !allocate(centers(2,neach))

  k1 = 2
  call xtri_uvnodes(k1, kout, uvs)
  ntot = 1

  !call prin2('quadratic uvs = *', uvs, 2*6)
  !stop

  !
  ! if necessary, recursively subdivide the triangle - first construct
  ! all the uv points
  !
  if (kover .gt. 0) then

    do i = 1,kover

      ltot = ntot

      do j = 1,ltot

        do k = 1,6
          uv0(1,k) = uvs(1,k,j)
          uv0(2,k) = uvs(2,k,j)
        end do

        ! create 1st sub-triangle
        uvs(1,1,j) = uv0(1,1)
        uvs(2,1,j) = uv0(2,1)
        uvs(1,2,j) = uv0(1,4)
        uvs(2,2,j) = uv0(2,4)
        uvs(1,3,j) = uv0(1,6)
        uvs(2,3,j) = uv0(2,6)
        uvs(1,4,j) = (uv0(1,1)+uv0(1,4))/2
        uvs(2,4,j) = (uv0(2,1)+uv0(2,4))/2
        uvs(1,5,j) = (uv0(1,4)+uv0(1,6))/2
        uvs(2,5,j) = (uv0(2,4)+uv0(2,6))/2
        uvs(1,6,j) = (uv0(1,1)+uv0(1,6))/2
        uvs(2,6,j) = (uv0(2,1)+uv0(2,6))/2

        ! create 2nd sub-triangle
        ntot = ntot + 1
        uvs(1,1,ntot) = (uv0(1,4) + uv0(1,4))/2
        uvs(2,1,ntot) = (uv0(2,4) + uv0(2,4))/2
        uvs(1,2,ntot) = (uv0(1,2) + uv0(1,2))/2
        uvs(2,2,ntot) = (uv0(2,2) + uv0(2,2))/2
        uvs(1,3,ntot) = (uv0(1,5) + uv0(1,5))/2
        uvs(2,3,ntot) = (uv0(2,5) + uv0(2,5))/2
        uvs(1,4,ntot) = (uv0(1,4) + uv0(1,2))/2
        uvs(2,4,ntot) = (uv0(2,4) + uv0(2,2))/2
        uvs(1,5,ntot) = (uv0(1,2) + uv0(1,5))/2
        uvs(2,5,ntot) = (uv0(2,2) + uv0(2,5))/2
        uvs(1,6,ntot) = (uv0(1,4) + uv0(1,5))/2
        uvs(2,6,ntot) = (uv0(2,4) + uv0(2,5))/2

        ! create 3rd sub-triangle
        ntot = ntot + 1
        uvs(1,1,ntot) = (uv0(1,6) + uv0(1,6))/2
        uvs(2,1,ntot) = (uv0(2,6) + uv0(2,6))/2
        uvs(1,2,ntot) = (uv0(1,5) + uv0(1,5))/2
        uvs(2,2,ntot) = (uv0(2,5) + uv0(2,5))/2
        uvs(1,3,ntot) = (uv0(1,3) + uv0(1,3))/2
        uvs(2,3,ntot) = (uv0(2,3) + uv0(2,3))/2
        uvs(1,4,ntot) = (uv0(1,6) + uv0(1,5))/2
        uvs(2,4,ntot) = (uv0(2,6) + uv0(2,5))/2
        uvs(1,5,ntot) = (uv0(1,3) + uv0(1,5))/2
        uvs(2,5,ntot) = (uv0(2,3) + uv0(2,5))/2
        uvs(1,6,ntot) = (uv0(1,6) + uv0(1,3))/2
        uvs(2,6,ntot) = (uv0(2,6) + uv0(2,3))/2




        ! create 4th sub-triangle
        ntot = ntot + 1
        uvs(1,1,ntot) = (uv0(1,5) + uv0(1,5))/2
        uvs(2,1,ntot) = (uv0(2,5) + uv0(2,5))/2
        uvs(1,2,ntot) = (uv0(1,6) + uv0(1,6))/2
        uvs(2,2,ntot) = (uv0(2,6) + uv0(2,6))/2
        uvs(1,3,ntot) = (uv0(1,4) + uv0(1,4))/2
        uvs(2,3,ntot) = (uv0(2,4) + uv0(2,4))/2
        uvs(1,4,ntot) = (uv0(1,5) + uv0(1,6))/2
        uvs(2,4,ntot) = (uv0(2,5) + uv0(2,6))/2
        uvs(1,5,ntot) = (uv0(1,4) + uv0(1,6))/2
        uvs(2,5,ntot) = (uv0(2,4) + uv0(2,6))/2
        uvs(1,6,ntot) = (uv0(1,4) + uv0(1,5))/2
        uvs(2,6,ntot) = (uv0(2,4) + uv0(2,5))/2

      end do
    end do

  end if

  ifplot = 0
  if (ifplot .eq. 1) then
    call prinf('kover = *', kover, 1)
    itype1 = 1
    iw = 66
    !call zpyplot(iw, uvs, 6*ntot, itype1, &
    !    'oversampled quadratic nodes*')
    stop
  end if




  ! !
  ! ! calculate the centers
  ! !
  ! do i = 1,neach
  !   x = (uvs(1,1,i) + uvs(1,2,i) + uvs(1,3,i))/3
  !   y = (uvs(2,1,i) + uvs(2,2,i) + uvs(2,3,i))/3
  !   centers(1,i) = x
  !   centers(2,i) = y
  ! end do


  !
  ! Evaluate the triangles and the interpolated sigma
  !
  ntriout = 0
  do itri = 1,ntri
    do j = 1,neach
      ntriout = ntriout + 1
      do l = 1,6
        u = uvs(1,l,j)
        v = uvs(2,l,j)
        call xeval(itri, u, v, xyz, dxyzduv, par1, par2, &
            par3, par4)
        xtriout(1,l,ntriout) = xyz(1)
        xtriout(2,l,ntriout) = xyz(2)
        xtriout(3,l,ntriout) = xyz(3)
      end do
    end do
  end do


  !call prinf('npols = *', npols, 1)

  ijk = 0
  !call prinf('ntri = *', ntri, 1)

  do itri = 1,ntri
    !call prinf('itri = *', itri, 1)
    ind = (itri-1)*npols + 1
    call dmatvec(npols, npols, umatr, sigma(ind), coefs)
    !call prin2('coefs = *', coefs, npols)

    do j = 1,neach
      ijk = ijk + 1

      do k = 1,6
        !u = centers(1,j)
        !v = centers(2,j)
        u = uvs(1,k,j)
        v = uvs(2,k,j)
        call ortho2sipols(u, v, norder, pols)
        !call prin2('pols = *', pols, npols)

        sigmaout(k,ijk) = 0
        do l = 1,npols
          sigmaout(k,ijk) = sigmaout(k,ijk) + coefs(l)*pols(l)
        end do
      end do

    end do
  end do

  return
end subroutine xtri_quadratic





subroutine dmatvec(m, n, a, x, y)
  implicit double precision (a-h,o-z)
  double precision :: a(m,n), x(n), y(m)

  do i = 1,m
    dd = 0
    do j = 1,n
      dd = dd + a(i,j)*x(j)
    end do
    y(i) = dd
  end do
  

  return
end subroutine dmatvec











subroutine xtri_vtk_flat_scalars(iw, ntri, xtri1s, m, sigma, title)
  implicit real *8 (a-h,o-z)
  real *8 :: xtri1s(3,3,ntri), sigma(m,3,ntri)
  character(len=*) :: title

  character(len=1024) :: filename, dataname, valsname, imgname
  character(len=1024) :: trisname, vecsname, centname
  character(len=12) :: fmt, fmt3, fmt4
  character(len=25) :: fmt2

  !
  ! This routine plots a sequence of FLAT TRIANGLES with surface
  ! color vals.
  !
  ! Input:
  !   iw - plot number, controls the filenames
  !   ntri - number of flat triangles
  !   xtri1s - full triangle information
  !   sigma - function to plot, tabulated at corners of triangles
  !
  ! Output:
  !   a vtk file plotIW.vtk that can be plotted in paraview
  !
  !

  if (iw .lt. 10) then
    fmt = "(A4,I1,A4)"
  elseif (iw .lt. 100) then
    fmt = "(A4,I2,A4)"
  elseif (iw .lt. 1000) then
    fmt = "(A4,I3,A4)"
  elseif (iw .lt. 10000) then
    fmt = "(A4,I4,A4)"
  end if

  write(filename, fmt) 'plot', iw, '.vtk'

  !
  ! write the vtk plotting script
  !
  iunit1 = 877
  open(unit = iunit1, file=trim(filename), status='replace')

  write(iunit1,'(a)') "# vtk DataFile Version 3.0"
  write(iunit1,'(a)') trim(title)
  write(iunit1,'(a)') "ASCII"
  write(iunit1,'(a)') "DATASET UNSTRUCTURED_GRID"
  write(iunit1,'(a,i8,a)') "POINTS ", ntri*3, " float"

  fmt2 = "(E11.5,2X,E11.5,2X,E11.5)"
  do i = 1,ntri
    write(iunit1,fmt2) xtri1s(1,1,i), xtri1s(2,1,i), xtri1s(3,1,i)
    write(iunit1,fmt2) xtri1s(1,2,i), xtri1s(2,2,i), xtri1s(3,2,i)
    write(iunit1,fmt2) xtri1s(1,3,i), xtri1s(2,3,i), xtri1s(3,3,i)
  end do


  write(iunit1,'(a,i8,i8)') "CELLS ", ntri, ntri*4

  do i = 1,ntri
    i1 = 3*(i-1) + 1
    write(iunit1,'(a,i8,i8,i8)') "3 ", i1-1, i1, i1+1
  end do

  write(iunit1,'(a,i8)') "CELL_TYPES ", ntri
  do i = 1,ntri
    write(iunit1,'(a)') "5"
  end do

  write(iunit1,'(a)') ""
  write(iunit1,'(a,i8)') "POINT_DATA ", ntri*3
  write(iunit1,'(a,i4)') "SCALARS scalar_function float ", m
  write(iunit1,'(a)') "LOOKUP_TABLE default"
  do i = 1,ntri
    do j = 1,3
      do k = 1,m
        write(iunit1,'(E11.5)') sigma(k,j,i)
      end do
    end do
  end do

  close(iunit1)

  return
end subroutine xtri_vtk_flat_scalars





subroutine xtri_vtk_quadratic_scalars(iw, ntri, xtri1s, m, sigma, title)
  implicit real *8 (a-h,o-z)
  real *8 :: xtri1s(3,6,ntri), sigma(m,6,ntri)
  character(len=*) :: title

  character(len=1024) :: filename, dataname, valsname, imgname
  character(len=1024) :: trisname, vecsname, centname
  character(len=12) :: fmt, fmt3, fmt4
  character(len=25) :: fmt2

  !
  ! This routine plots a sequence of FLAT TRIANGLES with surface
  ! color vals.
  !
  ! Input:
  !   iw - plot number, controls the filenames
  !   ntri - number of flat triangles
  !   xtri1s - full triangle information
  !   sigma - function to plot, tabulated at corners of triangles
  !
  ! Output:
  !   a vtk file plotIW.vtk that can be plotted in paraview
  !
  !

  if (iw .lt. 10) then
    fmt = "(A4,I1,A4)"
  elseif (iw .lt. 100) then
    fmt = "(A4,I2,A4)"
  elseif (iw .lt. 1000) then
    fmt = "(A4,I3,A4)"
  elseif (iw .lt. 10000) then
    fmt = "(A4,I4,A4)"
  end if

  write(filename, fmt) 'plot', iw, '.vtk'

  !
  ! write the vtk plotting script
  !
  iunit1 = 877
  open(unit = iunit1, file=trim(filename), status='replace')

  write(iunit1,'(a)') "# vtk DataFile Version 3.0"
  write(iunit1,'(a)') trim(title)
  write(iunit1,'(a)') "ASCII"
  write(iunit1,'(a)') "DATASET UNSTRUCTURED_GRID"
  write(iunit1,'(a,i8,a)') "POINTS ", ntri*6, " float"

  fmt2 = "(E11.5,2X,E11.5,2X,E11.5)"
  do i = 1,ntri
    write(iunit1,fmt2) xtri1s(1,1,i), xtri1s(2,1,i), xtri1s(3,1,i)
    write(iunit1,fmt2) xtri1s(1,2,i), xtri1s(2,2,i), xtri1s(3,2,i)
    write(iunit1,fmt2) xtri1s(1,3,i), xtri1s(2,3,i), xtri1s(3,3,i)
    write(iunit1,fmt2) xtri1s(1,4,i), xtri1s(2,4,i), xtri1s(3,4,i)
    write(iunit1,fmt2) xtri1s(1,5,i), xtri1s(2,5,i), xtri1s(3,5,i)
    write(iunit1,fmt2) xtri1s(1,6,i), xtri1s(2,6,i), xtri1s(3,6,i)
  end do


  write(iunit1,'(a,i8,i8)') "CELLS ", ntri, ntri*7

  do i = 1,ntri
    i1 = 6*(i-1) + 1
    write(iunit1,'(a,i8,i8,i8,i8,i8,i8)') "6 ", i1-1, i1, i1+1, &
        i1+2, i1+3, i1+4
  end do

  write(iunit1,'(a,i8)') "CELL_TYPES ", ntri
  do i = 1,ntri
    write(iunit1,'(a)') "22"
  end do

  write(iunit1,'(a)') ""
  write(iunit1,'(a,i8)') "POINT_DATA ", ntri*6
  write(iunit1,'(a,i4)') "SCALARS scalar_function float ", m
  write(iunit1,'(a)') "LOOKUP_TABLE default"
  do i = 1,ntri
    do j = 1,6
      do k = 1,m
        write(iunit1,'(E11.5)') sigma(k,j,i)
      end do
    end do
  end do

  close(iunit1)

  return
end subroutine xtri_vtk_quadratic_scalars





subroutine xtri_vtk_flat(iw, ntri, xtri1s, title)
  implicit real *8 (a-h,o-z)
  real *8 :: xtri1s(3,3,ntri)
  character(len=*) :: title

  character(len=1024) :: filename, dataname, valsname, imgname
  character(len=1024) :: trisname, vecsname, centname
  character(len=12) :: fmt, fmt3, fmt4
  character(len=25) :: fmt2

  !
  ! This routine plots a sequence of FLAT TRIANGLES with surface
  ! color vals.
  !
  ! Input:
  !   iw - plot number, controls the filenames
  !   ntri - number of flat triangles
  !   xtri1s - full triangle information
  !
  ! Output:
  !   files which can be executed in matlab to plot the surface
  !
  !

  if (iw .lt. 10) then
    fmt = "(A4,I1,A4)"
    fmt3 = "(A8,I1,A4)"
    fmt4 = "(A5,I1,A4)"
  elseif (iw .lt. 100) then
    fmt = "(A4,I2,A4)"
    fmt3 = "(A8,I2,A4)"
    fmt4 = "(A5,I2,A4)"
  elseif (iw .lt. 1000) then
    fmt = "(A4,I3,A4)"
    fmt3 = "(A8,I3,A4)"
    fmt4 = "(A5,I3,A4)"
  elseif (iw .lt. 10000) then
    fmt = "(A4,I4,A4)"
    fmt3 = "(A8,I4,A4)"
    fmt4 = "(A5,I4,A4)"
  end if

  write(filename, fmt) 'plot', iw, '.vtk'

  !
  ! write the vtk plotting script
  !
  iunit1 = 877
  open(unit = iunit1, file=trim(filename), status='replace')

  write(iunit1,'(a)') "# vtk DataFile Version 3.0"
  write(iunit1,'(a)') "vtk output"
  write(iunit1,'(a)') "ASCII"
  !write(iunit1,'(a)') "DATASET POLYDATA"
  write(iunit1,'(a)') "DATASET UNSTRUCTURED_GRID"
  write(iunit1,'(a,i8,a)') "POINTS ", ntri*3, " float"

  fmt2 = "(E11.5,2X,E11.5,2X,E11.5)"
  do i = 1,ntri
    write(iunit1,fmt2) xtri1s(1,1,i), xtri1s(2,1,i), xtri1s(3,1,i)
    write(iunit1,fmt2) xtri1s(1,2,i), xtri1s(2,2,i), xtri1s(3,2,i)
    write(iunit1,fmt2) xtri1s(1,3,i), xtri1s(2,3,i), xtri1s(3,3,i)
  end do


  write(iunit1,'(a,i8,i8)') "CELLS ", ntri, ntri*4

  do i = 1,ntri
    i1 = 3*(i-1) + 1
    write(iunit1,'(a,i8,i8,i8)') "3 ", i1-1, i1, i1+1
  end do

  write(iunit1,'(a,i8)') "CELL_TYPES ", ntri
  do i = 1,ntri
    write(iunit1,'(a)') "5"
  end do

  write(iunit1,'(a,i8)') "POINT_DATA ", ntri*3
  write(iunit1,'(a)') "SCALARS scalars float 1"
  write(iunit1,'(a)') "LOOKUP_TABLE default"
  do i = 1,ntri
    do j = 1,3
      write(iunit1,'(E11.5)') xtri1s(3,j,i)
    end do
  end do



  write(iunit1,'(a)') ""
  write(iunit1,'(a)') ""
  write(iunit1,'(a,i8)') "CELL_DATA ", ntri
  write(iunit1,'(a)') "SCALARS scalars float 1"
  write(iunit1,'(a)') "LOOKUP_TABLE default"
  do i = 1,ntri
    write(iunit1,'(E13.5)') (xtri1s(3,1,i) + &
        xtri1s(3,2,i) + xtri1s(3,3,i))/3
  end do







  close(iunit1)




  return
end subroutine xtri_vtk_flat








subroutine xtri_flatten(ntri, xeval, par1, par2, &
    par3, par4, norder, sigma, npols, umatr, &
    kover, ntriout, xtriout, sigmaout)
  implicit real *8 (a-h,o-z)
  real *8 :: sigma(*), xtriout(3,3,1)
  real *8 :: sigmaout(3,*), umatr(npols,npols)
  external xeval

  real *8 :: uv1(2), uv2(2), uv3(2), xyz(3), dxyzduv(3,20)
  real *8 :: coefs(10000), pols(10000)
  real *8, allocatable :: uvs(:,:,:), centers(:,:)

  !
  ! Oversample each curvilinear triangle to FLAT ones, korder times,
  ! and evaluate sigma at the corners of each triangle - this routine
  ! is mainly for plotting purposes.
  !
  ! Input:
  !   ntri -
  !   xeval - subroutine evaluating the triangles, must have the
  !       calling sequence:
  !
  !          xeval(itri, u, v, xyz, dxyzduv, par1, par2, par3, par4)
  !   norder - order of sigma on each triangle
  !   npols - number of point at which sigma is sampled on each
  !       triangle
  !   umatr - matrix converting values of sigma to expansion coefs
  !   kover - number of times to oversample each triangle
  !
  ! Output:
  !   xtriout - output flat triangles
  !   sigmaout - values of sigma at the corners of each triangle in
  !       xtriout
  !
  !

  neach = 4**kover
  allocate(uvs(2,3,neach))
  allocate(centers(2,neach))

  k1 = 1
  call xtri_uvnodes(k1, kout, uvs)

  !
  ! if necessary, recursively subdivide the triangle - first construct
  ! all the uv points
  !
  if (kover .gt. 0) then

    ntot = 1
    do i = 1,kover

      ltot = ntot

      do j = 1,ltot
        uv1(1) = uvs(1,1,j)
        uv1(2) = uvs(2,1,j)
        uv2(1) = uvs(1,2,j)
        uv2(2) = uvs(2,2,j)
        uv3(1) = uvs(1,3,j)
        uv3(2) = uvs(2,3,j)

        uvs(1,1,j) = uv1(1)
        uvs(2,1,j) = uv1(2)
        uvs(1,2,j) = (uv1(1) + uv2(1))/2
        uvs(2,2,j) = (uv1(2) + uv2(2))/2
        uvs(1,3,j) = (uv1(1) + uv3(1))/2
        uvs(2,3,j) = (uv1(2) + uv3(2))/2

        ntot = ntot + 1
        uvs(1,1,ntot) = (uv1(1) + uv2(1))/2
        uvs(2,1,ntot) = (uv1(2) + uv2(2))/2
        uvs(1,2,ntot) = uv2(1)
        uvs(2,2,ntot) = uv2(2)
        uvs(1,3,ntot) = (uv2(1) + uv3(1))/2
        uvs(2,3,ntot) = (uv2(2) + uv3(2))/2

        ntot = ntot + 1
        uvs(1,1,ntot) = (uv1(1) + uv3(1))/2
        uvs(2,1,ntot) = (uv1(2) + uv3(2))/2
        uvs(1,2,ntot) = (uv2(1) + uv3(1))/2
        uvs(2,2,ntot) = (uv2(2) + uv3(2))/2
        uvs(1,3,ntot) = uv3(1)
        uvs(2,3,ntot) = uv3(2)

        ntot = ntot + 1
        uvs(1,1,ntot) = (uv2(1) + uv3(1))/2
        uvs(2,1,ntot) = (uv2(2) + uv3(2))/2
        uvs(1,2,ntot) = (uv1(1) + uv3(1))/2
        uvs(2,2,ntot) = (uv1(2) + uv3(2))/2
        uvs(1,3,ntot) = (uv1(1) + uv2(1))/2
        uvs(2,3,ntot) = (uv1(2) + uv2(2))/2

      end do
    end do

  end if


  !
  ! calculate the centers
  !
  do i = 1,neach
    x = (uvs(1,1,i) + uvs(1,2,i) + uvs(1,3,i))/3
    y = (uvs(2,1,i) + uvs(2,2,i) + uvs(2,3,i))/3
    centers(1,i) = x
    centers(2,i) = y
  end do


  !
  ! evaluate the triangles and the interpolated sigma
  !
  ntriout = 0
  do itri = 1,ntri
    do j = 1,neach
      ntriout = ntriout + 1
      do l = 1,3
        u = uvs(1,l,j)
        v = uvs(2,l,j)
        call xeval(itri, u, v, xyz, dxyzduv, par1, par2, &
            par3, par4)
        xtriout(1,l,ntriout) = xyz(1)
        xtriout(2,l,ntriout) = xyz(2)
        xtriout(3,l,ntriout) = xyz(3)
      end do
    end do
  end do


  !call prinf('npols = *', npols, 1)

  ijk = 0
  !call prinf('ntri = *', ntri, 1)

  do itri = 1,ntri
    !call prinf('itri = *', itri, 1)
    ind = (itri-1)*npols + 1
    call dmatvec(npols, npols, umatr, sigma(ind), coefs)
    !call prin2('coefs = *', coefs, npols)

    do j = 1,neach
      ijk = ijk + 1

      do k = 1,3
        !u = centers(1,j)
        !v = centers(2,j)
        u = uvs(1,k,j)
        v = uvs(2,k,j)
        call ortho2sipols(u, v, norder, pols)
        !call prin2('pols = *', pols, npols)

        sigmaout(k,ijk) = 0
        do l = 1,npols
          sigmaout(k,ijk) = sigmaout(k,ijk) + coefs(l)*pols(l)
        end do
      end do

    end do
  end do

  return
end subroutine xtri_flatten





subroutine xtri_uvnodes(korder, kout, uvs)
  implicit real *8 (a-h,o-z)
  real *8 :: uvs(2,1)

  real *8 :: uvtemp(2,100)

  !
  ! if low order, grab the nodes exactly
  !
  done = 1
  if (korder .eq. 1) call xtri_uvnodes1(kout, uvs)
  if (korder .eq. 2) call xtri_uvnodes2(kout, uvs)
  if (korder .eq. 3) call xtri_uvnodes3(kout, uvs)

  !
  ! if higher order, they are computed recursively as per gmsh
  ! definitions
  !
  if ((korder .ge. 4) .and. (korder .le. 6)) then

    call xtri_uvnodes1(kout1,uvs)
    nedge = korder - 1
    h = done/(nedge+1)
    nk = kout1
    
    do i = 1,nedge
      nk = nk + 1
      uvs(1,nk) = i*h
      uvs(2,nk) = 0
    end do



    do i = 1,nedge
      nk = nk + 1
      uvs(1,nk) = 1 - i*h
      uvs(2,nk) = i*h
    end do


    do i = 1,nedge
      nk = nk + 1
      uvs(1,nk) = 0
      uvs(2,nk) = 1 - i*h
    end do
    

    if (korder .eq. 4) call xtri_uvnodes1(kout2, uvtemp)
    if (korder .eq. 5) call xtri_uvnodes2(kout2, uvtemp)
    if (korder .eq. 6) call xtri_uvnodes3(kout2, uvtemp)

    sc = 1-3*h
    do i = 1,kout2
      uvtemp(1,i) = uvtemp(1,i)*sc
      uvtemp(2,i) = uvtemp(2,i)*sc
    end do

    do i = 1,kout2
      nk = nk + 1
      uvs(1,nk) = uvtemp(1,i) + h      
      uvs(2,nk) = uvtemp(2,i) + h
    end do

    kout = nk
    
  end if

  return
end subroutine xtri_uvnodes





subroutine xtri_uvnodes1(kout, uvs)
  implicit real *8 (a-h,o-z)
  real *8 :: uvs(2,*)
  
  !
  ! Output:
  !   kout - the number of nodes on the triangle
  !   uvs - nodes describing the triangle, compatible with gmsh, see
  !         the document below for more information
  !              http://gmsh.info/doc/texinfo/gmsh.html
  !

  uvs(1,1) = 0
  uvs(2,1) = 0
  uvs(1,2) = 1
  uvs(2,2) = 0
  uvs(1,3) = 0
  uvs(2,3) = 1
  kout = 3

    return
end subroutine xtri_uvnodes1





subroutine xtri_uvnodes2(kout, uvs)
  implicit real *8 (a-h,o-z)
  real *8 :: uvs(2,*)
  
  !
  ! This routine returns the gmsh compatible nodes for the simplex
  !
  ! Output:
  !   kout - the number of nodes on the triangle
  !   uvs - nodes describing the triangle, compatible with gmsh, see
  !         the document below for more information
  !              http://gmsh.info/doc/texinfo/gmsh.html
  !

  uvs(1,1) = 0
  uvs(2,1) = 0
  uvs(1,2) = 1
  uvs(2,2) = 0
  uvs(1,3) = 0
  uvs(2,3) = 1
  uvs(1,4) = .5d0
  uvs(2,4) = 0
  uvs(1,5) = .5d0
  uvs(2,5) = .5d0
  uvs(1,6) = 0
  uvs(2,6) = .5d0
  kout = 6

  return
end subroutine xtri_uvnodes2





subroutine xtri_uvnodes3(kout, uvs)
  implicit real *8 (a-h,o-z)
  real *8 :: uvs(2,*)
  
  !
  ! This routine returns the gmsh compatible nodes for the simplex
  ! Output:
  !   kout - the number of nodes on the triangle
  !   uvs - nodes describing the triangle, compatible with gmsh, see
  !         the document below for more information
  !              http://gmsh.info/doc/texinfo/gmsh.html
  !

  done = 1.0d0
  d13 = done/3
  d23 = 2*done/3
  uvs(1,1) = 0
  uvs(2,1) = 0
  uvs(1,2) = 1
  uvs(2,2) = 0
  uvs(1,3) = 0
  uvs(2,3) = 1
  uvs(1,4) = d13
  uvs(2,4) = 0
  uvs(1,5) = d23
  uvs(2,5) = 0
  uvs(1,6) = d23
  uvs(2,6) = d13
  uvs(1,7) = d13
  uvs(2,7) = d23
  uvs(1,8) = 0
  uvs(2,8) = d23
  uvs(1,9) = 0
  uvs(2,9) = d13
  uvs(1,10) = d13
  uvs(2,10) = d13
  kout = 10

  return
end subroutine xtri_uvnodes3
