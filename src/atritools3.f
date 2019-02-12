c
c
c
c***********************************************************************
        subroutine atri3init(ier,ir,nverts,nfaces)
c***********************************************************************
        implicit real *8 (a-h,o-z)
        dimension verts(3,1),ifaces(3,1)
c
c       read geometry description file information in cart3d (.tri) format. 
c
c       Input parameters: 
c
c       ir : integer                      input file unit
c     
c       Output parameters:
c
c       ier : integer                     error code
c       nverts : integer                  number of vertices
c       nfaces : integer                  number of triangles
c
        ier=0
c
        open(unit=ir, status="OLD", iostat=istat)
        if (istat .ne. 0) then
cccc          write(*,*)"in atriread, OPEN statement iostat = ", istat
          ier=1
          return
        endif
c
        read(ir,*) nVerts, nFaces
c
        return
        end
c
c
c
c
c
c***********************************************************************
        subroutine atriread3(ier,ir,verts,nverts,ifaces,nfaces)
c***********************************************************************
        implicit real *8 (a-h,o-z)
        dimension verts(3,1),ifaces(3,1)
c
c       read geometry description file in cart3d (.tri) format. 
c
c       Input parameters: 
c
c       ir : integer                      input file unit
c     
c       Output parameters:
c
c       ier : integer                     error code
c       nverts : integer                  number of vertices
c       nfaces : integer                  number of triangles
c       verts  : real*8(3,nverts)         array of vertices
c       ifaces : real*8(3,nfaces)         indices of triangle vertices
c
        ier=0
c
        open(unit=ir, status="OLD", iostat=istat)
        if (istat .ne. 0) then
cccc          write(*,*)"in atriread, OPEN statement iostat = ", istat
          ier=1
          return
        endif
c
        read(ir,*) nVerts, nFaces
c
        read(ir,*) (verts(1,j),verts(2,j),verts(3,j),j=1,nVerts) 
        read(ir,*) (ifaces(1,j),ifaces(2,j),ifaces(3,j),j=1,nFaces)
c
        return
        end
c
c
c
c
c
c***********************************************************************
        subroutine atriwrite3(ier,iw,verts,nverts,ifaces,nfaces)
c***********************************************************************
        implicit real *8 (a-h,o-z)
        dimension verts(3,1),ifaces(3,1)
c
c       write geometry description file in cart3d (.tri) format. 
c
c       Input parameters: 
c
c       iw : integer                      output file unit
c       nverts : integer                  number of vertices
c       nfaces : integer                  number of triangles
c       verts  : real*8(3,nverts)         array of vertices
c       ifaces : real*8(3,nfaces)         indices of triangle vertices
c     
c       Output parameters:
c
c       ier : integer                     error code
c
c
        ier=0
c
        write(iw,*) nverts, nfaces
        write(iw,1800) (verts(1,j),verts(2,j),verts(3,j),j=1,nverts) 
        write(iw,1900) (ifaces(1,j),ifaces(2,j),ifaces(3,j),j=1,nfaces) 
 1800   format(3(1x,e22.16))
 1900   format(3(1x,i7))
c
        return
        end
c
c
c
c
c
c***********************************************************************
        subroutine qtri3init(ier,ir,nverts,nfaces)
c***********************************************************************
        implicit real *8 (a-h,o-z)
        dimension verts(3,1),ifaces(6,1)
c
c       read geometry description file information in cart3d (q.tri) format. 
c
c       Input parameters: 
c
c       ir : integer                      input file unit
c     
c       Output parameters:
c
c       ier : integer                     error code
c       nverts : integer                  number of vertices
c       nfaces : integer                  number of triangles
c
        ier=0
c
        open(unit=ir, status="OLD", iostat=istat)
        if (istat .ne. 0) then
cccc          write(*,*)"in atriread, OPEN statement iostat = ", istat
          ier=1
          return
        endif
c
        read(ir,*) nVerts, nFaces
c
        return
        end
c
c
c
c***********************************************************************
        subroutine qtriread3(ier,ir,verts,nverts,ifaces,nfaces)
c***********************************************************************
        implicit real *8 (a-h,o-z)
        dimension verts(3,1),ifaces(6,1)
c
c       read geometry description file in cart3d (q.tri) format. 
c
c       Input parameters: 
c
c       ir : integer                      input file unit
c     
c       Output parameters:
c
c       ier : integer                     error code
c       nverts : integer                  number of vertices
c       nfaces : integer                  number of triangles
c       verts  : real*8(3,nverts)         array of vertices
c       ifaces : real*8(6,nfaces)         indices of triangle vertices
c
        ier=0
c
        open(unit=ir, status="OLD", iostat=istat)
        if (istat .ne. 0) then
cccc          write(*,*)"in atriread, OPEN statement iostat = ", istat
          ier=1
          return
        endif
c
        read(ir,*) nVerts, nFaces
c
        read(ir,*) (verts(1,j),verts(2,j),verts(3,j),j=1,nVerts) 
        read(ir,*) (ifaces(1,j),ifaces(2,j),ifaces(3,j),
     $     ifaces(4,j),ifaces(5,j),ifaces(6,j),j=1,nFaces)
c
        return
        end
c
c
c
c***********************************************************************
        subroutine qtriwrite3(ier,iw,verts,nverts,ifaces,nfaces)
c***********************************************************************
        implicit real *8 (a-h,o-z)
        dimension verts(3,1),ifaces(6,1)
c
c       write geometry description file in cart3d (q.tri) format. 
c
c       Input parameters: 
c
c       iw : integer                      output file unit
c       nverts : integer                  number of vertices
c       nfaces : integer                  number of triangles
c       verts  : real*8(3,nverts)         array of vertices
c       ifaces : real*8(6,nfaces)         indices of triangle vertices
c     
c       Output parameters:
c
c       ier : integer                     error code
c
c
        ier=0
c
        write(iw,*) nverts, nfaces
        write(iw,1800) (verts(1,j),verts(2,j),verts(3,j),j=1,nverts) 
        write(iw,1900) (ifaces(1,j),ifaces(2,j),ifaces(3,j),
     $     ifaces(4,j),ifaces(5,j),ifaces(6,j),j=1,nfaces) 
 1800   format(3(1x,e22.16))
 1900   format(6(1x,i7))
c
        return
        end
c
c
c
c
c
c***********************************************************************
        subroutine ctri3init(ier,ir,nverts,nfaces)
c***********************************************************************
        implicit real *8 (a-h,o-z)
        dimension verts(3,1),ifaces(10,1)
c
c       read geometry description file information in cart3d (c.tri) format. 
c
c       Input parameters: 
c
c       ir : integer                      input file unit
c     
c       Output parameters:
c
c       ier : integer                     error code
c       nverts : integer                  number of vertices
c       nfaces : integer                  number of triangles
c
        ier=0
c
        open(unit=ir, status="OLD", iostat=istat)
        if (istat .ne. 0) then
cccc          write(*,*)"in atriread, OPEN statement iostat = ", istat
          ier=1
          return
        endif
c
        read(ir,*) nVerts, nFaces
c
        return
        end
c
c
c
c***********************************************************************
        subroutine ctriread3(ier,ir,verts,nverts,ifaces,nfaces)
c***********************************************************************
        implicit real *8 (a-h,o-z)
        dimension verts(3,1),ifaces(10,1)
c
c       read geometry description file in cart3d (c.tri) format. 
c
c       Input parameters: 
c
c       ir : integer                      input file unit
c     
c       Output parameters:
c
c       ier : integer                     error code
c       nverts : integer                  number of vertices
c       nfaces : integer                  number of triangles
c       verts  : real*8(3,nverts)         array of vertices
c       ifaces : real*8(10,nfaces)        indices of triangle vertices
c
        ier=0
c
        open(unit=ir, status="OLD", iostat=istat)
        if (istat .ne. 0) then
cccc          write(*,*)"in atriread, OPEN statement iostat = ", istat
          ier=1
          return
        endif
c
        read(ir,*) nVerts, nFaces
c
        read(ir,*) (verts(1,j),verts(2,j),verts(3,j),j=1,nVerts) 
        read(ir,*) (ifaces(1,j),ifaces(2,j),ifaces(3,j),
     $     ifaces(4,j),ifaces(5,j),ifaces(6,j),ifaces(7,j),
     $     ifaces(8,j),ifaces(9,j),ifaces(10,j),j=1,nFaces)
c
        return
        end
c
c
c
c***********************************************************************
        subroutine ctriwrite3(ier,iw,verts,nverts,ifaces,nfaces)
c***********************************************************************
        implicit real *8 (a-h,o-z)
        dimension verts(3,1),ifaces(10,1)
c
c       write geometry description file in cart3d (c.tri) format. 
c
c       Input parameters: 
c
c       iw : integer                      output file unit
c       nverts : integer                  number of vertices
c       nfaces : integer                  number of triangles
c       verts  : real*8(3,nverts)         array of vertices
c       ifaces : real*8(10,nfaces)        indices of triangle vertices
c     
c       Output parameters:
c
c       ier : integer                     error code
c
c
        ier=0
c
        write(iw,*) nverts, nfaces
        write(iw,1800) (verts(1,j),verts(2,j),verts(3,j),j=1,nverts) 
        write(iw,1900) (ifaces(1,j),ifaces(2,j),ifaces(3,j),
     $     ifaces(4,j),ifaces(5,j),ifaces(6,j),ifaces(7,j),
     $     ifaces(8,j),ifaces(9,j),ifaces(10,j),j=1,nfaces) 
 1800   format(3(1x,e22.16))
 1900   format(10(1x,i7))
c
        return
        end
c
c
c
c
c
