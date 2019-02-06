c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging routines and the beginning of
c       the auxiliary routines for triangulation
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine gentriainfo(verts,nverts,ifaces,nfaces,triainfo)
        implicit real *8 (a-h,o-z)
        dimension verts(3,1),ifaces(3,1),triainfo(3,3,1)
c
        do 1200 i=1,nfaces
c
        triainfo(1,1,i)=verts(1,ifaces(1,i))
        triainfo(2,1,i)=verts(2,ifaces(1,i))
        triainfo(3,1,i)=verts(3,ifaces(1,i))
c
        triainfo(1,2,i)=verts(1,ifaces(2,i))
        triainfo(2,2,i)=verts(2,ifaces(2,i))
        triainfo(3,2,i)=verts(3,ifaces(2,i))
c
        triainfo(1,3,i)=verts(1,ifaces(3,i))
        triainfo(2,3,i)=verts(2,ifaces(3,i))
        triainfo(3,3,i)=verts(3,ifaces(3,i))
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging routines and the beginning of
c       the auxiliary routines for triangulation (quadratic triangles)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine genqtriainfo(verts,nverts,iqfaces,nfaces,qtriainfo)
        implicit real *8 (a-h,o-z)
        dimension verts(3,1),iqfaces(6,1),qtriainfo(3,6,1)
c
        do 1200 i=1,nfaces
c
        qtriainfo(1,1,i)=verts(1,iqfaces(1,i))
        qtriainfo(2,1,i)=verts(2,iqfaces(1,i))
        qtriainfo(3,1,i)=verts(3,iqfaces(1,i))
c
        qtriainfo(1,2,i)=verts(1,iqfaces(2,i))
        qtriainfo(2,2,i)=verts(2,iqfaces(2,i))
        qtriainfo(3,2,i)=verts(3,iqfaces(2,i))
c
        qtriainfo(1,3,i)=verts(1,iqfaces(3,i))
        qtriainfo(2,3,i)=verts(2,iqfaces(3,i))
        qtriainfo(3,3,i)=verts(3,iqfaces(3,i))
c
        qtriainfo(1,4,i)=verts(1,iqfaces(4,i))
        qtriainfo(2,4,i)=verts(2,iqfaces(4,i))
        qtriainfo(3,4,i)=verts(3,iqfaces(4,i))
c
        qtriainfo(1,5,i)=verts(1,iqfaces(5,i))
        qtriainfo(2,5,i)=verts(2,iqfaces(5,i))
        qtriainfo(3,5,i)=verts(3,iqfaces(5,i))
c
        qtriainfo(1,6,i)=verts(1,iqfaces(6,i))
        qtriainfo(2,6,i)=verts(2,iqfaces(6,i))
        qtriainfo(3,6,i)=verts(3,iqfaces(6,i))
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
        subroutine genqtriainfo_flat
     $     (verts,nverts,ifaces,nfaces,qtriainfo)
        implicit real *8 (a-h,o-z)
        dimension verts(3,1),ifaces(3,1),qtriainfo(3,6,1)
c
        do 1200 i=1,nfaces
c
        qtriainfo(1,1,i)=verts(1,ifaces(1,i))
        qtriainfo(2,1,i)=verts(2,ifaces(1,i))
        qtriainfo(3,1,i)=verts(3,ifaces(1,i))
c
        qtriainfo(1,2,i)=verts(1,ifaces(2,i))
        qtriainfo(2,2,i)=verts(2,ifaces(2,i))
        qtriainfo(3,2,i)=verts(3,ifaces(2,i))
c
        qtriainfo(1,3,i)=verts(1,ifaces(3,i))
        qtriainfo(2,3,i)=verts(2,ifaces(3,i))
        qtriainfo(3,3,i)=verts(3,ifaces(3,i))
c
        qtriainfo(1,4,i)=(verts(1,ifaces(1,i))+verts(1,ifaces(2,i)))/2
        qtriainfo(2,4,i)=(verts(2,ifaces(1,i))+verts(2,ifaces(2,i)))/2
        qtriainfo(3,4,i)=(verts(3,ifaces(1,i))+verts(3,ifaces(2,i)))/2
c
        qtriainfo(1,5,i)=(verts(1,ifaces(2,i))+verts(1,ifaces(3,i)))/2
        qtriainfo(2,5,i)=(verts(2,ifaces(2,i))+verts(2,ifaces(3,i)))/2
        qtriainfo(3,5,i)=(verts(3,ifaces(2,i))+verts(3,ifaces(3,i)))/2
c
        qtriainfo(1,6,i)=(verts(1,ifaces(3,i))+verts(1,ifaces(1,i)))/2
        qtriainfo(2,6,i)=(verts(2,ifaces(3,i))+verts(2,ifaces(1,i)))/2
        qtriainfo(3,6,i)=(verts(3,ifaces(3,i))+verts(3,ifaces(1,i)))/2
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging routines and the beginning of
c       the auxiliary routines for triangulation (cubic triangles)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine genctriainfo(verts,nverts,icfaces,nfaces,ctriainfo)
        implicit real *8 (a-h,o-z)
        dimension verts(3,1),icfaces(10,1),ctriainfo(3,10,1)
c
        do 1200 i=1,nfaces
c
        do k=1,10
        ctriainfo(1,k,i)=verts(1,icfaces(k,i))
        ctriainfo(2,k,i)=verts(2,icfaces(k,i))
        ctriainfo(3,k,i)=verts(3,icfaces(k,i))
        enddo
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
