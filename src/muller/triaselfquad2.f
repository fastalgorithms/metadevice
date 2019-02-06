c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       this is the end of the debugging code and the beginning 
c       of the triangle self interaction quadrature routines
c
c       integrate single and double layer singularities on a triangle
c
c       Quadratures for self interaction of nodes in ortho2exps2.f
c       Orders 0, 1, 2, and 3
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine triaselfquad
     $     (ier,itype,norder,inode,xs,ys,ws,npts,x0,y0)
        implicit real *8 (a-h,o-z)
        dimension xs(1),ys(1),ws(1)
c       
        ier=1
c
        if( norder .eq. 0 ) then
c
        ier=0
        call tria_ord0_ret(xs,ys,ws,npts,x0,y0)
c
        endif
c
c
        if( norder .eq. 1 ) then
c
        ier=0
        call tria_ord1_ret(xs,ys,ws,npts,x0,y0)
c
        if(inode. eq. 2) then
        call xyrotate120(npts,xs,ys)
        call xyrotate120(1,x0,y0)
        endif
c        
        if(inode. eq. 3) then
        call xyrotate120(npts,xs,ys)
        call xyrotate120(1,x0,y0)
        call xyrotate120(npts,xs,ys)
        call xyrotate120(1,x0,y0)
        endif
        
        if(inode .eq. 1) then
c       ... do nothing
        endif
c
        endif
c
c
        if( norder .eq. 2 ) then
c
        ier=0

        if(inode. eq. 1 .or. inode .eq. 2 .or. inode .eq. 5 ) then
        call tria_ord2_ret2(xs,ys,ws,npts,x0,y0)
c
        if(inode. eq. 1) then
        call xyrotate120(npts,xs,ys)
        call xyrotate120(1,x0,y0)
        endif
c        
        if(inode. eq. 2) then
        call xyrotate120(npts,xs,ys)
        call xyrotate120(1,x0,y0)
        call xyrotate120(npts,xs,ys)
        call xyrotate120(1,x0,y0)
        endif
        
        if(inode .eq. 5) then
c       ... do nothing
        endif
c
        endif
c
        if(inode. eq. 3 .or. inode .eq. 4 .or. inode .eq. 6 ) then
        call tria_ord2_ret1(xs,ys,ws,npts,x0,y0)
c
        if(inode. eq. 6) then
        call xyrotate120(npts,xs,ys)
        call xyrotate120(1,x0,y0)
        endif
c        
        if(inode. eq. 4) then
        call xyrotate120(npts,xs,ys)
        call xyrotate120(1,x0,y0)
        call xyrotate120(npts,xs,ys)
        call xyrotate120(1,x0,y0)
        endif
       
        if(inode .eq. 3) then
c       ... do nothing
        endif
c
        endif
c
        endif
c
c
c
c
        if( norder .eq. 3 ) then
c
        ier=0

        if(inode. eq. 1 ) then
        call tria_ord0_ret(xs,ys,ws,npts,x0,y0)
        endif
c
        if(inode. eq. 2 .or. inode .eq. 3 .or. inode .eq. 5 ) then
        call tria_ord3_ret1(xs,ys,ws,npts,x0,y0)
c
        if(inode. eq. 3) then
        call xyrotate120(npts,xs,ys)
        call xyrotate120(1,x0,y0)
        endif
c        
        if(inode. eq. 5) then
        call xyrotate120(npts,xs,ys)
        call xyrotate120(1,x0,y0)
        call xyrotate120(npts,xs,ys)
        call xyrotate120(1,x0,y0)
        endif
c        
        if(inode .eq. 2) then
c       ... do nothing
        endif
c
        endif
c
c
        if(inode. eq. 4 .or. inode .eq. 6 .or. inode .eq. 10 ) then
        call tria_ord3_ret2(xs,ys,ws,npts,x0,y0)
c
        if(inode. eq. 6) then
        call xyrotate120(npts,xs,ys)
        call xyrotate120(1,x0,y0)
        endif
c        
        if(inode. eq. 10) then
        call xyrotate120(npts,xs,ys)
        call xyrotate120(1,x0,y0)
        call xyrotate120(npts,xs,ys)
        call xyrotate120(1,x0,y0)
        endif
        
        if(inode .eq. 4) then
c       ... do nothing
        endif
c
        endif
c
c
        if(inode. eq. 8 .or. inode .eq. 9 .or. inode .eq. 7 ) then
        call tria_ord3_ret2(xs,ys,ws,npts,x0,y0)
c
        call xyflipx(npts,xs,ys)
        call xyflipx(1,x0,y0)
c
        if(inode. eq. 7) then
        call xyrotate120(npts,xs,ys)
        call xyrotate120(1,x0,y0)
        endif
c        
        if(inode. eq. 8) then
        call xyrotate120(npts,xs,ys)
        call xyrotate120(1,x0,y0)
        call xyrotate120(npts,xs,ys)
        call xyrotate120(1,x0,y0)
        endif
        
        if(inode .eq. 9) then
c       ... do nothing
        endif
c
        endif
c
        endif
c
c
c
c
        if( itype .eq. 2 ) then
c
        do i=1,npts
        call xywmaptosimplex(xs(i),ys(i),ws(i),u,v,w)
        xs(i)=u
        ys(i)=v
        ws(i)=w
        enddo
c
        call xywmaptosimplex(x0,y0,w,u,v,w)
        x0=u
        y0=v
c
        endif
c
        return
        end
c
c
c
c
c
