cc Copyright (C) 2009-2010: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c       
c     This file contains interaction kernels for the Laplace,
c     Helmholtz, and Maxwell equations in R^3.
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging routines and the beginning of
c       the interaction routines in R^3
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       The calling sequence is:
c
c       subroutine interact(srcinfo,targinfo,cout,par1,par2,par3,par4)
c
c       srcinfo - geometry information at the source
c       targinfo - geometry information at the target
c       cout - the complex *16 value of interaction kernel
c       par1, par2, par3, par4 - extra parameters for future extensions
c
c       geometry information is encoded as a linear array containing
c          location, normal, and two tangent vectors, for a total of 12
c          real *8 elements, e.g. srcinfo 
c
c          src(1)=srcinfo(1)
c          src(2)=srcinfo(2)
c          src(3)=srcinfo(3)
c          srcnorm(1)=srcinfo(4)
c          srcnorm(2)=srcinfo(5)
c          srcnorm(3)=srcinfo(6)
c          srctang1(1)=srcinfo(7)
c          srctang1(2)=srcinfo(8)
c          srctang1(3)=srcinfo(9)
c          srctang2(1)=srcinfo(10)
c          srctang2(2)=srcinfo(11)
c          srctang2(3)=srcinfo(12)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       ... Laplace kernels
c
c
        subroutine lfinter1(srcinfo,targinfo,cout,par1,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension srcinfo(*),targinfo(*)
        dimension src(3),targ(3)
        complex *16 cout
c
c       ... single layer potential, S_0
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        src(3)=srcinfo(3)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targ(3)=targinfo(3)
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)
c
        r=sqrt(dx**2+dy**2+dz**2)
c
        cout=1/r
c
        return
        end
c
c
c
c
c
        subroutine lfinter2(srcinfo,targinfo,cout,par1,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension srcinfo(*),targinfo(*)
        dimension src(3),targ(3)
        dimension srcnorm(3),targnorm(3)
        complex *16 cout
c
c       ... double layer potential, D_0
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        src(3)=srcinfo(3)
        srcnorm(1)=srcinfo(4)
        srcnorm(2)=srcinfo(5)
        srcnorm(3)=srcinfo(6)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targ(3)=targinfo(3)
        targnorm(1)=targinfo(4)
        targnorm(2)=targinfo(5)
        targnorm(3)=targinfo(6)
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)
c
        r=sqrt(dx**2+dy**2+dz**2)
        d=dx*srcnorm(1)+dy*srcnorm(2)+dz*srcnorm(3)
c
        cout=d/r**3
c
        return
        end
c
c
c
c
c
        subroutine lfinter3(srcinfo,targinfo,cout,par1,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension srcinfo(*),targinfo(*)
        dimension src(3),targ(3)
        dimension srcnorm(3),targnorm(3)
        complex *16 cout
c
c       ... derivative of single layer potential at the target
c       warning, this routine actually computes -S_0'
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        src(3)=srcinfo(3)
        srcnorm(1)=srcinfo(4)
        srcnorm(2)=srcinfo(5)
        srcnorm(3)=srcinfo(6)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targ(3)=targinfo(3)
        targnorm(1)=targinfo(4)
        targnorm(2)=targinfo(5)
        targnorm(3)=targinfo(6)
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)
c
        r=sqrt(dx**2+dy**2+dz**2)
        d=dx*targnorm(1)+dy*targnorm(2)+dz*targnorm(3)
c
        cout=d/r**3
c
        return
        end
c
c
c
c
c
        subroutine lfinter4(srcinfo,targinfo,cout,par1,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        dimension srcinfo(*),targinfo(*)
        dimension src(3),targ(3)
        dimension srcnorm(3),targnorm(3)
        complex *16 cout
c
c       ... derivative of double layer potential at the target
c       warning, this routine actually computes -D_0'
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        src(3)=srcinfo(3)
        srcnorm(1)=srcinfo(4)
        srcnorm(2)=srcinfo(5)
        srcnorm(3)=srcinfo(6)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targ(3)=targinfo(3)
        targnorm(1)=targinfo(4)
        targnorm(2)=targinfo(5)
        targnorm(3)=targinfo(6)
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)
c
        r=sqrt(dx**2+dy**2+dz**2)
c
        d=0
        d=d+(2*dx**2-dy**2-dz**2)*targnorm(1)*srcnorm(1)
        d=d+(2*dy**2-dz**2-dx**2)*targnorm(2)*srcnorm(2)
        d=d+(2*dz**2-dx**2-dy**2)*targnorm(3)*srcnorm(3)
c
        d=d+3*(dx*dy)*targnorm(1)*srcnorm(2)
        d=d+3*(dy*dz)*targnorm(2)*srcnorm(3)
        d=d+3*(dz*dx)*targnorm(3)*srcnorm(1)
c
        d=d+3*(dy*dx)*targnorm(2)*srcnorm(1)
        d=d+3*(dz*dy)*targnorm(3)*srcnorm(2)
        d=d+3*(dx*dz)*targnorm(1)*srcnorm(3)
c
        cout=d/r**5
c
        return
        end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       ... Helmholtz kernels
c
c
        subroutine hfinter1(srcinfo,targinfo,cout,zk,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        complex *16 zk,ima
        dimension srcinfo(*),targinfo(*)
        dimension src(3),targ(3)
        complex *16 cout
        data ima/(0.0d0,1.0d0)/
c       
c       ... single layer potential, S_k
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        src(3)=srcinfo(3)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targ(3)=targinfo(3)
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)
c
        r=sqrt(dx**2+dy**2+dz**2)
c
        cout=exp(ima*zk*r)/r
c
        return
        end
c
c
c
c
c
        subroutine hfinter2(srcinfo,targinfo,cout,zk,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        complex *16 zk,ima
        dimension srcinfo(*),targinfo(*)
        dimension src(3),targ(3)
        dimension srcnorm(3),targnorm(3)
        complex *16 cout,cd
        data ima/(0.0d0,1.0d0)/
c
c       ... double layer potential, D_k
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        src(3)=srcinfo(3)
        srcnorm(1)=srcinfo(4)
        srcnorm(2)=srcinfo(5)
        srcnorm(3)=srcinfo(6)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targ(3)=targinfo(3)
        targnorm(1)=targinfo(4)
        targnorm(2)=targinfo(5)
        targnorm(3)=targinfo(6)
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)
c
        r=sqrt(dx**2+dy**2+dz**2)
        cd=dx*srcnorm(1)+dy*srcnorm(2)+dz*srcnorm(3)
        cd=cd*(1-ima*zk*r)
c
        cout=cd*exp(ima*zk*r)/r**3
c
        return
        end
c
c
c
c
c
        subroutine hfinter3(srcinfo,targinfo,cout,zk,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        complex *16 zk,ima
        dimension srcinfo(*),targinfo(*)
        dimension src(3),targ(3)
        dimension srcnorm(3),targnorm(3)
        complex *16 cout,cd
        data ima/(0.0d0,1.0d0)/
c
c       ... derivative of single layer potential at the target
c       warning, this routine actually computes -S_k'
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        src(3)=srcinfo(3)
        srcnorm(1)=srcinfo(4)
        srcnorm(2)=srcinfo(5)
        srcnorm(3)=srcinfo(6)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targ(3)=targinfo(3)
        targnorm(1)=targinfo(4)
        targnorm(2)=targinfo(5)
        targnorm(3)=targinfo(6)
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)
c
        r=sqrt(dx**2+dy**2+dz**2)
        cd=dx*targnorm(1)+dy*targnorm(2)+dz*targnorm(3)
        cd=cd*(1-ima*zk*r)
c
        cout=cd*exp(ima*zk*r)/r**3
c
        return
        end
c
c
c
c
c
        subroutine hfinter4(srcinfo,targinfo,cout,zk,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        complex *16 zk,ima
        dimension srcinfo(*),targinfo(*)
        dimension src(3),targ(3)
        dimension srcnorm(3),targnorm(3)
        complex *16 cout,cd
        data ima/(0.0d0,1.0d0)/
c
c       ... derivative of double layer potential at the target
c       warning, this routine actually computes -D_k'
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        src(3)=srcinfo(3)
        srcnorm(1)=srcinfo(4)
        srcnorm(2)=srcinfo(5)
        srcnorm(3)=srcinfo(6)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targ(3)=targinfo(3)
        targnorm(1)=targinfo(4)
        targnorm(2)=targinfo(5)
        targnorm(3)=targinfo(6)
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)
c
        r=sqrt(dx**2+dy**2+dz**2)
c
        cd=0
        cd=cd+(2*dx**2-dy**2-dz**2)*(1-ima*zk*r)*targnorm(1)*srcnorm(1)
        cd=cd+(2*dy**2-dz**2-dx**2)*(1-ima*zk*r)*targnorm(2)*srcnorm(2)
        cd=cd+(2*dz**2-dx**2-dy**2)*(1-ima*zk*r)*targnorm(3)*srcnorm(3)
c
        cd=cd+(-zk**2*r**2*dx**2)*targnorm(1)*srcnorm(1)
        cd=cd+(-zk**2*r**2*dy**2)*targnorm(2)*srcnorm(2)
        cd=cd+(-zk**2*r**2*dz**2)*targnorm(3)*srcnorm(3)
c
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dx*dy)*targnorm(1)*srcnorm(2)
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dy*dz)*targnorm(2)*srcnorm(3)
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dz*dx)*targnorm(3)*srcnorm(1)
c
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dy*dx)*targnorm(2)*srcnorm(1)
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dz*dy)*targnorm(3)*srcnorm(2)
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dx*dz)*targnorm(1)*srcnorm(3)
c
        cout=cd*exp(ima*zk*r)/r**5
c
        return
        end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       ... Helmholtz kernels (subtract the singularity)
c
c
c
        subroutine hfinter1m(srcinfo,targinfo,cout,zk,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        complex *16 zk,ima
        dimension srcinfo(*),targinfo(*)
        dimension src(3),targ(3)
        complex *16 cout
        data ima/(0.0d0,1.0d0)/
c       
c       ... single layer potential, S_k - S_0
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        src(3)=srcinfo(3)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targ(3)=targinfo(3)
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)
c
        r=sqrt(dx**2+dy**2+dz**2)
c
        cout=exp(ima*zk*r)/r
c
c
c       ... subtract the singularity
c
        cout=cout-1/r
c
        return
        end
c
c
c
c
c
        subroutine hfinter2m(srcinfo,targinfo,cout,zk,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        complex *16 zk,ima
        dimension srcinfo(*),targinfo(*)
        dimension src(3),targ(3)
        dimension srcnorm(3),targnorm(3)
        complex *16 cout,cd
        data ima/(0.0d0,1.0d0)/
c
c       ... double layer potential, D_k - D_0
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        src(3)=srcinfo(3)
        srcnorm(1)=srcinfo(4)
        srcnorm(2)=srcinfo(5)
        srcnorm(3)=srcinfo(6)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targ(3)=targinfo(3)
        targnorm(1)=targinfo(4)
        targnorm(2)=targinfo(5)
        targnorm(3)=targinfo(6)
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)
c
        r=sqrt(dx**2+dy**2+dz**2)
        cd=dx*srcnorm(1)+dy*srcnorm(2)+dz*srcnorm(3)
        cd=cd*(1-ima*zk*r)
c
        cout=cd*exp(ima*zk*r)/r**3
c
c
c       ... subtract the singularity
c
        d=dx*srcnorm(1)+dy*srcnorm(2)+dz*srcnorm(3)
        cout=cout-d/r**3
c
        return
        end
c
c
c
c
c
        subroutine hfinter3m(srcinfo,targinfo,cout,zk,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        complex *16 zk,ima
        dimension srcinfo(*),targinfo(*)
        dimension src(3),targ(3)
        dimension srcnorm(3),targnorm(3)
        complex *16 cout,cd
        data ima/(0.0d0,1.0d0)/
c
c       ... derivative of single layer potential at the target
c       warning, this routine actually computes - (S_k' - S_0')
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        src(3)=srcinfo(3)
        srcnorm(1)=srcinfo(4)
        srcnorm(2)=srcinfo(5)
        srcnorm(3)=srcinfo(6)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targ(3)=targinfo(3)
        targnorm(1)=targinfo(4)
        targnorm(2)=targinfo(5)
        targnorm(3)=targinfo(6)
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)
c
        r=sqrt(dx**2+dy**2+dz**2)
        cd=dx*targnorm(1)+dy*targnorm(2)+dz*targnorm(3)
        cd=cd*(1-ima*zk*r)
c
        cout=cd*exp(ima*zk*r)/r**3
c
c
c       ... subtract the singularity
c
        d=dx*targnorm(1)+dy*targnorm(2)+dz*targnorm(3)
        cout=cout-d/r**3
c
        return
        end
c
c
c
c
c
        subroutine hfinter4m(srcinfo,targinfo,cout,zk,par2,par3,par4)
        implicit real *8 (a-h,o-z)
        complex *16 zk,ima
        dimension srcinfo(*),targinfo(*)
        dimension src(3),targ(3)
        dimension srcnorm(3),targnorm(3)
        complex *16 cout,cd
        data ima/(0.0d0,1.0d0)/
c
c       ... derivative of double layer potential at the target
c       warning, this routine actually computes - (D_k' - D_0')
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        src(3)=srcinfo(3)
        srcnorm(1)=srcinfo(4)
        srcnorm(2)=srcinfo(5)
        srcnorm(3)=srcinfo(6)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targ(3)=targinfo(3)
        targnorm(1)=targinfo(4)
        targnorm(2)=targinfo(5)
        targnorm(3)=targinfo(6)
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)
c
        r=sqrt(dx**2+dy**2+dz**2)
c
        cd=0
        cd=cd+(2*dx**2-dy**2-dz**2)*(1-ima*zk*r)*targnorm(1)*srcnorm(1)
        cd=cd+(2*dy**2-dz**2-dx**2)*(1-ima*zk*r)*targnorm(2)*srcnorm(2)
        cd=cd+(2*dz**2-dx**2-dy**2)*(1-ima*zk*r)*targnorm(3)*srcnorm(3)
c
        cd=cd+(-zk**2*r**2*dx**2)*targnorm(1)*srcnorm(1)
        cd=cd+(-zk**2*r**2*dy**2)*targnorm(2)*srcnorm(2)
        cd=cd+(-zk**2*r**2*dz**2)*targnorm(3)*srcnorm(3)
c
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dx*dy)*targnorm(1)*srcnorm(2)
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dy*dz)*targnorm(2)*srcnorm(3)
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dz*dx)*targnorm(3)*srcnorm(1)
c
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dy*dx)*targnorm(2)*srcnorm(1)
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dz*dy)*targnorm(3)*srcnorm(2)
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dx*dz)*targnorm(1)*srcnorm(3)
c
        cout=cd*exp(ima*zk*r)/r**5
c
c
c       ... subtract the singularity
c
        d=0
        d=d+(2*dx**2-dy**2-dz**2)*targnorm(1)*srcnorm(1)
        d=d+(2*dy**2-dz**2-dx**2)*targnorm(2)*srcnorm(2)
        d=d+(2*dz**2-dx**2-dy**2)*targnorm(3)*srcnorm(3)
c
        d=d+3*(dx*dy)*targnorm(1)*srcnorm(2)
        d=d+3*(dy*dz)*targnorm(2)*srcnorm(3)
        d=d+3*(dz*dx)*targnorm(3)*srcnorm(1)
c
        d=d+3*(dy*dx)*targnorm(2)*srcnorm(1)
        d=d+3*(dz*dy)*targnorm(3)*srcnorm(2)
        d=d+3*(dx*dz)*targnorm(1)*srcnorm(3)
c
        cout=cout-d/r**5
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
c
c       ... Maxwell kernels
c
c
        subroutine eminter1(srcinfo,targinfo,cout,zk,info,par1,par2)
        implicit real *8 (a-h,o-z)
        complex *16 zk,ima
        dimension srcinfo(*),targinfo(*),info(*)
        dimension src(3),targ(3)
        dimension srcnorm(3),targnorm(3)
        dimension srctang1(3),srctang2(3)
        dimension targtang1(3),targtang2(3)
        dimension srcvec(3),targvec(3)
        complex *16 cout,cd
        data ima/(0.0d0,1.0d0)/
c
c
c       EFIE, part 1, vector potential A, single layer
c       n x \int S_k J
c
c
        ipar=info(1)
        jpar=info(2)
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        src(3)=srcinfo(3)
        srcnorm(1)=srcinfo(4)
        srcnorm(2)=srcinfo(5)
        srcnorm(3)=srcinfo(6)
        srctang1(1)=srcinfo(7)
        srctang1(2)=srcinfo(8)
        srctang1(3)=srcinfo(9)
        srctang2(1)=srcinfo(10)
        srctang2(2)=srcinfo(11)
        srctang2(3)=srcinfo(12)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targ(3)=targinfo(3)
        targnorm(1)=targinfo(4)
        targnorm(2)=targinfo(5)
        targnorm(3)=targinfo(6)
        targtang1(1)=targinfo(7)
        targtang1(2)=targinfo(8)
        targtang1(3)=targinfo(9)
        targtang2(1)=targinfo(10)
        targtang2(2)=targinfo(11)
        targtang2(3)=targinfo(12)
c
c
        if( ipar .eq. 1 ) then
        srcvec(1)=srctang1(1)
        srcvec(2)=srctang1(2)
        srcvec(3)=srctang1(3)
        endif
c
        if( ipar .eq. 2 ) then
        srcvec(1)=srctang2(1)
        srcvec(2)=srctang2(2)
        srcvec(3)=srctang2(3)
        endif
c
        if( ipar .eq. 3 ) then
        srcvec(1)=srcnorm(1)
        srcvec(2)=srcnorm(2)
        srcvec(3)=srcnorm(3)
        endif
c
c
        if( jpar .eq. 1 ) then
        targvec(1)=targtang1(1)
        targvec(2)=targtang1(2)
        targvec(3)=targtang1(3)
        endif
c
        if( jpar .eq. 2 ) then
        targvec(1)=targtang2(1)
        targvec(2)=targtang2(2)
        targvec(3)=targtang2(3)
        endif
c
        if( jpar .eq. 3 ) then
        targvec(1)=targnorm(1)
        targvec(2)=targnorm(2)
        targvec(3)=targnorm(3)
        endif
c
c       a.(bxc) = b.(cxa) = c.(axb) 
ccc        call dot_cross_prod3d(srcvec,targnorm,targvec,d)
        call dot_cross_prod3d(targvec,targnorm,srcvec,d)
        d=d
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)
c
        r=sqrt(dx**2+dy**2+dz**2)
c
        cout=d*exp(ima*zk*r)/r
c
        return
        end
c
c
c
c
c
        subroutine eminter4(srcinfo,targinfo,cout,zk,info,par1,par2)
        implicit real *8 (a-h,o-z)
        complex *16 zk,ima
        dimension srcinfo(*),targinfo(*),info(*)
        dimension src(3),targ(3)
        dimension srcnorm(3),targnorm(3)
        dimension srctang1(3),srctang2(3)
        dimension targtang1(3),targtang2(3)
        dimension srcvec(3),targvec(3)
        complex *16 cout,cd
        data ima/(0.0d0,1.0d0)/
c
c
c       EFIE, part 2, hypersingular, with singularity extraction
c       \int (grad grad S_k) J  - \int (grad grad S_0) J
c
c       Apply tangential projection only, without n \cross operator.
c
c
        ipar=info(1)
        jpar=info(2)
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        src(3)=srcinfo(3)
        srcnorm(1)=srcinfo(4)
        srcnorm(2)=srcinfo(5)
        srcnorm(3)=srcinfo(6)
        srctang1(1)=srcinfo(7)
        srctang1(2)=srcinfo(8)
        srctang1(3)=srcinfo(9)
        srctang2(1)=srcinfo(10)
        srctang2(2)=srcinfo(11)
        srctang2(3)=srcinfo(12)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targ(3)=targinfo(3)
        targnorm(1)=targinfo(4)
        targnorm(2)=targinfo(5)
        targnorm(3)=targinfo(6)
        targtang1(1)=targinfo(7)
        targtang1(2)=targinfo(8)
        targtang1(3)=targinfo(9)
        targtang2(1)=targinfo(10)
        targtang2(2)=targinfo(11)
        targtang2(3)=targinfo(12)
c
c
        if( ipar .eq. 1 ) then
        srcvec(1)=srctang1(1)
        srcvec(2)=srctang1(2)
        srcvec(3)=srctang1(3)
        endif
c
        if( ipar .eq. 2 ) then
        srcvec(1)=srctang2(1)
        srcvec(2)=srctang2(2)
        srcvec(3)=srctang2(3)
        endif
c
        if( ipar .eq. 3 ) then
        srcvec(1)=srcnorm(1)
        srcvec(2)=srcnorm(2)
        srcvec(3)=srcnorm(3)
        endif
c
c
        if( jpar .eq. 1 ) then
        targvec(1)=targtang1(1)
        targvec(2)=targtang1(2)
        targvec(3)=targtang1(3)
        endif
c
        if( jpar .eq. 2 ) then
        targvec(1)=targtang2(1)
        targvec(2)=targtang2(2)
        targvec(3)=targtang2(3)
        endif
c
        if( jpar .eq. 3 ) then
        targvec(1)=targnorm(1)
        targvec(2)=targnorm(2)
        targvec(3)=targnorm(3)
        endif
c
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)
c
        r=sqrt(dx**2+dy**2+dz**2)
c
c
c       ... construct the hypersingular part
c
        cd=0
        cd=cd+(2*dx**2-dy**2-dz**2)*(1-ima*zk*r)*srcvec(1)*targvec(1)
        cd=cd+(2*dy**2-dz**2-dx**2)*(1-ima*zk*r)*srcvec(2)*targvec(2)
        cd=cd+(2*dz**2-dx**2-dy**2)*(1-ima*zk*r)*srcvec(3)*targvec(3)
c
        cd=cd+(-zk**2*r**2*dx**2)*srcvec(1)*targvec(1)
        cd=cd+(-zk**2*r**2*dy**2)*srcvec(2)*targvec(2)
        cd=cd+(-zk**2*r**2*dz**2)*srcvec(3)*targvec(3)
c
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dx*dy)*srcvec(1)*targvec(2)
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dy*dz)*srcvec(2)*targvec(3)
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dz*dx)*srcvec(3)*targvec(1)
c
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dy*dx)*srcvec(2)*targvec(1)
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dz*dy)*srcvec(3)*targvec(2)
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dx*dz)*srcvec(1)*targvec(3)
c
        cout=cd*exp(ima*zk*r)/r**5
c
c
c       ... and extract the singularity
c
        d=0
        d=d+(2*dx**2-dy**2-dz**2)*srcvec(1)*targvec(1)
        d=d+(2*dy**2-dz**2-dx**2)*srcvec(2)*targvec(2)
        d=d+(2*dz**2-dx**2-dy**2)*srcvec(3)*targvec(3)
c
        d=d+3*(dx*dy)*srcvec(1)*targvec(2)
        d=d+3*(dy*dz)*srcvec(2)*targvec(3)
        d=d+3*(dz*dx)*srcvec(3)*targvec(1)
c
        d=d+3*(dy*dx)*srcvec(2)*targvec(1)
        d=d+3*(dz*dy)*srcvec(3)*targvec(2)
        d=d+3*(dx*dz)*srcvec(1)*targvec(3)
c
        cout=cout-d/r**5
c
        return
        end
c
c
c
c
c
        subroutine eminter3(srcinfo,targinfo,cout,zk,info,par1,par2)
        implicit real *8 (a-h,o-z)
        complex *16 zk,ima
        dimension srcinfo(*),targinfo(*),info(*)
        dimension src(3),targ(3)
        dimension srcnorm(3),targnorm(3)
        dimension srctang1(3),srctang2(3)
        dimension targtang1(3),targtang2(3)
        dimension srcvec(3),targvec(3)
        dimension gradvec(3)
        complex *16 cout,cd
        data ima/(0.0d0,1.0d0)/
c
c
c       MFIE, curl of the vector potential
c       - n x (-grad) x \int S_k J 
c
c
        ipar=info(1)
        jpar=info(2)
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        src(3)=srcinfo(3)
        srcnorm(1)=srcinfo(4)
        srcnorm(2)=srcinfo(5)
        srcnorm(3)=srcinfo(6)
        srctang1(1)=srcinfo(7)
        srctang1(2)=srcinfo(8)
        srctang1(3)=srcinfo(9)
        srctang2(1)=srcinfo(10)
        srctang2(2)=srcinfo(11)
        srctang2(3)=srcinfo(12)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targ(3)=targinfo(3)
        targnorm(1)=targinfo(4)
        targnorm(2)=targinfo(5)
        targnorm(3)=targinfo(6)
        targtang1(1)=targinfo(7)
        targtang1(2)=targinfo(8)
        targtang1(3)=targinfo(9)
        targtang2(1)=targinfo(10)
        targtang2(2)=targinfo(11)
        targtang2(3)=targinfo(12)
c
c
        if( ipar .eq. 1 ) then
        srcvec(1)=srctang1(1)
        srcvec(2)=srctang1(2)
        srcvec(3)=srctang1(3)
        endif
c
        if( ipar .eq. 2 ) then
        srcvec(1)=srctang2(1)
        srcvec(2)=srctang2(2)
        srcvec(3)=srctang2(3)
        endif
c
        if( ipar .eq. 3 ) then
        srcvec(1)=srcnorm(1)
        srcvec(2)=srcnorm(2)
        srcvec(3)=srcnorm(3)
        endif
c
c
        if( jpar .eq. 1 ) then
        targvec(1)=targtang1(1)
        targvec(2)=targtang1(2)
        targvec(3)=targtang1(3)
        endif
c
        if( jpar .eq. 2 ) then
        targvec(1)=targtang2(1)
        targvec(2)=targtang2(2)
        targvec(3)=targtang2(3)
        endif
c
        if( jpar .eq. 3 ) then
        targvec(1)=targnorm(1)
        targvec(2)=targnorm(2)
        targvec(3)=targnorm(3)
        endif
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)
c
        gradvec(1)=dx
        gradvec(2)=dy
        gradvec(3)=dz
c
c
ccc        call dot_cross_prod3d(targvec,gradvec,srcvec,d)
ccc        dd=d
c
        call dot_cross_cross_prod3d(targvec,targnorm,gradvec,srcvec,d)
        d=-d
c
        r=sqrt(dx**2+dy**2+dz**2)
c
        cd=d*(1-ima*zk*r)
        cout=cd*exp(ima*zk*r)/r**3
c
        return
        end
c
c
c
c
c
        subroutine eminter3s(srcinfo,targinfo,cout,zk,info,par1,par2)
        implicit real *8 (a-h,o-z)
        complex *16 zk,ima
        dimension srcinfo(*),targinfo(*),info(*)
        dimension src(3),targ(3)
        dimension srcnorm(3),targnorm(3)
        dimension srctang1(3),srctang2(3)
        dimension targtang1(3),targtang2(3)
        dimension srcvec(3),targvec(3)
        dimension gradvec(3)
        complex *16 cout,cd
        data ima/(0.0d0,1.0d0)/
c
c
c       MFIE, curl of the vector potential with singularity extraction
c       - n x (-grad) x \int S_k J + n x (-grad) x \int S_0 J 
c
c
        ipar=info(1)
        jpar=info(2)
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        src(3)=srcinfo(3)
        srcnorm(1)=srcinfo(4)
        srcnorm(2)=srcinfo(5)
        srcnorm(3)=srcinfo(6)
        srctang1(1)=srcinfo(7)
        srctang1(2)=srcinfo(8)
        srctang1(3)=srcinfo(9)
        srctang2(1)=srcinfo(10)
        srctang2(2)=srcinfo(11)
        srctang2(3)=srcinfo(12)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targ(3)=targinfo(3)
        targnorm(1)=targinfo(4)
        targnorm(2)=targinfo(5)
        targnorm(3)=targinfo(6)
        targtang1(1)=targinfo(7)
        targtang1(2)=targinfo(8)
        targtang1(3)=targinfo(9)
        targtang2(1)=targinfo(10)
        targtang2(2)=targinfo(11)
        targtang2(3)=targinfo(12)
c
c
        if( ipar .eq. 1 ) then
        srcvec(1)=srctang1(1)
        srcvec(2)=srctang1(2)
        srcvec(3)=srctang1(3)
        endif
c
        if( ipar .eq. 2 ) then
        srcvec(1)=srctang2(1)
        srcvec(2)=srctang2(2)
        srcvec(3)=srctang2(3)
        endif
c
        if( ipar .eq. 3 ) then
        srcvec(1)=srcnorm(1)
        srcvec(2)=srcnorm(2)
        srcvec(3)=srcnorm(3)
        endif
c
c
        if( jpar .eq. 1 ) then
        targvec(1)=targtang1(1)
        targvec(2)=targtang1(2)
        targvec(3)=targtang1(3)
        endif
c
        if( jpar .eq. 2 ) then
        targvec(1)=targtang2(1)
        targvec(2)=targtang2(2)
        targvec(3)=targtang2(3)
        endif
c
        if( jpar .eq. 3 ) then
        targvec(1)=targnorm(1)
        targvec(2)=targnorm(2)
        targvec(3)=targnorm(3)
        endif
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)
c
        gradvec(1)=dx
        gradvec(2)=dy
        gradvec(3)=dz
c
c
ccc        call dot_cross_prod3d(targvec,gradvec,srcvec,d)
ccc        dd=d
c
        call dot_cross_cross_prod3d(targvec,targnorm,gradvec,srcvec,d)
        d=-d
c
        r=sqrt(dx**2+dy**2+dz**2)
c
        cd=d*(1-ima*zk*r)
        cout=cd*exp(ima*zk*r)/r**3
c
c
c       ... and extract the singularity
c
        cd=d
        cout=cout-cd*1/r**3
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
c
c       ... Maxwell kernels, modified for the Muller equation, v2
c
c
        subroutine eminter1h(srcinfo,targinfo,cout,zk,info,par1,par2)
        implicit real *8 (a-h,o-z)
        complex *16 zk,ima
        dimension srcinfo(*),targinfo(*),info(*)
        dimension src(3),targ(3)
        dimension srcnorm(3),targnorm(3)
        dimension srctang1(3),srctang2(3)
        dimension targtang1(3),targtang2(3)
        dimension srcvec(3),targvec(3)
        complex *16 cout,cd
        data ima/(0.0d0,1.0d0)/
c
c
c       EFIE with singularity extraction for the hypersingular part
c
c       EFIE, part 1, vector potential A, single layer
c       (\int S_k J)  * (ima*zk)
c
c       EFIE, part 2, hypersingular, with singularity extraction
c       - (\int (grad grad S_k) J  - \int (grad grad S_0) J) / (ima*zk)
c
c       Apply tangential projection only, without n \cross operator.
c
c
c
        ipar=info(1)
        jpar=info(2)
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        src(3)=srcinfo(3)
        srcnorm(1)=srcinfo(4)
        srcnorm(2)=srcinfo(5)
        srcnorm(3)=srcinfo(6)
        srctang1(1)=srcinfo(7)
        srctang1(2)=srcinfo(8)
        srctang1(3)=srcinfo(9)
        srctang2(1)=srcinfo(10)
        srctang2(2)=srcinfo(11)
        srctang2(3)=srcinfo(12)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targ(3)=targinfo(3)
        targnorm(1)=targinfo(4)
        targnorm(2)=targinfo(5)
        targnorm(3)=targinfo(6)
        targtang1(1)=targinfo(7)
        targtang1(2)=targinfo(8)
        targtang1(3)=targinfo(9)
        targtang2(1)=targinfo(10)
        targtang2(2)=targinfo(11)
        targtang2(3)=targinfo(12)
c
c
        if( ipar .eq. 1 ) then
        srcvec(1)=srctang1(1)
        srcvec(2)=srctang1(2)
        srcvec(3)=srctang1(3)
        endif
c
        if( ipar .eq. 2 ) then
        srcvec(1)=srctang2(1)
        srcvec(2)=srctang2(2)
        srcvec(3)=srctang2(3)
        endif
c
        if( ipar .eq. 3 ) then
        srcvec(1)=srcnorm(1)
        srcvec(2)=srcnorm(2)
        srcvec(3)=srcnorm(3)
        endif
c
c
        if( jpar .eq. 1 ) then
        targvec(1)=targtang1(1)
        targvec(2)=targtang1(2)
        targvec(3)=targtang1(3)
        endif
c
        if( jpar .eq. 2 ) then
        targvec(1)=targtang2(1)
        targvec(2)=targtang2(2)
        targvec(3)=targtang2(3)
        endif
c
        if( jpar .eq. 3 ) then
        targvec(1)=targnorm(1)
        targvec(2)=targnorm(2)
        targvec(3)=targnorm(3)
        endif
c
        call dot_prod3d(targvec,srcvec,d)
        dd=d
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)
c
        r=sqrt(dx**2+dy**2+dz**2)
c
        cout=dd*exp(ima*zk*r)/r *(ima*zk)
c
c
c       ... construct the hypersingular part
c
        cd=0
        cd=cd+(2*dx**2-dy**2-dz**2)*(1-ima*zk*r)*srcvec(1)*targvec(1)
        cd=cd+(2*dy**2-dz**2-dx**2)*(1-ima*zk*r)*srcvec(2)*targvec(2)
        cd=cd+(2*dz**2-dx**2-dy**2)*(1-ima*zk*r)*srcvec(3)*targvec(3)
c
        cd=cd+(-zk**2*r**2*dx**2)*srcvec(1)*targvec(1)
        cd=cd+(-zk**2*r**2*dy**2)*srcvec(2)*targvec(2)
        cd=cd+(-zk**2*r**2*dz**2)*srcvec(3)*targvec(3)
c
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dx*dy)*srcvec(1)*targvec(2)
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dy*dz)*srcvec(2)*targvec(3)
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dz*dx)*srcvec(3)*targvec(1)
c
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dy*dx)*srcvec(2)*targvec(1)
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dz*dy)*srcvec(3)*targvec(2)
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dx*dz)*srcvec(1)*targvec(3)
c
        cout=cout + cd*exp(ima*zk*r)/r**5 *(1/(-ima*zk))
c
c
c       ... and extract the singularity
c
        d=0
        d=d+(2*dx**2-dy**2-dz**2)*srcvec(1)*targvec(1)
        d=d+(2*dy**2-dz**2-dx**2)*srcvec(2)*targvec(2)
        d=d+(2*dz**2-dx**2-dy**2)*srcvec(3)*targvec(3)
c
        d=d+3*(dx*dy)*srcvec(1)*targvec(2)
        d=d+3*(dy*dz)*srcvec(2)*targvec(3)
        d=d+3*(dz*dx)*srcvec(3)*targvec(1)
c
        d=d+3*(dy*dx)*srcvec(2)*targvec(1)
        d=d+3*(dz*dy)*srcvec(3)*targvec(2)
        d=d+3*(dx*dz)*srcvec(1)*targvec(3)
c
        cout=cout - d/r**5 *(1/(-ima*zk))
c
        return
        end
c
c
c
c
c
        subroutine eminter3h(srcinfo,targinfo,cout,zk,info,par1,par2)
        implicit real *8 (a-h,o-z)
        complex *16 zk,ima
        dimension srcinfo(*),targinfo(*),info(*)
        dimension src(3),targ(3)
        dimension srcnorm(3),targnorm(3)
        dimension srctang1(3),srctang2(3)
        dimension targtang1(3),targtang2(3)
        dimension srcvec(3),targvec(3)
        dimension gradvec(3)
        complex *16 cout,cd
        data ima/(0.0d0,1.0d0)/
c
c
c       MFIE, curl of the vector potential
c       - (-grad) x \int S_k J 
c
c       Apply tangential projection only, without n \cross operator.
c
c
        ipar=info(1)
        jpar=info(2)
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        src(3)=srcinfo(3)
        srcnorm(1)=srcinfo(4)
        srcnorm(2)=srcinfo(5)
        srcnorm(3)=srcinfo(6)
        srctang1(1)=srcinfo(7)
        srctang1(2)=srcinfo(8)
        srctang1(3)=srcinfo(9)
        srctang2(1)=srcinfo(10)
        srctang2(2)=srcinfo(11)
        srctang2(3)=srcinfo(12)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targ(3)=targinfo(3)
        targnorm(1)=targinfo(4)
        targnorm(2)=targinfo(5)
        targnorm(3)=targinfo(6)
        targtang1(1)=targinfo(7)
        targtang1(2)=targinfo(8)
        targtang1(3)=targinfo(9)
        targtang2(1)=targinfo(10)
        targtang2(2)=targinfo(11)
        targtang2(3)=targinfo(12)
c
c
        if( ipar .eq. 1 ) then
        srcvec(1)=srctang1(1)
        srcvec(2)=srctang1(2)
        srcvec(3)=srctang1(3)
        endif
c
        if( ipar .eq. 2 ) then
        srcvec(1)=srctang2(1)
        srcvec(2)=srctang2(2)
        srcvec(3)=srctang2(3)
        endif
c
        if( ipar .eq. 3 ) then
        srcvec(1)=srcnorm(1)
        srcvec(2)=srcnorm(2)
        srcvec(3)=srcnorm(3)
        endif
c
c
        if( jpar .eq. 1 ) then
        targvec(1)=targtang1(1)
        targvec(2)=targtang1(2)
        targvec(3)=targtang1(3)
        endif
c
        if( jpar .eq. 2 ) then
        targvec(1)=targtang2(1)
        targvec(2)=targtang2(2)
        targvec(3)=targtang2(3)
        endif
c
        if( jpar .eq. 3 ) then
        targvec(1)=targnorm(1)
        targvec(2)=targnorm(2)
        targvec(3)=targnorm(3)
        endif
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)
c
        gradvec(1)=dx
        gradvec(2)=dy
        gradvec(3)=dz
c
c
        call dot_cross_prod3d(targvec,gradvec,srcvec,d)
        d=-d
c
        r=sqrt(dx**2+dy**2+dz**2)
c
        cd=d*(1-ima*zk*r)
        cout=cd*exp(ima*zk*r)/r**3
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
c
c       ... Maxwell kernels, modified for the Muller equation, v3
c
c
        subroutine eminter1n(srcinfo,targinfo,cout,zk,info,par1,par2)
        implicit real *8 (a-h,o-z)
        complex *16 zk,ima
        dimension srcinfo(*),targinfo(*),info(*)
        dimension src(3),targ(3)
        dimension srcnorm(3),targnorm(3)
        dimension srctang1(3),srctang2(3)
        dimension targtang1(3),targtang2(3)
        dimension srcvec(3),targvec(3),targvec0(3)
        complex *16 cout,cd
        data ima/(0.0d0,1.0d0)/
c
c
c       EFIE with singularity extraction for the hypersingular part
c
c       EFIE, part 1, vector potential A, single layer
c       n x (\int S_k J)  * (ima*zk)
c
c       EFIE, part 2, hypersingular, with singularity extraction
c       - n x (\int (grad grad S_k) J  - \int (grad grad S_0) J) / (ima*zk)
c
c
c
        ipar=info(1)
        jpar=info(2)
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        src(3)=srcinfo(3)
        srcnorm(1)=srcinfo(4)
        srcnorm(2)=srcinfo(5)
        srcnorm(3)=srcinfo(6)
        srctang1(1)=srcinfo(7)
        srctang1(2)=srcinfo(8)
        srctang1(3)=srcinfo(9)
        srctang2(1)=srcinfo(10)
        srctang2(2)=srcinfo(11)
        srctang2(3)=srcinfo(12)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targ(3)=targinfo(3)
        targnorm(1)=targinfo(4)
        targnorm(2)=targinfo(5)
        targnorm(3)=targinfo(6)
        targtang1(1)=targinfo(7)
        targtang1(2)=targinfo(8)
        targtang1(3)=targinfo(9)
        targtang2(1)=targinfo(10)
        targtang2(2)=targinfo(11)
        targtang2(3)=targinfo(12)
c
c
        if( ipar .eq. 1 ) then
        srcvec(1)=srctang1(1)
        srcvec(2)=srctang1(2)
        srcvec(3)=srctang1(3)
        endif
c
        if( ipar .eq. 2 ) then
        srcvec(1)=srctang2(1)
        srcvec(2)=srctang2(2)
        srcvec(3)=srctang2(3)
        endif
c
        if( ipar .eq. 3 ) then
        srcvec(1)=srcnorm(1)
        srcvec(2)=srcnorm(2)
        srcvec(3)=srcnorm(3)
        endif
c
c
        if( jpar .eq. 1 ) then
        targvec(1)=targtang1(1)
        targvec(2)=targtang1(2)
        targvec(3)=targtang1(3)
        endif
c
        if( jpar .eq. 2 ) then
        targvec(1)=targtang2(1)
        targvec(2)=targtang2(2)
        targvec(3)=targtang2(3)
        endif
c
        if( jpar .eq. 3 ) then
        targvec(1)=targnorm(1)
        targvec(2)=targnorm(2)
        targvec(3)=targnorm(3)
        endif
c
        call dot_cross_prod3d(targvec,targnorm,srcvec,d)
        dd=d
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)
c
        r=sqrt(dx**2+dy**2+dz**2)
c
        cout=dd*exp(ima*zk*r)/r *(ima*zk)
c
c
c       a.(bxc) = b.(cxa) = c.(axb) 
        call cross_prod3d(targvec,targnorm,targvec0)
        targvec(1)=targvec0(1)
        targvec(2)=targvec0(2)
        targvec(3)=targvec0(3)
c
c       ... construct the hypersingular part
c
        cd=0
        cd=cd+(2*dx**2-dy**2-dz**2)*(1-ima*zk*r)*srcvec(1)*targvec(1)
        cd=cd+(2*dy**2-dz**2-dx**2)*(1-ima*zk*r)*srcvec(2)*targvec(2)
        cd=cd+(2*dz**2-dx**2-dy**2)*(1-ima*zk*r)*srcvec(3)*targvec(3)
c
        cd=cd+(-zk**2*r**2*dx**2)*srcvec(1)*targvec(1)
        cd=cd+(-zk**2*r**2*dy**2)*srcvec(2)*targvec(2)
        cd=cd+(-zk**2*r**2*dz**2)*srcvec(3)*targvec(3)
c
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dx*dy)*srcvec(1)*targvec(2)
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dy*dz)*srcvec(2)*targvec(3)
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dz*dx)*srcvec(3)*targvec(1)
c
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dy*dx)*srcvec(2)*targvec(1)
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dz*dy)*srcvec(3)*targvec(2)
        cd=cd+(3-zk**2*r**2-3*ima*zk*r)*(dx*dz)*srcvec(1)*targvec(3)
c
        cout=cout + cd*exp(ima*zk*r)/r**5 *(1/(-ima*zk))
c
c
c       ... and extract the singularity
c
        d=0
        d=d+(2*dx**2-dy**2-dz**2)*srcvec(1)*targvec(1)
        d=d+(2*dy**2-dz**2-dx**2)*srcvec(2)*targvec(2)
        d=d+(2*dz**2-dx**2-dy**2)*srcvec(3)*targvec(3)
c
        d=d+3*(dx*dy)*srcvec(1)*targvec(2)
        d=d+3*(dy*dz)*srcvec(2)*targvec(3)
        d=d+3*(dz*dx)*srcvec(3)*targvec(1)
c
        d=d+3*(dy*dx)*srcvec(2)*targvec(1)
        d=d+3*(dz*dy)*srcvec(3)*targvec(2)
        d=d+3*(dx*dz)*srcvec(1)*targvec(3)
c
        cout=cout - d/r**5 *(1/(-ima*zk))
c
        return
        end
c
c
c
c
c
        subroutine eminter3n(srcinfo,targinfo,cout,zk,info,par1,par2)
        implicit real *8 (a-h,o-z)
        complex *16 zk,ima
        dimension srcinfo(*),targinfo(*),info(*)
        dimension src(3),targ(3)
        dimension srcnorm(3),targnorm(3)
        dimension srctang1(3),srctang2(3)
        dimension targtang1(3),targtang2(3)
        dimension srcvec(3),targvec(3)
        dimension gradvec(3)
        complex *16 cout,cd
        data ima/(0.0d0,1.0d0)/
c
c
c       MFIE, curl of the vector potential
c       - n x (-grad) x \int S_k J 
c
c
        ipar=info(1)
        jpar=info(2)
c
        src(1)=srcinfo(1)
        src(2)=srcinfo(2)
        src(3)=srcinfo(3)
        srcnorm(1)=srcinfo(4)
        srcnorm(2)=srcinfo(5)
        srcnorm(3)=srcinfo(6)
        srctang1(1)=srcinfo(7)
        srctang1(2)=srcinfo(8)
        srctang1(3)=srcinfo(9)
        srctang2(1)=srcinfo(10)
        srctang2(2)=srcinfo(11)
        srctang2(3)=srcinfo(12)
c
        targ(1)=targinfo(1)
        targ(2)=targinfo(2)
        targ(3)=targinfo(3)
        targnorm(1)=targinfo(4)
        targnorm(2)=targinfo(5)
        targnorm(3)=targinfo(6)
        targtang1(1)=targinfo(7)
        targtang1(2)=targinfo(8)
        targtang1(3)=targinfo(9)
        targtang2(1)=targinfo(10)
        targtang2(2)=targinfo(11)
        targtang2(3)=targinfo(12)
c
c
        if( ipar .eq. 1 ) then
        srcvec(1)=srctang1(1)
        srcvec(2)=srctang1(2)
        srcvec(3)=srctang1(3)
        endif
c
        if( ipar .eq. 2 ) then
        srcvec(1)=srctang2(1)
        srcvec(2)=srctang2(2)
        srcvec(3)=srctang2(3)
        endif
c
        if( ipar .eq. 3 ) then
        srcvec(1)=srcnorm(1)
        srcvec(2)=srcnorm(2)
        srcvec(3)=srcnorm(3)
        endif
c
c
        if( jpar .eq. 1 ) then
        targvec(1)=targtang1(1)
        targvec(2)=targtang1(2)
        targvec(3)=targtang1(3)
        endif
c
        if( jpar .eq. 2 ) then
        targvec(1)=targtang2(1)
        targvec(2)=targtang2(2)
        targvec(3)=targtang2(3)
        endif
c
        if( jpar .eq. 3 ) then
        targvec(1)=targnorm(1)
        targvec(2)=targnorm(2)
        targvec(3)=targnorm(3)
        endif
c
        dx=targ(1)-src(1)
        dy=targ(2)-src(2)
        dz=targ(3)-src(3)
c
        gradvec(1)=dx
        gradvec(2)=dy
        gradvec(3)=dz
c
c
ccc        call dot_cross_prod3d(targvec,gradvec,srcvec,d)
ccc        dd=d
c
        call dot_cross_cross_prod3d(targvec,targnorm,gradvec,srcvec,d)
        d=-d
c
        r=sqrt(dx**2+dy**2+dz**2)
c
        cd=d*(1-ima*zk*r)
        cout=cd*exp(ima*zk*r)/r**3
c
        return
        end
c
c
c
c
c
