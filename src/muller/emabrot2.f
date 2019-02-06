cc Copyright (C) 2010: Leslie Greengard and Zydrunas Gimbutas
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
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
c       This is the end of the debugging code and the beginning 
c       of the electromagnetic multipole expansion rotation routines
c       
c       This version is compatible with the notation _without_ (-1)**m
c       scaling factor
c
c       Maxwell's equations in R^3
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       This file contains a set of subroutines for the handling of
c       electromagnetic field multipole expansions. It contains 6
c       use-callable subroutines. Following is a brief description of
c       these subroutines.
c
c
c   emabrotyzf - rotate the z-axis oriented em-multipole
c       expansion to point into the direction dir
c
c   emabrotyzb - rotates the dir-oriented em-multipole
c       expansion to point into the direction z
c
c   emabrotaf - perform the forward Euler rotation of the em-multipole 
c
c   emabrotab - perform the backward Euler rotation of the em-multipole 
c
c
c
        subroutine emabrotyzf(nterms,dir,coefs,couts)
        implicit real *8 (a-h,o-z)
c
c       This subroutine rotates the z-axis oriented em-multipole
c       expansion to point into the direction dir. Such a transformation
c       has one degree of freedom and therefore is intended to be used
c       together with emabrotyzb
c
c       Input parameters:
c
c       nterms - the number of terms in em-multipole expansion
c       dir - the new orientation direction, real *8 (3)
c       coefs - the coefficients of em-multipole expansion
c                complex *16(0:nterms,-nterms:nterms)
c
c       Output parameters:
c
c       couts - the coefficients of rotated em-multipole expansion
c                complex *16(0:nterms,-nterms:nterms)
c
c
        dimension dir(3)
        complex *16 coefs(0:nterms,-nterms:nterms)
        complex *16 couts(0:nterms,-nterms:nterms)
        complex *16 mpole2(1 000 000)
        dimension w(1 000 000)
c
c       em-fields
c       rotate em-multipole expansion 
c
c
        x=dir(1)
        y=dir(2)
        z=dir(3)
c
        r=sqrt(x*x+y*y+z*z)
c
        x=x/r
        y=y/r
        z=z/r
c       
        if( r .eq. 0 ) then
          phi=0
          theta=0
        else
          phi=atan2(y,x)
          theta=acos(z)
        endif
c
        lw=1 000 000
c
        lmp = (nterms+1)*(2*nterms+1)        
        if( lmp .gt. lw) then
        call prinf('in emabrotyzf, lmp=*',lmp,1)
        endif
c
        lused = 2*(nterms+1)*(2*nterms+1) + 2*(2*nterms+1) + 7
        if( lmp .gt. lw) then
        call prinf('in emabrotyzf, lused=*',lused,1)
        endif
c
        call mpolerotz(nterms,coefs,phi,mpole2)
        call mpoleroty(nterms,mpole2,theta,couts,w,lw,lused)
c       
        return
        end
c
c
c
c
c
        subroutine emabrotyzb(nterms,dir,coefs,couts)
        implicit real *8 (a-h,o-z)
c
c       This subroutine rotates the dir-oriented em-multipole
c       expansion to point into the direction z. Such a transformation has
c       one degree of freedom and therefore is intended to be used
c       together with emabrotyzf
c
c       Input parameters:
c
c       nterms - the number of terms in em-multipole expansion
c       dir - the new orientation direction, real *8 (3)
c       coefs - the coefficients of em-multipole expansion
c                complex *16(0:nterms,-nterms:nterms)
c
c       Output parameters:
c
c       couts - the coefficients of rotated em-multipole expansion
c                complex *16(0:nterms,-nterms:nterms)
c
        dimension dir(3)
        complex *16 coefs(0:nterms,-nterms:nterms)
        complex *16 couts(0:nterms,-nterms:nterms)
        complex *16 mpole2(1 000 000)
        dimension w(1 000 000)
c
c       em-fields
c       rotate back
c       
        x=dir(1)
        y=dir(2)
        z=dir(3)
c
        r=sqrt(x*x+y*y+z*z)
c
        x=x/r
        y=y/r
        z=z/r
c        
        if( r .eq. 0 ) then
        phi=0
        theta=0
        else
        phi=atan2(y,x)
        theta=acos(z)
        endif
c
        lw=1 000 000
c
        lmp = (nterms+1)*(2*nterms+1)        
        if( lmp .gt. lw) then
        call prinf('in emabrotyzb, lmp=*',lmp,1)
        endif
c
        lused = 2*(nterms+1)*(2*nterms+1) + 2*(2*nterms+1) + 7
        if( lmp .gt. lw) then
        call prinf('in emabrotyzb, lused=*',lused,1)
        endif
c        
        call mpoleroty(nterms,coefs,-theta,mpole2,w,lw,lused)
        call mpolerotz(nterms,mpole2,-phi,couts)
c       
        return
        end
c
c
c
c
c
        subroutine emabrotab(nterms,rota,coefs,couts)
        implicit real *8 (a-h,o-z)
c
c       This subroutine performs the forward 
c       Euler rotations the em-multipole expansion.
c
c       Input parameters:
c
c       nterms - the number of terms in em-multipole expansion
c       rota - the Euler angles, real *8 (3)
c       coefs - the coefficients of em-multipole expansion
c                complex *16(0:nterms,-nterms:nterms)
c
c       Output parameters:
c
c       couts - the coefficients of rotated em-multipole expansion
c                complex *16(0:nterms,-nterms:nterms)
c
        dimension rota(3)
        complex *16 coefs(0:nterms,-nterms:nterms)
        complex *16 couts(0:nterms,-nterms:nterms)
        complex *16 mpole1(1 000 000)
        complex *16 mpole2(1 000 000)
        dimension w(1 000 000)
c
c       em-fields
c       rotate em-multipole expansion 
c
        lw=1 000 000
c
        lmp = (nterms+1)*(2*nterms+1)        
        if( lmp .gt. lw) then
        call prinf('in emabrotaf, lmp=*',lmp,1)
        endif
c
        lused = 2*(nterms+1)*(2*nterms+1) + 2*(2*nterms+1) + 7
        if( lmp .gt. lw) then
        call prinf('in emabrotaf, lused=*',lused,1)
        endif
c
        call mpolerotz(nterms,coefs,rota(3),mpole2)
        call mpoleroty(nterms,mpole2,rota(2),mpole1,w,lw,lused)
        call mpolerotz(nterms,mpole1,rota(1),couts)
c       
        return
        end
c
c
c
c
c
        subroutine emabrotaf(nterms,rota,coefs,couts)
        implicit real *8 (a-h,o-z)
c
c       This subroutine performs the backward
c       Euler rotations the em-multipole expansion.
c
c       Input parameters:
c
c       nterms - the number of terms in em-multipole expansion
c       rota - the Euler angles, real *8 (3)
c       coefs - the coefficients of em-multipole expansion
c                complex *16(0:nterms,-nterms:nterms)
c
c       Output parameters:
c
c       couts - the coefficients of rotated em-multipole expansion
c                complex *16(0:nterms,-nterms:nterms)
c
        dimension rota(3)
        complex *16 coefs(0:nterms,-nterms:nterms)
        complex *16 couts(0:nterms,-nterms:nterms)
        complex *16 mpole1(1 000 000)
        complex *16 mpole2(1 000 000)
        dimension w(1 000 000)
c
c       em-fields
c       rotate em-multipole expansion 
c
        lw=1 000 000
c
        lmp = (nterms+1)*(2*nterms+1)        
        if( lmp .gt. lw) then
        call prinf('in emabrotab, lmp=*',lmp,1)
        endif
c
        lused = 2*(nterms+1)*(2*nterms+1) + 2*(2*nterms+1) + 7
        if( lmp .gt. lw) then
        call prinf('in emabrotab, lused=*',lused,1)
        endif
c
        call mpolerotz(nterms,coefs,-rota(1),mpole2)
        call mpoleroty(nterms,mpole2,-rota(2),mpole1,w,lw,lused)
        call mpolerotz(nterms,mpole1,-rota(3),couts)
c       
        return
        end
c
c
c
c
c
        subroutine mpolerotz(nterms,mpole,phi,mpout)
        implicit real *8 (a-h,o-z)
c
c       This subroutine performs the rotation of the multipole expansion
c       about the z-axis
c
c     Input parameters:
c
c     nterms  (integer *4)  : number of terms in the multipole expansion
c     mpole  (complex *16)  : the multipole expansion
c     phi        (real *8)  : the angle of rotation
c
c     Output parameters:
c 
c     mpout  (complex *16)  : the rotated multipole expansion
c
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 mpout(0:nterms,-nterms:nterms)
        complex *16 ephi,ima
        data ima/(0.0d0,1.0d0)/
c
c
c----- a rotation of PHI radians about the Z-axis 
c
c
        do l=0,nterms
        do m=-l,l
c
        ephi=exp(ima*m*phi)
        mpout(l,m)=ephi*mpole(l,m)
c
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
        subroutine mpoleroty(nterms,mpole,theta,mpout,w,lw,lused)
        implicit real *8 (a-h,o-z)
c
c       This subroutine performs the rotation of the multipole expansion
c       about the y-axis
c
c     Input parameters:
c
c     nterms  (integer *4)  : number of terms in the multipole expansion
c     mpole  (complex *16)  : the multipole expansion
c     theta
c
c     Output parameters:
c 
c     mpout  (complex *16)  : the rotated multipole expansion
c
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 mpout(0:nterms,-nterms:nterms)
c
c
c----- a rotation of THETA radians about the Y-axis 
c
        call rotviarecur3(theta,nterms,nterms,nterms,mpole,nterms,
     $     mpout,nterms,w,lw,lused)
c
        return
        end
c
c
c
c
c       
