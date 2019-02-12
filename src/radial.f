cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code and beginning of the 
c       quadrature code.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains subroutines for generating quadrature rules
c       for integrating a class of radially singular functions over
c       more or less arbitrary triangles containing the origin.
c
c       The quadrature rules are designed to integrate functions
c       which admit an expansion of the form
c
c                                N    (j-1)
c       f(r\cos(t),r\sin(t)) = \sum  r      p  (t) ,                       (1)
c                               j=0          3j+2
c
c       where p_{3j+2} denotes a trigonometric polynomial of order 3j+2,
c       over more or less arbitrary triangles containing the origin.
c       We say that the order of the quadrature rule (1) is N.
c
c       These rules are constructed on the fly using precomputed tables 
c       stored on the disk.  The possible orders of the quadrature 
c       rules are fixed at the time of precomputation.
c
c       Note also that the routines will fail for sufficiently 
c       degenerate triangles; the precise tolerances depend on the
c       limits of the parameters for the precomputed quadrature rules.
c
c       The quadrature rules are accurate to roughly 25 digits.
c
c       The following subroutines are user-callable:
c
c   radial_init - read a precomputed table of quadraturs from a 
c       text file on the disk.
c
c   raddiag - return a quadrature for integrating a function of the 
c        form (1) over a triangle containing the origin.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



        subroutine radial_init(ier,norder,rad,lrad,lkeep)
        implicit double precision (a-h,o-z)
        dimension rad(1)
        dimension rs(2,1000),ts(2,1000)
        character*2 str
        character*27 filename
c
c       Read a table of precomputed quadrature rules from a file
c       on the disk into the user-supplied array.  The table is
c       stored in the file ``radquads??.txt'' where ?? is the order
c       of the quadrature rules.
c
c                            Input Parameters:
c
c   norder - the order of the quadrature rules; at the present time,
c       norder can be 4, 8, 12, 16 or 20
c   lrad - length of the user-supplied array rad which will hold the
c       quadrature table
c
c                            Output Parameters:
c
c   ier - an error return code;
c       ier = 0    indicates successful execution
c       ier = 4    mean that the user-supplied array rad is of 
c                  insufficient length
c       ier = 128  means that the input file could not be opened
c       ier = 1024 means that a formatting error was encountered while
c                  attempting to read the input file
c
c   rad - on return, this user-supplied array will contain a quadrature
c       table and a structure header describing the quadratures
c
c
 0050 format ("../data/quad/radquads",I2.2,".txt")
 0100 format (I3.3)
 0200 format (D44.36)
c
        ier   = 0
        lkeep = 0
c
        max   = 1000
c
        write(filename,0050) norder
c
        iw = 101
        open(iw,FILE=filename,STATUS='old',ERR=1000)
c
c       Grab the header data.        
c
        read(iw,0100) nr
        read(iw,0100) nt
c
c       Read the parameter intervals.
c
        read (iw,0200) (rs(1,j),rs(2,j),j=1,nr)
        read (iw,0200) (ts(1,j),ts(2,j),j=1,nt)
c
cccc        call prina("in radial_init, filename = *",filename,16)
cccc        call prin2("in radial_init, rs = *",rs,2*nr)
cccc        call prin2("in radial_init, ts = *",ts,2*nt)
c
c       Fill out the header of the rad structure.
c
        nlege = (norder)/2+1
c
        irs = 100
        lrs = 2*nr
c
        its = irs+lrs
        lts = 2*nt
c
        ixslege = its+lts
        lxslege = nlege
c
        iwhtslege = ixslege+lxslege
        lwhtslege = nlege
c
        ins = iwhtslege+lwhtslege
        lns = nr*nt

        ixs = ins+lns
        lxs = max*nt*nr
c
        iwhts = ixs+lxs
        lwhts = max*nr*nt
c
        lkeep = iwhts+lwhts
c
        if (lkeep .gt. lrad) then
        ier = 4
        return
        endif
c
        rad(1)  = nr
        rad(2)  = nt
        rad(3)  = irs
        rad(4)  = its
        rad(5)  = ins
        rad(6)  = ixs
        rad(7)  = iwhts
        rad(8)  = max
c
        rad(10) = nlege
        rad(11) = ixslege
        rad(12) = iwhtslege
c
c       Construct the Legendre quadrature.
c
        call legequad(nlege,rad(ixslege),rad(iwhtslege))
c
c       Copy the as and bs into the array.
c
        call radmove(nr*2,rs,rad(irs))
        call radmove(nt*2,ts,rad(its))
c
c       Call the auxillary routine to read the quadrature rules.
c
        call radinit0(ier,iw,max,nr,nt,rad(ins),rad(ixs),rad(iwhts),
     -    rad(irs),rad(its),norder)
        if (ier .ne. 0) return
c
        close(iw)
        return
c
c       We arrive here on IO errors.
c
 1000 continue
        ier = 128
        return
        end
c
c
c
        subroutine radinit0(ier,iw,max,nr,nt,ns,xs,whts,rs,ts,norder)
        implicit double precision (a-h,o-z)
        dimension ns(nt,nr),rs(2,1),ts(2,1)
        dimension xs(max,nt,nr),whts(max,nt,nr)
        character*2 str
c
        ier = 0
 0100 format (I3.3)
 0200 format (D44.36)
c
c       Read each of the quadrature rules from the text file.
c        
        do 1100 ir=1,nr
        do 1200 it=1,nt
c
c       Read the quadrature.
c
        read (iw,0100) nn
        read (iw,0200) (xs(i,it,ir),i=1,nn)
        read (iw,0200) (whts(i,it,ir),i=1,nn)
c

c$$$ 3000 format("$",I2,"$ &","$",E13.6,"$ & $",E13.6,"$ & $",E13.6,"$ & $",
c$$$     -    E13.6," $ & $",I3,"$ \\")
c$$$
c$$$        call corrand3(1,dd)
c$$$        if (dd .gt. .93d0) then
c$$$        write (*,3000) norder,rs(1,ir),rs(2,ir),ts(1,it),ts(2,it),nn
c$$$        endif
c
        ns(it,ir)=nn
c
 1200 continue
 1100 continue
c

        end



        subroutine raddiag(ier,rad,verts,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension rad(1),xs(1),ys(1),whts(1),verts(2,3)
c
c       Return a quadrature for the evaluation of radially singular
c       functions of the form (1) over an arbitrary triangle
c       containing the origin.
c
c       Note that this routine will *fail* if the origin is too close 
c       to the boundary; the precise tolerances depend on the
c       precomputed quadrature table.
c
c                            Input Parameters:
c
c   rad - the radial quadrature structure generated by radial_init
c   verts - a (2,3) array each column of which gives the coordinates
c       of one vertex of the user-supplied triangle
c
c                           Output Parameters:
c
c   ier - an error return code;
c       ier = 0     indicates succesful execution
c
c       ier = 64    means that one of the necessary quadrature rules is 
c                   missing from the table in the rad structure; typically, 
c                   this  indicates that the triangle is very close to 
c                   degenerate
c
c       ier = 128   means that one of more of the vertices is *very*
c                   close to the origin
c
c       ier = 256   means the origin is *very* close to the boundary
c                   of the triangle
c
c       ier = 512   means that the vertices are nearly colinear
c
c
c   nquad - the size of the quadrature rule
c   (xs,ys) - the coordinates of the quadrature
c                   
c      
        ier   = 0
        nquad = 0
c
c       Extract the coordinates of the vertices.
c
        x1 = verts(1,1)
        y1 = verts(2,1)
c
        x2 = verts(1,2)
        y2 = verts(2,2)
c
        x3 = verts(1,3)
        y3 = verts(2,3)
c
c       Build quadratures for each subtriangle.
c
        call raddiag0(ier,rad,x1,y1,x2,y2,nquad1,xs(1),ys(1),whts(1),
     1    w,lw)
        if (ier .ne. 0) return
c
        call raddiag0(ier,rad,x1,y1,x3,y3,nquad2,xs(nquad1+1),
     1    ys(nquad1+1),whts(nquad1+1),w,lw)
        if (ier .ne. 0) return
c
        call raddiag0(ier,rad,x2,y2,x3,y3,nquad3,xs(nquad1+nquad2+1),
     1    ys(nquad1+nquad2+1),whts(nquad1+nquad2+1),w,lw)
        if (ier .ne. 0) return
c
        nquad=nquad1+nquad2+nquad3
c

        end
c
c
c
        subroutine raddiag0(ier,rad,x1_,y1_,x2_,y2_,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension rad(1),xs(1),ys(1),whts(1),verts(2,3)
        dimension amatr(2,2),ts(1000),twhts(1000)
c
c       Return a quadrature for evaluating radial singular functions
c       of the form (1) over a triangle one of whose vertices is
c       the origin.
c
c       Triangles of this type can be transformed via rotation, scaling,
c       and, if necessary, reflection into the canonical form shown
c       below:
c
c                              z3 = a e^{ib}
c                              0
c                             / \          
c                            /   \         0<a<1
c                           /     \        0<b<Pi
c                          /       \         
c                         /         \        
c                        /           \ 
c                       /             \
c                      0---------------0
c                   (0,0)            (1,0)
c
c
c       Note: the user is responsible for ensure that the ouput arrays
c       xs,ys, and whts are sufficiently large to hold the resulting
c       quadrature formula.
c
c                            Input Parameters:
c
c   rad - the radial structure as returned by radial_init
c   (x1,y1) - the coordinates of one of the vertices which is not the
c       origin
c   (x2,y2) - the coordinates of one of the other vertex which is not
c       the origin
c
c                           Output Parameters:
c
c   ier - an error return code;
c       ier = 0     indicates successful execution
c
c       ier = 64    means that one of the necessary quadrature rules is 
c                   missing from the table in the rad structure; typically, 
c                   this  indicates that the triangle is very close to 
c                   degenerate
c
c       ier = 128   means that one of the user-supplied vertices is
c                   the origin
c
c   nquad - the number of nodes in the resulting quadrature
c   (xs,ys) - the coordinates of the quadrature nodes
c   whts - the quadrature weights
c
        ier   = 0
        nquad = 0
c        
c
c       Fetch parameters from the rad structure.
c
        nlege     = rad(10) 
        ixslege   = rad(11) 
        iwhtslege = rad(12) 
c
c       Compute the radii of the vertices.
c
        rr1 = sqrt(x1_**2+y1_**2)
        rr2 = sqrt(x2_**2+y2_**2)
c      
c       Perform a sanity check ... make sure neither of the specified
c       vertices is 0.
c      
        if (rr1 .eq. 0 .OR. rr2 .eq. 0) then
        ier = 128
        return
        endif
c
c       Swap the two vertices in order to ensure (x2,y2) is the most
c       distant from the origin.
c      
        if (rr1 .gt. rr2) then
        x2 = x1_
        y2 = y1_
        x1 = x2_
        y1 = y2_
        else
        x1 = x1_
        y1 = y1_
        x2 = x2_
        y2 = y2_
        endif
c
c       Build a transformation matrix taking the vertex of largest
c       radius to (1,0).
c
        rr = (x2**2+y2**2)
c
        amatr(1,1) = x2/rr
        amatr(2,2) = x2/rr
        amatr(1,2) = y2/rr
        amatr(2,1) = -y2/rr
c
        u = amatr(1,1)*x1+amatr(1,2)*y1
        v = amatr(2,1)*x1+amatr(2,2)*y1
c
        if (v .lt. 0) then
        amatr(2,1)=-amatr(2,1)
        amatr(2,2)=-amatr(2,2)
        v=-v
        endif
c
c       Fetch the theta quadrature.
c
        r = sqrt(u**2+v**2)
        t = atan2(v,u)
c
        lw = 10 000 000
c
        call radfetch(ier,rad,r,t,nt,its,itwhts)
c
        if (ier .ne. 0) then
        call prin2("radfetch failed with r = *",r,1)
        call prin2("and t = *",t,1)
c
        print *,"r = ",r
        print *,"t = ",t
        stop
c
        ier = 64
        return
        endif
c
c       Call an auxillary routine to build the product quadrature.
c
        call raddiag1(nt,rad(its),rad(itwhts),
     1    nlege,rad(ixslege),rad(iwhtslege),
     2    amatr,nquad,xs,ys,whts,r,t)
c
        end

c
        subroutine raddiag1(nt,ts,twhts,nlege,xslege,whtslege,amatr,
     1     nquad,xs,ys,whts,a,b)
        implicit double precision (a-h,o-z)
        dimension ts(1),twhts(1),xslege(1),whtslege(1)
        dimension amatr(2,2),ainv(2,2)
        dimension xs(1),ys(1),whts(1)
c
        nquad=0
c
c        print *,"raddiag: ",a,b
c        print *,""
c
c       Build the product quadrature.
c
        do 1000 j=1,nt
        t    = ts(j)*b
        twht = twhts(j)*b
c
        r1 = 0
        r2 = a*sin(b)/(a*sin(b-t)+sin(t))
c
        alpha = (r2-r1)/2
        beta  = (r2+r1)/2
c
        do 1100 i=1,nlege
        r = xslege(i)*alpha+beta
        rwht = whtslege(i)*alpha
c
        x = r*cos(t)
        y = r*sin(t)
        wht = rwht*twht*r
c
        nquad=nquad+1
        xs(nquad)=x
        ys(nquad)=y
        whts(nquad)=wht
c
 1100 continue
 1000 continue
c
c       Transform the quadrature.
c
        call rad2x2inv(amatr,ainv,det)
        det=abs(1/det)
c
        do 2000 j=1,nquad
        x = xs(j)
        y = ys(j)
        wht = whts(j)
c
        u = ainv(1,1)*x+ainv(1,2)*y
        v = ainv(2,1)*x+ainv(2,2)*y
        wht = wht*det
c
        xs(j)=u
        ys(j)=v
        whts(j)=wht
 2000 continue
c
        end



        subroutine radfetch(ier,rad,r,t,n,ixs0,iwhts0)
        implicit double precision (a-h,o-z)
        dimension rad(1)
c
c       Return pointers to the nodes and weights of one of the 
c       precomputed quadrature formulae residing in the rad structure.
c
c                          Input Parameters:
c
c   rad - the structure returned by radial_init
c   r - the value of the radial parameter
c   t - the value of the angular parameter
c
c                         Output Parameters:
c
c   ier - an error return code;
c       ier = 0   indicates succcessful execution
c       ier = 64  means that the requested quadrature rule is mising
c                 from the table stored in the rad structure
c
c   n - the number of quadrature nodes
c   ixs0 - a pointer into the rad structure to the quadrature nodes
c   iwhts0 - a pointer into the rad structure to the quadrature weights
c 
c
        ier = 0
        n   = 0 
c
c       Fetch data from the structure's header.
c
        nr     =  rad(1) 
        nt     =  rad(2) 
        irs    =  rad(3) 
        its    =  rad(4) 
        ins    =  rad(5) 
        ixs    =  rad(6)
        iwhts  =  rad(7) 
        max    =  rad(8)
c
c       Call an auxillary routine to find the quadrature.
c
        call radfetch0(ier,max,nr,nt,rad(irs),rad(its),rad(ins),
     1     rad(ixs),rad(iwhts),n,ii,r,t)
c
        ixs0   = ixs+ii-1
        iwhts0 = iwhts+ii-1
c
        end
c
c
c
        subroutine radfetch0(ier,max,nr,nt,rs,ts,ns,xs,xwhts,
     1    n,ii,r,t)
        implicit double precision (a-h,o-z)
        dimension rs(2,nr),ts(2,nt)
        dimension ns(nt,nr),xs(max,nt,nr),whts(max,nt,nr)
c
        data eps / 1.0d-12 /
        ier = 0
c
c       Figure out which quadrature rule to return.
c
        do 1000 ir=1,nr
        if (rs(1,ir) .le. r .AND. r .le. rs(2,ir)+eps) goto 1100
 1000 continue
        ier = 64
        return
 1100 continue
c
        do 2000 it=1,nt       
        if (ts(1,it) .le. t .AND. t .le. ts(2,it)+eps) goto 2100
 2000 continue
        ier = 64
        return
 2100 continue
c
c       Get the number of nodes and set the pointer.
c
        n   = ns(it,ir)
        ii  = 1+max*(it-1)+max*nt*(ir-1)
c
        end


        subroutine rayintersect(ier,x1,y1,x2,y2,th,x,y)
        implicit double precision (a-h,o-z)
c
c       Find the point of intersection of the ray \theta=th
c       and the line segment between two points --- presuming
c       there is such an intersection.  
c
c                          Input Parameters:
c
c   (x1,y1) - coordinates of one point forming the line segment
c   (x2,y2) - coordinates of the second point forming the line segment
c   th - real number in the interval [-pi,pi] specifiying the angle
c       of the ray
c               
c                         Output Parameters:
c
c   ier - an error return code;
c       ier = 0   indicates successful execution
c       ier = 16  means that there were no intersections
c       ier = 64  means that there are an infinite number of 
c                 intersections
c
c  (x,y) - coordinates of the intersection
c
        ier = 0
c        eps0 = 1.0d-15
c
c       Compute some initial parameters.
c
        ct = cos(th)
        st = sin(th)
        dx = x2-x1
        dy = y2-y1
c
c       If the line segments are parallel, there are either no
c       intersections or an infinite number, depending on the 
c       direction of the ray.
c
        det = -dx*st+dy*ct
c
        if (abs(det) .le. 0) then
        ier=16
        th1 = atan2(y1,x1)        
        if (th1*th .gt. 0) ier = 64
        return
        endif
c
c       If not, then solve the linear system to check for 
c       intersection.
c
        t = (st*x1-ct*y1)/det
        s = (dy*x1-dx*y1)/det
c
        ier = 0
        if (0 .le. t .AND. t. le. 1 .AND. s .gt. 0) then
        x=(1-t)*x1+t*x2
        y=(1-t)*y1+t*y2
        return
        endif
c
        ier = 16
c
        end
c
c
c
        subroutine rad2x2inv(amatr,ainv,det)
        implicit double precision (a-h,o-z)
        dimension amatr(2,2),ainv(2,2)
c
        det = amatr(1,1)*amatr(2,2)-amatr(1,2)*amatr(2,1)
        ainv(1,1) =  amatr(2,2)/det
        ainv(2,2) =  amatr(1,1)/det
        ainv(1,2) = -amatr(1,2)/det
        ainv(2,1) = -amatr(2,1)/det
        end
c
c
c
        subroutine radinsort2(k,a,b)
        implicit double precision (a-h,o-z)
        dimension a(1),b(1)
c
        do 1000 i=2,k
        val=a(i)
        val2=b(i)
        j=i-1
        do 1100 while (j .ge. 1 .AND. a(j) .gt. val) 
        a(j+1)=a(j)
        b(j+1)=b(j)
        j=j-1
 1100 continue
        a(j+1)=val
        b(j+1)=val2
 1000 continue
        end
c
c
c
        subroutine radmove(k,a,b)
        implicit double precision (a-h,o-z)
        dimension a(1),b(1)
        do 1000 j=1,k
        b(j)=a(j)
 1000 continue
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code and beginning of the
c       code proper.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains code that adaptively generates a quadrature
c       for the evaluation of the integral of a user-specified radially 
c       singular function on a more or less arbitrary triangle 
c       containing the origin.  By radially singular function, we mean 
c       a function f(x,y) which admits representation as 
c
c                                f_{-1}(t)
c         f(r\cos(t),r\sin(t)) = ---------  + f_0(t) + f_1(t) r + ...     (1)
c                                    r       
c
c       with the f_j(t) analytic on a strip containing the real line.
c
c       In order to construct the quadrature, the integral is first
c       written in polar coordinates as
c
c             2\pi    R(t)
c         \int     \int    f(r\cos(t),r\sin(t)) r dr dt.
c             0        0
c
c       Here,  R(t) gives a parameterization of the boundary of the 
c       user-specified triangle in polar coordinates.  The inner
c       integral is approximated using a Legendre quadrature 
c       rules, while ``double exponential'' quadrature rules are used 
c       to evaluate the outer integral.  The orders of these quadratures
c       are adaptively increased until the approximations converge.
c
c       This code is a fairly reliable mechanism for the evaluation of
c       the singular integrals which arise in the discretization of 
c       the weakly singular boundary integral operators arising 
c       from the reformulation of elliptic boundary value problems
c       as integral equations.  
c
c       The following subroutines are user-callable:
c
c   radadap - adaptively construct a quadrature for the evaluation of
c        a user-specified radially singular function given over an
c        arbitrary triangle containing the origin
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



        subroutine radadap(ier,eps,verts,maxquad,nquad,xs,ys,
     -    whts,funuser,par1,par2,par3,par4,par5,par6,par7,par8)
        implicit double precision (a-h,o-z)
        dimension verts(2,3),xs(1),ys(1),whts(1)
        external funuser
c
c        Adaptively construct a quadrature for the evaluation of
c        a user-specified radially singular function given over an
c        arbitrary triangle containing the origin.
c
c                             Input Parameters:
c
c   eps - precision for the quadrature rule
c   verts - vertices of the triangle which comprises the integration 
c        domain
c   maxquad - the maximum number of quadrature nodes which the
c       user-specified arrays xs, ys and whts can accomdate
c
c   funuser - an external subroutine which returns the values of the
c       integrand
c
c       subroutine funuser(x,y,par1,par2,par3,par4,par5,par6,par7,par8,
c           val)
c
c       Return the value of the integrand f(x,y) in the parameter val.
c       The parameters par? are arbitrarily-typed and supplied by the 
c       user.
c
c                             Output Parameters:
c
c   ier - an error return code;
c
c       ier = 0   indicates successful execution
c       ier = 16  means that the user-specified maximum number of 
c                 quadrature nodes was exceeded before the procedure
c                 converged
c
c   nquad - the length of the adaptive 
c   (xs,ys) - these user-supplied arrays will contain the coordinates 
c       of the quadrature nodes
c   whts - this user-supplied array will contain the weights of the 
c       quadrature rule
c
        ier = 0
c
        x1 = verts(1,1)
        y1 = verts(2,1)
        x2 = verts(1,2)
        y2 = verts(2,2)
        x3 = verts(1,3)
        y3 = verts(2,3)
c
        nquad    = 0
        maxquad0 = maxquad
c
        call radadap_0(ier,eps,x1,y1,x2,y2,maxquad0,nquad0,
     -    xs(nquad+1),ys(nquad+1),whts(nquad+1),
     -    funuser,par1,par2,par3,par4,par5,par6,par7,par8)
c
        if (ier .ne. 0) return
        nquad    = nquad+nquad0
        maxquad0 = maxquad0-nquad0
c
        call radadap_0(ier,eps,x2,y2,x3,y3,maxquad0,nquad0,
     -    xs(nquad+1),ys(nquad+1),whts(nquad+1),
     -    funuser,par1,par2,par3,par4,par5,par6,par7,par8)
        if (ier .ne. 0) return
        nquad    = nquad+nquad0
        maxquad0 = maxquad0-nquad0
c
        call radadap_0(ier,eps,x3,y3,x1,y1,maxquad0,nquad0,
     -    xs(nquad+1),ys(nquad+1),whts(nquad+1),
     -    funuser,par1,par2,par3,par4,par5,par6,par7,par8)
        if (ier .ne. 0) return
        nquad    = nquad+nquad0
        maxquad0 = maxquad0-nquad0
c
        return
        end



        subroutine radadap_0(ier,eps,x1,y1,x2,y2,maxquad,nquad,xs,ys,
     -   whts,funuser,par1,par2,par3,par4,par5,par6,par7,par8)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1)
        double complex sum,sum0,val
        external funuser
c
        nlege   = 5
        n       = 10
c
        call radadap_00(ier,eps,nlege,n,x1,y1,x2,y2,maxquad,nquad,
     -   xs,ys,whts,funuser,par1,par2,par3,par4,par5,par6,par7,par8,
     -   sum0)
        if (ier .ne. 0) return
c
 1000 continue
        nlege = nlege + 1
        call radadap_00(ier,eps,nlege,n,x1,y1,x2,y2,maxquad,nquad,
     -   xs,ys,whts,funuser,par1,par2,par3,par4,par5,par6,par7,par8,
     -   sum)
        if (ier .ne. 0) return
c
        errabs = abs(sum-sum0)
c
c        print *,"     ",nlege,n,errabs
c
c$$$        print *,sum
c$$$        print *,sum0
c$$$        print *,""
c
        sum0   = sum
c
        if (errabs .gt. eps) goto 1000
c
        return
        end


        subroutine radadap_00(ier,eps,nlege,n,x1,y1,x2,y2,maxquad,nquad,
     -   xs,ys,whts,funuser,par1,par2,par3,par4,par5,par6,par7,par8,
     -   sum)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1)
        double complex sum,sum0,val
        external funuser
c
        if (maxquad .lt. 2*n*nlege) then
        ier = 16
        return
        endif
c
        call radadap_1(n,nlege,x1,y1,x2,y2,nquad,xs,ys,whts)
c
        sum0 = 0
        do 1000 i=1,nquad
        x = xs(i)
        y = ys(i)
        wht = whts(i)
c
        call funuser(x,y,par1,par2,par3,par4,par5,par6,par7,par8,val)
        sum0 = sum0 + val*wht
 1000 continue
c
 1100 continue
        n     = n*2
c
        if (maxquad .lt. 2*n*nlege) then
        ier = 16
        return
        endif
c
        call radadap_1(n,nlege,x1,y1,x2,y2,nquad,xs,ys,whts)
c
        sum = 0
        do 1400 i=1,nquad
        x = xs(i)
        y = ys(i)
        wht = whts(i)
c
        call funuser(x,y,par1,par2,par3,par4,par5,par6,par7,par8,val)
        sum = sum + val*wht
 1400 continue
c
        errabs = abs(sum-sum0)
        sum0   = sum
c
        print *,"   ",nlege,n,errabs
c
        if (errabs .gt. eps) goto 1100
c
        n = n/2
c
        end



        subroutine radadap_1(n,nlege,x1_,y1_,x2_,y2_,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1)
        dimension amatr(2,2),ainv(2,2)
c
c       Return a quadrature for the triangle with vertices (0,0),
c       (x1,y1) and (x2,y2).
c
c                           Input Parameters:
c
c   n - the number of nodes for the double exponential quadrature
c   nlege - the number of nodes for the Legendre quadrature
c   (x1,y1) - the coordinates of one vertex of the triangle which
c       comprises the integration domain
c   (x2,y2) - the coordinates of another of the vertices
c
c                            Output Parameters:
c
c   nquad - the number of nodes in the resulting quadrature
c   (xs,ys) - the location of the nodes in the resulting quadrature
c       rule
c   whts - the weights of the resulting quadrature rule
c
c
c        eps0  = 1.0d-15
        eps0 = 0
c
c       Compute the radii of the vertices.
c
        rr1 = sqrt(x1_**2+y1_**2)
        rr2 = sqrt(x2_**2+y2_**2)
c
c       Perform a sanity check ... make sure neither of the specified
c       vertices is 0.
c      
        if (rr1 .le. eps0 .OR. rr2 .le. eps0) then
        ier = 128
        return
        endif
c
c       Swap the two vertices in order to ensure (x2,y2) is the most
c       distant from the origin.
c      
        if (rr1 .ge. rr2) then
        x2 = x1_
        y2 = y1_
        x1 = x2_
        y1 = y2_
        else
        x1 = x1_
        y1 = y1_
        x2 = x2_
        y2 = y2_
        endif
c
c       Build a transformation matrix taking the vertex of largest
c       radius to (1,0).
c
        rr = (x2**2+y2**2)
c
        amatr(1,1) = x2/rr
        amatr(2,2) = x2/rr
        amatr(1,2) = y2/rr
        amatr(2,1) = -y2/rr
c
        u = amatr(1,1)*x1+amatr(1,2)*y1
        v = amatr(2,1)*x1+amatr(2,2)*y1
c
        if (v .lt. 0) then
        amatr(2,1)=-amatr(2,1)
        amatr(2,2)=-amatr(2,2)
        v=-v
        endif
c
c       Build a quadrature for the modified triangle.
c
        r0 = sqrt(u**2+v**2)
        t0 = atan2(v,u)
c
        call radadap_2(n,nlege,r0,t0,nquad,xs,ys,whts)
c
c       Transform the quadrature.
c
        call radadap_2x2inv(amatr,ainv,det)
        det=abs(1/det)
c
        do 2000 j=1,nquad
        x = xs(j)
        y = ys(j)
        wht = whts(j)
c
        u = ainv(1,1)*x+ainv(1,2)*y
        v = ainv(2,1)*x+ainv(2,2)*y
        wht = wht*det
c
        xs(j)=u
        ys(j)=v
        whts(j)=wht
 2000 continue
c        
        end



        subroutine radadap_2(n,nlege,r0,t0,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1),xslege(10 000),whtslege(10 000)
        data pi    / 3.14159265358979323846264338327950288d0 /
        data nlege0 / -1 /
        data isinit / -1 /
        save
c
c       Return a quadrature for radially singular functions given on a 
c       triangle with vertices
c 
c       (0,0), (1,0) and (r0 \cos(t0), r0 \sin(t0) ).
c
c       The quadraturee is constructed by first passing to polar
c       coordinates --- that is, rewriting the integral as 
c
c           2\pi    R(t)
c       \int    \int     f(r\cos(t),r\sin(t)) r drdt,                     (2)
c           0       0
c
c       and then using double exponential quadratures in the 
c       t variable and a Gaussian quadrature in the r variable.
c
c
c                             Input Parameters:
c
c   n - the number of points for the quadrature rule used to
c       evaluate the outer integral in (2)
c   nlege - the number of Legendre nodes used to evaluate the inner
c       integral in (2)
c
c   (r0,t0) - the location, in polar coordinates, of the third vertex
c       of the triangle comprising the integration domain
c
c                            Output Parameters:
c
c   nquad - the number of nodes in the resulting quadrature
c   (xs,ys) - the location of the nodes in the resulting quadrature
c       rule
c   whts - the weights of the resulting quadrature rule
c
        if (isinit .eq. -1) then
        isinit = 1
        call mach_zero(eps)
        dd = 2.9d0
        if (eps .lt. 1.0d-16) dd = 3.9d0
        endif
c
        if (nlege .ne. nlege0) then
        call legequad(nlege,xslege,whtslege)
        nlege0 = nlege
        endif
c
        nquad = 0
c
        h = n
        h = dd/n
c     
        do 1000 k=-n,n
        t    = (tanh(pi/2*sinh(k*h))+1)*t0/2
        twht = (h*pi/2*cosh(k*h)*1.0d0/cosh(pi/2*sinh(k*h))**2)*t0/2
        df   = r0*sin(t0)/(r0*sin(t0-t)+sin(t))
c
        do 1100 j=1,nlege
        r    = (xslege(j)+1)*df/2 
        rwht = whtslege(j)*df/2
        x   = r*cos(t)
        y   = r*sin(t)
        wht = r*rwht*twht
        nquad       = nquad+1
        xs(nquad)   = x
        ys(nquad)   = y
        whts(nquad) = wht
 1100 continue
 1000 continue
c
c
        end



        subroutine radadap_2x2inv(amatr,ainv,det)
        implicit double precision (a-h,o-z)
        dimension amatr(2,2),ainv(2,2)
c
        det = amatr(1,1)*amatr(2,2)-amatr(1,2)*amatr(2,1)
        ainv(1,1) =  amatr(2,2)/det
        ainv(2,2) =  amatr(1,1)/det
        ainv(1,2) = -amatr(1,2)/det
        ainv(2,1) = -amatr(2,1)/det
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code and beginning of the
c       code proper.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains code for the adaptive integration of complex-
c       valued radially singular functions given on triangles.  More 
c       specifically, it contains a subroutine which returns a 
c       quadrature for the evaluation of integrals of the form
c
c         \int   f(x,y) dxdy                                              (1)
c              T
c
c       where T is a triangle containing the origin and f(x,y)
c       is a user-specified function which admits representation as
c
c                                f_{-1}(t)
c         f(r\cos(t),r\sin(t)) = ---------  + f_0(t) + f_1(t) r + ...
c                                    r
c
c       with the f_j(t) analytic on the real line.
c
c       The quadrature is formed by first changing to polar coordinates;
c       that is, rewriting the integral as
c
c              2\pi       R(t)
c         \int       \int       f(r\cos(t), r\sin(t)) r drdt,
c              0          0
c
c       where R(t) gives a parameterization of the boundary of T in
c       polar coordinates, and then using tensor products of Legendre
c       quadrature rules.  The domain of integral in the t variable is, 
c       of course, partitioned so as to ensure that the discontinuities
c       in R(t) appear at the endpoints of the intervals.
c
c       The following subroutines are user-callable:
c
c    radadap2 - construct a quadrature for the evaluation of an integral
c       of the form (1) with f and T specified by the user
c
c    radgauss - construct a quadrature for the evaluation of an integral
c       of the form (1) using Gaussian rules of a user-specified order
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



        subroutine radgauss(ngauss,verts,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension verts(2,3),xs(1),ys(1),whts(1)
        dimension xslege(ngauss),whtslege(ngauss)
        external funuser
c
c       Return a quadratue for the evaluation of integrals of the 
c       form
c
c         \int   f(x,y) dxdy                                              (1)
c              T
c
c       where T is a triangle containing the origin and f(x,y) is a
c       complex-valued function which admits representation as
c
c                                f_{-1}(t)
c         f(r\cos(t),r\sin(t)) = ---------  + f_0(t) + f_1(t) r + ...
c                                   r
c
c       The quadrature is formed using a Gaussian rule of order
c       specified by the user.
c
c                               Input Parameters:
c
c   verts - vertices of the triangle
c   ngauss - the length of the Legendre rule to use
c
c                              Output Parameters:
c
c   nquad - the length of the adaptive 
c   (xs,ys) - these user-supplied arrays will contain the coordinates 
c       of the quadrature nodes
c   whts - this user-supplied array will contain the weights of the 
c       quadrature rule
c
c                                    
        nlege = ngauss
        call legequad(nlege,xslege,whtslege)
c
        x1 = verts(1,1)
        y1 = verts(2,1)
        x2 = verts(1,2)
        y2 = verts(2,2)
        x3 = verts(1,3)
        y3 = verts(2,3)
c
        nquad    = 0
c
        call radadap2_1(nlege,x1,y1,x2,y2,nquad0,xs(nquad+1),
     -     ys(nquad+1),whts(nquad+1))
        nquad    = nquad+nquad0
c
        call radadap2_1(nlege,x2,y2,x3,y3,nquad0,xs(nquad+1),
     -     ys(nquad+1),whts(nquad+1))
        nquad    = nquad+nquad0
c
        call radadap2_1(nlege,x3,y3,x1,y1,nquad0,xs(nquad+1),
     -     ys(nquad+1),whts(nquad+1))
        nquad    = nquad+nquad0
c        
        end



        subroutine radadap2(ier,eps,verts,maxquad,nquad,xs,ys,
     -    whts,funuser,par1,par2,par3,par4,par5,par6,par7,par8)
        implicit double precision (a-h,o-z)
        dimension verts(2,3),xs(1),ys(1),whts(1)
        external funuser
c
c       Return a quadratue rule for the evaluation of an integral
c       of the form
c
c         \int   f(x,y) dxdy                                              (1)
c              T
c
c       where T is a triangle containing the origin and f(x,y) is a 
c       user-specified complex-valued function which admits 
c       representation as
c
c                                f_{-1}(t)
c         f(r\cos(t),r\sin(t)) = ---------  + f_0(t) + f_1(t) r + ...
c                                    r       
c
c                             Input Parameters:
c
c   eps - precision for the quadrature rule
c   verts - vertices of the triangle
c   maxquad - the maximum number of quadrature nodes which the
c       user-specified arrays xs, ys and whts can accomdate
c
c   funuser - an external subroutine which returns the values of the
c       integrand
c
c       subroutine funuser(x,y,par1,par2,par3,par4,par5,par6,par7,par8,
c           val)
c
c       Return the value of the integrand f(x,y) in the parameter val.
c       The parameters par? are arbitrarily-typed and supplied by the 
c       user.
c
c                             Output Parameters:
c
c   ier - an error return code;
c
c       ier = 0   indicates successful execution
c       ier = 16  means that the user-specified maximum number of 
c                 quadrature nodes was exceeded before the procedure
c                 converged
c
c   nquad - the length of the adaptive 
c   (xs,ys) - these user-supplied arrays will contain the coordinates 
c       of the quadrature nodes
c   whts - this user-supplied array will contain the weights of the 
c       quadrature rule
c
        ier = 0
c
        x1 = verts(1,1)
        y1 = verts(2,1)
        x2 = verts(1,2)
        y2 = verts(2,2)
        x3 = verts(1,3)
        y3 = verts(2,3)
c
        nquad    = 0
        maxquad0 = maxquad
c
        call radadap2_0(ier,eps,x1,y1,x2,y2,maxquad0,nquad0,
     -    xs(nquad+1),ys(nquad+1),whts(nquad+1),
     -    funuser,par1,par2,par3,par4,par5,par6,par7,par8)
c
        if (ier .ne. 0) return
        nquad    = nquad+nquad0
        maxquad0 = maxquad0-nquad0
c
        call radadap2_0(ier,eps,x2,y2,x3,y3,maxquad0,nquad0,
     -    xs(nquad+1),ys(nquad+1),whts(nquad+1),
     -    funuser,par1,par2,par3,par4,par5,par6,par7,par8)
        if (ier .ne. 0) return
        nquad    = nquad+nquad0
        maxquad0 = maxquad0-nquad0
c
        call radadap2_0(ier,eps,x3,y3,x1,y1,maxquad0,nquad0,
     -    xs(nquad+1),ys(nquad+1),whts(nquad+1),
     -    funuser,par1,par2,par3,par4,par5,par6,par7,par8)
        if (ier .ne. 0) return
        nquad    = nquad+nquad0
        maxquad0 = maxquad0-nquad0
c
        end



        subroutine radadap2_0(ier,eps,x1,y1,x2,y2,maxquad,nquad,
     -    xs,ys,whts,funuser,par1,par2,par3,par4,par5,par6,par7,par8)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1)
        double complex sum,sum0,val
        external funuser
c
        ier    = 0
        nlege1 = 5
        nlege2 = 5
c
        call radadap2_00(ier,nlege1,eps,x1,y1,x2,y2,maxquad,
     -    nquad,xs,ys,whts,funuser,par1,par2,par3,par4,par5,par6,par7,
     -    par8,nlege2)
        if (ier .ne. 0) return
c
        sum0 = 0
        do 1000 i=1,nquad
        x   = xs(i)
        y   = ys(i)
        wht = whts(i)
c
        call funuser(x,y,par1,par2,par3,par4,par5,par6,par7,par8,val)
        sum0 = sum0 + val*wht
 1000 continue
c
 1100 continue
c
c       Increase the number of points.
c
        nlege1 = nlege1+1
c
c       Compute the integral using the new quadrature.
c
        call radadap2_00(ier,nlege1,eps,x1,y1,x2,y2,maxquad,
     -    nquad,xs,ys,whts,funuser,par1,par2,par3,par4,par5,par6,par7,
     -    par8,nlege2)
        if (ier .ne. 0) return
c
        sum = 0
        do 1200 i=1,nquad
        x   = xs(i)
        y   = ys(i)
        wht = whts(i)
c
        call funuser(x,y,par1,par2,par3,par4,par5,par6,par7,par8,val)
c
        sum = sum + val*wht
 1200 continue
c
        errabs = abs(sum-sum0)
        sum0   = sum
c
c        print *,"   ",nlege1,errabs
c
c       If the error is too large, try again.
c
        if (errabs .gt. eps) goto 1100
c
c        nlege1=nlege1+8
c        nlege2=nlege2+8
c
c        call radadap2_1(nlege1,nlege2,x1,y1,x2,y2,
c     -    nquad,xs,ys,whts)
c        
        end


        subroutine radadap2_00(ier,nlege1,eps,x1,y1,x2,y2,maxquad,
     -    nquad,xs,ys,whts,funuser,par1,par2,par3,par4,par5,par6,par7,
     -    par8,nlege2)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1)
        double complex sum,sum0,val
        external funuser
c
        ier    = 0
        nlege2 = nlege2-1
c
c       Compute the value of the integral using nlege points.
c
        if (nlege1*nlege2 .gt. maxquad) then
        ier = 16
        return
        endif
c
        call radadap2_1(nlege1,nlege2,x1,y1,x2,y2,nquad,xs,ys,whts)
c
        sum0 = 0
        do 1000 i=1,nquad
        x   = xs(i)
        y   = ys(i)
        wht = whts(i)
c
        call funuser(x,y,par1,par2,par3,par4,par5,par6,par7,par8,val)
        sum0 = sum0 + val*wht
 1000 continue
c
 1100 continue
c
c       Increase the number of points.
c
        nlege2 = nlege2+1
        if (nlege1*nlege2 .gt. maxquad) then
        ier = 16
        return
        endif
c
c       Compute the integral using the new quadrature.
c
        call radadap2_1(nlege1,nlege2,x1,y1,x2,y2,nquad,xs,ys,whts)
c
        sum = 0
        do 1200 i=1,nquad
        x   = xs(i)
        y   = ys(i)
        wht = whts(i)
c
        call funuser(x,y,par1,par2,par3,par4,par5,par6,par7,par8,val)
c
        sum = sum + val*wht
 1200 continue
c
        errabs = abs(sum-sum0)
        sum0   = sum
c
c        print *,"   ",nlege1,nlege2,errabs
c
c       If the error is too large, try again.
c
        if (errabs .gt. eps) goto 1100
c
        end



        subroutine radadap2_1(nlege1,nlege2,x1_,y1_,x2_,y2_,
     -    nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1)
        dimension amatr(2,2),ainv(2,2)
c
c       Return a quadrature for the triangle with vertices (0,0),
c       (x1,y1) and (x2,y2).
c
c                           Input Parameters:
c
c
c   nlege - order of the Legendre rule to use
c   xslege - the nodes of the nlege-point Legendre quadrature
c   whtslege - the weights of the nlege-point Legendre quadrature
c
c   (x1,y1) - the coordinates
c
c   (x2,y2) - 
c
c                            Output Parameters:
c
c   nquad - the number of nodes in the resulting quadrature
c   (xs,ys) - the location of the nodes in the resulting quadrature
c       rule
c   whts - the weights of the resulting quadrature rule
c
c
        eps0  = 1.0d-15
        eps0  = 0
c
c       Compute the radii of the vertices.
c
        rr1 = sqrt(x1_**2+y1_**2)
        rr2 = sqrt(x2_**2+y2_**2)
c
c       Perform a sanity check ... make sure neither of the specified
c       vertices is 0.
c      
        if (rr1 .le. eps0 .OR. rr2 .le. eps0) then
        ier = 128
        return
        endif
c
c       Swap the two vertices in order to ensure (x2,y2) is the most
c       distant from the origin.
c      
        if (rr1 .ge. rr2) then
        x2 = x1_
        y2 = y1_
        x1 = x2_
        y1 = y2_
        else
        x1 = x1_
        y1 = y1_
        x2 = x2_
        y2 = y2_
        endif
c
c       Build a transformation matrix taking the vertex of largest
c       radius to (1,0).
c
        rr = (x2**2+y2**2)
c
        amatr(1,1) = x2/rr
        amatr(2,2) = x2/rr
        amatr(1,2) = y2/rr
        amatr(2,1) = -y2/rr
c
        u = amatr(1,1)*x1+amatr(1,2)*y1
        v = amatr(2,1)*x1+amatr(2,2)*y1
c
        if (v .lt. 0) then
        amatr(2,1)=-amatr(2,1)
        amatr(2,2)=-amatr(2,2)
        v=-v
        endif
c
c       Build a quadrature for the modified triangle.
c
        r0 = sqrt(u**2+v**2)
        t0 = atan2(v,u)
c
        call radadap2_2(nlege1,nlege2,r0,t0,nquad,xs,ys,whts)
c
c       Transform the quadrature.
c
        call radadap2_2x2inv(amatr,ainv,det)
        det=abs(1/det)
c
        do 2000 j=1,nquad
        x = xs(j)
        y = ys(j)
        wht = whts(j)
c
        u = ainv(1,1)*x+ainv(1,2)*y
        v = ainv(2,1)*x+ainv(2,2)*y
        wht = wht*det
c
        xs(j)=u
        ys(j)=v
        whts(j)=wht
 2000 continue
c        
        end


        subroutine radadap2_2(nlege1,nlege2,r0,t0,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1)
        dimension xslege1(nlege1),whtslege1(nlege1)
        dimension xslege2(nlege2),whtslege2(nlege2)
c
c       Return a quadrature for radially singular functions given on a 
c       triangle with vertices
c 
c       (0,0), (1,0) and (r0 \cos(t0), r0 \sin(t0) ).
c
c       The quadrature is constructed using product Legendre rules
c       of a user-specified order.
c
c                             Input Parameters:
c
c   nlege - order of the Legendre rule to use
c   xslege - the nodes of the nlege-point Legendre quadrature
c   whtslege - the weights of the nlege-point Legendre quadrature
c
c   (r0,t0) - the location, in polar coordinates, of the third vertex
c       of the triangle comprising the integration domain
c
c                            Output Parameters:
c
c   nquad - the number of nodes in the resulting quadrature
c   (xs,ys) - the location of the nodes in the resulting quadrature
c       rule
c   whts - the weights of the resulting quadrature rule
c
        call legequad1(nlege1,xslege1,whtslege1)
        call legequad2(nlege2,xslege2,whtslege2)
c
        nquad = 0
c
        do 1000 i=1,nlege2
        t    = (xslege2(i)+1)*t0/2
        twht = whtslege2(i)*t0/2
c
        do 1100 j=1,nlege1
c
        df = r0*sin(t0)/(r0*sin(t0-t)+sin(t))
c
        r    = (xslege1(j)+1)*df/2 
        rwht = whtslege1(j)*df/2
c
        x   = r*cos(t)
        y   = r*sin(t)
        wht = r*rwht*twht
c
        nquad       = nquad+1
        xs(nquad)   = x
        ys(nquad)   = y
        whts(nquad) = wht
 1100 continue
 1000 continue
c
c
        end



        subroutine radadap2_2x2inv(amatr,ainv,det)
        implicit double precision (a-h,o-z)
        dimension amatr(2,2),ainv(2,2)
c
        det = amatr(1,1)*amatr(2,2)-amatr(1,2)*amatr(2,1)
        ainv(1,1) =  amatr(2,2)/det
        ainv(2,2) =  amatr(1,1)/det
        ainv(1,2) = -amatr(1,2)/det
        ainv(2,1) = -amatr(2,1)/det
        end

        

