cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code and beginning of the
c       code proper.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains code for generating quadratures for the
c       evaluation of integrals of the form
c
c           \int    K(x,y) \sigma(y) dS(y)                                (1)
c               S
c       where:
c
c           S is a surface element and dS is the corresponding
c           surface measure,
c
c           \sigma is a smooth function.
c           
c           K(x,y) is the double or single layer potential on S, and
c 
c           x is a specified point in S.
c
c       It is assumed that the surface element S is the image under
c       a parameterization
c
c          p: T \to \mathbb{R}^3
c
c       given over a triangle T.  
c
c       The behavior of the Jacobian dp of the parameterization
c       p at the point x has a large influence on the form of the 
c       integrand of (1).  This code proceeds by first composing the 
c       given parameterization p with an appropriate linear mapping 
c
c                 A: \mathbb{R}^2 \to \mathbb{R}^2
c
c       in order to  form a new parameterization p' = p \ocirc A such 
c       that the Jacobian of p' at the point x is conformal.  Then
c       the quadrature rules from radial.f are used to evaluate
c       the resulting integral.
c
c       The quadratures returned by raddiag are formed using precomputed
c       quadrature tables stored on the disk.  These precomputed
c       tables determine the possible orders for the polynomials p
c       and q in (2).  See radial_init for a list of possible orders.
c
c       The following subroutines are user-callable:
c
c
c   self_quadrature_new - return a quadrature for evaluating an integral of
c       the form (1), with no entries or save statements...
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



        subroutine self_quadrature_new(ier, rad, verts0, x0, y0, dr,
     1    nquad, xs, ys, whts)
        implicit double precision (a-h,o-z)
        real *8 rad(1)
        dimension verts(2,3),dr(3,2),a(2,2),ainv(2,2),verts0(2,3)
        dimension xs(1),ys(1),whts(1),b(3,2)
c
c        dimension rad(500 000)
c        double precision, allocatable :: rad(:)
c        save rad
c
c       Return a quadrature formula for evaluating integrals of the
c       form (1).
c
c                            Input Parameters:
c
c    verts - a (2,3) matrix which specifys the vertices of the triangle
c       over which the surface element is parameterized
c    (x0,y0) - the coordinates (in the parameterization variables)
c       of the target point x
c    dr - the (3,2) Jacobian of the parameterization at (x0,y0)
c
c                           Output Parameters:
c
c    ier - an error return code;
c       ier = 0     indicates successful execution
c       ier = 128   means that the quadrature formula could not be
c                   constructed; this usually means 
c
c    nquad - the number of nodes in the resulting quadrature formula
c    xs - this user-supplied array will contain the x coordinates of 
c       the quadrature nodes upon return
c    ys - this user-supplied array will contain the y coordinates of 
c       the quadrature nodes upone return
c    whts - this user-supplied array will contain the quadrature weights
c       upon return
c
c
        ier = 0
c
c
c       Find a linear map A which makes the Jacobian at the target node
c       conformal.  Also, compute the inverse of A.
c
        call self_findmap(dr,a,ainv)
c
c
c       Apply the inverse of A to the vertices of the triangle T in
c       order to find the vertices of the triangle T_0.
c
        do 2000 i=1,3
        x = verts0(1,i)-x0
        y = verts0(2,i)-y0
c
        xx = ainv(1,1)*x+ainv(1,2)*y
        yy = ainv(2,1)*x+ainv(2,2)*y
c
        verts(1,i) = xx
        verts(2,i) = yy
 2000 continue
c
c       Fetch a quadrature on T_0 for radially singular functions.
c
        call raddiag(jer,rad,verts,nquad,xs,ys,whts)
c
        if (jer .ne. 0) then
        ier = 128
        return
        endif
c
c       Apply the mapping A to the quadrature formula.
c
        det = abs(a(1,1)*a(2,2)-a(1,2)*a(2,1))
c
        sum=0
c
        do 3000 i=1,nquad
        s   = xs(i)
        t   = ys(i)
        wht = whts(i)
c       
        u = a(1,1)*s + a(1,2)*t
        v = a(2,1)*s + a(2,2)*t
        wht = wht*det
c
        sum=sum+wht
c
        xs(i)   = u+x0
        ys(i)   = v+y0
        whts(i) = wht
 3000 continue
c
        return
        end
