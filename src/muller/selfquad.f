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
c   self_quadrature_init - this routine must be called before 
c       self_quadrature
c
c   self_quadrature - return a quadrature for evaluating an integral of
c       the form (1)
c
c   selfadap - construct an adaptive quadrature for evaluating an 
c       integral of the form (1) using double exponential quadrature
c       rules
c
c   selfadap2 - construct an adaptive quadrature for evaluating an 
c       integral of the form (1) using Gaussian quadrature rules
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine selfadap2(ier,eps,ifconformal,verts0,x0,y0,dr,
     -    maxquad,nquad,xs,ys,whts,funuser,par1,par2,par3,par4,par5)
        implicit double precision (a-h,o-z)
        dimension verts(2,3),dr(3,2),xs(1),ys(1),whts(1)
        dimension verts0(2,3),a(2,2),ainv(2,2)
        external funuser,funself
c       
c       Adaptively construct a quadrature for the evaluation of an 
c       integral of the form (1).
c
c                             Input Parameters:
c
c   eps - precision for the quadrature
c   ifconformal - an integer parameter indicating whether or not a 
c       mapping which makes the user-specified parameterization 
c       conformal at the targe point 
c
c       ifconformal = 0  do not apply a mapping
c       ifconformal = 1  apply a linear mapping which makes the 
c
c   verts0 - a (2,3) array each of column gives the coordinates of a 
c       vertex of the triangle T comprising the integration domain
c   (x0,y0) - the coordinates of the target point, which must be 
c       contained in the triangle T
c   dr - the Jacobian of the user-specified parameterization at the
c       target point
c   maxquad - the size of the largest quadrature which the user-
c       supplied arrays xs, ys and whts can accomodate
c
c   funuser - an external subroutine which supplies the values of the
c       integrand in (1)
c       
c
c       subroutine funuser(x,y,x0,y0,par1,par2,par3,par4,par5,val)
c
c       Return the value of the integrand in (1) at the point (x,y) in
c       the parameter val.  The coordinates of the target point are
c       (x0,y0).
c
c
c                             Output Parameters:
c
c   ier - an error return code;
c       ier = 0   indicates successful execution
c       ier = 16  means that the maximum number of quadrature nodes
c
c   nquad - the size of the quadrature
c   (xs,ys) - these user-supplied arrays will contain the coordinates of
c       the nodes of the quadrature formula
c   whts - these user-supplied arrays will contain the weights of the
c       resulting quadrature formula
c     
c
        ier = 0
c
c       Find a linear map A which makes the Jacobian at the target node
c       conformal.  Also, compute the inverse of A.
c
        if (ifconformal .eq. 1) then
        call self_findmap(dr,a,ainv)
        else
        a(1,1) = 1
        a(2,2) = 1
        a(2,1) = 0
        a(1,2) = 0
        ainv(1,1) = 1
        ainv(2,2) = 1
        ainv(2,1) = 0
        ainv(1,2) = 0
        endif
c
c       Apply the inverse of A to the vertices of the triangle T in
c       order to find the vertices of the triangle T_0.
c
        x00 = ainv(1,1)*x0 + ainv(1,2)*y0
        y00 = ainv(2,1)*x0 + ainv(2,2)*y0
c
        do 2000 i=1,3
        x = verts0(1,i)
        y = verts0(2,i)
c
        xx = ainv(1,1)*x+ainv(1,2)*y
        yy = ainv(2,1)*x+ainv(2,2)*y
c
        verts(1,i) = xx-x00
        verts(2,i) = yy-y00
 2000 continue                
c
c       Build a quadrature.
c
        call radadap2(ier,eps,verts,maxquad,nquad,xs,ys,
     -    whts,funself,funuser,x00,y00,par1,par2,par3,par4,par5)
c
        if (ier .ne. 0) return
c
c       Map the quadrature.
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
        end



        subroutine selfadap(ier,eps,ifconformal,verts0,x0,y0,dr,
     -    maxquad,nquad,xs,ys,whts,funuser,par1,par2,par3,par4,par5)
        implicit double precision (a-h,o-z)
        dimension verts(2,3),dr(3,2),xs(1),ys(1),whts(1)
        dimension verts0(2,3),a(2,2),ainv(2,2)
        external funuser,funself
c       
c       Adaptively construct a quadrature for the evaluation of an 
c       integral of the form (1).
c
c                             Input Parameters:
c
c   eps - precision for the quadrature
c   ifconformal - an integer parameter indicating whether or not a 
c       mapping which makes the user-specified parameterization 
c       conformal at the targe point 
c
c       ifconformal = 0  do not apply a mapping
c       ifconformal = 1  apply a linear mapping which makes the 
c
c   verts0 - a (2,3) array each of column gives the coordinates of a 
c       vertex of the triangle T comprising the integration domain
c   (x0,y0) - the coordinates of the target point, which must be 
c       contained in the triangle T
c   dr - the Jacobian of the user-specified parameterization at the
c       target point
c   maxquad - the size of the largest quadrature which the user-
c       supplied arrays xs, ys and whts can accomodate
c
c   funuser - an external subroutine which supplies the values of the
c       integrand in (1)
c       
c
c       subroutine funuser(x,y,x0,y0,par1,par2,par3,par4,par5,val)
c
c       Return the value of the integrand in (1) at the point (x,y) in
c       the parameter val.  The coordinates of the target point are
c       (x0,y0).
c
c
c                             Output Parameters:
c
c   ier - an error return code;
c       ier = 0   indicates successful execution
c       ier = 16  means that the maximum number of quadrature nodes was
c                 exceeded
c
c   nquad - the size of the quadrature
c   (xs,ys) - these user-supplied arrays will contain the coordinates of
c       the nodes of the quadrature formula
c   whts - these user-supplied arrays will contain the weights of the
c       resulting quadrature formula
c     
c
        ier = 0
c
c       Find a linear map A which makes the Jacobian at the target node
c       conformal.  Also, compute the inverse of A.
c
        if (ifconformal .eq. 1) then
        call self_findmap(dr,a,ainv)
        else
        a(1,1) = 1
        a(2,2) = 1
        a(2,1) = 0
        a(1,2) = 0
        ainv(1,1) = 1
        ainv(2,2) = 1
        ainv(2,1) = 0
        ainv(1,2) = 0
        endif
c
c       Apply the inverse of A to the vertices of the triangle T in
c       order to find the vertices of the triangle T_0.
c
        x00 = ainv(1,1)*x0 + ainv(1,2)*y0
        y00 = ainv(2,1)*x0 + ainv(2,2)*y0
c
        do 2000 i=1,3
        x = verts0(1,i)
        y = verts0(2,i)
c
        xx = ainv(1,1)*x+ainv(1,2)*y
        yy = ainv(2,1)*x+ainv(2,2)*y
c
        verts(1,i) = xx-x00
        verts(2,i) = yy-y00
 2000 continue                
c
c       Build a quadrature.
c
        call radadap(ier,eps,verts,maxquad,nquad,xs,ys,
     -    whts,funself,funuser,x00,y00,par1,par2,par3,par4,par5)
c
        if (ier .ne. 0) return
c
c       Map the quadrature.
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
        end


        subroutine funself(x,y,funuser,x0,y0,par1,par2,par3,par4,par5,
     -    val)
        implicit double precision (a-h,o-z)
        double complex val
        external funuser
c
        call funuser(x0+x,y0+y,x0,y0,par1,par2,par3,par4,par5,val)
c
        end



        subroutine self_quadrature(ier,verts0,x0,y0,dr,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension verts(2,3),dr(3,2),a(2,2),ainv(2,2),verts0(2,3)
        dimension xs(1),ys(1),whts(1),b(3,2)
c
c        dimension rad(500 000)
        double precision, allocatable :: rad(:)
        save rad
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
c
        entry self_quadrature_init(ier0,norder)
c
c       This initialization routine must be called before the
c       self_quadrature subroutine.  The user must specify the
c       order for the 
c
c                          Input Parameters:
c
c   norder - order for the self-interaction quadratures;
c       right now, the possible orders are 4, 8, 12, 16 and 20
c
c                          Output Parameters:
c
c   ier0 - an error return code;
c       ier0 = 0     indicates successful executionc
c       ier0 = 128   means that the quadrature table for the specified
c                    order could not be found
c       ier0 = 1024  means that their was an I/O error reading the
c                    table from the disk; this usually means that the
c                    data file containing it was improperly formatted
c
        ier0 = 0
        lrad = 500 000
        if(allocated(rad)) deallocate(rad)
        allocate(rad(lrad))
c
        call radial_init(jer,norder,rad,lrad,lkeep)
        if (jer .ne. 0) then
        ier0 = 1024
        return
        endif
c
        return
        end


        subroutine self_findmap(a,b,binv)
        implicit double precision (a-h,o-z)
        dimension a(3,2),b(2,2),binv(2,2),q(3,2)
c
c       Given a (3,2) matrix A, return a (2,2) matrix B such that
c       A*B has orthonormal columns.  Also, return the inverse
c       BINV of the mapping B.
c
c
c                             --- WARNING ---
c       THIS ROUTINE IS A BIT ROUGH --- IT IS EXPECTED TO FAIL
c       IN THE EVENT THAT THE MATRIX A IS NEARLY SINGULAR.  BUT
c       WHAT ARE YOU DOING USING PARAMETERIZATION WITH JACOBIANS
c       WHICH ARE NEARLY SINGULAR ANYWAY?
c                              --------------
c
c
c                          Input Parameters:
c
c   a - the (3,2) input matrix 
c
c                         Output Parameters:
c
c   b - a (2,2) matrix such that a*b has orthonormal columns
c   binv - the (2,2) inverse of b
c
c
c       Orthonormalize the columns of the matrix a.
c
        sum1 = sqrt(a(1,1)**2 + a(2,1)** 2 + a(3,1)**2)
        q(1,1) = a(1,1) / sum1
        q(2,1) = a(2,1) / sum1
        q(3,1) = a(3,1) / sum1
c
        dip = a(1,2)*q(1,1)+a(2,2)*q(2,1)+a(3,2)*q(3,1)
c
        q(1,2) = a(1,2) - dip * q(1,1)
        q(2,2) = a(2,2) - dip * q(2,1)
        q(3,2) = a(3,2) - dip * q(3,1)
c
        sum2 = sqrt(q(1,2)**2 + q(2,2)** 2 + q(3,2)**2)
        q(1,2) = q(1,2) / sum2
        q(2,2) = q(2,2) / sum2
        q(3,2) = q(3,2) / sum2
c
c       Compute BINV = Q'*A.
c     
        binv(1,1) = q(1,1)*a(1,1) + q(2,1)*a(2,1) + q(3,1)*a(3,1)
        binv(1,2) = q(1,1)*a(1,2) + q(2,1)*a(2,2) + q(3,1)*a(3,2)
        binv(2,1) = q(1,2)*a(1,1) + q(2,2)*a(2,1) + q(3,2)*a(3,1)
        binv(2,2) = q(1,2)*a(1,2) + q(2,2)*a(2,2) + q(3,2)*a(3,2)
c
c       Compute B = (BINV)^(-1).
c
        det = binv(1,1)*binv(2,2) - binv(1,2)*binv(2,1)
        b(1,1) = binv(2,2)/det
        b(2,2) = binv(1,1)/det
        b(1,2) = -binv(1,2)/det
        b(2,1) = -binv(2,1)/det
c
        end

