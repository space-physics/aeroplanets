*
*----------------------------------------------------------------------|
*
      subroutine dgphiv( n, m, t, u, v, w, tol, anorm,
     .                   wsp,lwsp, iwsp,liwsp, itrace,iflag,
     .                   valcsr1, ipr1, jpr1) 

*----------------------------------------------------------------------|

      implicit none
      integer n, m, lwsp, liwsp, itrace, iflag, iwsp(liwsp)
!      integer ipr1(n+1), jpr1(*)
      integer ipr1(*), jpr1(*)
      double precision  valcsr1(*)
      double precision t, tol, anorm, u(n), v(n), w(n), wsp(lwsp)
*
*-----Purpose----------------------------------------------------------|
*
*---  DGPHIV computes w = exp(t*A)v + t*phi(tA)u which is the solution 
*     of the nonhomegenous linear ODE problem w' = Aw + u, w(0) = v.
*     phi(z) = (exp(z)-1)/z and A is a General matrix.
*
*     The techniques used are based on Krylov subspace projection 
*     methods (see reference below). The matrix under consideration 
*     interacts (implicity) via the external routine `matvec' perfoming
*     the matrix-vector product.
*
*------------------------MODIFICATION
* THE MATVEC, WHICH WAS INITIALLY AN ARGUMENT, HAS BEEN MODIFIED,
* THE ORDER ARGUMENT, THE IPR1 AND JPR1 PARAMETERS HAVE BEEN ADDED 
*
*
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*                      
*     m      : (input) maximum size for the Krylov basis.
*                      
*     t      : (input) time at wich the solution is needed (can be < 0).
*   
*     u(n)   : (input) operand vector with respect to the phi function
*              (forcing term of the ODE problem).
*
*     v(n)   : (input) operand vector with respect to the exp function
*              (initial condition of the ODE problem).
*  
*     w(n)   : (output) computed approximation of exp(t*A)v + t*phi(tA)u 
* 
*     tol    : (input/output) the requested accuracy tolerance on w. 
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+3)^2+4*(m+3)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*                   
* iwsp(liwsp): (workspace) liwsp .ge. m+3
*
*     itrace : (input) running mode. 0=silent, 1=print step-by-step info
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments 
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the 
*     computations. They are located in the workspace arrays wsp and 
*     iwsp as indicated below: 
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = istep, number of integration steps used up to completion 
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = dummy
*     wsp(4)  = dummy
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxstep, mxreject, ideg
      double precision delta, gamma
      parameter( mxstep   = 500, 
     .           mxreject = 0,
     .           ideg     = 6, 
     .           delta    = 1.2d0,
     .           gamma    = 0.9d0 )

*     mxstep  : maximum allowable number of integration steps.
*               The value 0 means an infinite number of steps.
* 
*     mxreject: maximum allowable number of rejections at each step. 
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as 
*               an approximation to exp(H).
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*-----Author-----------------------------------------------------------|
*
*     Roger B. Sidje (rbs@maths.uq.oz.au)     
*     Last update: April 3, 1996
*
*-----Reference--------------------------------------------------------|
*
*     R.B. Sidje, EXPOKIT. Software Package for Computing Matrix
*     Exponentials. Technical Report, Department of Mathematics 
*     University of Queensland - Brisbane QLD 4072, Australia.
*     April, 1996.
*
*----------------------------------------------------------------------|
*
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        istep, iphih
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 avnorm, hj1j, hij, SQR1

      intrinsic AINT,ABS,DBLE,LOG10,MAX,MIN,NINT,SIGN,SQRT
      double precision DDOT, DNRM2
      integer stoppons
*---  check restrictions on input parameters ...
      iflag = 0
      if ( lwsp.lt.n*(m+3)+5*(m+3)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+3 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) stop 'bad sizes (in input of DGPHIV)'
*
*---  initialisations ...
*
      k1 = 3
      mh = m + 3
      iv = 1 
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      sgn      = SIGN( 1.0d0,t )
      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      istep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm
 
      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

*
*---  step-by-step integration ...
*
      SQR1 = SQRT( 0.1d0 )
      call DCOPY( n, v,1, w,1 )

 100  if ( t_now.ge.t_out ) goto 500

      nmult =  nmult + 1
      call matvec( w, wsp(iv), n, valcsr1, ipr1, jpr1 )
      call DAXPY( n, 1.0d0, u,1, wsp(iv),1 )
      beta = DNRM2( n, wsp(iv),1 )
!      do i = 1, n
!        print*, beta, wsp(iv + i)
!        enddo
      if ( beta.eq.0.0d0 ) goto 500
      call DSCAL( n, 1.0d0/beta, wsp(iv),1 )
      do i = 1,mh*mh
         wsp(ih+i-1) = 0.0d0
      enddo
      if ( istep.eq.0 ) then
*---     obtain the very first stepsize ...
         xm = 1.0d0/DBLE( m )
         p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
c         print*, "The norm:", anorm,"Beta:", beta
         t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
         p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
         t_new = AINT( t_new/p1 + 0.55d0 ) * p1
      endif
      istep = istep + 1
      t_step = MIN( t_out-t_now, t_new )
*
*---  Arnoldi loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v), n, valcsr1, ipr1, jpr1 )
         do i = 1,j
            hij = DDOT( n, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            call DAXPY( n, -hij, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            wsp(ih+(j-1)*mh+i-1) = hij
         enddo
         hj1j = DNRM2( n, wsp(j1v),1 )
         wsp(ih+(j-1)*mh+j) = hj1j
         if ( hj1j.le.break_tol ) then
            print*,'happy breakdown: ','mbrkdwn =',j,' h =',hj1j
            mbrkdwn = j
            k1 = 0
            ibrkflag = 1
            goto 300
         endif
         call DSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
 300  continue
*
*---  set 1's for the 3-extended scheme ...
*
      wsp(ih+(m-1)*mh+m) = 0.0d0
      wsp(ih+mh*mbrkdwn) = 1.0d0
      do i = 1,k1-1
         wsp(ih+(mbrkdwn+i)*mh+mbrkdwn+i-1) = 1.0d0
      enddo
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 
 400  continue
*---  if `happy breakdown' go straightforward at the end ... 
      if ( k1.eq.0 ) then
         t_step = t_out-t_now
         tbrkdwn = t_now
       endif

 401  continue
*
*---  compute w = beta*t_step*V*phi(t_step*H)*e1 + w
*
      nexph = nexph + 1
      mx = mbrkdwn + k1

*---  irreducible rational Pade approximation ...
      call DGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .              wsp(ifree),lfree, iwsp, iexph, ns, iflag )
      iexph = ifree + iexph - 1
      iphih = iexph + mbrkdwn*mx
      nscale = nscale + ns
      wsp(iphih+mbrkdwn)   = hj1j*wsp(iphih+mx+mbrkdwn-1)
      wsp(iphih+mbrkdwn+1) = hj1j*wsp(iphih+2*mx+mbrkdwn-1)
      mx = mbrkdwn + MAX( 0,k1-2 )
      call DGEMV( 'n', n,mx,beta,wsp(iv),n,wsp(iphih),1,1.0d0,w,1 )
 
 402  continue
* 
*---  error estimate ...
*
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         if ( ireject.eq.0 ) then
            nmult = nmult + 1
            call matvec( wsp(j1v-n), wsp(j1v), n, valcsr1, ipr1, jpr1 )
            avnorm = DNRM2( n, wsp(j1v),1 )
         endif
         p1 = ABS( wsp(iphih+m) )   * beta
         p2 = ABS( wsp(iphih+m+1) ) * beta * avnorm 
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
*old            xm = 1.0d0/DBLE( m )
                xm = 1.0d0/DBLE( m+1 )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
*old            xm = 1.0d0/DBLE( m )
                xm = 1.0d0/DBLE( m+1 )
         else
            err_loc = p1
*old            xm = 1.0d0/DBLE( m-1 )
                xm = 1.0d0/DBLE( m )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*   
      if ( (err_loc.gt.delta*t_step*tol) .and. 
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ne.0 ) then
            print*,'t_step =',t_old
            print*,'err_loc =',err_loc
            print*,'err_required =',delta*t_old*tol
            print*,'stepsize rejected, stepping down to:',t_step
         endif 
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            print*,"Failure in expogk: ---"
            print*,"The requested tolerance is too high."
            Print*,"Rerun with a smaller value."
            iflag = 2
            stop
            return
         endif
         goto 401
      endif
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1 

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step 
*
*---  display and keep some information ...
*
      if ( itrace.ne.0 ) then
         print*,'integration',istep,'---------------------------------'
         print*,'scale-square =',ns
         print*,'step_size =',t_step
         print*,'err_loc   =',err_loc
         print*,'next_step =',t_new
      endif
 
      step_min = MIN( step_min, t_step ) 
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )
 
      if ( mxstep.eq.0 .or. istep.lt.mxstep ) goto 100
      iflag = 1
 
 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = istep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = step_min
      wsp(2)  = step_max
      wsp(3)  = 0.0d0
      wsp(4)  = 0.0d0
      wsp(5)  = x_error
      wsp(6)  = s_error
      wsp(7)  = tbrkdwn
      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
*
      subroutine DGPADM( ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag )

*----------------------------------------------------------------------|

      implicit none
      integer ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
      double precision t, H(ldh,m), wsp(lwsp)
*
*-----Purpose----------------------------------------------------------|
*
*     Computes exp(t*H) = r(t*H) = (+/-)( I + 2*( q(t*H)\p(t*H) ) ) 
*     where r(x) is the irreducible rational Pade approximation to
*     exp(x). This routine computes the matrix exponential of a 
*     general matrix in full.
*
*-----Arguments--------------------------------------------------------|
*
*     ideg      : (input) the degre of the diagonal pade used.
*                 a value of 6 is generally satisfactory.
*
*     m         : (input) order of H.
*
*     H(ldh,m)  : (input) argument matrix.
*
*     t         : (input) time-scale (can be < 0).
*                  
*     wsp(lwsp): (workspace/output) lwsp .ge. 4*m*m+ideg+1.
*
*     ipiv(m)  : (workspace)
*
*>>>> iexph    : (output) integer such that wsp(iexph) points on exp(H)
*                i.e., exp(H) is located at wsp(iexph ... iexph+m*m-1)
*                      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*     ns       : (output) number of scaling-squaring used.
*
*     iflag    : (output) exit flag.
*                      0 - no problem
*                     <0 - problem
*
*-----Author-----------------------------------------------------------|
*
*     Roger B. Sidje (rbs@maths.uq.oz.au)     
*     Last update: April 3, 1996
*
*-----Reference--------------------------------------------------------|
*
*     R.B. Sidje, EXPOKIT. Software Package for Computing Matrix
*     Exponentials. Technical Report, Department of Mathematics 
*     University of Queensland - Brisbane QLD 4072, Australia.
*     April, 1996.
*
*----------------------------------------------------------------------|
*
      integer mm,i,j,k,i1,i2,ih2,ip,iq,iused,ifree,iodd,icoef,iput,iget
      double precision hnorm,scale,scale2,cp,cq

      intrinsic INT,AINT,ABS,DBLE,LOG,MAX

*---  check restrictions on input parameters ...
      mm = m*m
      iflag = 0
      if ( ldh.lt.m ) iflag = -1
      if ( lwsp.lt.4*mm+ideg+1 ) iflag = -2
      if ( iflag.ne.0 ) stop 'bad sizes (in input of DGPADM)'
*
*---  initialise pointers ...
*
      icoef = 1
      ih2 = icoef + (ideg+1)
      ip  = ih2 + mm
      iq  = ip + mm
      ifree = iq + mm
*
*---  scaling: seek ns such that ||t*H/2^ns|| < 1/2; 
*     and set scale = t/2^ns ...
*
      do i = 1,m
         wsp(i) = 0.0d0
      enddo
      do j = 1,m
         do i = 1,m
            wsp(i) = wsp(i) + ABS( H(i,j) )
         enddo
      enddo
      hnorm = 0.0d0
      do i = 1,m
         hnorm = MAX( hnorm,wsp(i) )
      enddo
      hnorm = ABS( t*hnorm )
      if ( hnorm.eq.0.0d0 ) stop 'Error - null H in input of DGPADM.'
      ns = MAX( 0,INT(LOG(hnorm)/LOG(2.0d0))+2 )
      scale = t / DBLE(2**ns)
      scale2 = scale*scale
*
*---  compute pade coefficients ...
*
      i1 = ideg+1
      i2 = 2*ideg+1
      wsp(icoef) = 1.0d0
      do k = 1,ideg
         wsp(icoef+k) = (wsp(icoef+k-1)*dble( i1-k ))/dble( k*(i2-k) )
      enddo
*
*---  H2 = scale2*H*H ...
*
      call DGEMM( 'n','n',m,m,m,scale2,H,ldh,H,ldh,0.0d0,wsp(ih2),m )
*
*---  initialize p and q
*
      cp = wsp(icoef+ideg-1)
      cq = wsp(icoef+ideg)
      do j = 1,m
         do i = 1,m
            wsp(ip + (j-1)*m + i-1) = 0.0d0
            wsp(iq + (j-1)*m + i-1) = 0.0d0
         enddo
         wsp(ip + (j-1)*(m+1)) = cp
         wsp(iq + (j-1)*(m+1)) = cq
      enddo
*
*---  Apply Horner rule ...
*
      iodd = 1
      k = ideg - 1
 100  continue
      iused = iodd*iq + (1-iodd)*ip
      call DGEMM( 'n','n',m,m,m, 1.0d0,wsp(iused),m,
     .             wsp(ih2),m, 0.0d0,wsp(ifree),m )
      do j = 1,m
         wsp(ifree+(j-1)*(m+1)) = wsp(ifree+(j-1)*(m+1))+wsp(icoef+k-1)
      enddo
      ip = (1-iodd)*ifree + iodd*ip
      iq = iodd*ifree + (1-iodd)*iq
      ifree = iused
      iodd = 1-iodd
      k = k-1
      if ( k.gt.0 )  goto 100
*
*---  Obtain (+/-)(I + 2*(p\q)) ...
*
      if ( iodd .eq. 1 ) then
         call DGEMM( 'n','n',m,m,m, scale,wsp(iq),m,
     .                H,ldh, 0.0d0,wsp(ifree),m )
         iq = ifree
      else
         call DGEMM( 'n','n',m,m,m, scale,wsp(ip),m,
     .                H,ldh, 0.0d0,wsp(ifree),m )
         ip = ifree
      endif
      call DAXPY( mm, -1.0d0,wsp(ip),1, wsp(iq),1 )
      call DGESV( m,m, wsp(iq),m, ipiv, wsp(ip),m, iflag )
      if ( iflag.ne.0 ) stop 'Problem in DGPADM/DGESV'
      call DSCAL( mm, 2.0d0, wsp(ip), 1 )
      do j = 1,m
         wsp(ip+(j-1)*(m+1)) = wsp(ip+(j-1)*(m+1)) + 1.0d0
      enddo
      iput = ip
      if ( ns.eq.0 .and. iodd.eq.1 ) then
         call DSCAL( mm, -1.0d0, wsp(ip), 1 )
         goto 200
      endif
*
*--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...
*
      iodd = 1
      do k = 1,ns
         iget = iodd*ip + (1-iodd)*iq
         iput = (1-iodd)*ip + iodd*iq
         call DGEMM( 'n','n',m,m,m, 1.0d0,wsp(iget),m, wsp(iget),m,
     .                0.0d0,wsp(iput),m )
         iodd = 1-iodd
      enddo
 200  continue
      iexph = iput
      END
*----------------------------------------------------------------------|
*
*
* Modified by gg for the c++ interface
*
*
      subroutine matvec (x, y, n, valcsr1, ipr1, jpr1) 
      integer ipr1(*), jpr1(*)
      double precision  valcsr1(*)
      integer n
c        integer nnz_b
c      parameter(nnz_b=max(nnz_b1,nnz_b2,nnz_b3,nnz_b4))
c      common /MAT1/ valcsr1,ipr1,jc1
c      integer ipr1(nm/2+1),jc1(nnz_b)
c      real*8  valcsr1(nnz_b)
      double precision  x(*), y(*) 
c-----------------------------------------------------------------------
c         A times a vector
c----------------------------------------------------------------------- 
c multiplies the matrix A by a vector using the dot product form
c Matrix A is stored in compressed sparse row storage.
c
c on entry:
c----------
c x     = real array of length equal to the column dimension of
c         the A matrix.
c
c on return:
c-----------
c y     = real array of length n, containing the product y=Ax
c
c-----------------------------------------------------------------------
c local variables
c
      double precision t
      integer i, k
c-----------------------------------------------------------------------
      do 100 i = 1,n
c
c     compute the inner product of row i with vector x
c 
         t = 0.0d0
         do 99 k=ipr1(i), ipr1(i+1)-1 
            t = t + valcsr1(k)*x(jpr1(k))
 99      continue
c
c     store result in y(i) 
c
         y(i) = t
 100  continue
c
      return
      end
c For the test      
c----------------------------------------------------------------------- 
      subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
c----------------------------------------------------------------------- 
      double precision a(*),ao(*),x
      integer ir(*),jc(*),jao(*),iao(*)
c-----------------------------------------------------------------------
c  Coordinate     to   Compressed Sparse Row 
c----------------------------------------------------------------------- 
c converts a matrix that is stored in coordinate format
c  a, ir, jc into a row general sparse ao, jao, iao format.
c
c on entry:
c--------- 
c nrow	= dimension of the matrix 
c nnz	= number of nonzero elements in matrix
c a,
c ir, 
c jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
c         nonzero elements of the matrix with a(k) = actual real value of
c 	  the elements, ir(k) = its row number and jc(k) = its column 
c	  number. The order of the elements is arbitrary. 
c
c on return:
c----------- 
c ir 	is destroyed
c
c ao, jao, iao = matrix in general sparse matrix format with ao 
c 	continung the real values, jao containing the column indices, 
c	and iao being the pointer to the beginning of the row, 
c	in arrays ao, jao.
c
c Notes:
c------ This routine is NOT in place.  See coicsr
c
c------------------------------------------------------------------------
      if(nnz.lt.1) return 
      do 1 k=1,nrow+1
         iao(k) = 0
 1    continue
c determine row-lengths.
      do 2 k=1, nnz
         iao(ir(k)) = iao(ir(k))+1
 2    continue
c starting position of each row..
      k = 1
      do 3 j=1,nrow+1
         k0 = iao(j)
         iao(j) = k
         k = k+k0
 3    continue
c go through the structure  once more. Fill in output matrix.
      do 4 k=1, nnz
         i = ir(k)
         j = jc(k)
         x = a(k)
         iad = iao(i)
         ao(iad) =  x
         jao(iad) = j
         iao(i) = iad+1
 4    continue
c shift back iao
      do 5 j=nrow,1,-1
         iao(j+1) = iao(j)
 5    continue
      return
      end
c------------- end of coocsr ------------------------------------------- 
c----------------------------------------------------------------------- 
