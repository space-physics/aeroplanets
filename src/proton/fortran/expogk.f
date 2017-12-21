*
*----------------------------------------------------------------------|-----
*
      subroutine expogk( n, m, t, v, w, tol, anorm,
     .                   wsp,lwsp, iwsp,liwsp, matvec, iverbose,iflag ) 

*----------------------------------------------------------------------|-----

      implicit double precision( a-h, o-z )
      integer n, m, lwsp, liwsp, iverbose, iflag, iwsp(liwsp)
      double precision t, tol, anorm, v(n), w(n), wsp(lwsp)
      external matvec
*
*-----Purpose----------------------------------------------------------|-----
*
*---  EXPO-GK computes w = exp(t*A)*v - for a General matrix A ...
*
*     exp(t*A) is *not* computed in isolation before being applied to v. 
*     On the contrary, the action of the matrix exponential operator on the 
*     operand vector is evaluated directly.
*
*     In this way, large sparse problems can be addressed efficiently.
*     The techniques used are based on Krylov subspace projection methods.
*     The matrix under consideration interacts (implicity) via the external 
*     routine `matvec' perfoming the matrix X vector product.
*
*-----Arguments--------------------------------------------------------|-----
*
*     n          : (input) the order of the principal matrix A.
*
*     m          : (input) the maximum size for the Krylov basis.
*
*     t          : (input) the time at wich the solution is needed.
*
*     v(n)       : (input) the given operand vector.
*
*     w(n)       : (output) the computed approximation of exp(t*A)*v.
*
*     tol        : (input) the requested tolerance on the approximation w.
*                  (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol.)
*
*     anorm      : (input) an approximation of some norm of A.
*
*     wsp(lwsp)  : (workspace) lwsp .ge. n*(m+1)+n + (m+2)^2 +4*(m+2)^2+ideg+1
*                                       +----------+---------+---------------+
*                  (actually, ideg=6)        V          H     wsp for EXPPADE
*                   
*     iwsp(liwsp): (workspace) liwsp .ge. m+2
*
*     matvec     : external subroutine for matrix-vector multiplication.
*                  synopsis: matvec( x, y )
*                            double precision x(*), y(*)
*                  computes: y(1:n) <- A*x(1:n)
*                            where A is the principal matrix.
*
*     iverbose   : (input) running mode. 0=silent, 1=verbose.
*
*     iflag      : (output) exit flag.
*                    <0 - bad input arguments 
*                     0 - no problem
*                     1 - maximum number of steps reached without convergence
*
*-----Accounts on the computation--------------------------------------|-----
*     Upon exit, an interested user may retrieve accounts on the computations.
*     They are located in the workspace arrays wsp and iwsp as indicated below: 
*
*     location  mnemonic  description
*     -----------------------------------------------------------------------
*     iwsp(1) = nmult     number of matrix X vector multiplications used
*     iwsp(2) = nexph     number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale    number of repeated squaring involved in exppade
*     iwsp(4) = istep     number of integration steps used up to completion 
*     iwsp(5) = nreject   number of rejected step-sizes
*     iwsp(6) = ibrkflag  set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn   if `happy breakdown', basis-size when it occured
*     -----------------------------------------------------------------------
*     wsp(1)  = step_min  minimum step-size used during integration
*     wsp(2)  = step_max  maximum step-size used during integration
*     wsp(3)  = dummy
*     wsp(4)  = dummy
*     wsp(5)  = x_error   maximum among all local truncation errors
*     wsp(6)  = s_error   global sum of local truncation errors
*     wsp(7)  = tbrkdwn   if `happy breakdown', time when it occured
*
*----------------------------------------------------------------------|-----
*-----The following parameters may be also adjusted herein-------------|-----
*
      parameter( mxstep   = 500, 
     .           mxreject = 0,
     .           ideg     = 6, 
     .           delta    = 1.2d0,
     .           gamma    = 0.9d0 )

*     mxstep   : maximum allowable number of integration steps.
*                The value 0 corresponds to an infinite number of steps.
* 
*     mxreject : maximum allowable number of rejections at each step. 
*                The value 0 corresponds to an infinite number of rejections.
*
*     ideg     : the Pade approximation of type (ideg,ideg) is used as an
*                approximation to exp(H). The value 0 switches to the
*                uniform rational Chebyshev approximation of type (14,14).
*
*     delta    : local truncation error `safety factor'
*
*     gamma    : stepsize `shrinking factor'
*
*-----Author-----------------------------------------------------------|-----
*
*     Roger B. Sidje (rbs@maths.uq.oz.au), Department of Mathematics 
*     University of Queensland - Brisbane QLD 4072, Australia.
*
*     Last update: March 28, 1995.
*
*-----References-------------------------------------------------------|-----
*
*>>   B. Philippe and R.B. Sidje, Transient Solutions of Markov processes by 
*     Krylov subspaces. In W. J. Stewart, editor, 2nd International Workshop
*     on the Numerical Solution of Markov Chains, Raleigh, NC, USA, 1995.
*
*>>   R.B. Sidje, Parallel algorithms for large sparse matrix exponentials.
*     Application to numerical transient analysis of Markov processes.
*     PhD thesis, University of Rennes 1, France, July 1994.
*
*>>   E. Gallopoulos and Y. Saad, Efficient solution of parabolic equations
*     by Krylov approximation methods. SIAM J. Sci. Stat. Comput., 
*     13(5):1236--1264, Sept 1992.
*
*>>   Y. Saad, Analysis of some Krylov subspace approximations to the matrix
*     exponential operator. SIAM J. Numer. Anal., 29(1):208--227, Feb 1992.
*
*----------------------------------------------------------------------|-----
*
*---  check restrictions on input parameters ...
      iflag = 0
      if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+2 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) stop 'bad sizes (in input of expogk)'
*
*---  initialisations ...
*
      k1 = 2
      mh = m + 2
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

      tbrkdwn  = 0.0d0
      step_min = t
      step_max = 0.0d0
      istep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

*>>>  break_tol = tol
      break_tol = anorm*tol

      a1 = 4.0d0/3.0d0
 1    a2 = a1 - 1.0d0
      a3 = a2 + a2 + a2
      eps = DABS( a3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      rndoff = eps*anorm
 
      call DCOPY( n, v,1, w,1 )
      beta = DNRM2( n, w,1 )
*
*---  obtain the very first stepsize ...
*
      xm = 1.0d0/DBLE( m )
      fact = tol*(((m+1)/2.72D0)**(m+1))*DSQRT(2.0D0*3.14D0*(m+1))
      t_new = (1.0d0/anorm)*(fact/(4.0d0*beta*anorm))**xm
      p = 10.0d0**(NINT( DLOG10( t_new )-DSQRT( 0.1d0 ) )-1)
      t_new = DINT( t_new/p + 0.55d0 ) * p
*
*---  step-by-step integration ...
*
 100  if ( t_now.ge.t ) goto 500

      istep = istep + 1
      t_step = DMIN1( t-t_now, t_new )
 
      beta = DNRM2( n, w,1 )
      do i = 1,n
         wsp(iv + i-1) = w(i)/beta
      enddo
      do i = 1,mh*mh
         wsp(ih+i-1) = 0.0d0
      enddo
*
*---  Arnoldi loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v) )
         do i = 1,j
            hij = DDOT( n, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            call DAXPY( n, -hij, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            wsp(ih+(j-1)*mh+i-1) = hij
         enddo
         hj1j = DNRM2( n, wsp(j1v),1 )
         wsp(ih+(j-1)*mh+j) = hj1j
         if ( hj1j.le.break_tol ) then
            if ( iverbose.ne.0 ) print*,'happy breakdown: ',
     .                                  'mbrkdwn =',j,' h =',hj1j
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
*---  set 1 for the 2-corrected scheme ...
*
      wsp(ih+m*mh+m+1) = 1.0d0
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 
 400  continue
*---  if `happy breakdown' go straightforward at the end ... 
      if ( k1.eq.0 ) then
         t_step = t-t_now
         tbrkdwn = t_now
       endif

 401  continue
*---  compute w = beta*V*exp(t_step*H)*e1
      nexph = nexph + 1
      mx = mbrkdwn + k1
      if ( ideg.ne.0 ) then
*---     rational Pade approximation ...
         call exppade( ideg, mx, t_step, wsp(ih),mh,
     .                 wsp(ifree),lfree, iwsp, iexph, ns, iflag )
         iexph = ifree + iexph - 1
         nscale = nscale + ns
      else
*---     uniform rational Chebyshev approximation ...
         iexph = ifree
         do i = 1,mx
            wsp(iexph+i-1) = 0.0d0
         enddo
         wsp(iexph) = 1.0d0
         call expcheb( mx, t_step,wsp(ih),mh, wsp(iexph) )
      endif
      mx = mbrkdwn + MAX0( 0,k1-1 )
      call DGEMV( 'n', n,mx,beta,wsp(iv),n,wsp(iexph),1,0.0d0,w,1 )
 
 402  continue
* 
*---  error estimate ...
* 
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         if ( ireject.eq.0 ) then
            nmult = nmult + 1
            call matvec( wsp(j1v-n), wsp(j1v) )
            avnorm = DNRM2( n, wsp(j1v),1 )
         endif
         phi1 = DABS( beta*wsp(iexph+m) )
         phi2 = DABS( beta*wsp(iexph+m+1) * t_step*avnorm )
         if ( phi1.gt.10.0d0*phi2 ) then
            err_loc = phi2
            xm = 1.0d0/DBLE( m )
         elseif ( phi1.gt.phi2 ) then
            err_loc = (phi1*phi2)/(phi1-phi2)
            xm = 1.0d0/DBLE( m )
         else
            err_loc = phi1
            xm = 1.0d0/DBLE( m-1 )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*   
      if ( (err_loc.gt.delta*t_step*tol) .and. 
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p = 10.0d0**(NINT( DLOG10( t_step )-DSQRT( 0.1d0 ) )-1)
         t_step = DINT( t_step/p + 0.55d0 ) * p
         if ( iverbose.ne.0 ) then
            print*,'t_step =',t_old
            print*,'err_loc =',err_loc
            print*,'err_required =',delta*t_old*tol
            print*,'stepsize rejected, stepping down to:',t_step
         endif 
         ireject = ireject + 1
         nreject = nreject + 1
         goto 401
      endif
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p = 10.0d0**(NINT( DLOG10( t_new )-DSQRT( 0.1d0 ) )-1)
      t_new = DINT( t_new/p + 0.55d0 ) * p 

      err_loc = DMAX1( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step 
*
*---  display and keep some information ...
*
      if ( iverbose.ne.0 ) then
         print*,'integration',istep,'---------------------------------'
         print*,'scale-square =',ns
         print*,'step_size =',t_step
         print*,'err_loc   =',err_loc
         print*,'next_step =',t_new
      endif
 
      step_min = DMIN1( step_min, t_step ) 
      step_max = DMAX1( step_max, t_step )
      s_error = s_error + err_loc
      x_error = DMAX1( x_error, err_loc )
 
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
*----------------------------------------------------------------------|------
*----------------------------------------------------------------------|-----
*
      subroutine exppade( ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag )
 
*----------------------------------------------------------------------|-----

      implicit double precision( a-h, o-z )
      dimension  H(ldh,m), wsp(lwsp), ipiv(m)
*
*-----Purpose----------------------------------------------------------|-----
*
*     Computes exp(t*H) = r(t*H) = (+/-)( I + 2*( q(t*H)\p(t*H) ) ) 
*     where r(x) is the irreducible rational Pade approximation to exp(x).
*     There is no assumption on H, it may be a general square matrix.
*     Therefore, this routine can be used in its own right to compute
*     the full matrix exponential.
*
*-----Author-----------------------------------------------------------|-----
*
*     Roger B. Sidje (rbs@maths.uq.oz.au),Department of Mathematics 
*     University of Queensland - Brisbane QLD 4072, Australia.
*
*     Last update: March 28, 1995.
*
*-----Arguments--------------------------------------------------------|-----
*
*     ideg      : (input) the degre of the diagonal pade used.
*                 a value of 6 is generally satisfactory.
*
*     m         : (input) order of H.
*
*     H(ldh,m)  : (input) argument matrix.
*
*     t         : (input) time-scale.
*                  
*     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.
*
*     ipiv(m)   : (workspace)
*
*>>>> iexph     : (output) an integer such that wsp(iexph) points on exp(H)
*                 i.e., exp(H) is located at wsp(iexph ... iexph+m*m-1).
*                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*     ns        : (output) number of scaling-squaring used.
*
*     iflag     : (output) exit flag.
*                       0 - no problem
*                      <0 - problem
*
*----------------------------------------------------------------------|-----

*---  check restrictions on input parameters ...
      mm = m*m
      iflag = 0
      if ( ldh.lt.m ) iflag = -1
      if ( lwsp.lt.4*mm+ideg+1 ) iflag = -2
      if ( iflag.ne.0 ) stop 'bad sizes (in input of exppade)'
*
*---  initialise pointers ...
*
      icoef = 1
      ih2 = icoef + (ideg+1)
      ip  = ih2 + mm
      iq  = ip + mm
      ifree = iq + mm
*
*---  scaling: seek ns such that ||t*H/2^ns|| < 1/2; and set scale = t/2^ns ...
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
         hnorm = DMAX1( hnorm,wsp(i) )
      enddo
      hnorm = t*hnorm
      if ( hnorm.eq.0.0d0 ) stop 'Error - null H in input of exppade.' 
* // MODIF Marina le 23/08  //
*      ns = DMAX1( 0,IDINT(DLOG(hnorm)/DLOG(2.0d0))+2 )
      ns = DMAX1( 0.0d0,IDINT(DLOG(hnorm)/DLOG(2.0d0))+2.0d0 )
* // --------------------- //

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
*----------------------------------------------------------------------|-----
*----------------------------------------------------------------------|-----
*
      subroutine expcheb( m, t, H,ldh, y )
 
*----------------------------------------------------------------------|-----
 
      implicit double precision( a-h, o-z )
      dimension  H(ldh,m), y(m)
 
*-----Purpose----------------------------------------------------------|-----
*
*---  Computes y = exp(t*H)*y using the partial fraction expansion of the
*     uniform rational Chebyshev approximation to exp(-x) of type (14,14).
*     The matrix H is assumed to be upper-Hessenberg.
*
*-----Author-----------------------------------------------------------|-----
*
*     Roger B. Sidje (rbs@maths.uq.oz.au), Department of Mathematics 
*     University of Queensland - Brisbane QLD 4072, Australia.
*
*     Last update: March 28, 1995.
*
*-----Arguments--------------------------------------------------------|-----
*
*     m       : (input) order of the Hessenberg matrix H
*
*     t       : (input) time-scaling factor.
*
*     H(ldh,m): (input) upper Hessenberg matrix.
*
*     y(m)    : (input/output) on input the operand vector,
*               on output the resulting vector exp(t*H)*y.
*
*----------------------------------------------------------------------|-----
*
      parameter ( mmax=70, ndeg=7 ) 
      complex*16 Hc(mmax,mmax), yc(mmax), tmpc
      complex*16 alpha(ndeg), theta(ndeg)
      double precision ysave(mmax)

*---  Check restrictions on input parameters ...
      if ( m.gt.mmax ) stop 'm is too big in expcheb, try exppade.'

*---  Coefficients and poles of the partial fraction expansion
*
*     exp(-x) = alpha0 + sum   Real [ alpha(i) / (x - theta(i) ]
*                     i=1,ndeg
*  
      alpha0  =  0.183216998528140087E-11
      alpha(1)=( 0.557503973136501826E+02,-0.204295038779771857E+03)
      alpha(2)=(-0.938666838877006739E+02, 0.912874896775456363E+02)
      alpha(3)=( 0.469965415550370835E+02,-0.116167609985818103E+02)
      alpha(4)=(-0.961424200626061065E+01,-0.264195613880262669E+01)
      alpha(5)=( 0.752722063978321642E+00, 0.670367365566377770E+00)
      alpha(6)=(-0.188781253158648576E-01,-0.343696176445802414E-01)
      alpha(7)=( 0.143086431411801849E-03, 0.287221133228814096E-03)

      theta(1)=(-0.562314417475317895E+01, 0.119406921611247440E+01)
      theta(2)=(-0.508934679728216110E+01, 0.358882439228376881E+01)
      theta(3)=(-0.399337136365302569E+01, 0.600483209099604664E+01)
      theta(4)=(-0.226978543095856366E+01, 0.846173881758693369E+01)
      theta(5)=( 0.208756929753827868E+00, 0.109912615662209418E+02)
      theta(6)=( 0.370327340957595652E+01, 0.136563731924991884E+02)
      theta(7)=( 0.889777151877331107E+01, 0.166309842834712071E+02)
*     
*---  Accumulation of the contribution of each pole ...
*
      do j = 1,m
         ysave(j) = y(j)
         y(j) = y(j)*alpha0
      enddo
      do ip = 1, ndeg     
*---     Solve each fraction using Gaussian elimination with pivoting ...
         do j = 1,m
            do i = 1,MIN0( j+1,m )
               Hc(i,j) = DCMPLX( -t*H(i,j) )
            enddo
            Hc(j,j) = Hc(j,j) - theta(ip) 
            yc(j) = DCMPLX( ysave(j) ) 
         enddo 
         do i = 1,m-1
            i1 = i + 1
*---        Get pivot and exchange rows ...
            if ( ABS( Hc(i,i) ).lt.ABS( Hc(i1,i) ) ) then
               call ZSWAP( m-i+1, Hc(i,i),mmax, Hc(i1,i),mmax )
               call ZSWAP( 1, yc(i),1, yc(i1),1 )
            endif
*---        Forward eliminiation ... 
            tmpc = Hc(i1,i) / Hc(i,i)
            call ZAXPY( m-i1+1, -tmpc, Hc(i,i1),mmax, Hc(i1,i1),mmax )
            yc(i1) = yc(i1) - tmpc*yc(i)
         enddo
*---     Backward substitution ...    
         do i = m,1,-1 
            tmpc = yc(i)
            do j = i+1,m
               tmpc = tmpc - Hc(i,j)*yc(j)
            enddo
            yc(i) = tmpc / Hc(i,i)
         enddo
*---     Accumulate the partial result in y ...     
         do i = 1,m
            y(i) = y(i) + DBLE( alpha(ip)*yc(i) ) 
         enddo
      enddo
      END
*----------------------------------------------------------------------|-----


*----------------------------------------------------------------------|-----
*     External LAPACK routines used by EXPO-MK/GK.
*     This file is supplied in case where LAPACK is not installed in your
*     environment.
*----------------------------------------------------------------------|-----
      subroutine DGESV( N, M, A,LDA, IPIV, B,LDB, IFLAG )
      integer N, M, LDA, LDB, IPIV(M), IFLAG
      double precision A(LDA,N), B(LDB,N)
      call LUFACT( A,LDA, N, ipiv, iflag )
      if ( IFLAG.ne.0 ) stop 'Error in LU factorisation'
      do j = 1,M
         call LUSOLV( A,LDA, N, IPIV,B(1,j), 0 )
      enddo
      end
*----------------------------------------------------------------------|-----
      subroutine LUFACT(a,lda,n,ipvt,info)
      integer lda,n,ipvt(n),info
      double precision a(lda,n)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
*----------------------------------------------------------------------|-----
      subroutine LUSOLV(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(n),job
      double precision a(lda,n),b(n)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
