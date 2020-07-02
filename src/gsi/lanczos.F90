module lanczos
!$$$ module documentation block
!           .      .    .                                       .
! module:   lanczos
!   prgmmr: tremolet
!
! abstract: Contains variables and routines for preconditioned
!           lanczos minimizer following mike fisher's algorithm.
!
! program history log:
!   2007-05-16  tremolet
!   2007-07-11  tremolet - increment sensitivity to obs
!   2007-11-23  todling  - add timers
!   2009-01-18  todling  - minimal changes to interface w/ quad-based evaljgrad
!                          Note: no attempt made to reproduce across pe's yet
!   2009-08-18  lueken   - update documentation
!   2010-03-10  treadon  - add essl interface
!   2010-03-17  todling  - add analysis error estimate (congrad_siga)
!   2010-08-19  lueken   - add only to module use
!   2010-09-06  todling  - add ltcost parameter to allow writing out true cost
!   2010-09-20  el akkraoui - properly repositioned call related to ltcost  
!   2011-04-07  todling  - rename precond to lanczos_precond
!   2011-04-19  el akkraoui - avoid convert precond vectors for the next outer loop;
!                             and avoid lanczos decomposition at each iteration
!   2011-07-04  todling  - determine precision based on kinds
!
! Subroutines Included:
!   congrad       - Main minimization routine
!   setup_precond - Prepare the preconditioner
!   save_precond  - Save eigenvectors for constructing the next preconditioner
!   lanczos_precond - Preconditioner itself (called from congrad, internal)
!
! Variable Definitions:
!   lmpcgl  : .t. ====> precondition conjugate-gradient minimization
!   r_max_cnum_pc : Maximum allowed condition number for the preconditioner
!   npcvecs : number of vectors which make up the preconditioner
!
!   yvcglpc: eigenvectors (from an earlier minimization)
!            that are used to construct the preconditioner.
!   rcglpc : eigenvalues (from an earlier minimization)
!            that are used to construct the preconditioner.
!   nvcglpc: the number of eigenpairs used to form the preconditioner.
!
!   yvcglev: eigenvectors for the current minimization.
!   rcglev : eigenvalues for the current minimization.
!   nvcglev: the number of eigenpairs for the current minimization.
!   ltcost : .t. to calculate true cost function (unscalled), this adds
!           considerable computation cost; only used in test mode
!
!   yvcglwk: work array of eigenvectors
!
! attributes:
!   language: f90
!   machine:
!
! ------------------------------------------------------------------------------
use kinds, only: r_kind,i_kind,r_quad,r_single,r_double
use constants, only: zero, one, two, one_tenth
use jfunc, only: iter
use control_vectors, only: control_vector,allocate_cv,inquire_cv,deallocate_cv, &
    write_cv,read_cv,dot_product,assignment(=)
use file_utility, only : get_lun
use timermod, only: timer_ini, timer_fnl
! ------------------------------------------------------------------------------

implicit none
save
private
public congrad, setup_congrad, save_precond, congrad_ad, read_lanczos, &
       congrad_siga

! ------------------------------------------------------------------------------

logical :: ltcost_= .false.
logical :: lmpcgl = .false.
logical :: lconvert = .false. !if true, convert the preconditioner vectors for the next outer loop.
logical :: ldecomp  = .false. !if true, carry lanczos decomposition at each iteration
real(r_kind) :: r_max_cnum_pc = 10.0_r_kind
real(r_kind) :: xmin_ritz = one
real(r_kind) :: pkappa = one_tenth

integer(i_kind) :: npcvecs, nvcglpc, nvcglev, nwrvecs
real(r_kind), allocatable :: rcglpc(:)
real(r_kind), allocatable :: rcglev(:)

integer(i_kind) :: mype,nprt,jiter,maxiter
logical :: l4dvar,lanczosave
real(r_kind), allocatable :: zlancs(:,:)

type(control_vector), allocatable, dimension(:) :: yvcglpc
type(control_vector), allocatable, dimension(:) :: yvcglev
type(control_vector), allocatable, dimension(:) :: yvcglwk
type(control_vector), allocatable, dimension(:) :: cglwork

! --------------------------------------
integer(i_kind), parameter :: n_default_real_kind = r_single
integer(i_kind), parameter :: n_double_kind       = r_double
! --------------------------------------

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! congrad
! ------------------------------------------------------------------------------
subroutine setup_congrad(kpe,kprt,kiter,kiterstart,kmaxit,kwrvecs, &
                         ld4dvar,ldsave,ltcost)
!$$$  subprogram documentation block
!                .      .    .                                         .
! subprogram:    setup_congrad
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-08-05  lueken - added subprogram doc block
!
!   input argument list:
!    kpe,kprt,kiter,kiterstart,kmaxit,kwrvecs
!    ld4dvar,ldsave,ltcost
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

implicit none

integer(i_kind), intent(in   ) :: kpe,kprt,kiter,kiterstart,kmaxit,kwrvecs
logical        , intent(in   ) :: ld4dvar,ldsave,ltcost

integer(i_kind) :: ii

mype=kpe
nprt=kprt
jiter=kiter
maxiter=kmaxit
nwrvecs=kwrvecs
l4dvar=ld4dvar
lanczosave=ldsave
ltcost_=ltcost

if (allocated(zlancs)) deallocate(zlancs)
allocate(zlancs(maxiter+1,4))
zlancs=zero

allocate(cglwork(maxiter+1))
do ii=1,kmaxit+1
   call allocate_cv(cglwork(ii))
   cglwork(ii)=zero
enddo

if (jiter==kiterstart) then
   npcvecs=0
   nvcglpc=0
   nvcglev=0
endif

if (jiter>1) call setup_precond()

if (mype==0) write(6,*)'setup_congrad end'
call inquire_cv

end subroutine setup_congrad
! ------------------------------------------------------------------------------
subroutine congrad(xhat,pcost,gradx,preduc,kmaxit,iobsconv,lsavevecs)
!$$$  subprogram documentation block
!                .      .    .                                         .
! subprogram:    congrad
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-08-05  lueken - added subprogram doc block
!
!   input argument list:
!    xhat
!    gradx
!    preduc
!    kmaxit
!    iobsconv
!
!   output argument list:
!    xhat
!    pcost
!    gradx
!    preduc
!    kmaxit
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

implicit none

type(control_vector), intent(inout) :: xhat
real(r_kind)        , intent(  out) :: pcost
type(control_vector), intent(inout) :: gradx
real(r_kind)        , intent(inout) :: preduc
integer(i_kind)     , intent(inout) :: kmaxit
integer(i_kind)     , intent(in   ) :: iobsconv

logical             , intent(in   ) :: lsavevecs

character(len=*), parameter :: myname='congrad'
type(control_vector)        :: grad0,zww
type(control_vector)        :: gradf
type(control_vector)        :: xiter,xsens
real(r_quad)                :: pcostq
real(r_kind)                :: zbeta(2:kmaxit+1),zdelta(kmaxit),zv(kmaxit+1,kmaxit+1),&
   zbnds(kmaxit),zritz(kmaxit+1),zsave(kmaxit+1,4),&
   zqg0(kmaxit+1),zsstwrk(2*kmaxit)
real(r_kind)                :: zdla, zeta, preduc_norm
real(r_kind)                :: zbndlm, zgnorm, znorm2l1, zreqrd, ztheta1
integer(i_kind)             :: ingood,itheta1,jm,imaxevecs,ii,jj,jk,isize
integer(i_kind)             :: kminit, kmaxevecs,iunit,iprt
logical                     :: lsavinc, lldone
character(len=17)           :: clfile

! --------------------------------------

!--- initialize timer
call timer_ini('congrad')

iprt=nprt
if(ltcost_) iprt=0
kminit = kmaxit
kmaxevecs = kmaxit
imaxevecs = 0
lldone=.false.
if (kmaxit>maxiter) then
   write(6,*)'setup_congrad: kmaxit>maxiter',kmaxit,maxiter
   call stop2(138)
end if

if (mype==0) write(6,*) '---- Lanczos Solver ----'

!--- allocate distributed storage

call allocate_cv(grad0)
call allocate_cv(zww)

!--- 'zeta' is an upper bound on the relative error of the gradient.

zeta  = 1.0e-4_r_kind

zreqrd = preduc

!--- change of variable to account for preconditioning

if (lmpcgl) call lanczos_precond(xhat,+2)

zgnorm = sqrt( dot_product (gradx,gradx))

if (mype==0) write (6,*)'grepmin Starting point: Estimated gradient norm=',zgnorm

if (lmpcgl) call lanczos_precond(gradx,-2)

cglwork(1) = gradx
znorm2l1 = dot_product(cglwork(1),cglwork(1))
cglwork(1)%values = cglwork(1)%values / sqrt(znorm2l1)

!--- save initial control vector and gradient

grad0 = gradx

zqg0(1) = dot_product(cglwork(1),grad0)

if(nprt>=1.and.ltcost_) then
   if (mype==0) then
      write(6,*)' True cost function (at considerable additional computational cost): '
      write(6,*)' ----------------------------------------------------------------- '
   endif
endif

!--- Lanczos iteration starts here

ingood = 0
iter   = 1
lanczos_loop : do

!--- evaluate the hessian applied to the latest lanczos vector

   do jj=1,zww%lencv
      zww%values(jj) = xhat%values(jj) + cglwork(iter)%values(jj)
   enddo

   if (lmpcgl) call lanczos_precond(zww,-2)

   lsavinc=.false.
   call evaljgrad(zww,pcostq,gradx,lsavinc,iprt,myname)
   pcost=pcostq

   if (lmpcgl) call lanczos_precond(gradx,-2)

   do jj=1,gradx%lencv
      gradx%values(jj) = gradx%values(jj) - grad0%values(jj)
   enddo

!--- calculate zdelta

   zdelta(iter) = dot_product(cglwork(iter),gradx)

   if (zdelta(iter)<=zero) then
      if (mype==0) write(6,*)'congrad stopping: J" not positive definite',zdelta(iter)
      iter = iter-1
      exit lanczos_loop
   endif

!--- Calculate the new lanczos vector (this is the lanczos recurrence)

   do jj=1,gradx%lencv
      gradx%values(jj) = gradx%values(jj) - zdelta(iter) * cglwork(iter)%values(jj)
   enddo
   if (iter>1) then
      do jj=1,gradx%lencv
         gradx%values(jj) = gradx%values(jj) - zbeta(iter) * cglwork(iter-1)%values(jj)
      enddo
   endif

!--- orthonormalize gradient against previous gradients

   do jm=iter,1,-1
      zdla = dot_product(gradx,cglwork(jm))
      do jj=1,gradx%lencv
         gradx%values(jj) = gradx%values(jj) - zdla*cglwork(jm)%values(jj)
      enddo
   enddo

   zbeta(iter+1) = sqrt(dot_product(gradx,gradx))

   do jj=1,gradx%lencv
      cglwork(iter+1)%values(jj) = gradx%values(jj) / zbeta(iter+1)
   enddo

   zqg0(iter+1) = dot_product(cglwork(iter+1),grad0)

!--- calculate the reduction in the gradient norm and cost

   zlancs(1:iter,1) =  zdelta(1:iter)
   zlancs(2:iter,2) =  zbeta (2:iter)
   zlancs(1:iter,3) = -zqg0  (1:iter)
 
   call ptsv

   do jj=1,zww%lencv
      zww%values(jj) = grad0%values(jj) &
         + (zbeta(iter+1)*zlancs(iter,3))*cglwork(iter+1)%values(jj)
   enddo

   do jj=1,iter
      do ii=1,zww%lencv
         zww%values(ii)  = zww%values(ii)  - cglwork(jj)%values(ii)*zqg0(jj)
      enddo
   enddo

   if (lmpcgl) call lanczos_precond(zww,+2)

   preduc_norm = sqrt(dot_product(zww,zww))
   preduc = preduc_norm/zgnorm
   if (mype==0) write (6,'(2(1X,A,ES25.18))') &
      'Estimated gradient norm=',preduc_norm,' reduction = ',preduc


!--- determine eigenvalues and eigenvectors of the tri-diagonal problem
   if((ldecomp .or. (iter==kmaxit)) .and. lsavevecs) then 
      zlancs(1:iter  ,4) = zdelta(1:iter)
      zlancs(1:iter-1,1) = zbeta (2:iter)

      if (iter /= 1) then
         call steqr
      else
         zv(1,1) = one
      endif

      zritz(1:iter) = zlancs(1:iter,4)

      if (mype==0) write(6,*)'congrad: ritz values are: ',zritz(1:iter)

!--- estimate error bounds

      zbndlm = zeta*zritz(iter)
 
      zbnds(1:iter) = abs(zbeta(iter+1)*zv(iter,1:iter))
      if (mype==0) write (6,*)'congrad: error bounds are: ',zbnds(1:iter)

!--- Check for exploding or negative ritz values

      if (any(zritz(1:iter)<zero)) then
         if (mype==0) write(6,*)'congrad stopping: negative ritz value'
         iter = iter-1
         zlancs(1:iter  ,4) = zdelta(1:iter)
         zlancs(1:iter-1,1) = zbeta(2:iter)
 
         if (iter > 1) then
            call steqr
         else
            zv(1,1) = one
         endif

         zritz(1:iter) = zlancs(1:iter,4)
 
         zbnds(1:iter) = abs(zbeta(iter+1)*zv(iter,1:iter))
         exit lanczos_loop
      endif

      if (ingood>0) then
         if (zritz(itheta1)>1.01_r_kind*ztheta1) then
            if (mype==0) write(6,*)'congrad stopping: ritz values exploding'
            if (mype==0) write(6,*)'leading ritz value=',zritz(itheta1)
            if (mype==0) write(6,*)'leading converged eigenvalue=',ztheta1
         endif
      endif

!--- Count the converged eigenvectors

      ingood = 0
      do jm=1,iter
         if (zbnds(jm)<=zbndlm) then
            ingood = ingood+1
            if (mype==0) write(6,*)'congrad: converged eigenvalue ',zritz(jm)
         endif
      enddo

!--- save leading converged eigenvalue for explosion test

      if (ingood > 0) then
         do jm=iter,1,-1
            if (zbnds(jm) <= zbndlm) then
               ztheta1 = zritz(jm)
               itheta1 = jm
               exit
            endif
         enddo
      endif

      if (mype==0) write(6,*)'congrad: End of iteration: ',iter
      if (mype==0) write(6,'(/)')
 
!     count how many eigenpairs have converged to pkappa precision and have
!     eigenvalue > xmin_ritz (which is 1 by default)
!     (For the analysis, all eigenvalues should be >1. For the singular vector calculation,
!     we are not interested in decaying modes.)
!     However, when svs are computed with projection operators, 1 may not
!     be an appropriate choice for xmin_ritz

      imaxevecs = count(zbnds(1:iter)/zritz(1:iter)<=pkappa .and. zritz(1:iter)>xmin_ritz)

      if (imaxevecs >= kmaxevecs) then
         if (mype==0) write(6,*)imaxevecs,' eigenpairs converged to precision ',pkappa
         if (mype==0) write(6,'(/)')
         exit lanczos_loop
      endif

   end if

!  Tests for end of iterations
   if (iter >= kmaxit .or. (preduc <= zreqrd .and. iter >= kminit)) &
      exit lanczos_loop

 

   if (nprt>=1.and.ltcost_.and.iobsconv==0) then

!     Compute actual increment
      zsave=zlancs
      zlancs(1:iter,1) =  zdelta(1:iter)
      zlancs(2:iter,2) =  zbeta (2:iter)
      zlancs(1:iter,3) = -zqg0  (1:iter)
      call ptsv
 
      call allocate_cv(xiter)
      xiter=zero
      do jj=1,iter
         do ii=1,xiter%lencv
            xiter%values(ii) = xiter%values(ii)  + cglwork(jj)%values(ii)*zlancs(jj,3)
         enddo
      enddo

      call allocate_cv(gradf)
      gradf=zero
      if (lmpcgl) then
         call lanczos_precond(xiter,-2)
      endif
      call evaljgrad(xiter,pcostq,gradf,lsavinc,nprt,myname)

      call deallocate_cv(gradf)
      call deallocate_cv(xiter)
      zlancs=zsave

   endif

!  Test convergence in observation space
   if (iobsconv>=2) then

!     Compute actual increment
      zsave=zlancs
      zlancs(1:iter,1) =  zdelta(1:iter)
      zlancs(2:iter,2) =  zbeta (2:iter)
      zlancs(1:iter,3) = -zqg0  (1:iter)
      call ptsv
 
      call allocate_cv(xiter)
      call allocate_cv(xsens)
      xiter=zero
      do jj=1,iter
         do ii=1,xiter%lencv
            xiter%values(ii) = xiter%values(ii)  + cglwork(jj)%values(ii)*zlancs(jj,3)
         enddo
      enddo
      if (lmpcgl) call lanczos_precond(xiter,-2)
      xsens=xiter

!     Compute observation impact
      call congrad_ad(xsens,iter)
      call test_obsens(xiter,xsens)

!     Clean-up
      call deallocate_cv(xiter)
      call deallocate_cv(xsens)
      zlancs=zsave
   endif

!--- Increment the iteration counter

   iter = iter+1
   if (ingood>0) itheta1 = itheta1+1

enddo lanczos_loop

!--- end of Lanczos iteration

lldone=.true.

if (mype==0) then
   write(6,*)'Summary of Lanczos iteration:'
   write(6,*)'   Number of iterations performed: ',iter
   write(6,*)'   Maximum allowed number of iterations: ',kmaxit
   write(6,*)'   Minimum allowed number of iterations: ',kminit
   write(6,*)'   Required reduction in norm of gradient: ',zreqrd
   write(6,*)'   Achieved reduction in norm of gradient: ',preduc
   if (preduc > zreqrd) then
      write(6,*)'   *** Failed to meet convergence criterion ***'
   endif
   write(6,*)'   Number of sufficiently-converged eigenpairs: ',imaxevecs
endif

!--- Calculate the solution vector and gradient

zlancs(1:iter,1) =  zdelta(1:iter)
zlancs(2:iter,2) =  zbeta (2:iter)
zlancs(1:iter,3) = -zqg0  (1:iter)

call ptsv

do jj=1,gradx%lencv
   gradx%values(jj) = grad0%values(jj) &
      + zbeta(iter+1)*cglwork(iter+1)%values(jj)*zlancs(iter,3)
enddo

do jj=1,iter
   do ii=1,xhat%lencv
      xhat%values(ii)  = xhat%values(ii)  + cglwork(jj)%values(ii)*zlancs(jj,3)
      gradx%values(ii) = gradx%values(ii) - cglwork(jj)%values(ii)*zqg0(jj)
   enddo
enddo

!--- transform control variable and gradient back to unpreconditioned space

if (lmpcgl) then
   call lanczos_precond(xhat ,-2)
   call lanczos_precond(gradx,+2)
endif

!--- Compute observation impact
if (iobsconv>=1) then
   call allocate_cv(xsens)
   xsens=xhat

   call congrad_ad(xsens,iter)
   call test_obsens(xhat,xsens)

   call deallocate_cv(xsens)
endif

!--- Save lanczos vectors (if required for adjoint)

if (lanczosave) then
   do jj=1,iter
      clfile='lanczvec.XXX.YYYY'
      write(clfile(10:12),'(I3.3)') jiter
      write(clfile(14:17),'(I4.4)') jj
      call write_cv(cglwork(jj),clfile)
   enddo

   if (mype==0) then
      iunit=get_lun()
      clfile='zlanczos.XXX'
      write(clfile(10:12),'(I3.3)') jiter
      write(6,*)'Writing Lanczos coef. to file ',clfile
 
      open(iunit,file=trim(clfile),form='unformatted')
      write(iunit)maxiter
      write(iunit)zlancs(1:maxiter+1,1:4)
      close(iunit)
   endif
endif

!--- Calculate sufficiently converged eigenvectors of the preconditioned Hessian

if (l4dvar .and. lsavevecs) then
   zbnds(1:iter) = zbnds(1:iter)/zritz(1:iter)

   nvcglev = 0
   do jk=iter,1,-1
      if (zbnds(jk) <= pkappa .AND. zritz(jk) > xmin_ritz) then
         nvcglev=nvcglev+1
      endif
   enddo
   if (mype==0) write(6,*) &
      'Number of eigenpairs converged to requested accuracy NVCGLEV=',nvcglev

   nvcglev= min(nvcglev,nwrvecs)
   if (mype==0) write(6,*) &
      'Number of eigenevctors to be calculated is NVCGLEV=',nvcglev

   allocate(rcglev(nvcglev))
   allocate(yvcglev(nvcglev))
   do ii=1,nvcglev
      call allocate_cv(yvcglev(ii))
   enddo

   ii=0
   do jk=iter,1,-1
      if (zbnds(jk) <= pkappa .and. zritz(jk) > xmin_ritz .and. ii<nwrvecs) then
         ii = ii+1
         rcglev(ii) = zritz(jk)
         yvcglev(ii) = zero
         isize=size(yvcglev(ii)%values)
         do jm=1,iter
            do jj=1,isize
               yvcglev(ii)%values(jj)=yvcglev(ii)%values(jj) + cglwork(jm)%values(jj)*zv(jm,jk)
            enddo
         enddo

         do jm=1,ii-1
            zdla=dot_product (yvcglev(jm),yvcglev(ii))
            do jj=1,isize
               yvcglev(ii)%values(jj) = yvcglev(ii)%values(jj) - zdla*yvcglev(jm)%values(jj)
            enddo
         enddo

         zdla=dot_product (yvcglev(ii),yvcglev(ii))
         yvcglev(ii)%values = yvcglev(ii)%values / sqrt(zdla)
      endif
   enddo

   if (mype==0.and.nvcglev>0) then
      write(6,'(/)')
      write(6,*)'Calculated eigenvectors for the following eigenvalues:'
      write(6,*)'RCGLEV=',rcglev(1:nvcglev)
      write(6,'(/)')
   endif

endif

!--- release memory, etc.

call deallocate_cv(grad0)
call deallocate_cv(zww)

!--- return the number of iterations actually performed

kmaxit=iter

!--- finalize timer
call timer_fnl('congrad')

return

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------
!   steqr - Simplified interface to lapack routines ssteqr/dsteqr
!           and to essl routines sgeev/dgeev
!-----------------------------------------------------------------------
  subroutine steqr
!$$$  subprogram documentation block
!                .      .    .                                         .
! subprogram:    steqr
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-08-05  lueken - added subprogram doc block
!   2010-03-10  treadon - add essl interface
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block
    implicit none

    integer(i_kind) :: info

#ifdef ibm_sp
    logical,allocatable:: select(:)
    integer(i_kind):: n,i,j,jj,ldz
    integer(i_kind):: iopt,lda,naux
    integer(i_kind),allocatable:: indx(:)
    real(r_kind),allocatable:: a(:,:),aux(:),w_order(:)
    complex(r_kind),allocatable:: w(:),z(:,:)

!   Use essl
    iopt=1
    n=iter
    lda=n
    ldz=kmaxit+1
    naux=2*n
    allocate(select(n),indx(n),a(lda,n),w(n),z(ldz,n),aux(naux),w_order(n))

!   Load matrix a.
    a=zero
    do i=1,n-1
       a(i,  i)=zlancs(i,4)  ! load diagonal
       a(i+1,i)=zlancs(i,1)  ! load sub-diagnonal
       a(i,i+1)=zlancs(i,1)  ! load super-diagonal
    end do
    a(n,n)=zlancs(n,4)       ! load diagonal

!   Additional initializations
    select=.false.    ! select not used for iopt=1
    w=zero
    z=zero
    aux=zero

!   Call essl routines
    if (r_kind == n_default_real_kind) then
       call sgeev(iopt, a, lda, w, z, ldz, select, n, aux, naux)
       do i=1,n
          w_order(i)=real(w(i),r_kind)
       end do
       call ssortx(w_order,1,n,indx)  ! sort eigenvalues into ascending order
    elseif (r_kind == n_double_kind) then
       call dgeev(iopt, a, lda, w, z, ldz, select, n, aux, naux)
       do i=1,n
          w_order(i)=real(w(i),r_kind)
       end do
       call dsortx(w_order,1,n,indx)  ! sort eigenvalues into ascending order
    else
       write(6,*)'STEQR: r_kind is neither default real nor double precision'
       call stop2(319)
    endif

!   Load essl eigenvalues and eigenvectors into output arrays
    do j=1,n
       zlancs(j,4)=w_order(j)          ! eigenvalues
       jj=indx(j)
       do i=1,ldz
          zv(i,j)=real(z(i,jj),r_kind) ! eigenvectors
       end do
    end do

!   Deallocate work arrays    
    deallocate(select,indx,a,w,z,aux,w_order)
    
#else

!   Use lapack
    if (r_kind == n_default_real_kind) then
       call ssteqr ('I',iter,zlancs(1,4),zlancs,zv,kmaxit+1,zsstwrk,info)
    elseif (r_kind == n_double_kind) then
       call dsteqr ('I',iter,zlancs(1,4),zlancs,zv,kmaxit+1,zsstwrk,info)       
    else
       write(6,*)'STEQR: r_kind is neither default real nor double precision'
       call stop2(319)
    endif
    if (info /= 0) then
       write (6,*)'Error in congrad: SSTEQR/DSTEQR returned info=',info
       write(6,*) 'STEQR: SSTEQR/DSTEQR returned non-zero info'
       call stop2(320)
    endif
#endif


end subroutine steqr

!-----------------------------------------------------------------------
!   ptsv - Simplified interface to lapack routines sptsv/dptsv
!          and to essl routines sptf,s/dptf,s
!-----------------------------------------------------------------------
subroutine ptsv
!$$$  subprogram documentation block
!                .      .    .                                         .
! subprogram:    ptsv
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-08-05  lueken - added subprogram doc block
!   2010-03-10  treadon- add essl interface
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block
  implicit none
  
  integer(i_kind) :: info

#ifdef ibm_sp
  integer(i_kind) :: i,n,iopt
  real(r_kind),allocatable:: c(:),d(:),bx(:)

! Use essl
  iopt=0
  n=iter
  allocate(c(n),d(n),bx(n))

! Load matrices
  c=zero
  do i=1,n-1
     c(i+1) = zlancs(i+1,2) ! lower subdiagonal of a
     d(i) = zlancs(i,1)     ! main diagonal of a
     bx(i)=zlancs(i,3)      ! right hand side b
  end do
  d(n) =zlancs(n,1)
  bx(n)=zlancs(n,3)          

! Factorize and solve system of equations using essl routines
  if (r_kind == n_default_real_kind) then
     call sptf(n,c,d,iopt)   ! factorize
     call spts(n,c,d,bx)     ! solve
  elseif (r_kind == n_double_kind) then
     call dptf(n,c,d,iopt)   ! factorize
     call dpts(n,c,d,bx)     ! solve
  else
     write(6,*) 'r_kind is neither default real nor double precision'
     call stop2(321)
  endif
  
! Load essl result into output arrays
  do i=1,n
     zlancs(i,3)=bx(i)
  end do
  
! Deallocate work arrays
  deallocate(c,d,bx)

#else

! Use lapack
  if (r_kind == n_default_real_kind) then
     call sptsv (iter,1,zlancs(1,1),zlancs(2,2),zlancs(1,3),kmaxit+1,info)
  elseif (r_kind == n_double_kind) then
     call dptsv (iter,1,zlancs(1,1),zlancs(2,2),zlancs(1,3),kmaxit+1,info)
  else
     write(6,*) 'r_kind is neither default real nor double precision'
     call stop2(321)
  endif
  if (info /= 0) then
     write (6,*) 'Error in congrad: SPTSV/DPTSV returned ',info
     write(6,*)'CONGRAD: SPTSV/DPTSV returned non-zero info'
     call stop2(322)
  endif

#endif

end subroutine ptsv

! ------------------------------------------------------------------------------
end subroutine congrad
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
subroutine congrad_ad(xsens,kiter)
!$$$  subprogram documentation block
!                .      .    .                                         .
! subprogram:    congrad_ad
!   prgmmr:
!
! abstract: Apply product of adjoint of estimated hessian to a vector.
!
! program history log:
!   2009-08-05  lueken - added subprogram doc block
!
!   input argument list:
!    xsens
!    kiter
!
!   output argument list:
!    xsens
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

implicit none

type(control_vector), intent(inout) :: xsens
integer(i_kind)     , intent(in   ) :: kiter

real(r_kind) :: zaa(kiter),zzz
integer(i_kind) :: ii,jj

!--- initialize timer
call timer_ini('congrad_ad')

zzz=dot_product(xsens,xsens)
if (mype==0) write(6,888)'congrad_ad: Norm  input=',sqrt(zzz)

if (lmpcgl) call lanczos_precond(xsens,-2)

zaa=zero
do jj=1,kiter
   zaa(jj)=dot_product(xsens,cglwork(jj))
enddo
do jj=2,kiter
   zaa(jj)=zaa(jj)-zlancs(jj,2)*zaa(jj-1)
enddo
zaa(kiter)=zaa(kiter)/zlancs(kiter,1)
do jj=kiter-1,1,-1
   zaa(jj)=zaa(jj)/zlancs(jj,1) - zaa(jj+1)*zlancs(jj+1,2)
enddo
xsens=zero
do jj=1,kiter
   do ii=1,xsens%lencv
      xsens%values(ii) = xsens%values(ii) + zaa(jj) * cglwork(jj)%values(ii)
   enddo
enddo

if (lmpcgl) call lanczos_precond(xsens,-2)

zzz=dot_product(xsens,xsens)
if (mype==0) write(6,888)'congrad_ad: Norm output=',sqrt(zzz)
888 format(A,3(1X,ES25.18))

!--- finalize timer
call timer_fnl('congrad_ad')

return
end subroutine congrad_ad
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
subroutine congrad_siga(siga,ivecs,rc)
!$$$  subprogram documentation block
!                .      .    .                                         .
! subprogram:    congrad_siga
!   prgmmr: todling
!
! abstract: Calculate estimate of analysis error
!
! program history log:
!  2010-03-17  todling  - initia code
!  2010-05-16  todling  - update to use gsi_bundle
!  2013-01-26  parrish  - wcoss debug compile flagged type mismatch error for
!                          "call bkg_stddev(aux,mval(ii))".
!                         I changed to 
!                          "call bkg_stddev(aux%step(ii),mval(ii))".
!                         Don't know if this is the correct modification.
!
!   input argument list:
!    siga
!
!   output argument list:
!    siga
!    ivecs
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block
use gsi_4dvar, only : nsubwin
use bias_predictors, only: predictors,allocate_preds,deallocate_preds
use state_vectors, only: allocate_state,deallocate_state
use gsi_bundlemod, only: gsi_bundlehadamard
use gsi_bundlemod, only: gsi_bundle
use gsi_bundlemod, only: assignment(=)
implicit none
type(gsi_bundle),intent(inout) :: siga  ! analysis errors
integer(i_kind), intent(  out) :: ivecs ! 
integer(i_kind), intent(  out) :: rc    ! error return code
! local variables
type(control_vector) :: aux
type(gsi_bundle)     :: mval(nsubwin)
type(predictors)     :: sbias
integer(i_kind)      :: ii,jj
real(r_kind)         :: zz

rc=0
npcvecs = nvcglev
ivecs=min(npcvecs,nwrvecs)
if (ivecs<1) then
   if (mype==0) write(6,*)'save_precond: cannot get siga, ivecs=', ivecs
   rc=1
   return
endif

call allocate_preds(sbias)
do ii=1,nsubwin
   call allocate_state(mval(ii))
end do
call allocate_cv(aux)

!-- calculate increment on analysis error covariance diag(delta p)
siga=zero
do jj=1,ivecs
   zz=sqrt(one-one/sqrt(rcglev(jj)))
   aux%values = zz * yvcglev(jj)%values
   call control2model(aux,mval,sbias)
   do ii=1,nsubwin
      call gsi_bundlehadamard(siga,mval(ii),mval(ii))
   enddo
enddo

do ii=1,nsubwin
!-- get b standard deviations
   call bkg_stddev(aux%step(ii),mval(ii))
!-- calculate diag(pa) = diag(b) - diag(delta p)
!   i.e., add diag(b) as rank-1 update to diag(delta p)
   call gsi_bundlehadamard(siga,mval(ii),mval(ii))
enddo

call deallocate_cv(aux)
do ii=1,nsubwin
   call deallocate_state(mval(ii))
end do
call deallocate_preds(sbias)

end subroutine congrad_siga
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!   save_precond - Save eigenvectors from congrad for next minimization
! ------------------------------------------------------------------------------
subroutine save_precond(ldsave)
!$$$  subprogram documentation block
!                .      .    .                                         .
! subprogram:    save_precond
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-08-05  lueken - added subprogram doc block
!
!   input argument list:
!    ldsave
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

implicit none

logical, intent(in   ) :: ldsave

real(r_kind), allocatable :: zmat(:,:)
integer(i_kind) :: ii,jj, info, iunit, ivecs
real(r_kind) :: zz
character(len=13) :: clfile

if (ldsave) then

!--- read eigenvalues of the preconditioner

   npcvecs = nvcglev+nvcglpc
   if (mype==0) write(6,*)'save_precond: NVCGLEV,NVCGLPC,NPCVECS=', &
                                         nvcglev,nvcglpc,npcvecs

   allocate(yvcglwk(npcvecs))
   ii=0
   if(.not.lconvert) then
      do jj=1,nvcglev
         ii=ii+1
         !  zz=sqrt(rcglpc(jj)-one)
         call allocate_cv(yvcglwk(ii))
         yvcglwk(ii)%values = yvcglev(jj)%values
         call deallocate_cv(yvcglev(jj))
      enddo
      if (allocated(yvcglev)) deallocate(yvcglev)
      
      nvcglev=0

   else
!--- copy preconditioner vectors to work file

      if (mype==0.and.nvcglpc>0) write(6,*)'save_precond: RCGLPC=',rcglpc
      do jj=1,nvcglpc
         ii=ii+1
         zz=sqrt(rcglpc(jj)-one)
         call allocate_cv(yvcglwk(ii))
         yvcglwk(ii)%values = zz * yvcglpc(jj)%values
         call deallocate_cv(yvcglpc(jj))
      enddo
      if (allocated(yvcglpc)) deallocate(yvcglpc)
!     if (allocated( rcglpc)) deallocate( rcglpc)
      nvcglpc=0

!--- copy and transform eigenvectors of preconditioned hessian

      if (mype==0.and.nvcglev>0) write(6,*)'save_precond: RCGLEV=',rcglev
      do jj=1,nvcglev
         ii=ii+1
         zz=sqrt(rcglev(jj)-one)
         call allocate_cv(yvcglwk(ii))
         yvcglwk(ii)%values = zz * yvcglev(jj)%values
         call deallocate_cv(yvcglev(jj))
      enddo
      if (allocated(yvcglev)) deallocate(yvcglev)
      if (allocated( rcglev)) deallocate( rcglev)
      nvcglev=0

      if (mype==0) write(6,*)'save_precond: NVCGLPC,NVCGLEV,npcvecs,ii=', &
                      nvcglpc,nvcglev,npcvecs,ii
      if (ii/=npcvecs) then
         write(6,*)'save_precond: error number of vectors',ii,npcvecs
         call stop2(139)
      end if

!---  form the inner matrix for the shermann-morrison-woodbury inversion

      allocate(zmat(npcvecs,npcvecs))
      do jj=1,npcvecs
         do ii=jj,npcvecs
            zmat(ii,jj) = dot_product (yvcglwk(jj),yvcglwk(ii))
         enddo
         zmat(jj,jj) = zmat(jj,jj) + one
      enddo

!--- Cholesky decompose

      if (mype==0) write(6,*)'save_precond: call dpotrf npcvecs=',npcvecs
      if (r_kind==n_default_real_kind) then
         call spotrf('L',npcvecs,zmat,npcvecs,info)
      elseif (r_kind==n_double_kind) then
         call dpotrf('L',npcvecs,zmat,npcvecs,info)
      else
         write(6,*)'save_precond: r_kind is neither default real nor double precision'
         call stop2(323)
      endif

      if (info/=0) then
         write(6,*)'save_precond: error computing Cholesky decomposition'
         write(6,*)'SPOTRF/DPOTRF returns info=',info
         call stop2(324)
      endif

!--- transform vectors

      do jj=1,npcvecs
         do ii=1,jj-1
            yvcglwk(jj)%values = yvcglwk(jj)%values - zmat(jj,ii)*yvcglwk(ii)%values
         enddo
         yvcglwk(jj)%values = yvcglwk(jj)%values / zmat(jj,jj)
      enddo

   endif

!--- Save the eigenvectors

   if (l4dvar) then
      ivecs=min(npcvecs,nwrvecs)
      do jj=1,ivecs
         clfile='evec.XXX.YYYY'
         write(clfile(6:8) ,'(I3.3)') jiter
         write(clfile(10:13),'(I4.4)') jj
         call write_cv(yvcglwk(jj),clfile)
      enddo

      if (mype==0) then
         iunit=78
         clfile='eval.XXX'
         write(clfile(6:8),'(I3.3)') jiter
         open(iunit,file=clfile)
         write(iunit,*)rcglev
         close(iunit)
      endif

      if (mype==0) then
         iunit=78
         clfile='numpcvecs.XXX'
         write(clfile(11:13),'(I3.3)') jiter
         open(iunit,file=clfile)
         write(iunit,*)ivecs
         close(iunit)
      endif

      do ii=1,npcvecs
         call deallocate_cv(yvcglwk(ii))
      enddo
      deallocate(yvcglwk)
   else
      do ii=nwrvecs+1,npcvecs
         call deallocate_cv(yvcglwk(ii))
      enddo
      npcvecs=min(npcvecs,nwrvecs)
   endif

   if (allocated(zmat)) deallocate(zmat)
endif

do ii=1,maxiter+1
   call deallocate_cv(cglwork(ii))
enddo
deallocate(cglwork)

return
end subroutine save_precond

! ------------------------------------------------------------------------------
!   setup_precond - Calculates the preconditioner for congrad
! ------------------------------------------------------------------------------
subroutine setup_precond()
!$$$  subprogram documentation block
!                .      .    .                                         .
! subprogram:    setup_precond
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-08-05  lueken - added subprogram doc block
!   2010-03-10  treadon - add essl interface
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

implicit none

integer(i_kind), allocatable :: indarr(:)
real(r_kind), allocatable :: zq(:),zlam(:),zu(:,:),zuut(:,:),zwork(:),zzz(:)
integer(i_kind) :: info,ik,inpcv,ji,jj,jk,ii,iunit
real(r_kind) :: za, zps
character(len=13) :: clfile

#ifdef ibm_sp
! Declare variables and work arrays used by essl
integer(i_kind):: iopt, ldz, n, naux
real(r_kind), allocatable ::  w(:),z(:,:),ap(:),aux(:)
#endif

!--- read vectors, apply change of variable and copy to work file

if (l4dvar) then
   iunit=78
   clfile='numpcvecs.XXX'
   write(clfile(11:13),'(I3.3)') jiter-1
   open(iunit,file=clfile)
   read(iunit,*)npcvecs
   close(iunit)
  
   if (npcvecs<1) then
      write(6,*)'SETUP_PRECOND: no vectors for preconditioner',npcvecs
      call stop2(140)
   end if
     
   allocate(yvcglwk(npcvecs))
   do ii=1,npcvecs
      call allocate_cv(yvcglwk(ii))
   enddo
     
   do jj=1,npcvecs
      clfile='evec.XXX.YYYY'
      write(clfile(6:8) ,'(I3.3)') jiter-1
      write(clfile(10:13),'(I4.4)') jj
      call read_cv(yvcglwk(jj),clfile)
   enddo
endif
if(.not. lconvert) then 
   nvcglpc=npcvecs   
   if(nvcglpc > 0) then   
      allocate (rcglpc(nvcglpc))

      iunit=78
      clfile='eval.XXX'
      write(clfile(6:8),'(I3.3)') jiter-1
      open(iunit,file=clfile)
      read(iunit,*)rcglpc
      close(iunit)
      do ii=1,nvcglpc
         rcglpc(ii) = min(r_max_cnum_pc,rcglpc(ii))
      end do
  

      allocate (yvcglpc(nvcglpc))
      do jj=1,nvcglpc
         call allocate_cv(yvcglpc(jj))
      enddo
      do jj=1,nvcglpc
         yvcglpc(jj) = zero
         do jk=1,yvcglpc(jj)%lencv
            yvcglpc(jj)%values(jk) = yvcglwk(jj)%values(jk)
         enddo
      enddo
 
      lmpcgl = .true.
   else
      nvcglpc = 0
      lmpcgl = .false.
   endif
else
   if (mype==0) write(6,*)'allocate arrays with npcvecs=',npcvecs
   allocate(indarr(npcvecs))
   allocate(zq(npcvecs),zlam(npcvecs),zu(npcvecs,npcvecs))
   allocate(zuut(npcvecs,npcvecs),zwork(3*npcvecs),zzz(npcvecs))

!--- Perform householder transformations to reduce the matrix of vectors
!--- to upper triangular

   do jj=1,npcvecs
      call allgather_cvsection(yvcglwk(jj),zq(1:jj),1,jj)
     
      zps = dot_product(yvcglwk(jj),yvcglwk(jj)) - dot_product(zq(1:jj),zq(1:jj))
     
      if (zq(jj) < zero) then
         zu(jj,jj) = -sqrt(zps+zq(jj)*zq(jj))
      else
         zu(jj,jj) =  sqrt(zps+zq(jj)*zq(jj))
      endif
     
      zq(jj) = zq(jj) - zu(jj,jj)

      do jk=1,jj-1
         zu(jk,jj) = zq(jk)
      enddo
     
      zps = zps + zq(jj)*zq(jj)
     
      zzz(1:jj-1)=zero
      zzz(jj)=zq(jj)
      call set_cvsection(zzz(1:jj),yvcglwk(jj),1,jj)

      do jk=1,yvcglwk(jj)%lencv
         yvcglwk(jj)%values(jk) = yvcglwk(jj)%values(jk) * sqrt(two/zps)
      enddo

!--- we now have the householder vector in yvcglwk(jj), and the non-zero
!--- elements of the transformed vector in zu. Now apply the householder
!--- transformations to the remaining vectors.
     
      do ji=jj+1,npcvecs
         zps = dot_product (yvcglwk(jj),yvcglwk(ji))
         do jk=1,yvcglwk(ji)%lencv
            yvcglwk(ji)%values(jk) = yvcglwk(ji)%values(jk) - zps*yvcglwk(jj)%values(jk)
         enddo
      enddo
   enddo

!--- Multiply the upper triangle by its transpose and find eigenvectors
!--- and eigenvalues

   do jj=1,npcvecs
      do ji=jj+1,npcvecs
         zu(ji,jj) = zero
      enddo
   enddo

   do jj=1,npcvecs
      do ji=jj,npcvecs
         zuut(ji,jj) = zero
         do jk=ji,npcvecs
            zuut(ji,jj) = zuut(ji,jj) + zu(ji,jk)*zu(jj,jk)
         enddo
      enddo
   enddo


#ifdef ibm_sp

!  Use essl
   iopt=1
   n=npcvecs
   ldz=n
   naux=3*n

   allocate(w(n),z(n,n),aux(naux),ap(n*n))
   w=zero
   z=zero
   aux=zero

!  Load zuut in essl lower-packed storage mode
   ap=zero
   jk=0
   do jj=1,n
      do ii=jj,n
         jk=jk+1
         ap(jk)=zuut(ii,jj)
      end do
   end do

!  Call essl routines
   if (r_kind==n_default_real_kind) then
      call sspev(iopt,ap,w,z,ldz,n,aux,naux)
   elseif (r_kind==n_double_kind) then
      call dspev(iopt,ap,w,z,ldz,n,aux,naux)
   else
      write(6,*)'SETUP_PRECOND: r_kind is neither default real nor double precision'
      call stop2(325)
   endif

!  Load essl results into output arrays
   do jj=1,n
      zlam(jj)=w(jj)
      do ii=1,n
         zuut(ii,jj)=z(ii,jj)
      end do
   end do
 
!  Deallocate work arrays
   deallocate(w,z,aux,ap)
 

#else
!  Use lapack routines
   if (r_kind==n_default_real_kind) then
      call ssyev('V','L',npcvecs,zuut,npcvecs,zlam,zwork,size(zwork),info)
   elseif (r_kind==n_double_kind) then
      call dsyev('V','L',npcvecs,zuut,npcvecs,zlam,zwork,size(zwork),info)
   else
      write(6,*)'SETUP_PRECOND: r_kind is neither default real nor double precision'
      call stop2(325)
   endif
   if (info/=0) then
      write(6,*)'SETUP_PRECOND: SSYEV/DSYEV returned with info=',info
      write(6,*)'SETUP_PRECOND: SSYEV/DSYEV returned non-zero return code'
      call stop2(326)
   endif
 
#endif

!--- convert to eigenvalues of the preconditioner

   do jk=1,npcvecs
      zlam(jk) = one / (one - zlam(jk))
   enddo

   if (mype==0) write(6,*)'SETUP_PRECOND: eigenvalues found are: ',(zlam(ji),ji=1,npcvecs)

!--- sort eigenvalues with eigenvalues larger than 1 after eigenvalues
!--- smaller than 1 and with eigenvalues larger than 1 sorted in decreasing
!--- order

   do ji=1,npcvecs
      indarr(ji) = ji
   enddo

!--- straight insertion sort courtesy of numerical recipies

   do jj=2,npcvecs
      za = zlam(jj)
      ik = indarr(jj)
      do ji=jj-1,1,-1
         if (zlam(ji)>one .and. (zlam(ji)>=za .or. za<=one)) then
            ii=ji
            exit
         else
            ii=0
         endif
         zlam(ji+1) = zlam(ji)
         indarr(ji+1) = indarr(ji)
      enddo
      zlam(ii+1) = za
      indarr(ii+1) = ik
   enddo

   inpcv = npcvecs

   do while (zlam(inpcv) <= zero)
      if (mype==0) write(6,*)'Warning - eigenvalue less than 1: ',zlam(inpcv)
      inpcv = inpcv-1
      if (inpcv == 0) then
         if (mype==0) write(6,*)'SETUP_PRECOND: cannot form preconditioner - '//&
            'no positive eigenvalues.'
         if (mype==0) write(6,*)'SETUP_PRECOND: minimisation will not be preconditioned.'
         exit
      endif
   enddo

   if (inpcv>0) then
      if (mype==0) write(6,*)'Number of preconditioning vectors selected is ',inpcv
      if (mype==0) write(6,*)'SETUP_PRECOND: selected eigenvalues are: ',(zlam(ji),ji=1,inpcv)
 
      if (allocated(yvcglpc)) then
         do jj=1,nvcglpc
            call deallocate_cv(yvcglpc(jj))
         enddo
         deallocate(yvcglpc)
         nvcglpc=0
      endif
      if (allocated(rcglpc)) deallocate(rcglpc)

      !--- Save eigenvalues
      nvcglpc = inpcv
      allocate (rcglpc(nvcglpc))
      rcglpc(:) = min(r_max_cnum_pc,zlam(1:nvcglpc))
 
      allocate (yvcglpc(nvcglpc))
      do jj=1,nvcglpc
         call allocate_cv(yvcglpc(jj))
      enddo

!--- apply householder transformations to the eigenvectors to get the
!--- eigenvectors of the preconditioner

      do jj=1,nvcglpc
         yvcglpc(jj) = zero
         call set_cvsection(zuut(1:npcvecs,indarr(jj)),yvcglpc(jj),1,npcvecs)
 
         do ji=npcvecs,1,-1
            zps = dot_product (yvcglwk(ji),yvcglpc(jj))
            do jk=1,yvcglpc(jj)%lencv
               yvcglpc(jj)%values(jk) = yvcglpc(jj)%values(jk) - zps*yvcglwk(ji)%values(jk)
            enddo
         enddo
      enddo
      lmpcgl = .true.
   else
      nvcglpc = 0
      lmpcgl = .false.
   endif

   deallocate(indarr)
   deallocate(zq,zlam,zu,zuut,zwork,zzz)
endif 
npcvecs = 0
do jj=1,npcvecs
   call deallocate_cv(yvcglwk(jj))
enddo
deallocate(yvcglwk)
  
  
return
end subroutine setup_precond

! ------------------------------------------------------------------------------
!   precond - Preconditioner for minimization
! ------------------------------------------------------------------------------
subroutine lanczos_precond(ycvx,kmat)
!$$$  subprogram documentation block
!                .      .    .                                         .
! subprogram:    precond
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-08-05  lueken - added subprogram doc block
!   2011-04-07  todling - renamed to from old precond name
!
!   input argument list:
!    ycvx
!    kmat
!
!   output argument list:
!    ycvx
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

implicit none

type(control_vector),intent(inout) :: ycvx
integer(i_kind)     ,intent(in   ) :: kmat

real(r_kind) :: zevals(nvcglpc),zdp(nvcglpc)
integer(i_kind) :: jk, ji

if     (kmat== 1    ) then
   zevals(:) = rcglpc(:)
elseif (kmat==-1    ) then
   zevals(:) = one/rcglpc(:)
elseif (kmat== 2) then
   zevals(1:nvcglpc) = sqrt(rcglpc(:))
elseif (kmat==-2) then
   zevals(1:nvcglpc) = one/sqrt(rcglpc(:))
else
   write(6,*)'Error: invalid value for kmat in precond: ',kmat
   write(6,*)'PRECOND: invalid value for kmat' 
   call stop2(327)
endif

do jk=1,nvcglpc
   zdp(jk) = (zevals(jk)-one)*dot_product(ycvx,yvcglpc(jk))
enddo

do jk=1,nvcglpc
   do ji=1,ycvx%lencv
      ycvx%values(ji) = ycvx%values(ji) + yvcglpc(jk)%values(ji) * zdp(jk)
   enddo
enddo

return
end subroutine lanczos_precond
! ------------------------------------------------------------------------------
subroutine read_lanczos(kmaxit)
!$$$  subprogram documentation block
!                .      .    .                                         .
! subprogram:    read_lanczos
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-08-05  lueken - added subprogram doc block
!
!   input argument list:
!    kmaxit
!
!   output argument list:
!    kmaxit
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

implicit none

integer(i_kind) , intent(inout) :: kmaxit

integer(i_kind) :: jj, iunit, kiter, ilen
character(len=17) :: clfile

if (kmaxit>maxiter) then
   write(6,*) 'read_lanczos: kmaxit>maxiter',kmaxit,maxiter
   call stop2(141)
end if

do jj=1,kmaxit
   clfile='lanczvec.XXX.YYYY'
   write(clfile(10:12),'(I3.3)') jiter
   write(clfile(14:17),'(I4.4)') jj
   call read_cv(cglwork(jj),clfile)
enddo

if (mype==0) then
   iunit=get_lun()
   clfile='zlanczos.XXX'
   write(clfile(10:12),'(I3.3)') jiter
   write(6,*)'Reading Lanczos coef. from file ',clfile

   open(iunit,file=trim(clfile),form='unformatted')
   read(iunit)kiter
   if (kiter>maxiter) then
      write(6,*)'read_laczos: kiter>maxiter',kiter,maxiter
      call stop2(142)
   end if
   read(iunit)zlancs(1:kiter+1,1:4)
   close(iunit)
endif
ilen=(kmaxit+1)*4
call mpl_bcast(0,ilen,zlancs)

end subroutine read_lanczos
! ------------------------------------------------------------------------------
end module lanczos
