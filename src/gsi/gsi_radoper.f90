module gsi_radoper
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module gsi_radoper
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2018-08-10
!
! abstract: an oboper extension for radnode type
!
! program history log:
!   2018-08-10  j guo   - added this document block for initial implementation
!
!   input argument list: see Fortran 90 style document below
!
!   output argument list: see Fortran 90 style document below
!
! attributes:
!   language: Fortran 90 and/or above
!   machine:
!
!$$$  end subprogram documentation block

! module interface:

  use gsi_oboper, only: oboper
  use m_radnode , only: radnode
  implicit none
  public:: radoper      ! data stracture

  type,extends(oboper):: radoper
  contains
     procedure,nopass:: mytype
     procedure,nopass:: nodemold
     procedure:: setup_
     procedure:: intjo1_
     procedure:: stpjo1_
  end type radoper

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname='gsi_radoper'
  type(radnode),save,target:: mynodemold_

contains
  function mytype(nodetype)
    implicit none
    character(len=:),allocatable:: mytype
    logical,optional, intent(in):: nodetype
    mytype="[radoper]"
    if(present(nodetype)) then
       if(nodetype) mytype=mynodemold_%mytype()
    endif
  end function mytype

  function nodemold()
  !> %nodemold() returns a mold of its corresponding obsnode
    use m_obsnode, only: obsnode
    implicit none
    class(obsnode),pointer:: nodemold
    nodemold => mynodemold_
  end function nodemold

  subroutine setup_(self, lunin, mype, is, nobs, init_pass,last_pass)
    use rad_setup, only: setup
    use kinds, only: i_kind
    use gsi_oboper, only: len_obstype
    use gsi_oboper, only: len_isis

    use m_rhs  , only: aivals => rhs_aivals
    use m_rhs  , only: stats  => rhs_stats

    use obsmod , only: write_diag
    use radinfo, only: diag_rad
    use jfunc  , only: jiter

    use mpeu_util, only: die
    implicit none
    class(radoper ), intent(inout):: self
    integer(i_kind), intent(in):: lunin
    integer(i_kind), intent(in):: mype
    integer(i_kind), intent(in):: is
    integer(i_kind), intent(in):: nobs
    logical        , intent(in):: init_pass     ! supporting multi-pass setup()
    logical        , intent(in):: last_pass     ! with incremental backgrounds.

    character(len=*),parameter:: myname_=myname//"::setup_"

    character(len=len_obstype):: obstype
    character(len=len_isis   ):: isis
    integer(i_kind):: nreal,nchanl,ier
    logical:: diagsave

    if(nobs == 0) return

    read(lunin,iostat=ier) obstype,isis,nreal,nchanl
    if(ier/=0) call die(myname_,'read(obstype,...), iostat =',ier)

    diagsave = write_diag(jiter) .and. diag_rad

    call setup(self%obsll(:), self%odiagll(:), lunin, mype,   &
        aivals,stats,nchanl,nreal,nobs,obstype,isis,is,diagsave, &
        init_pass,last_pass)

  end subroutine setup_

  subroutine intjo1_(self, ibin, rval,sval, qpred,sbias)
  !-- a single bin interface, with %pred resolved
    use intradmod, only: intjo => intrad
    use gsi_bundlemod  , only: gsi_bundle
    use bias_predictors, only: predictors
    use bias_predictors, only: predictors_getdim
    use m_obsnode , only: obsnode
    use m_obsllist, only: obsllist_headnode
    use kinds     , only: i_kind, r_quad
    implicit none
    class(radoper  ),intent(in   ):: self
    integer(i_kind ),intent(in   ):: ibin
    type(gsi_bundle),intent(inout):: rval   ! (ibin)
    type(gsi_bundle),intent(in   ):: sval   ! (ibin)
    real(r_quad    ),target,dimension(:),intent(inout):: qpred  ! (ibin)
    type(predictors),target,intent(in   ):: sbias

        !!$omp ...
        ! do ib=1,nobs_bins
        !    do it=1,nobs_type
        !       iop => oboper_create(mold=oboper_typemold(it))
        !       call iop%intjo(ib, rval(ib),sval(ib), qpred(:,ib),sbias)
        !       iop => null()
        !    enddo
        ! enddo

    character(len=*),parameter:: myname_=myname//"::intjo1_"
    integer(i_kind):: i,l
    class(obsnode),pointer:: headnode

    call predictors_getdim(lbnd_s=i,ubnd_s=l)
    headnode => obsllist_headnode(self%obsll(ibin))
    call intjo(headnode, rval,sval, qpred(i:l),sbias%predr)
    headnode => null()

  end subroutine intjo1_

  subroutine stpjo1_(self, ibin, dval,xval,pbcjo,sges,nstep,dbias,xbias)
  !-- a single bin interface, with %pred resolved
    use stpradmod, only: stpjo => stprad
    use kinds, only: r_quad,r_kind,i_kind
    use gsi_bundlemod, only: gsi_bundle
    use bias_predictors, only: predictors,predictors_getdim
    use radinfo, only: npred,jpch_rad
    use m_obsnode , only: obsnode
    use m_obsllist, only: obsllist_headnode
    implicit none
    class(radoper  ),intent(in):: self
    integer(i_kind ),intent(in):: ibin
    type(gsi_bundle),intent(in):: dval
    type(gsi_bundle),intent(in):: xval
    real(r_quad    ),dimension(:),intent(inout):: pbcjo ! (1:4)
    real(r_kind    ),dimension(:),intent(in   ):: sges
    integer(i_kind),intent(in):: nstep

    type(predictors),target, intent(in):: dbias
    type(predictors),target, intent(in):: xbias

    character(len=*),parameter:: myname_=myname//"::stpjo1_"
    class(obsnode),pointer:: headnode
    real(r_kind),pointer,dimension(:,:):: dpred,xpred
    integer(i_kind):: n

    headnode => obsllist_headnode(self%obsll(ibin))

    call predictors_getdim(size_s=n)
    dpred(1:npred,1:jpch_rad) => dbias%predr(1:n)
    xpred(1:npred,1:jpch_rad) => xbias%predr(1:n)

    call stpjo(headnode,dval,xval, dpred,xpred,pbcjo(:),sges,nstep)

    dpred => null()
    xpred => null()

    headnode => null()
  end subroutine stpjo1_

end module gsi_radoper
