module gsi_toper
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module gsi_toper
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2018-08-10
!
! abstract: an oboper extension for tnode type
!
! program history log:
!   2018-08-10  j guo   - added this document block
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
  use m_tnode   , only: tnode
  implicit none
  public:: toper      ! data stracture

  type,extends(oboper):: toper
  contains
     procedure,nopass:: mytype
     procedure,nopass:: nodemold
     procedure:: setup_
     procedure:: intjo1_
     procedure:: stpjo1_
  end type toper

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname='gsi_toper'
  type(tnode),save,target:: mynodemold_

contains
  function mytype(nodetype)
    implicit none
    character(len=:),allocatable:: mytype
    logical,optional, intent(in):: nodetype
    mytype="[toper]"
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
    use t_setup, only: setup
    use kinds, only: i_kind
    use gsi_oboper, only: len_obstype
    use gsi_oboper, only: len_isis

    use m_rhs  , only: awork => rhs_awork
    use m_rhs  , only: bwork => rhs_bwork
    use m_rhs  , only: iwork => i_t

    use obsmod  , only: write_diag
    use convinfo, only: diag_conv
    use jfunc   , only: jiter

    use mpeu_util, only: die
    implicit none
    class(toper ), intent(inout):: self
    integer(i_kind), intent(in):: lunin
    integer(i_kind), intent(in):: mype
    integer(i_kind), intent(in):: is
    integer(i_kind), intent(in):: nobs
    logical        , intent(in):: init_pass     ! supporting multi-pass setup()
    logical        , intent(in):: last_pass     ! with incremental backgrounds.

    !----------------------------------------
    character(len=*),parameter:: myname_=myname//"::setup_"

    character(len=len_obstype):: obstype
    character(len=len_isis   ):: isis
    integer(i_kind):: nreal,nchanl,ier,nele
    logical:: diagsave

    if(nobs == 0) return

    read(lunin,iostat=ier) obstype,isis,nreal,nchanl
    if(ier/=0) call die(myname_,'read(obstype,...), iostat =',ier)
    nele = nreal+nchanl

    diagsave  = write_diag(jiter) .and. diag_conv

    call setup(self%obsll(:), self%odiagll(:), &
        lunin,mype,bwork,awork(:,iwork),nele,nobs,is,diagsave)

  end subroutine setup_

  subroutine intjo1_(self, ibin, rval,sval, qpred,sbias)
    use inttmod, only: intjo => intt
    use gsi_bundlemod  , only: gsi_bundle
    use bias_predictors, only: predictors
    use bias_predictors, only: predictors_getdim
    use m_obsnode , only: obsnode
    use m_obsllist, only: obsllist_headnode
    use kinds     , only: i_kind, r_quad
    implicit none
    class(toper  ),intent(in   ):: self
    integer(i_kind ),intent(in   ):: ibin
    type(gsi_bundle),intent(inout):: rval   ! (ibin)
    type(gsi_bundle),intent(in   ):: sval   ! (ibin)
    real(r_quad    ),target,dimension(:),intent(inout):: qpred  ! (ibin)
    type(predictors),target,             intent(in   ):: sbias

    !----------------------------------------
    character(len=*),parameter:: myname_=myname//"::intjo1_"
    integer(i_kind):: i,l,n
    class(obsnode),pointer:: headnode

! Are the different calls to intt() with optional arguments realy needed? 
! There is no checking of present(rpred) or present(spred) inside intt_()
! anyway.  Other logic is used to avoid accessing non-present rpred(:) and
! spred(:).

    call predictors_getdim(lbnd_t=i,ubnd_t=l,size_t=n)
    headnode => obsllist_headnode(self%obsll(ibin))
    if(n>0) then
       call intjo(headnode, rval,sval, qpred(i:l),sbias%predt)
    else
       call intjo(headnode, rval,sval)
    endif
    headnode => null()

  end subroutine intjo1_

  subroutine stpjo1_(self, ibin, dval,xval,pbcjo,sges,nstep,dbias,xbias)
    use stptmod, only: stpjo => stpt
    use gsi_bundlemod, only: gsi_bundle
    use bias_predictors, only: predictors, predictors_getdim
    use aircraftinfo, only: npredt,ntail,aircraft_t_bc_pof,aircraft_t_bc
    use m_obsnode , only: obsnode
    use m_obsllist, only: obsllist_headnode
    use kinds, only: r_quad,r_kind,i_kind
    implicit none
    class(toper  ),intent(in):: self
    integer(i_kind ),intent(in):: ibin
    type(gsi_bundle),intent(in):: dval
    type(gsi_bundle),intent(in):: xval
    real(r_quad    ),dimension(:),intent(inout):: pbcjo ! (1:4)
    real(r_kind    ),dimension(:),intent(in   ):: sges
    integer(i_kind),intent(in):: nstep

    type(predictors),target, intent(in):: dbias
    type(predictors),target, intent(in):: xbias

    !----------------------------------------
    character(len=*),parameter:: myname_=myname//"::stpjo1_"
    class(obsnode),pointer:: headnode
    real(r_kind),pointer,dimension(:,:) :: dpred,xpred
    integer(i_kind):: n

! Are the different calls to stpt() with optional arguments realy needed? 
! There is no checking of present(rpred) or present(spred) inside intt_()
! anyway.  Other logic is used to avoid accessing non-present rpred(:) and
! spred(:).

    headnode => obsllist_headnode(self%obsll(ibin))
    call predictors_getdim(size_t=n)
    if(n<=0 .or. .not. (aircraft_t_bc_pof .or. aircraft_t_bc)) then
       call stpjo(headnode,dval,xval,pbcjo(:),sges,nstep)
    else
       dpred(1:npredt,1:ntail) => dbias%predt(1:n)
       xpred(1:npredt,1:ntail) => xbias%predt(1:n)
       call stpjo(headnode,dval,xval,pbcjo(:),sges,nstep,dpred,xpred)
       dpred => null()
       xpred => null()
    endif
    headnode => null()
  end subroutine stpjo1_

end module gsi_toper
