module gsi_lcbasoper
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module gsi_lcbasoper
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2018-08-10
!
! abstract: an oboper extension for lcbasnode type
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

  use gsi_oboper , only: oboper
  use m_lcbasnode, only: lcbasnode
  implicit none
  public:: lcbasoper      ! data stracture

  type,extends(oboper):: lcbasoper
  contains
     procedure,nopass:: mytype
     procedure,nopass:: nodemold
     procedure:: setup_
     procedure:: intjo1_
     procedure:: stpjo1_
  end type lcbasoper

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname='gsi_lcbasoper'
  type(lcbasnode),save,target:: mynodemold_

contains
  function mytype(nodetype)
    implicit none
    character(len=:),allocatable:: mytype
    logical,optional, intent(in):: nodetype
    mytype="[lcbasoper]"
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
    use lcbas_setup, only: setup
    use kinds, only: i_kind
    use gsi_oboper, only: len_obstype
    use gsi_oboper, only: len_isis

    use m_rhs  , only: awork => rhs_awork
    use m_rhs  , only: bwork => rhs_bwork
    use m_rhs  , only: iwork => i_lcbas

    use obsmod  , only: write_diag
    use convinfo, only: diag_conv
    use jfunc   , only: jiter

    use mpeu_util, only: die
    implicit none
    class(lcbasoper ), intent(inout):: self
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
    use intlcbasmod, only: intjo => intlcbas
    use gsi_bundlemod  , only: gsi_bundle
    use bias_predictors, only: predictors
    use m_obsnode , only: obsnode
    use m_obsllist, only: obsllist_headnode
    use kinds     , only: i_kind, r_quad
    implicit none
    class(lcbasoper  ),intent(in   ):: self
    integer(i_kind ),intent(in   ):: ibin
    type(gsi_bundle),intent(inout):: rval   ! (ibin)
    type(gsi_bundle),intent(in   ):: sval   ! (ibin)
    real(r_quad    ),target,dimension(:),intent(inout):: qpred  ! (ibin)
    type(predictors),target,             intent(in   ):: sbias

    !----------------------------------------
    character(len=*),parameter:: myname_=myname//"::intjo1_"
    class(obsnode),pointer:: headnode

    headnode => obsllist_headnode(self%obsll(ibin))
    call intjo(headnode, rval,sval)
    headnode => null()

  end subroutine intjo1_

  subroutine stpjo1_(self, ibin, dval,xval,pbcjo,sges,nstep,dbias,xbias)
    use stplcbasmod, only: stpjo => stplcbas
    use gsi_bundlemod, only: gsi_bundle
    use bias_predictors, only: predictors
    use m_obsnode , only: obsnode
    use m_obsllist, only: obsllist_headnode
    use kinds, only: r_quad,r_kind,i_kind
    implicit none
    class(lcbasoper  ),intent(in):: self
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

    headnode => obsllist_headnode(self%obsll(ibin))
    call stpjo(headnode,dval,xval,pbcjo(:),sges,nstep)
    headnode => null()
  end subroutine stpjo1_

end module gsi_lcbasoper
