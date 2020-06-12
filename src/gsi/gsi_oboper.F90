module gsi_oboper
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module gsi_oboper
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2018-06-26
!
! abstract: gsi observation operator, bundling obs_diags and obsllist objects
!
! program history log:
!   2018-06-26  j guo   - a new module for abstract GSI oboper.
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

  use m_obsdiagnode, only: obs_diags
  use m_obsllist   , only: obsllist

  use kinds, only: i_kind
  use mpeu_util, only: assert_
  implicit none
  private       ! except
  public :: oboper              ! data structure
  public :: len_obstype
  public :: len_isis

  integer(i_kind),parameter:: len_obstype=10
  integer(i_kind),parameter:: len_isis   =20

        ! oboper is a bundle of observation operator arrays (or lists), such as
        ! linked-lists of obs_diag (obs_diags) and obsnode (obsllist), plus type
        ! specific parameters.
        !
        ! In this implementation, an oboper, with pointers _associated_ to
        ! rank-1 arrays of obs_diags and obsllist, where both targets are
        ! instantiated separately with own fixed dimensions in nobs_type and
        ! nobs_bins.
        !
        ! It is planned in the future, to implement an oboper _contains_ dynamic
        ! components of these rank-1 arrays.

  type,abstract:: oboper
    !private
        ! In the first oboper implementation, %obsll(:) and odiagll(:) are
        ! treated as aliases to the instances of m_obsdiags::obsdiags(:,:) and
        ! m_obsdiagss::obsllists(:,:).  Both linked-lists are dimensioned for
        ! 1:nobs_type in the current implementation, and accesssed once per type
        ! and per bin, in intjo() and stpjo().
        !
        ! On the other hand, in the current setuprhsall() implementation, oboper
        ! objects are accessed for 1:ndat, or once per obs-stream, where each
        ! type is in general accessed in zero or multiple times.

     type(obs_diags),pointer,dimension(:):: odiagll    ! (1:nobs_bins)
     type(obsllist ),pointer,dimension(:)::   obsll    ! (1:nobs_bins)

     contains
       procedure(mytype  ),deferred,nopass:: mytype    ! type information
       procedure(nodemold),deferred,nopass:: nodemold  ! type information

       procedure, non_overridable:: init  => init_         ! initialize
       procedure, non_overridable:: clean => clean_        ! finalize

       generic:: setup => setup_
          procedure(setup_ ),deferred:: setup_      ! incremental object initialization
       generic:: intjo => intjo_, intjo1_
          procedure,  non_overridable:: intjo_      ! interface supporting intjo()
          procedure(intjo1_),deferred:: intjo1_     ! interface for 1-bin intjo()
       generic:: stpjo => stpjo_, stpjo1_
          procedure,  non_overridable:: stpjo_      ! interface supporting stpjo()
          procedure(stpjo1_),deferred:: stpjo1_     ! interface for 1-bin stpjo()

  end type oboper

        ! In setuprhsall(),
        !
        !   | use m_obsdiags, only: oboper_associate, oboper_dissociate
        !   | use gsi_oboper, only: oboper
        !   | use gsi_obopertypemanager, only: oboper_typemold
        !   | use obsmod, only: ndat,dtype
        !
        ! then in a loop of obs-streams
        !
        !   | class(oboper),pointer:: my_oboper
        !   | do is=1,ndat
        !   |    my_oboper => oboper_associate(oboper_typemold(dtype(is)))
        !   |    call my_oboper%setup(...)
        !   |    call oboper_dissociate(my_oboper)
        !   | enddo
        !

        ! In intjo() or stpjo(),
        !
        !   | use gsi_obopertypemanager, only: lbound_oboper
        !   | use gsi_obopertypemanager, only: ubound_oboper
        !   | use gsi_obopertypemanager, only: oboper_typemold
        ! 
        ! then in a loop of oboper
        !
        !   | class(oboper),pointer:: my_oboper
        !   | do iop=lbound_oboper,ubound_oboper
        !   |    my_oboper => oboper_associate(oboper_typemold(iOp))
        !   |    call my_oboper%intjo(...)
        !   |    call oboper_dissociate(my_oboper)
        !   | enddo

!--- Design considerations ---
! (1) Fully objectize oboper, meaning, capable of being instantiated where and
!     when it is needed.
!
!
! (2) For continuity, its instantiation is a type-indexed array of polymorphic
!     class(oboper), containing rank-1 pointers aliased to obsllist(1:nobs_bins)
!     and diagllist(1:nobs_bins).  This means its current instantiation is
!     declared based on a type-wrapper-array structure,
!
!       type,abstract:: oboper; ...
!       type:: oboper_element; class(oboper),pointer:: ptr; ...
!       type(oboper_element),dimension(nobs_type):: obopers
!
!     defined in a type-loop, (m_obsdiags?)
!
!       allocate(obopers(it)%ptr,mold=oboper_typemold(it))
!
!       | oboper_typemold(it) result(mold)
!       |   select case(it)
!       |      case(iobtype_rad); mold => radoper_mold()
!       |      case ...
!
!     followed by
!
!       associate(i_op => obopers(...)%ptr)
!          call i_op%init(...)   # type-bound init(), with a line of
!                                # self%nodetype=oboper%mytype(nodetype=.true.)
!       end associate
!
!
! (3) In future implementations, one might want to define oboper on a per-stream
!     base.  Then it would be instantiated in a stream-loop,
!
!       allocate(obopers(is)%ptr,mold=oboper_typemold(dtype(is)))
!
!       | oboper_typemold(dtype) result(mold)
!       |  select case(dtype)
!       |     case("rad","amsua",...); mold => radoper_mold()
!       |     case ...
!  
! (4) So types of obopers are now one-to-one mapped to obsnode types.  This means
!     that each oboper type must be hardwired to a known obsnode type, while
!     dtype(:) to obopers(:) types are not.
!     

!--------- interfaces
abstract interface
  function mytype(nodetype)
        ! %mytype() for self's typename
        ! %mytype(nodetype=.true.) for self's corresponding node type name
    implicit none
    character(len=:), allocatable:: mytype
    logical, optional, intent(in):: nodetype      ! if .true., return %mytype() of its obsnode

    ! logical:: nodetype_
    ! nodetype_=.false.
    ! if(present(nodetype)) nodetype_=nodetype
    ! if(nodetype_) then
    !    if(nodetype) mytype=mynodemold_%mytype()
    ! else
    !    mytype="[radoper]"
    ! endif

  end function mytype
end interface

abstract interface
  function nodemold()
  !> %nodemold() returns a mold of its corresponding obsnode
    use m_obsnode, only: obsnode
    implicit none
    class(obsnode),pointer:: nodemold

    !> For a given
    !>   type(someoper):: myoper
    !>
    !> then code
    !>
    !>   class(obsnode),pointer:: mynodemold_
    !>   mynodemold_ =>  myoper%nodemold()
    !>
    !> would return a mold of myoper's corresponding obsnode type

  end function nodemold
end interface

abstract interface
  subroutine setup_(self, lunin, mype, is, nobs, init_pass,last_pass)
    use kinds, only: i_kind
    import:: oboper
    implicit none
    class(oboper  ), intent(inout):: self
    integer(i_kind), intent(in):: lunin
    integer(i_kind), intent(in):: mype
    integer(i_kind), intent(in):: is
    integer(i_kind), intent(in):: nobs
    logical        , intent(in):: init_pass     ! supporting multi-pass setup()
    logical        , intent(in):: last_pass     ! with incremental backgrounds.

    ! An example in radoper%setup(),
    !
    !   if(nobs == 0) return
    !
    !   read(lunin,iostat=ier) obstype,isis,nreal,nchanl
    !   if(ier/=0) call die(myname_,'read(), iostat =',ier)
    !   nele=nreal+nchanl
    !
    !   call setuprad(self%obsll(:),self%odiagll(:), lunin, mype,   &
    !     aivals,stats,nchanl,nreal,nobs,obstype,isis,is,rad_diagsave,init_pass,last_pass)

  end subroutine setup_
end interface

abstract interface
  !>> This is the interface for single bin intjo().
  !>> call self%intjo(ib,rval(ib),sval(ib),qpred(:,ib),sbias)

  subroutine intjo1_(self, ibin, rval, sval, qpred, sbias)
    use gsi_bundlemod  , only: gsi_bundle
    use bias_predictors, only: predictors
    use kinds          , only: i_kind, r_quad
    import:: oboper
    implicit none
    class(oboper   ), intent(in   ):: self
    integer(i_kind ), intent(in   ):: ibin
    type(gsi_bundle), intent(inout):: rval
    type(gsi_bundle), intent(in   ):: sval
    real(r_quad    ), target, dimension(:),intent(inout):: qpred ! a buffer of rbias
    type(predictors), target,              intent(in   ):: sbias

        ! This implementation can be used both to an oboper instance with
        ! multiple bins, or a "slice" of oboper instance with a single bin,
        ! where the slice of self contains arrays (ibin:ibin) of components.

    !do ibin=lbound(self%obsll,1),ubound(self%obsll,1)
    !   call self%intjo(ibin, rval(ibin),sval(ibin), qpred(:,ibin),sbias)
    !enddo
  end subroutine intjo1_
end interface

abstract interface
  !>> This is the interface for single bin stpjo().
  !>> call self%stpjo(ib,dval(ib),xval(ib),pbcjo(:,it,ib),sges,nstep,dbias,xbias)

  subroutine stpjo1_(self, ibin,dval,xval,pbcjo,sges,nstep,dbias,xbias)
    use gsi_bundlemod  , only: gsi_bundle
    use bias_predictors, only: predictors
    use kinds          , only: r_quad,r_kind,i_kind
    import:: oboper
    implicit none
    class(oboper  ),intent(in):: self
    integer(i_kind),intent(in):: ibin

    type(gsi_bundle),intent(in   ):: dval
    type(gsi_bundle),intent(in   ):: xval
    real(r_quad    ),dimension(:),intent(inout):: pbcjo ! (1:4)
    real(r_kind    ),dimension(:),intent(in   ):: sges
    integer(i_kind ),             intent(in   ):: nstep
    type(predictors), target, intent(in):: dbias
    type(predictors), target, intent(in):: xbias

  end subroutine stpjo1_
end interface

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  character(len=*),parameter:: myname="gsi_oboper"

contains
#include "myassert.H"

subroutine init_(self,obsll,odiagll)
  implicit none
  class(oboper),intent(inout):: self
  type(obsllist ),target,dimension(:),intent(in):: obsll
  type(obs_diags),target,dimension(:),intent(in):: odiagll

  self%odiagll => odiagll(:)
  self%  obsll => obsll(:)
end subroutine init_

subroutine clean_(self)
  implicit none
  class(oboper),intent(inout):: self
  self%odiagll => null()
  self%  obsll => null()
end subroutine clean_

subroutine intjo_(self, rval,sval,qpred,sbias)
  use gsi_bundlemod  , only: gsi_bundle
  use bias_predictors, only: predictors
  use kinds, only: i_kind, r_quad
  implicit none
  class(oboper   ), intent(in):: self
  type(gsi_bundle), dimension(:  ),intent(inout):: rval
  type(gsi_bundle), dimension(:  ),intent(in   ):: sval
  real(r_quad    ), dimension(:,:),intent(inout):: qpred
  type(predictors)                ,intent(in   ):: sbias

        ! nb=nobs_bins
        ! do ityp=1,nobs_type
        !   iop => oboper_associate(mold=oboper_typemold(ityp))
        !   call iop%intjo(rval(:nb),sval(:nb), qpred(:,:nb),sbias)
        !   call oboper_dissociate(iop)
        ! enddo
        !
        ! This implementation can be used both to an oboper instance with
        ! multiple bins, or a "slice" of oboper instance with a single bin,
        ! where the slice of self contains arrays (ibin:ibin) of components.

  character(len=*),parameter:: myname_=myname//"::intjo_"
  integer(i_kind):: lbnd,ubnd,ibin

  lbnd = lbound(self%obsll,1)
  ubnd = ubound(self%obsll,1)
  ASSERT(lbnd == lbound( rval,1) .and. ubnd == ubound( rval,1))
  ASSERT(lbnd == lbound( sval,1) .and. ubnd == ubound( sval,1))
  ASSERT(lbnd == lbound(qpred,2) .and. ubnd == ubound(qpred,2))

  do ibin=lbnd,ubnd
     call self%intjo(ibin,rval(ibin),sval(ibin),qpred(:,ibin),sbias)
  enddo
end subroutine intjo_

subroutine stpjo_(self, dval,xval, pbcjo,sges,nstep, dbias,xbias)
  use kinds, only: r_quad,r_kind,i_kind
  use gsi_bundlemod, only: gsi_bundle
  use bias_predictors, only: predictors
  implicit none
  class(oboper   ),intent(in):: self
  type(gsi_bundle),dimension(  :),intent(in   ):: dval
  type(gsi_bundle),dimension(  :),intent(in   ):: xval
  real(r_quad    ),dimension(:,:),intent(inout):: pbcjo ! (1:4,1:nbin)
  real(r_kind    ),dimension(:  ),intent(in   ):: sges
  integer(i_kind),intent(in):: nstep

  type(predictors), intent(in):: dbias
  type(predictors), intent(in):: xbias

  integer(i_kind):: lbnd,ubnd,ibin

  lbnd = lbound(self%obsll,1)
  ubnd = ubound(self%obsll,1)
  ASSERT(lbnd == lbound( dval,1) .and. ubnd == ubound( dval,1))
  ASSERT(lbnd == lbound( xval,1) .and. ubnd == ubound( xval,1))
  ASSERT(lbnd == lbound(pbcjo,2) .and. ubnd == ubound(pbcjo,2))

  do ibin=lbnd,ubnd
     call self%stpjo(ibin,dval(ibin),xval(ibin),pbcjo(:,ibin),sges,nstep,dbias,xbias)
  enddo
end subroutine stpjo_
end module gsi_oboper
!.
