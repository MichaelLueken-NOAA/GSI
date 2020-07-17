module m_lagnode
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module m_lagnode
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2016-05-18
!
! abstract: class-module of obs-type lagnode (lagrangian data)
!
! program history log:
!   2016-05-18  j guo   - added this document block for the initial polymorphic
!                         implementation.
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
  use m_obsdiagnode, only: obs_diag
  use m_obsdiagnode, only: obs_diags
  use kinds , only: i_kind,r_kind
  use mpeu_util, only: assert_,die,perr,warn,tell
  use m_obsnode, only: obsnode
  implicit none
  private

  public:: lagnode

  type,extends(obsnode):: lagnode
     !type(lag_ob_type),pointer :: llpoint => null()
     type(obs_diag), pointer :: diag_lon => null()
     type(obs_diag), pointer :: diag_lat => null()
     real(r_kind)    :: res_lon       ! residual
     real(r_kind)    :: res_lat       ! residual
     real(r_kind)    :: err2_lon      ! error squared
     real(r_kind)    :: err2_lat      ! error squared
     real(r_kind)    :: raterr2       ! square of ratio of final obs error 
                                      !  to original obs error
     real(r_kind)    :: obslon        ! observed longitude (rad)
     real(r_kind)    :: obslat        ! observed latitude  (rad)
     real(r_kind)    :: geslon        ! guessed longitude (rad)
     real(r_kind)    :: geslat        ! guessed latitude  (rad)
     real(r_kind)   ,dimension(:),allocatable :: specr  ! tl parameter
     !real(r_kind)    :: time          ! observation time in sec     
     real(r_kind)    :: b             ! variational quality control parameter
     real(r_kind)    :: pg            ! variational quality control parameter
     integer(i_kind),dimension(:),allocatable :: speci  ! tl parameter
     integer(i_kind) :: intnum        ! internal number of balloon
     !logical         :: luse          ! flag indicating if ob is used in pen.

     !integer(i_kind) :: idv,iob	      ! device id and obs index for sorting
     !real   (r_kind) :: elat, elon      ! earth lat-lon for redistribution
     !real   (r_kind) :: dlat, dlon      ! earth lat-lon for redistribution
  contains
     procedure,nopass::  mytype
     procedure::  sethop => obsnode_sethop_
     procedure::   xread => obsnode_xread_
     procedure::  xwrite => obsnode_xwrite_
     procedure:: isvalid => obsnode_isvalid_
     procedure::  gettlddp => gettlddp_

     ! procedure, nopass:: headerread  => obsheader_read_
     ! procedure, nopass:: headerwrite => obsheader_write_
     procedure:: init  => obsnode_init_
     procedure:: clean => obsnode_clean_
  end type lagnode

  public:: lagnode_typecast
  public:: lagnode_nextcast
     interface lagnode_typecast; module procedure typecast_ ; end interface
     interface lagnode_nextcast; module procedure nextcast_ ; end interface

  public:: lagnode_appendto
     interface lagnode_appendto; module procedure appendto_ ; end interface

  character(len=*),parameter:: myname="m_lagnode"

#include "myassert.H"
#include "mytrace.H"
contains
function typecast_(anode) result(ptr_)
!-- cast a class(obsnode) to a type(lagnode)
  use m_obsnode, only: obsnode
  implicit none
  type(lagnode ),pointer:: ptr_
  class(obsnode),pointer,intent(in):: anode
  ptr_ => null()
  if(.not.associated(anode)) return
        ! logically, typecast of a null-reference is a null pointer.
  select type(anode)
     type is(lagnode)
        ptr_ => anode
  end select
  return
end function typecast_

function nextcast_(anode) result(ptr_)
!-- cast an obsnode_next(obsnode) to a type(lagnode)
  use m_obsnode, only: obsnode,obsnode_next
  implicit none
  type(lagnode ),pointer:: ptr_
  class(obsnode),target ,intent(in):: anode

  class(obsnode),pointer:: inode_
  inode_ => obsnode_next(anode)
  ptr_ => typecast_(inode_)
  return
end function nextcast_

subroutine appendto_(anode,oll)
!-- append anode to linked-list oll
  use m_obsnode , only: obsnode
  use m_obsllist, only: obsllist,obsllist_appendnode
  implicit none
  type(lagnode),pointer,intent(in):: anode
  type(obsllist),intent(inout):: oll

  class(obsnode),pointer:: inode_
  inode_ => anode
  call obsllist_appendnode(oll,inode_)
  inode_ => null()
end subroutine appendto_

! obsnode implementations

function mytype()
  implicit none
  character(len=:),allocatable:: mytype
  mytype="[lagnode]"
end function mytype

subroutine obsnode_init_(anode)
  use lag_traj, only: lag_rk2itenpara_i,lag_rk2itenpara_r
  implicit none
  class(lagnode),intent(out):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_init_'
_ENTRY_(myname_)
  allocate(anode%speci(lag_rk2itenpara_i), &
           anode%specr(lag_rk2itenpara_r)  )
  return
end subroutine obsnode_init_

subroutine obsnode_clean_(anode)
  implicit none
  class(lagnode),intent(inout):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_clean_'
_ENTRY_(myname_)
!_TRACEV_(myname_,'%mytype() =',anode%mytype())
  if(allocated(anode%speci)) deallocate(anode%speci)
  if(allocated(anode%specr)) deallocate(anode%specr)
_EXIT_(myname_)
  return
end subroutine obsnode_clean_

subroutine obsnode_xread_(anode,iunit,istat,diaglookup,skip)
  use m_obsdiagnode, only: obsdiaglookup_locate
  implicit none
  class(lagnode) , intent(inout):: anode
  integer(i_kind) , intent(in   ):: iunit
  integer(i_kind) , intent(  out):: istat
  type(obs_diags) , intent(in   ):: diaglookup
  logical,optional, intent(in   ):: skip

  character(len=*),parameter:: myname_=myname//'::obsnode_xread_'
  logical:: skip_
_ENTRY_(myname_)
  skip_=.false.
  if(present(skip)) skip_=skip

  istat=0
  if(skip_) then
     read(iunit,iostat=istat)
     if(istat/=0) then
        call perr(myname_,'skipping read(%(res,err2,...)), istat =',istat)
        _EXIT_(myname_)
        return
     endif

  else
     read(iunit,iostat=istat)    anode%res_lon , &
                                 anode%res_lat , &
                                 anode%err2_lon, &
                                 anode%err2_lat, &
                                 anode%raterr2 , &
                                 anode%b       , &
                                 anode%pg      , &
                                 anode%obslon  , &
                                 anode%obslat  , &
                                 anode%geslon  , &
                                 anode%geslat  , &
                                 anode%intnum  , &
                                 anode%speci   , & !(lag_rk2itenpara_i)
                                 anode%specr       !(lag_rk2itenpara_r)
     if (istat/=0) then
        call perr(myname_,'read(%(res,err2,...)), istat =',istat)
        _EXIT_(myname_)
        return
     end if

     anode%diag_lon => obsdiaglookup_locate(diaglookup,anode%idv,anode%iob,1)
     anode%diag_lat => obsdiaglookup_locate(diaglookup,anode%idv,anode%iob,2)

     if(.not.( associated(anode%diag_lon).and. &
               associated(anode%diag_lat)      ) ) then
        call perr(myname_,'obsdiaglookup_locate(lon,lat), %idv =',anode%idv)
        call perr(myname_,'                               %iob =',anode%iob)
        if(.not.associated(anode%diag_lon)) &
        call perr(myname_,'       can not locate diag_lon, ich =',1)
        if(.not.associated(anode%diag_lat)) &
        call perr(myname_,'       can not locate diag_lat, ich =',2)
        call  die(myname_)
     endif
  endif
_EXIT_(myname_)
  return
end subroutine obsnode_xread_

subroutine obsnode_xwrite_(anode,junit,jstat)
  implicit none
  class(lagnode),intent(in):: anode
  integer(i_kind),intent(in   ):: junit
  integer(i_kind),intent(  out):: jstat

  character(len=*),parameter:: myname_=myname//'::obsnode_xwrite_'
_ENTRY_(myname_)

  jstat=0
  write(junit,iostat=jstat)     anode%res_lon , &
                                anode%res_lat , &
                                anode%err2_lon, &
                                anode%err2_lat, &
                                anode%raterr2 , &
                                anode%b       , &
                                anode%pg      , &
                                anode%obslon  , &
                                anode%obslat  , &
                                anode%geslon  , &
                                anode%geslat  , &
                                anode%intnum  , &
                                anode%speci   , & !(lag_rk2itenpara_i)
                                anode%specr       !(lag_rk2itenpara_r)
  if (jstat/=0) then
     call perr(myname_,'write(%(res,err2,...)), jstat =',jstat)
     _EXIT_(myname_)
     return
  end if
_EXIT_(myname_)
  return
end subroutine obsnode_xwrite_

subroutine obsnode_sethop_(anode)
  use m_cvgridlookup, only: cvgridlookup_getiw
  implicit none
  class(lagnode),intent(inout):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_sethop_'
  real(r_kind) :: dum
_ENTRY_(myname_)
  !-- yet to be defined
  call perr(myname_,'nothing about sethop has been defined')
  call  die(myname_)
  !-- following is here to satisfy var-usage requirement
  dum=anode%elat
_EXIT_(myname_)
  return
end subroutine obsnode_sethop_

function obsnode_isvalid_(anode) result(isvalid_)
  implicit none
  logical:: isvalid_
  class(lagnode),intent(in):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_isvalid_'
_ENTRY_(myname_)
  isvalid_=associated(anode%diag_lat) .and. &
           associated(anode%diag_lat)
_EXIT_(myname_)
  return
end function obsnode_isvalid_

pure subroutine gettlddp_(anode,jiter,tlddp,nob)
  use kinds, only: r_kind
  implicit none
  class(lagnode), intent(in):: anode
  integer(kind=i_kind),intent(in):: jiter
  real(kind=r_kind),intent(inout):: tlddp
  integer(kind=i_kind),optional,intent(inout):: nob

  tlddp = tlddp + anode%diag_lat%tldepart(jiter)*anode%diag_lat%tldepart(jiter)
  tlddp = tlddp + anode%diag_lon%tldepart(jiter)*anode%diag_lon%tldepart(jiter)
  if(present(nob)) nob=nob+2
  return
end subroutine gettlddp_

end module m_lagnode
