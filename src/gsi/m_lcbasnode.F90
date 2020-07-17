module m_lcbasnode
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module m_lcbasnode
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2016-05-18
!
! abstract: class-module of obs-type lcbasnode (base height of lowest cloud seen)
!
! program history log:
!   2016-05-18  j guo   - added this document block for the initial polymorphic
!                         implementation.
!   2017-01-19  j guo   - fixed a null reference bug in nextcast_().
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

  public:: lcbasnode

  type,extends(obsnode):: lcbasnode
     ! private
     !type(lcbas_ob_type),pointer :: llpoint => null()
     type(obs_diag), pointer :: diags => null()
     real(r_kind)    :: res           !  lcbas residual
     real(r_kind)    :: err2          !  lcbas error squared
     real(r_kind)    :: raterr2       !  square of ratio of final obs error
                                      !  to original obs error
     !real(r_kind)    :: time          !  observation time in sec
     real(r_kind)    :: b             !  variational quality control parameter
     real(r_kind)    :: pg            !  variational quality control parameter
     real(r_kind)    :: wij(4)        !  horizontal interpolation weights
     integer(i_kind) :: ij(4)         !  horizontal locations
     !logical         :: luse          !  flag indicating if ob is used in pen.

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
     ! procedure:: init  => obsnode_init_
     ! procedure:: clean => obsnode_clean_
  end type lcbasnode

  public:: lcbasnode_typecast
  public:: lcbasnode_nextcast
     interface lcbasnode_typecast; module procedure typecast_ ; end interface
     interface lcbasnode_nextcast; module procedure nextcast_ ; end interface
 
  public:: lcbasnode_appendto
     interface lcbasnode_appendto; module procedure appendto_ ; end interface

  character(len=*),parameter:: myname="m_lcbasnode"

#include "myassert.H"
#include "mytrace.H"
contains
function typecast_(anode) result(ptr_)
!-- cast a class(obsnode) to a type(lcbasnode)
  use m_obsnode, only: obsnode
  implicit none
  type(lcbasnode),pointer:: ptr_
  class(obsnode ),pointer,intent(in):: anode
  ptr_ => null()
  if(.not.associated(anode)) return
        ! logically, typecast of a null-reference is a null pointer.
  select type(anode)
     type is(lcbasnode)
        ptr_ => anode
  end select
  return
end function typecast_

function nextcast_(anode) result(ptr_)
!-- cast an obsnode_next(obsnode) to a type(lcbasnode)
  use m_obsnode, only: obsnode,obsnode_next
  implicit none
  type(lcbasnode),pointer:: ptr_
  class(obsnode ),target ,intent(in):: anode

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
  type(lcbasnode),pointer,intent(in):: anode
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
  mytype="[lcbasnode]"
end function mytype

subroutine obsnode_xread_(anode,iunit,istat,diaglookup,skip)
  use m_obsdiagnode, only: obsdiaglookup_locate
  implicit none
  class(lcbasnode), intent(inout):: anode
  integer(i_kind), intent(in   ):: iunit
  integer(i_kind), intent(  out):: istat
  type(obs_diags), intent(in   ):: diaglookup
  logical,optional,intent(in   ):: skip

  character(len=*),parameter:: myname_=myname//'::obsnode_xread_'
  logical:: skip_
_ENTRY_(myname_)

  skip_=.false.
  if(present(skip)) skip_=skip

  if(skip_) then
     read(iunit,iostat=istat)
     if(istat/=0) then
        call perr(myname_,'skipping read(%(res,...)), iostat =',istat)
        _EXIT_(myname_)
        return
     endif

  else
     read(iunit,iostat=istat)    anode%res    , &
                                 anode%err2   , &
                                 anode%raterr2, &
                                 anode%b      , &
                                 anode%pg     , &
                                 anode%wij    , &
                                 anode%ij
     if (istat/=0) then
        call perr(myname_,'read(%(res,...)), iostat =',istat)
        _EXIT_(myname_)
        return
     end if

     anode%diags => obsdiaglookup_locate(diaglookup,anode%idv,anode%iob,1)
     if(.not.  associated(anode%diags)) then
        call perr(myname_,'obsdiaglookup_locate(), %idv =',anode%idv)
        call perr(myname_,'                        %iob =',anode%iob)
        call  die(myname_)
     endif
  endif ! if(skip_); else
_EXIT_(myname_)
  return
end subroutine obsnode_xread_

subroutine obsnode_xwrite_(anode,junit,jstat)
  implicit none
  class(lcbasnode),intent(in):: anode
  integer(i_kind),intent(in   ):: junit
  integer(i_kind),intent(  out):: jstat

  character(len=*),parameter:: myname_=myname//'::obsnode_xwrite_'
_ENTRY_(myname_)

  write(junit,iostat=jstat)     anode%res    , &
                                anode%err2   , &
                                anode%raterr2, &
                                anode%b      , &
                                anode%pg     , &
                                anode%wij    , &
                                anode%ij
  if (jstat/=0) then
     call perr(myname_,'write(%(res,...)), iostat =',jstat)
     _EXIT_(myname_)
     return
  end if
_EXIT_(myname_)
  return
end subroutine obsnode_xwrite_

subroutine obsnode_sethop_(anode)
  use m_cvgridlookup, only: cvgridlookup_getiw
  implicit none
  class(lcbasnode),intent(inout):: anode

!  character(len=*),parameter:: myname_=myname//'::obsnode_sethop_'
!_ENTRY_(myname_)
  call cvgridlookup_getiw(anode%elat,anode%elon,anode%ij,anode%wij)
!_EXIT_(myname_)
  return
end subroutine obsnode_sethop_

function obsnode_isvalid_(anode) result(isvalid_)
  implicit none
  logical:: isvalid_
  class(lcbasnode),intent(in):: anode

!  character(len=*),parameter:: myname_=myname//'::obsnode_isvalid_'
!_ENTRY_(myname_)
  isvalid_= associated(anode%diags)
!_EXIT_(myname_)
  return
end function obsnode_isvalid_

pure subroutine gettlddp_(anode,jiter,tlddp,nob)
  use kinds, only: r_kind
  implicit none
  class(lcbasnode), intent(in):: anode
  integer(kind=i_kind),intent(in):: jiter
  real(kind=r_kind),intent(inout):: tlddp
  integer(kind=i_kind),optional,intent(inout):: nob

  tlddp = tlddp + anode%diags%tldepart(jiter)*anode%diags%tldepart(jiter)
  if(present(nob)) nob=nob+1
  return
end subroutine gettlddp_

end module m_lcbasnode
