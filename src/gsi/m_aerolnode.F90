module m_aerolnode
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module m_aerolnode
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2016-05-18
!
! abstract: class-module of obs-type aerolnode (unfinished leveled aerosol data?)
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
  use kinds , only: i_kind,r_kind
  use mpeu_util, only: assert_,die,perr,warn,tell
  use m_obsnode, only: obsnode
  implicit none
  private

  public:: aerolnode

  type,extends(obsnode):: aerolnode
     !type(aerol_ob_type),pointer :: llpoint => null()
     type(obs_diag), pointer :: diags => null()
     real(r_kind)    :: res    =0._r_kind    !  aerosol residual
     real(r_kind)    :: err2   =0._r_kind    !  aerosol obs error squared
     real(r_kind)    :: raterr2=0._r_kind    !  square of ratio of final obs error
                                             !  to original obs error
     !real(r_kind)    :: time                !  observation time
     real(r_kind)    :: b      =0._r_kind    !  variational quality control parameter
     real(r_kind)    :: pg     =0._r_kind    !  variational quality control parameter
     real(r_kind)    :: wij(8) =0._r_kind    !  horizontal interpolation weights
     integer(i_kind) :: ij(8)  =0_i_kind     !  horizontal locations
     !logical         :: luse          !  flag indicating if ob is used in pen.

     !integer(i_kind) :: idv,iob         ! device id and obs index for sorting
     !real   (r_kind) :: elat, elon      ! earth lat-lon for redistribution
     !real   (r_kind) :: dlat, dlon      ! earth lat-lon for redistribution
     real   (r_kind) :: dlev   =0._r_kind      ! reference to the vertical grid
  contains
     procedure,nopass::  mytype
     procedure::  sethop => obsnode_sethop_
     procedure::   xread => obsnode_xread_
     procedure::  xwrite => obsnode_xwrite_
     procedure:: isvalid => obsnode_isvalid_
     procedure:: gettlddp => gettlddp_

    ! procedure, nopass:: headerread  => obsheader_read_
    ! procedure, nopass:: headerwrite => obsheader_write_
    ! procedure:: init  => obsnode_init_
    ! procedure:: clean => obsnode_clean_
  end type aerolnode

  public:: aerolnode_typecast
  public:: aerolnode_nextcast
     interface aerolnode_typecast; module procedure typecast_ ; end interface
     interface aerolnode_nextcast; module procedure nextcast_ ; end interface

  public:: aerolnode_appendto
     interface aerolnode_appendto; module procedure appendto_ ; end interface

  character(len=*),parameter:: myname="m_aerolnode"

#include "myassert.H"
#include "mytrace.H"
contains
function typecast_(anode) result(ptr_)
!-- cast a class(obsnode) to a type(aerolnode)
  use m_obsnode, only: obsnode
  implicit none
  type(aerolnode),pointer:: ptr_
  class(obsnode ),pointer,intent(in):: anode
  ptr_ => null()
  if(.not.associated(anode)) return
        ! logically, typecast of a null-reference is a null pointer.
  select type(anode)
     type is(aerolnode)
        ptr_ => anode
  end select
  return
end function typecast_

function nextcast_(anode) result(ptr_)
!-- cast an obsnode_next(obsnode) to a type(aerolnode)
  use m_obsnode, only: obsnode,obsnode_next
  implicit none
  type(aerolnode),pointer:: ptr_
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
  type(aerolnode),pointer,intent(in):: anode
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
  mytype="[aerolnode]"
end function mytype

subroutine obsnode_xread_(anode,iunit,istat,diaglookup,skip)
  use m_obsdiagnode, only: obsdiaglookup_locate
  use m_obsdiagnode, only: obs_diags
  implicit none
  class(aerolnode), intent(inout):: anode
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
        call perr(myname_,'skipping read(%(res,...)), istat =',istat)
        _EXIT_(myname_)
        return
     endif
  else
     read(iunit,iostat=istat)    anode%res    , &
                                 anode%err2   , &
                                 anode%raterr2, &
                                 anode%b      , &
                                 anode%pg     , &
                                 anode%dlev   , &
                                 anode%wij    , &
                                 anode%ij
     if(istat/=0) then
        call perr(myname_,'read(%(res,...)), istat =',istat)
        _EXIT_(myname_)
        return
     endif

     anode%diags => obsdiaglookup_locate(diaglookup,anode%idv,anode%iob,1)
     if(.not.  associated(anode%diags)) then
        call perr(myname_,'obsdiaglookup_locate(), %idv =',anode%idv)
        call perr(myname_,'                        %iob =',anode%iob)
        call  die(myname_)
     endif
  endif
_EXIT_(myname_)
  return
end subroutine obsnode_xread_

subroutine obsnode_xwrite_(anode,junit,jstat)
  implicit none
  class(aerolnode),intent(in):: anode
  integer(i_kind),intent(in   ):: junit
  integer(i_kind),intent(  out):: jstat

  character(len=*),parameter:: myname_=myname//'::obsnode_xwrite_'
_ENTRY_(myname_)

  write(junit,iostat=jstat)     anode%res    , &
                                anode%err2   , &
                                anode%raterr2, &
                                anode%b      , &
                                anode%pg     , &
                                anode%dlev   , &
                                anode%wij    , &
                                anode%ij
_EXIT_(myname_)
  return
end subroutine obsnode_xwrite_

subroutine obsnode_sethop_(anode)
  use m_cvgridlookup, only: cvgridlookup_getiw
  implicit none
  class(aerolnode),intent(inout):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_sethop_'
_ENTRY_(myname_)
  call cvgridlookup_getiw(anode%elat,anode%elon,anode%dlev,anode%ij,anode%wij)
_EXIT_(myname_)
  return
end subroutine obsnode_sethop_

function obsnode_isvalid_(anode) result(isvalid_)
  implicit none
  logical:: isvalid_
  class(aerolnode),intent(in):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_isvalid_'
_ENTRY_(myname_)
  isvalid_=associated(anode%diags)
_EXIT_(myname_)
  return
end function obsnode_isvalid_

pure subroutine gettlddp_(anode,jiter,tlddp,nob)
  use kinds, only: r_kind
  implicit none
  class(aerolnode), intent(in):: anode
  integer(kind=i_kind),intent(in):: jiter
  real(kind=r_kind),intent(inout):: tlddp
  integer(kind=i_kind),optional,intent(inout):: nob

  tlddp = tlddp + anode%diags%tldepart(jiter)*anode%diags%tldepart(jiter)
  if(present(nob)) nob=nob+1
return
end subroutine gettlddp_

end module m_aerolnode
