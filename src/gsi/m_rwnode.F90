module m_rwnode
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module m_rwnode
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2016-05-18
!
! abstract: class-module of obs-type rwnode (radar radial winds)
!
! program history log:
!   2016-05-18  j guo   - added this document block for the initial polymorphic
!                         implementation.
!   2016-06-23  lippi   - add costilt and sintilt for radial wind calculations. Also,
!                         add cosazm_costilt and sinazm_costilt as u, v
!                         factors respectively.
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

  public:: rwnode

  type,extends(obsnode):: rwnode
     !type(rw_ob_type),pointer :: llpoint => null()
     type(obs_diag), pointer :: diags => null()     
     real(r_kind)    :: res           !  radial wind residual
     real(r_kind)    :: err2          !  radial wind error squared
     real(r_kind)    :: raterr2       !  square of ratio of final obs error 
                                      !  to original obs error
     !real(r_kind)    :: time          !  observation time in sec     
     real(r_kind)    :: b             !  variational quality control parameter
     real(r_kind)    :: pg            !  variational quality control parameter
     real(r_kind)    :: cosazm        !  v factor
     real(r_kind)    :: sinazm        !  u factor
     real(r_kind)    :: costilt       !  u,v factor
     real(r_kind)    :: sintilt       !  w factor
     real(r_kind)    :: cosazm_costilt!  u factor
     real(r_kind)    :: sinazm_costilt!  v factor
     real(r_kind)    :: wij(8)        !  horizontal interpolation weights
     integer(i_kind) :: ij(8)         !  horizontal locations
     !logical         :: luse          !  flag indicating if ob is used in pen.

     !integer(i_kind) :: idv,iob	      ! device id and obs index for sorting
     !real   (r_kind) :: elat, elon      ! earth lat-lon for redistribution
     !real   (r_kind) :: dlat, dlon      ! earth lat-lon for redistribution
     real   (r_kind) :: dlev            ! reference to the vertical grid
     real   (r_kind) :: factw           ! factor of 10m wind
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
  end type rwnode

  public:: rwnode_typecast
  public:: rwnode_nextcast
     interface rwnode_typecast; module procedure typecast_ ; end interface
     interface rwnode_nextcast; module procedure nextcast_ ; end interface

  public:: rwnode_appendto
     interface rwnode_appendto; module procedure appendto_ ; end interface

  character(len=*),parameter:: myname="m_rwnode"

#include "myassert.H"
#include "mytrace.H"
contains
function typecast_(anode) result(ptr_)
!-- cast a class(obsnode) to a type(rwnode)
  use m_obsnode, only: obsnode
  implicit none
  type(rwnode  ),pointer:: ptr_
  class(obsnode),pointer,intent(in):: anode
  ptr_ => null()
  if(.not.associated(anode)) return
        ! logically, typecast of a null-reference is a null pointer.
  select type(anode)
     type is(rwnode)
        ptr_ => anode
  end select
  return
end function typecast_

function nextcast_(anode) result(ptr_)
!-- cast an obsnode_next(obsnode) to a type(rwnode)
  use m_obsnode, only: obsnode,obsnode_next
  implicit none
  type(rwnode  ),pointer:: ptr_
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
  type(rwnode),pointer,intent(in):: anode
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
  mytype="[rwnode]"
end function mytype

subroutine obsnode_xread_(anode,iunit,istat,diaglookup,skip)
  use m_obsdiagnode, only: obsdiaglookup_locate
  implicit none
  class(rwnode),intent(inout):: anode
  integer(i_kind),intent(in   ):: iunit
  integer(i_kind),intent(  out):: istat
  type(obs_diags),intent(in   ):: diaglookup
  logical,optional,intent(in   ):: skip

  character(len=*),parameter:: myname_=myname//'.obsnode_xread_'
  logical:: skip_
_ENTRY_(myname_)
  skip_=.false.
  if(present(skip)) skip_=skip

  istat=0
  if(skip_) then
     read(iunit,iostat=istat)
     if(istat/=0) then
        call perr(myname_,'skipping read(%(res,err2,...)), iostat =',istat)
        _EXIT_(myname_)
        return
     endif

  else
     read(iunit,iostat=istat)    anode%res    , &
                                 anode%err2   , &
                                 anode%raterr2, &
                                 anode%b      , &
                                 anode%pg     , &
                                 anode%cosazm , &
                                 anode%sinazm , &
                                 anode%sintilt , &
                                 anode%costilt , &
                                 anode%cosazm_costilt , &
                                 anode%sinazm_costilt , &
                                 anode%dlev   , &
                                 anode%factw  , &
                                 anode%wij    , &
                                 anode%ij
     if (istat/=0) then
        call perr(myname_,'read(%(res,err2,...)), iostat =',istat)
        _EXIT_(myname_)
        return
     end if

     anode%diags => obsdiaglookup_locate(diaglookup,anode%idv,anode%iob,1)
     if(.not.associated(anode%diags)) then
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
  class(rwnode),intent(in):: anode
  integer(i_kind),intent(in   ):: junit
  integer(i_kind),intent(  out):: jstat

  character(len=*),parameter:: myname_=myname//'.obsnode_xwrite_'
_ENTRY_(myname_)

  jstat=0
  write(junit,iostat=jstat)     anode%res    , &
                                anode%err2   , &
                                anode%raterr2, &
                                anode%b      , &
                                anode%pg     , &
                                anode%cosazm , &
                                anode%sinazm , &
                                anode%sintilt , &
                                anode%costilt , &
                                anode%cosazm_costilt , &
                                anode%sinazm_costilt , &
                                anode%dlev   , &
                                anode%factw  , &
                                anode%wij    , &
                                anode%ij
  if (jstat/=0) then
     call perr(myname_,'write(%(res,err2,...)), iostat =',jstat)
     _EXIT_(myname_)
     return
  end if
_EXIT_(myname_)
  return
end subroutine obsnode_xwrite_

subroutine obsnode_sethop_(anode)
  use m_cvgridlookup, only: cvgridlookup_getiw
  implicit none
  class(rwnode),intent(inout):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_sethop_'
_ENTRY_(myname_)
  call cvgridlookup_getiw(anode%elat,anode%elon,anode%dlev,anode%ij,anode%wij)
  anode%wij(1:8) = anode%wij(1:8)*anode%factw
_EXIT_(myname_)
  return
end subroutine obsnode_sethop_

function obsnode_isvalid_(anode) result(isvalid_)
  implicit none
  logical:: isvalid_
  class(rwnode),intent(in):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_isvalid_'
_ENTRY_(myname_)
  isvalid_=associated(anode%diags)
_EXIT_(myname_)
  return
end function obsnode_isvalid_

pure subroutine gettlddp_(anode,jiter,tlddp,nob)
  use kinds, only: r_kind
  implicit none
  class(rwnode), intent(in):: anode
  integer(kind=i_kind),intent(in):: jiter
  real(kind=r_kind),intent(inout):: tlddp
  integer(kind=i_kind),optional,intent(inout):: nob

  tlddp = tlddp + anode%diags%tldepart(jiter)*anode%diags%tldepart(jiter)
  if(present(nob)) nob=nob+1
  return
end subroutine gettlddp_

end module m_rwnode
