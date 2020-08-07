module m_qnode
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module m_qnode
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2016-05-18
!
! abstract: class-module of obs-type qnode (moisture)
!
! program history log:
!   2016-05-18  j guo   - added this document block for the initial polymorphic
!                         implementation.
!   2019-09-20  X.Su    - add new variational QC parameters
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

  public:: qnode

  type,extends(obsnode):: qnode
     !type(q_ob_type),pointer :: llpoint => null()
     type(obs_diag), pointer :: diags => null() 
     real(r_kind)    :: res           !  moisture residual
     real(r_kind)    :: err2          !  moisture error squared
     real(r_kind)    :: raterr2       !  square of ratio of final obs error 
                                      !  to original obs error
     !real(r_kind)    :: time          !  observation time in sec     
     real(r_kind)    :: b             !  variational quality control parameter
     real(r_kind)    :: pg            !  variational quality control parameter
     real(r_kind)    :: jb            !  variational quality control parameter
     integer(i_kind) :: ib            !  new variational quality control parameter
     integer(i_kind) :: ik            !  new variational quality control parameter
     real(r_kind)    :: wij(8)        !  horizontal interpolation weights
     real(r_kind)    :: qpertb        !  random number adding to the obs
     integer(i_kind) :: ij(8)         !  horizontal locations
     integer(i_kind) :: k1            !  level of errtable 1-33
     integer(i_kind) :: kx            !  ob type
     !logical         :: luse          !  flag indicating if ob is used in pen.

     !integer(i_kind) :: idv,iob	      ! device id and obs index for sorting
     !real   (r_kind) :: elat, elon      ! earth lat-lon for redistribution
     !real   (r_kind) :: dlat, dlon      ! earth lat-lon for redistribution
     real   (r_kind) :: dlev            ! reference to the vertical grid

     integer(i_kind) :: ich0=0  ! ich code to mark derived data.  See
                                ! qnode_ich0 and qnode_ich0_pbl_Pseudo below
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
  end type qnode

  public:: qnode_typecast
  public:: qnode_nextcast
     interface qnode_typecast; module procedure typecast_ ; end interface
     interface qnode_nextcast; module procedure nextcast_ ; end interface

  public:: qnode_appendto
     interface qnode_appendto; module procedure appendto_ ; end interface

  ! Because there are two components in qnode for an ordinary wind obs,
  ! ich values are set to (1,2).  Therefore, ich values for pbl_pseudo_surfobsuv
  ! are set to (3,4), and qnode_ich0_pbl_pseudo is set to 2.

  public:: qnode_ich0
  public:: qnode_ich0_pbl_pseudo

  integer(i_kind),parameter :: qnode_ich0            = 0             ! ich=0+1
  integer(i_kind),parameter :: qnode_ich0_pbl_pseudo = qnode_ich0+1  ! ich=1+1

  character(len=*),parameter:: myname="m_qnode"

#include "myassert.H"
#include "mytrace.H"
contains
function typecast_(anode) result(ptr_)
!-- cast a class(obsnode) to a type(qnode)
  use m_obsnode, only: obsnode
  implicit none
  type(qnode   ),pointer:: ptr_
  class(obsnode),pointer,intent(in):: anode
  ptr_ => null()
  if(.not.associated(anode)) return
        ! logically, typecast of a null-reference is a null pointer.
  select type(anode)
     type is(qnode)
        ptr_ => anode
  end select
  return
end function typecast_

function nextcast_(anode) result(ptr_)
!-- cast an obsnode_next(obsnode) to a type(qnode)
  use m_obsnode, only: obsnode,obsnode_next
  implicit none
  type(qnode   ),pointer:: ptr_
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
  type(qnode),pointer,intent(in):: anode
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
  mytype="[qnode]"
end function mytype

subroutine obsnode_xread_(anode,iunit,istat,diaglookup,skip)
  use m_obsdiagnode, only: obsdiaglookup_locate
  implicit none
  class(qnode) , intent(inout):: anode
  integer(i_kind) , intent(in   ):: iunit
  integer(i_kind) , intent(  out):: istat
  type(obs_diags) , intent(in   ):: diaglookup
  logical,optional, intent(in   ):: skip

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
                                 anode%jb     , &
                                 anode%ib     , &
                                 anode%ik     , &
                                 anode%qpertb , &
                                 anode%k1     , &
                                 anode%kx     , &
                                 anode%dlev   , &
                                 anode%ich0   , &
                                 anode%wij    , &
                                 anode%ij
     if (istat/=0) then
        call perr(myname_,'read(%(res,err2,...)), iostat =',istat)
        _EXIT_(myname_)
        return
     end if

     anode%diags => obsdiaglookup_locate(diaglookup,anode%idv,anode%iob,anode%ich0+1)
     if(.not.associated(anode%diags)) then
        call perr(myname_,'obsdiaglookup_locate(), %idv =',anode%idv)
        call perr(myname_,'                        %iob =',anode%iob)
        call perr(myname_,'                       %ich0 =',anode%ich0)
        call  die(myname_)
     endif
  endif
_EXIT_(myname_)
  return
end subroutine obsnode_xread_

subroutine obsnode_xwrite_(anode,junit,jstat)
  implicit none
  class(qnode),intent(in):: anode
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
                                anode%jb     , &
                                anode%ib     , &
                                anode%ik     , &
                                anode%qpertb , &
                                anode%k1     , &
                                anode%kx     , &
                                anode%dlev   , &
                                anode%ich0   , &
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
  class(qnode),intent(inout):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_sethop_'
_ENTRY_(myname_)
  call cvgridlookup_getiw(anode%elat,anode%elon,anode%dlev,anode%ij,anode%wij)
_EXIT_(myname_)
  return
end subroutine obsnode_sethop_

function obsnode_isvalid_(anode) result(isvalid_)
  implicit none
  logical:: isvalid_
  class(qnode),intent(in):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_isvalid_'
_ENTRY_(myname_)
  isvalid_=associated(anode%diags)
_EXIT_(myname_)
  return
end function obsnode_isvalid_

pure subroutine gettlddp_(anode,jiter,tlddp,nob)
  use kinds, only: r_kind
  implicit none
  class(qnode), intent(in):: anode
  integer(kind=i_kind),intent(in):: jiter
  real(kind=r_kind),intent(inout):: tlddp
  integer(kind=i_kind),optional,intent(inout):: nob

  tlddp = tlddp + anode%diags%tldepart(jiter)*anode%diags%tldepart(jiter)
  if(present(nob)) nob=nob+1
  return
end subroutine gettlddp_

end module m_qnode
