module m_pcpnode
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module m_pcpnode
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2016-05-18
!
! abstract: class-module of obs-type pcpnode (precipitation)
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

  public:: pcpnode

  type,extends(obsnode):: pcpnode
     !type(pcp_ob_type),pointer :: llpoint => null()
     type(obs_diag), pointer :: diags => null()
     real(r_kind)    :: obs           !  observed precipitation value 
     real(r_kind)    :: err2          !  error variances squared
     real(r_kind)    :: raterr2       !  ratio of error variances squared 
     !real(r_kind)    :: time          !  observation time in sec     
     real(r_kind)    :: ges           !  guess observation value
     real(r_kind)    :: wij(4)        !  horizontal interpolation weights
     real(r_kind),dimension(:),pointer :: predp => null()
                                      !  predictors (npredp)
     real(r_kind),dimension(:),pointer :: dpcp_dvar => null()
                                      !  error variances squared (nsig5)
     integer(i_kind) :: ij(4)         !  horizontal locations
     integer(i_kind) :: icxp          !  type of precipitation rate observation
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
 
     procedure, nopass:: headerread  => obsheader_read_
     procedure, nopass:: headerwrite => obsheader_write_
     procedure:: init  => obsnode_init_
     procedure:: clean => obsnode_clean_
  end type pcpnode

  public:: pcpnode_typecast
  public:: pcpnode_nextcast
     interface pcpnode_typecast; module procedure typecast_ ; end interface
     interface pcpnode_nextcast; module procedure nextcast_ ; end interface

  public:: pcpnode_appendto
     interface pcpnode_appendto; module procedure appendto_ ; end interface

  character(len=*),parameter:: myname="m_pcpnode"

#include "myassert.H"
#include "mytrace.H"
contains
function typecast_(anode) result(ptr_)
!-- cast a class(obsnode) to a type(pcpnode)
  use m_obsnode, only: obsnode
  implicit none
  type(pcpnode ),pointer:: ptr_
  class(obsnode),pointer,intent(in):: anode
  ptr_ => null()
  if(.not.associated(anode)) return
        ! logically, typecast of a null-reference is a null pointer.
  select type(anode)
     type is(pcpnode)
        ptr_ => anode
  end select
  return
end function typecast_

function nextcast_(anode) result(ptr_)
!-- cast an obsnode_next(obsnode) to a type(pcpnode)
  use m_obsnode, only: obsnode,obsnode_next
  implicit none
  type(pcpnode ),pointer:: ptr_
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
  type(pcpnode),pointer,intent(in):: anode
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
  mytype="[pcpnode]"
end function mytype

subroutine obsheader_read_(iunit,mobs,jread,istat)
  use gridmod, only: nsig5
  use pcpinfo, only: npredp
  implicit none
  integer(i_kind),intent(in ):: iunit
  integer(i_kind),intent(out):: mobs
  integer(i_kind),intent(out):: jread
  integer(i_kind),intent(out):: istat

  character(len=*),parameter:: myname_=myname//'.obsheader_read_'
  integer(i_kind):: mpredp,msig5
_ENTRY_(myname_)
  
  read(iunit,iostat=istat) mobs,jread, mpredp,msig5
  if(istat==0 .and. (npredp/=mpredp .or. nsig5/=msig5)) then
     call perr(myname_,'unmatched dimension information, npredp or nsig5')
     if(npredp/=mpredp) then
        call perr(myname_,'  expecting npredp =',npredp)
        call perr(myname_,'   but read mpredp =',mpredp)
     endif
     if(nsig5/=msig5) then
        call perr(myname_,'   expecting nsig5 =',nsig5)
        call perr(myname_,'    but read msig5 =',msig5)
     endif
     call die(myname_)
  endif
_EXIT_(myname_)
  return
end subroutine obsheader_read_

subroutine obsheader_write_(junit,mobs,jwrite,jstat)
  use gridmod, only: nsig5
  use pcpinfo, only: npredp
  implicit none
  integer(i_kind),intent(in ):: junit
  integer(i_kind),intent(in ):: mobs
  integer(i_kind),intent(in ):: jwrite
  integer(i_kind),intent(out):: jstat

  character(len=*),parameter:: myname_=myname//'.obsheader_write_'
_ENTRY_(myname_)
  
  write(junit,iostat=jstat) mobs,jwrite, npredp,nsig5
_EXIT_(myname_)
  return
end subroutine obsheader_write_

subroutine obsnode_init_(anode)
  use gridmod, only: nsig5
  use pcpinfo, only: npredp
  implicit none
  class(pcpnode),intent(out):: anode

  character(len=*),parameter:: myname_=myname//'.obsnode_init_'
_ENTRY_(myname_)
  !anode = _obsnode_()
  anode%llpoint => null()
  anode%luse = .false.
  anode%time = 0._r_kind
  anode%elat = 0._r_kind
  anode%elon = 0._r_kind
  anode%idv  =-1
  anode%iob  =-1
  !-anode%dlev = 0._r_kind
  !-anode%ich  =-1

  allocate(anode%predp(npredp), &
           anode%dpcp_dvar(1:nsig5) )
_EXIT_(myname_)
  return
end subroutine obsnode_init_

subroutine obsnode_clean_(anode)
  implicit none
  class(pcpnode),intent(inout):: anode

  character(len=*),parameter:: myname_=myname//'.obsnode_clean_'
_ENTRY_(myname_)
!_TRACEV_(myname_,'%mytype() =',anode%mytype())
  if(associated(anode%predp    )) deallocate(anode%predp)
  if(associated(anode%dpcp_dvar)) deallocate(anode%dpcp_dvar)
_EXIT_(myname_)
  return
end subroutine obsnode_clean_

subroutine obsnode_xread_(anode,iunit,istat,diaglookup,skip)
  use m_obsdiagnode, only: obsdiaglookup_locate
  implicit none
  class(pcpnode) , intent(inout):: anode
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
     read(iunit,iostat=istat)    anode%obs    , &
                                 anode%err2   , &
                                 anode%raterr2, &
                                 anode%ges    , &
                                 anode%icxp   , &
                                 anode%predp(:)    , &
                                 anode%dpcp_dvar(:), &
                                 anode%wij(:) , &
                                 anode%ij(:)
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
  class(pcpnode),intent(in):: anode
  integer(i_kind),intent(in   ):: junit
  integer(i_kind),intent(  out):: jstat

  character(len=*),parameter:: myname_=myname//'.obsnode_xwrite_'
_ENTRY_(myname_)

  jstat=0
  write(junit,iostat=jstat)     anode%obs    , &
                                anode%err2   , &
                                anode%raterr2, &
                                anode%ges    , &
                                anode%icxp   , &
                                anode%predp(:)    , &
                                anode%dpcp_dvar(:), &
                                anode%wij(:) , &
                                anode%ij(:)
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
  class(pcpnode),intent(inout):: anode

!  character(len=*),parameter:: myname_=myname//'::obsnode_sethop_'
!_ENTRY_(myname_)
  call cvgridlookup_getiw(anode%elat,anode%elon,anode%ij,anode%wij)
!_EXIT_(myname_)
  return
end subroutine obsnode_sethop_

function obsnode_isvalid_(anode) result(isvalid_)
  implicit none
  logical:: isvalid_
  class(pcpnode),intent(in):: anode

!  character(len=*),parameter:: myname_=myname//'::obsnode_isvalid_'
!_ENTRY_(myname_)
  isvalid_=associated(anode%diags)
!_EXIT_(myname_)
  return
end function obsnode_isvalid_

pure subroutine gettlddp_(anode,jiter,tlddp,nob)
  use kinds, only: r_kind
  implicit none
  class(pcpnode), intent(in):: anode
  integer(kind=i_kind),intent(in):: jiter
  real(kind=r_kind),intent(inout):: tlddp
  integer(kind=i_kind),optional,intent(inout):: nob

  tlddp = tlddp + anode%diags%tldepart(jiter)*anode%diags%tldepart(jiter)
  if(present(nob)) nob=nob+1
  return
end subroutine gettlddp_

end module m_pcpnode
