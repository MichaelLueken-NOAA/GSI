module m_swcpnode
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module m_swcpnode
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2016-05-18
!
! abstract: class-module of obs-type swcpnode (solid-water content path)
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

  public:: swcpnode

  type,extends(obsnode):: swcpnode
     !type(swcp_ob_type),pointer :: llpoint => null()
     type(obs_diag), pointer :: diags => null()
     real(r_kind)    :: res           !  solid-water content path residual
     real(r_kind)    :: err2          !  solid-water content path error squared
     real(r_kind)    :: raterr2       !  square of ratio of final obs error 
                                      !  to original obs error
     !real(r_kind)    :: time          !  observation time in sec     
     real(r_kind)    :: b             !  variational quality control parameter
     real(r_kind)    :: pg            !  variational quality control parameter
     real(r_kind)    :: wij(4)        !  horizontal interpolation weights
     real(r_kind),dimension(:),pointer :: jac_t => null()
                                      !  t jacobian 
     real(r_kind),dimension(:),pointer :: jac_p => null()
                                      !  p jacobian
     real(r_kind),dimension(:),pointer :: jac_q => null()
                                      !  q jacobian 
     real(r_kind),dimension(:),pointer :: jac_qi => null()
                                      !  qi jacobian 
     real(r_kind),dimension(:),pointer :: jac_qs => null()
                                      !  qs jacobian 
     real(r_kind),dimension(:),pointer :: jac_qg => null()
                                      !  qg jacobian 
     real(r_kind),dimension(:),pointer :: jac_qh => null()
                                      !  qh jacobian 
!     real(r_kind),dimension(:),pointer :: dp  => null()
!                                      !  delta pressure at mid layers at obs locations
     integer(i_kind),dimension(:,:),pointer :: ij  => null()
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
  end type swcpnode

  public:: swcpnode_typecast
  public:: swcpnode_nextcast
     interface swcpnode_typecast; module procedure typecast_ ; end interface
     interface swcpnode_nextcast; module procedure nextcast_ ; end interface

  public:: swcpnode_appendto
     interface swcpnode_appendto; module procedure appendto_ ; end interface

  character(len=*),parameter:: myname="m_swcpnode"

#include "myassert.H"
#include "mytrace.H"
contains
function typecast_(anode) result(ptr_)
!-- cast a class(obsnode) to a type(swcpnode)
  use m_obsnode, only: obsnode
  implicit none
  type(swcpnode),pointer:: ptr_
  class(obsnode),pointer,intent(in):: anode
  ptr_ => null()
  if(.not.associated(anode)) return
        ! logically, typecast of a null-reference is a null pointer.
  select type(anode)
     type is(swcpnode)
        ptr_ => anode
  end select
  return
end function typecast_

function nextcast_(anode) result(ptr_)
!-- cast an obsnode_next(obsnode) to a type(swcpnode)
  use m_obsnode, only: obsnode,obsnode_next
  implicit none
  type(swcpnode),pointer:: ptr_
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
  type(swcpnode),pointer,intent(in):: anode
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
  mytype="[swcpnode]"
end function mytype

subroutine obsheader_read_(iunit,mobs,jread,istat)
  use gridmod, only: nsig
  implicit none
  integer(i_kind),intent(in ):: iunit
  integer(i_kind),intent(out):: mobs
  integer(i_kind),intent(out):: jread
  integer(i_kind),intent(out):: istat

  character(len=*),parameter:: myname_=myname//".obsheader_read_"
  integer(i_kind):: msig
_ENTRY_(myname_)
  
  read(iunit,iostat=istat) mobs,jread, msig
  if(istat==0 .and. nsig/=msig) then
     call perr(myname_,'unexpected dimension information, nsig =',nsig)
     call perr(myname_,'                         but read msig =',msig)
     call  die(myname_)
  endif
_EXIT_(myname_)
  return
end subroutine obsheader_read_

subroutine obsheader_write_(junit,mobs,jwrite,jstat)
  use gridmod, only: nsig
  implicit none
  integer(i_kind),intent(in ):: junit
  integer(i_kind),intent(in ):: mobs
  integer(i_kind),intent(in ):: jwrite
  integer(i_kind),intent(out):: jstat
  
  character(len=*),parameter:: myname_=myname//".obsheader_write_"
_ENTRY_(myname_)
  write(junit,iostat=jstat) mobs,jwrite, nsig
_EXIT_(myname_)
  return
end subroutine obsheader_write_

subroutine obsnode_init_(anode)
  use gridmod, only: nsig
  implicit none
  class(swcpnode),intent(out):: anode

  character(len=*),parameter:: myname_=myname//'.obsnode_init_'
_ENTRY_(myname_)
  anode%llpoint => null()
  anode%luse = .false.
  anode%elat = 0._r_kind
  anode%elon = 0._r_kind
  anode%time = 0._r_kind
  anode%idv  =-1
  anode%iob  =-1
  allocate(anode%jac_t(nsig  ), &
           anode%jac_p(nsig+1), &
           anode%jac_q(nsig  ), &
           anode%jac_qi(nsig  ), &
           anode%jac_qs(nsig  ), &
           anode%jac_qg(nsig  ), &
           anode%jac_qh(nsig  ), &
           anode%ij(4, nsig  )  )
!  allocate(anode%dp(nsig))
_EXIT_(myname_)
  return
end subroutine obsnode_init_

subroutine obsnode_clean_(anode)
  implicit none
  class(swcpnode),intent(inout):: anode

  character(len=*),parameter:: myname_=myname//'.obsnode_clean_'
_ENTRY_(myname_)
!_TRACEV_(myname_,'%mytype() =',anode%mytype())
  if(associated(anode%jac_t )) deallocate(anode%jac_t )
  if(associated(anode%jac_p )) deallocate(anode%jac_p )
  if(associated(anode%jac_q )) deallocate(anode%jac_q )
  if(associated(anode%jac_qi)) deallocate(anode%jac_qi)
  if(associated(anode%jac_qs)) deallocate(anode%jac_qs)
  if(associated(anode%jac_qg)) deallocate(anode%jac_qg)
  if(associated(anode%jac_qh)) deallocate(anode%jac_qh)
!  if(associated(anode%dp    )) deallocate(anode%dp    )
  if(associated(anode%ij    )) deallocate(anode%ij    )
_EXIT_(myname_)
  return
end subroutine obsnode_clean_

subroutine obsnode_xread_(anode,iunit,istat,diaglookup,skip)
  use m_obsdiagnode, only: obsdiaglookup_locate
  implicit none
  class(swcpnode) , intent(inout):: anode
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
                                 anode%wij    , & !(4)
                                 anode%jac_t  , & !(  nsig)
                                 anode%jac_p  , & !(  nsig)
                                 anode%jac_q  , & !(  nsig)
                                 anode%jac_qi , & !(  nsig)
                                 anode%jac_qs , & !(  nsig)
                                 anode%jac_qg , & !(  nsig)
                                 anode%jac_qh , & !(4,nsig)
                                 anode%ij
!                                 anode%dp
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
  class(swcpnode),intent(in):: anode
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
                                anode%wij    , &
                                anode%jac_t  , &
                                anode%jac_p  , &
                                anode%jac_q  , &
                                anode%jac_qi , &
                                anode%jac_qs , &
                                anode%jac_qg , &
                                anode%jac_qh , &
                                anode%ij
!                                anode%dp
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
  use gridmod, only: nsig,latlon11
  implicit none
  class(swcpnode),intent(inout):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_sethop_'
  integer(i_kind):: k
_ENTRY_(myname_)

  ASSERT(size(anode%ij,2)==nsig)
  ASSERT(nsig>0)

  call cvgridlookup_getiw(anode%elat,anode%elon,anode%ij(:,1),anode%wij)
  do k=2,nsig
     anode%ij(:,k) = anode%ij(:,1)+(k-1)*latlon11
  enddo
_EXIT_(myname_)
  return
end subroutine obsnode_sethop_

function obsnode_isvalid_(anode) result(isvalid_)
  implicit none
  logical:: isvalid_
  class(swcpnode),intent(in):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_isvalid_'
_ENTRY_(myname_)
  isvalid_=associated(anode%diags)
_EXIT_(myname_)
  return
end function obsnode_isvalid_

pure subroutine gettlddp_(anode,jiter,tlddp,nob)
  use kinds, only: r_kind
  implicit none
  class(swcpnode), intent(in):: anode
  integer(kind=i_kind),intent(in):: jiter
  real(kind=r_kind),intent(inout):: tlddp
  integer(kind=i_kind),optional,intent(inout):: nob

  tlddp = tlddp + anode%diags%tldepart(jiter)*anode%diags%tldepart(jiter)
  if(present(nob)) nob=nob+1
  return
end subroutine gettlddp_

end module m_swcpnode
