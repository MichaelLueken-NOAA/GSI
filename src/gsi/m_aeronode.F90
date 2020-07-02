module m_aeronode
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module m_aeronode
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2016-05-18
!
! abstract: class-module of obs-type aeronode (aerosal aod)
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
  use m_obsdiagnode, only: obs_diag,aofp_obs_diag => fptr_obsdiagnode
  use m_obsdiagnode, only: obs_diags

  use kinds , only: i_kind,r_kind
  use mpeu_util, only: assert_,die,perr,warn,tell
  use m_obsnode, only: obsnode
  implicit none
  private

  public:: aeronode

  type,extends(obsnode):: aeronode
     !type(aero_ob_type),pointer :: llpoint => null()
     type(aofp_obs_diag), dimension(:), pointer :: diags => null()
     real(r_kind),dimension(:),pointer    :: res  => null()    !  aerosol property residual
     real(r_kind),dimension(:),pointer    :: err2 => null()    !  aerosol property error squared
     real(r_kind),dimension(:),pointer    :: raterr2 => null() !  square of ratio of final obs error
                                                               !  to original obs error
     !real(r_kind)                         :: time             !  observation time in sec
     real(r_kind)    :: wij(4)  =0._r_kind                     !  horizontal interpolation weights
     real(r_kind),dimension(:,:),pointer :: daod_dvar => null()! jacobians_aero (nsig*n_aerosols,nchan)
     real(r_kind),dimension(:),pointer    :: prs => null()     !  pressure levels
     integer(i_kind),dimension(:),pointer :: ipos  => null()
     integer(i_kind),dimension(:),pointer :: icx  => null()
     integer(i_kind) :: ij(4)   =0                             !  horizontal locations
     integer(i_kind) :: nlaero  =0                             !  number of channels
     !logical         :: luse                                   !  flag indicating if ob is used in pen.
     !integer(i_kind) :: idv,iob                                !  device id and obs index for sorting
     !real   (r_kind) :: elat, elon      ! earth lat-lon for redistribution
     !real   (r_kind) :: dlat, dlon      ! earth lat-lon for redistribution
     integer(i_kind),dimension(:),pointer :: ich => null()
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
     procedure:: clean => obsnode_clean_
  end type aeronode

  public:: aeronode_typecast
  public:: aeronode_nextcast
     interface aeronode_typecast; module procedure typecast_ ; end interface
     interface aeronode_nextcast; module procedure nextcast_ ; end interface

  public:: aeronode_appendto
     interface aeronode_appendto; module procedure appendto_ ; end interface

  character(len=*),parameter:: myname="m_aeronode"

#include "myassert.H"
#include "mytrace.H"
contains
function typecast_(anode) result(ptr_)
!-- cast a class(obsnode) to a type(aeronode)
  use m_obsnode, only: obsnode
  implicit none
  type(aeronode),pointer:: ptr_
  class(obsnode),pointer,intent(in):: anode
  ptr_ => null()
  if(.not.associated(anode)) return
        ! logically, typecast of a null-reference is a null pointer.
  select type(anode)
     type is(aeronode)
        ptr_ => anode
  end select
  return
end function typecast_

function nextcast_(anode) result(ptr_)
!-- cast an obsnode_next(obsnode) to a type(aeronode)
  use m_obsnode, only: obsnode,obsnode_next
  implicit none
  type(aeronode),pointer:: ptr_
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
  type(aeronode),pointer,intent(in):: anode
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
  mytype="[aeronode]"
end function mytype

subroutine obsnode_clean_(anode)
  implicit none
  class(aeronode),intent(inout):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_clean_'
_ENTRY_(myname_)
!_TRACEV_(myname_,'%mytype() =',anode%mytype())
  if(associated(anode%daod_dvar)) deallocate(anode%daod_dvar)
  if(associated(anode%diags  )) deallocate(anode%diags  )
  if(associated(anode%res    )) deallocate(anode%res    )
  if(associated(anode%err2   )) deallocate(anode%err2   )
  if(associated(anode%raterr2)) deallocate(anode%raterr2)
  if(associated(anode%prs    )) deallocate(anode%prs    )
  if(associated(anode%ipos   )) deallocate(anode%ipos   )
  if(associated(anode%icx    )) deallocate(anode%icx    )
_EXIT_(myname_)
  return
end subroutine obsnode_clean_

subroutine obsnode_xread_(anode,iunit,istat,diaglookup,skip)
  use aeroinfo, only: nsigaerojac
  use m_obsdiagnode, only: obsdiaglookup_locate
  implicit none
  class(aeronode), intent(inout):: anode
  integer(i_kind), intent(in   ):: iunit
  integer(i_kind), intent(  out):: istat
  type(obs_diags), intent(in   ):: diaglookup
  logical,optional,intent(in   ):: skip

  character(len=*),parameter:: myname_=myname//'::obsnode_xread_'
  integer(i_kind),allocatable,dimension(:):: ich
  integer(i_kind):: nlaero,nlevp,k
  logical:: skip_
_ENTRY_(myname_)

  skip_=.false.
  if(present(skip)) skip_=skip

  if(skip_) then
     read(iunit,iostat=istat)
     if(istat/=0) then
        call perr(myname_,'skipping read(%nlaero), iostat =',istat)
        _EXIT_(myname_)
        return
     endif

     read(iunit,iostat=istat)
     if(istat/=0) then
        call perr(myname_,'skipping read(%(ich,...)), iostat =',istat)
        _EXIT_(myname_)
        return
     endif

  else
     read(iunit,iostat=istat)    anode%nlaero
     if (istat/=0) then
        call perr(myname_,'read(%nlaero), iostat =',istat)
        _EXIT_(myname_)
        return
     end if

     if(associated(anode%daod_dvar)) deallocate(anode%daod_dvar)
     if(associated(anode%diags  )) deallocate(anode%diags  )
     if(associated(anode%res    )) deallocate(anode%res    )
     if(associated(anode%err2   )) deallocate(anode%err2   )
     if(associated(anode%raterr2)) deallocate(anode%raterr2)
     if(associated(anode%prs    )) deallocate(anode%prs    )
     if(associated(anode%ipos   )) deallocate(anode%ipos   )
     if(associated(anode%icx    )) deallocate(anode%icx    )

     nlaero=anode%nlaero
     nlevp = max(nlaero,1)

     allocate( anode%daod_dvar(nsigaerojac,nlaero), &
               anode%diags(nlaero), &
               anode%res  (nlaero), &
               anode%err2 (nlaero), &
               anode%raterr2(nlaero), &
               anode%prs  (nlevp ), &
               anode%ipos (nlaero), &
               anode%icx  (nlaero)  )

     allocate(       ich  (nlaero)  )

     read(iunit,iostat=istat)          ich    , &
                                 anode%res    , &
                                 anode%err2   , &
                                 anode%raterr2, &
                                 anode%prs    , &
                                 anode%ipos   , &
                                 anode%icx    , &
                                 anode%daod_dvar, &
                                 anode%wij    , &
                                 anode%ij
     if (istat/=0) then
        call perr(myname_,'read(ich,%(...)), iostat =',istat)
        _EXIT_(myname_)
        return
     end if

     do k=1,anode%nlaero
        anode%diags(k)%ptr => obsdiaglookup_locate(diaglookup,anode%idv,anode%iob,ich(k))
        if(.not. associated(anode%diags(k)%ptr)) then
           call perr(myname_,'obsdiaglookup_locate(k), k =',k)
           call perr(myname_,'                      %idv =',anode%idv)
           call perr(myname_,'                      %iob =',anode%iob)
           call perr(myname_,'                       ich =',ich(k))
           call  die(myname_)
        endif
     enddo
     deallocate(ich)
    
  endif ! if(skip_); else
_EXIT_(myname_)
  return
end subroutine obsnode_xread_

subroutine obsnode_xwrite_(anode,junit,jstat)
  implicit none
  class(aeronode),intent(in):: anode
  integer(i_kind),intent(in   ):: junit
  integer(i_kind),intent(  out):: jstat

  character(len=*),parameter:: myname_=myname//'::obsnode_xwrite_'
  integer(i_kind):: k
_ENTRY_(myname_)

  write(junit,iostat=jstat) anode%nlaero
  if(jstat/=0) then
     call perr(myname_,'write(%nlaero), iostat =',jstat)
     _EXIT_(myname_)
     return
  endif

  write(junit,iostat=jstat)  (/(anode%diags(k)%ptr%ich,k=1,anode%nlaero)/), & !(nlaero)
                                anode%res    , & !(nlaero)
                                anode%err2   , & !(nlaero)
                                anode%raterr2, & !(nlaero)
                                anode%prs    , & !(nlevp )
                                anode%ipos   , & !(nlaero)
                                anode%icx    , & !(nlaero)
                                anode%daod_dvar, & !(nsigaerojac,nlaero)
                                anode%wij    , &
                                anode%ij
  if (jstat/=0) then
     call perr(myname_,'write(ich,%(...)), iostat =',jstat)
     _EXIT_(myname_)
     return
  end if
_EXIT_(myname_)
  return
end subroutine obsnode_xwrite_

subroutine obsnode_sethop_(anode)
  use m_cvgridlookup, only: cvgridlookup_getiw
  implicit none
  class(aeronode),intent(inout):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_sethop_'
_ENTRY_(myname_)
  call cvgridlookup_getiw(anode%elat,anode%elon,anode%ij,anode%wij)
_EXIT_(myname_)
  return
end subroutine obsnode_sethop_

function obsnode_isvalid_(anode) result(isvalid_)
  implicit none
  logical:: isvalid_
  class(aeronode),intent(in):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_isvalid_'
  integer(i_kind):: k
_ENTRY_(myname_)
  isvalid_= all ( (/ (associated(anode%diags(k)%ptr), k=1,anode%nlaero) /) )
_EXIT_(myname_)
  return
end function obsnode_isvalid_

pure subroutine gettlddp_(anode,jiter,tlddp,nob)
  use kinds, only: r_kind
  implicit none
  class(aeronode), intent(in):: anode
  integer(kind=i_kind),intent(in):: jiter
  real(kind=r_kind),intent(inout):: tlddp
  integer(kind=i_kind),optional,intent(inout):: nob

  integer(kind=i_kind):: k
  do k=1,anode%nlaero
     tlddp = tlddp + anode%diags(k)%ptr%tldepart(jiter)*anode%diags(k)%ptr%tldepart(jiter)
  enddo
  if(present(nob)) nob=nob+anode%nlaero
  return
end subroutine gettlddp_

end module m_aeronode
