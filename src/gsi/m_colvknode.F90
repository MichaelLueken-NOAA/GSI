module m_colvknode
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module m_colvknode
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2016-05-18
!
! abstract: class-module of obs-type colvknode (sbuv co)
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

  public:: colvknode

  type,extends(obsnode):: colvknode
     !type(colvk_ob_type),pointer :: llpoint => null()
     type(aofp_obs_diag), dimension(:), pointer :: diags => null()
     real(r_kind),dimension(:),pointer :: res => null()
                                      !  co residual
     real(r_kind),dimension(:),pointer :: err2 => null()
                                      !  co error squared
     real(r_kind),dimension(:),pointer :: raterr2 => null()
                                      !  square of ratio of final obs error
                                      !  to original obs error
     !real(r_kind)    :: time          !  observation time in sec
     real(r_kind),dimension(  :),pointer :: wk  => null()
     real(r_kind),dimension(:,:),pointer :: wij => null()
                                      !  horizontal interpolation weights
     real(r_kind),dimension(:),pointer :: prs => null()
                                      !  pressure levels
     integer(i_kind),dimension(:),pointer :: ipos  => null()
     real(r_kind),dimension(:,:),pointer :: ak  => null()   
                                      ! mopitt vertical averaging kernel
     real(r_kind),dimension(:),pointer :: ap  => null()   
                                      ! mopitt a priori
     !real(r_kind),dimension(:),pointer   :: wkk1 => null()
     !real(r_kind),dimension(:),pointer   :: wkk2 => null()
                                      ! vertical intropolation weights for mopitt

     integer(i_kind) :: nlco =0       ! number of levels for this profile
     integer(i_kind) :: ij(4)=0       !  horizontal locations
     !logical         :: luse          !  flag indicating if ob is used in pen.

     !integer(i_kind) :: idv,iob         ! device id and obs index for sorting
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
    procedure:: clean => obsnode_clean_

  end type colvknode

  public:: colvknode_typecast
  public:: colvknode_nextcast
     interface colvknode_typecast; module procedure typecast_ ; end interface
     interface colvknode_nextcast; module procedure nextcast_ ; end interface

  public:: colvknode_appendto
     interface colvknode_appendto; module procedure appendto_ ; end interface

  character(len=*),parameter:: myname="m_colvknode"

#include "myassert.H"
#include "mytrace.H"
contains
function typecast_(anode) result(ptr_)
!-- cast a class(obsnode) to a type(colvknode)
  use m_obsnode, only: obsnode
  implicit none
  type(colvknode),pointer:: ptr_
  class(obsnode ),pointer,intent(in):: anode
  ptr_ => null()
  if(.not.associated(anode)) return
        ! logically, typecast of a null-reference is a null pointer.
  select type(anode)
     type is(colvknode)
        ptr_ => anode
  end select
  return
end function typecast_

function nextcast_(anode) result(ptr_)
!-- cast an obsnode_next(obsnode) to a type(colvknode)
  use m_obsnode, only: obsnode,obsnode_next
  implicit none
  type(colvknode),pointer:: ptr_
  class(obsnode ),target ,intent(in):: anode

  class(obsnode),pointer:: inode_
  inode_ => obsnode_next(anode)
  ptr_ => typecast_(inode_)
  return
end function nextcast_

subroutine appendto_(anode,oll)
  use m_obsnode , only: obsnode
  use m_obsllist, only: obsllist,obsllist_appendnode
  implicit none
  type(colvknode),pointer,intent(in):: anode
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
  mytype="[colvknode]"
end function mytype

subroutine obsnode_clean_(anode)
  implicit none
  class(colvknode),intent(inout):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_clean_'
_ENTRY_(myname_)
!_TRACEV_(myname_,'%mytype() =',anode%mytype())
  if(associated(anode%res    )) deallocate(anode%res    )
  if(associated(anode%diags  )) deallocate(anode%diags  )
  if(associated(anode%err2   )) deallocate(anode%err2   )
  if(associated(anode%raterr2)) deallocate(anode%raterr2)
  if(associated(anode%prs    )) deallocate(anode%prs    )
  if(associated(anode%ipos   )) deallocate(anode%ipos   )
  if(associated(anode%ak     )) deallocate(anode%ak     )
  if(associated(anode%ap     )) deallocate(anode%ap     )
  if(associated(anode%wk     )) deallocate(anode%wk     )
  if(associated(anode%wij    )) deallocate(anode%wij    )

_EXIT_(myname_)
  return
end subroutine obsnode_clean_

subroutine obsnode_xread_(anode,iunit,istat,diaglookup,skip)
  use gridmod, only: nsig
  use m_obsdiagnode, only: obsdiaglookup_locate
  implicit none
  class(colvknode) , intent(inout):: anode
  integer(i_kind) , intent(in   ):: iunit
  integer(i_kind) , intent(  out):: istat
  type(obs_diags) , intent(in   ):: diaglookup
  logical,optional, intent(in   ):: skip

  character(len=*),parameter:: myname_=myname//'::obsnode_xread_'
  integer(i_kind),allocatable,dimension(:):: ich
  integer(i_kind):: k,nlco,nlevp
  logical:: skip_
_ENTRY_(myname_)
  skip_=.false.
  if(present(skip)) skip_=skip

  istat=0
  if(skip_) then
     read(iunit,iostat=istat)
     if (istat/=0) then
        call perr(myname_,'skipping read(%nlco), istat =',istat)
        _EXIT_(myname_)
        return
     end if
     read(iunit,iostat=istat)
     if(istat/=0) then
        call perr(myname_,'skipping read(%(res,err2,...)), istat =',istat)
        _EXIT_(myname_)
        return
     endif

  else
     read(iunit,iostat=istat)    anode%nlco
     if (istat/=0) then
        call perr(myname_,'read(%nlco), istat =',istat)
        _EXIT_(myname_)
        return
     end if

     call die(myname_,'unfinished implementation here')
     ! The size of %wij(8,:) seems in a mess.  It is nlev, nlevs, or nsig?
     ! See setupco() for additional details.

                ! re-allocation is needed, because a node may have been read in
                ! but excluded later.
     if(associated(anode%res    )) deallocate(anode%res    )
     if(associated(anode%diags  )) deallocate(anode%diags  )
     if(associated(anode%err2   )) deallocate(anode%err2   )
     if(associated(anode%raterr2)) deallocate(anode%raterr2)
     if(associated(anode%prs    )) deallocate(anode%prs    )
     if(associated(anode%ipos   )) deallocate(anode%ipos   )
     if(associated(anode%ak     )) deallocate(anode%ak     )
     if(associated(anode%ap     )) deallocate(anode%ap     )
     if(associated(anode%wk     )) deallocate(anode%wk     )
     if(associated(anode%wij    )) deallocate(anode%wij    )

     nlco = anode%nlco
     nlevp= max(nlco,1)

     allocate( anode%res    (nlco) , &
               anode%diags  (nlco) , &
               anode%err2   (nlco) , &
               anode%raterr2(nlco) , &
               anode%prs    (nlevp), &
               anode%ipos   (nlco) , & 
               anode%ak(nlco,nlco) , & 
               anode%ap     (nlco) , & 
               anode%wk     (nsig) , &
               anode%wij  (8,nsig)   )

     allocate(       ich    (nlco)   )

     read(iunit,iostat=istat)          ich    , & !(nlco)
                                 anode%res    , & !(nlco)
                                 anode%err2   , & !(nlco)
                                 anode%raterr2, & !(nlco)
                                 anode%prs    , & !(nlevp)
                                 anode%ipos   , & !(nlco)
                                 anode%ak     , & !(nlco,nlco)
                                 anode%ap     , & !(nlco)
                                 anode%wk     , & !(  nsig)
                                 anode%wij    , & !(8,nsig)
                                 anode%ij         !(4)
     if (istat/=0) then
        call perr(myname_,'read(ich,%(...)), istat =',istat)
        _EXIT_(myname_)
        return
     end if

     do k=1,nlco
      ! Not sure if ich=k or ich=ich(k) for now.  More verification is needed.
        anode%diags(k)%ptr => obsdiaglookup_locate(diaglookup,anode%idv,anode%iob,ich(k))
        if(.not.associated(anode%diags(k)%ptr)) then
           call perr(myname_,'obsdiaglookup_locate(k), k =',anode%idv)
           call perr(myname_,'                      %idv =',anode%idv)
           call perr(myname_,'                      %iob =',anode%iob)
           call perr(myname_,'                       ich =',ich(k))
           call  die(myname_)
           istat=-1
           _EXIT_(myname_)
           return
        endif
     enddo
     deallocate(ich)
  endif
_EXIT_(myname_)
  return
end subroutine obsnode_xread_

subroutine obsnode_xwrite_(anode,junit,jstat)
  implicit none
  class(colvknode),intent(in):: anode
  integer(i_kind),intent(in   ):: junit
  integer(i_kind),intent(  out):: jstat

  character(len=*),parameter:: myname_=myname//'::obsnode_xwrite_'
  integer(i_kind):: k
_ENTRY_(myname_)

  jstat=0
  write(junit,iostat=jstat) anode%nlco
  if (jstat/=0) then
     call perr(myname_,'write(%nlco), jstat =',jstat)
     _EXIT_(myname_)
     return
  end if

  write(junit,iostat=jstat) (/ (anode%diags(k)%ptr%ich, k=1,anode%nlco) /), & !(nlco)
                                anode%res    , & !(nlco)
                                anode%err2   , & !(nlco)
                                anode%raterr2, & !(nlco)
                                anode%prs    , & !(nlevp)
                                anode%ipos   , & !(nlco)
                                anode%ak     , & !(nlco,nlco)
                                anode%ap     , & !(nlco)
                                anode%wk     , & !(  nsig)
                                anode%wij    , & !(8,nsig)
                                anode%ij         !(4)
  if (jstat/=0) then
     call perr(myname_,'write(%(res,err2,...)), istat =',jstat)
     _EXIT_(myname_)
     return
  end if
_EXIT_(myname_)
  return
end subroutine obsnode_xwrite_

subroutine obsnode_sethop_(anode)
  !use m_obsnode, only: obstype_gethop => obsnode_gethop
  use m_cvgridlookup, only: cvgridlookup_getiw
  use gridmod, only: nsig
  implicit none
  class(colvknode),intent(inout):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_sethop_'
  real(r_kind),dimension(4):: zwij
  integer(i_kind):: k
_ENTRY_(myname_)
  call cvgridlookup_getiw(anode%elat,anode%elon,anode%ij,zwij)
  do k=1,nsig
     anode%wij(1:4,k) = zwij(1:4)*anode%wk(k)
     anode%wij(5:8,k) = zwij(1:4)*(1._r_kind-anode%wk(k))
  enddo
_EXIT_(myname_)
  return
end subroutine obsnode_sethop_

function obsnode_isvalid_(anode) result(isvalid_)
  implicit none
  logical:: isvalid_
  class(colvknode),intent(in):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_isvalid_'
  integer(i_kind):: k
_ENTRY_(myname_)
  isvalid_ = all( (/ (associated(anode%diags(k)%ptr), k=1,anode%nlco) /) )
_EXIT_(myname_)
  return
end function obsnode_isvalid_

pure subroutine gettlddp_(anode,jiter,tlddp,nob)
  use kinds, only: r_kind
  implicit none
  class(colvknode), intent(in):: anode
  integer(kind=i_kind),intent(in):: jiter
  real(kind=r_kind),intent(inout):: tlddp
  integer(kind=i_kind),optional,intent(inout):: nob

  integer(kind=i_kind):: k
  do k=1,anode%nlco
     tlddp = tlddp + anode%diags(k)%ptr%tldepart(jiter)*anode%diags(k)%ptr%tldepart(jiter)
  enddo
  if(present(nob)) nob=nob+anode%nlco
  return
end subroutine gettlddp_

end module m_colvknode
