module m_oznode
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module m_oznode
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2016-05-18
!
! abstract: class-module of obs-type oznode (ozone layer amounts and total column)
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

  public:: oznode

  type,extends(obsnode):: oznode
     !type(oz_ob_type),pointer :: llpoint => null()
     type(aofp_obs_diag), dimension(:), pointer :: diags => null()
     real(r_kind),dimension(:),pointer :: res => null()
                                      !  ozone residual
     real(r_kind),dimension(:),pointer :: err2 => null()
                                      !  ozone error squared
     real(r_kind),dimension(:),pointer :: raterr2 => null()
                                      !  square of ratio of final obs error 
                                      !  to original obs error
     !real(r_kind)    :: time          !  observation time in sec     
     real(r_kind),dimension(:,:),pointer :: wij => null()
                                      !  horizontal interpolation weights
     real(r_kind),dimension(:),pointer :: prs => null()
                                      !  pressure levels
     integer(i_kind),dimension(:),pointer :: ipos  => null()
     integer(i_kind) :: nloz          ! number of levels for this profile
     integer(i_kind) :: ij(4)         !  horizontal locations
     !logical         :: luse          !  flag indicating if ob is used in pen.

     !integer(i_kind) :: idv,iob	      ! device id and obs index for sorting
     !real   (r_kind) :: elat, elon      ! earth lat-lon for redistribution
     !real   (r_kind) :: dlat, dlon      ! earth lat-lon for redistribution
     real   (r_kind) :: rozcon          ! see setupozlay() for about rozcon
     real   (r_kind),dimension(:),pointer :: dprsi           ! see setupozlay() for about prsitmp(k)-prsitmp(k+1)

     logical:: has_eff
     real(r_kind),dimension(:),pointer :: apriori    ! OMI retrieval first guess
     real(r_kind),dimension(:),pointer :: efficiency ! OMI efficiency factor
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
  end type oznode

  public:: oznode_typecast
  public:: oznode_nextcast
     interface oznode_typecast; module procedure typecast_ ; end interface
     interface oznode_nextcast; module procedure nextcast_ ; end interface

  public:: oznode_appendto
     interface oznode_appendto; module procedure appendto_ ; end interface

  character(len=*),parameter:: myname="m_oznode"

#include "myassert.H"
#include "mytrace.H"
contains
function typecast_(anode) result(ptr_)
!-- cast a class(obsnode) to a type(oznode)
  use m_obsnode, only: obsnode
  implicit none
  type (oznode ),pointer:: ptr_
  class(obsnode),pointer,intent(in):: anode
  ptr_ => null()
  if(.not.associated(anode)) return
        ! logically, typecast of a null-reference is a null pointer.
  select type(anode)
     type is(oznode)
        ptr_ => anode
  end select
  return
end function typecast_

function nextcast_(anode) result(ptr_)
!-- cast an obsnode_next(obsnode) to a type(oznode)
  use m_obsnode, only: obsnode,obsnode_next
  implicit none
  type (oznode ),pointer:: ptr_
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
  type(oznode),pointer,intent(in):: anode
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
  mytype="[oznode]"
end function mytype

subroutine obsnode_clean_(anode)
  implicit none
  class(oznode),intent(inout):: anode

  character(len=*),parameter:: myname_=myname//'.obsnode_clean_'
_ENTRY_(myname_)
!_TRACEV_(myname_,'%mytype() =',anode%mytype())
  if(associated(anode%diags  )) deallocate(anode%diags  )
  if(associated(anode%res    )) deallocate(anode%res    )
  if(associated(anode%err2   )) deallocate(anode%err2   )
  if(associated(anode%raterr2)) deallocate(anode%raterr2)
  if(associated(anode%prs    )) deallocate(anode%prs    )
  if(associated(anode%ipos   )) deallocate(anode%ipos   )
  if(associated(anode%apriori)) deallocate(anode%apriori)
  if(associated(anode%efficiency)) deallocate(anode%efficiency)
  if(associated(anode%dprsi  )) deallocate(anode%dprsi  )
  if(associated(anode%wij    )) deallocate(anode%wij    )
_EXIT_(myname_)
  return
end subroutine obsnode_clean_

subroutine obsnode_xread_(anode,iunit,istat,diaglookup,skip)
  use m_obsdiagnode, only: obsdiaglookup_locate
  use gridmod, only: nsig
  implicit none
  class(oznode) , intent(inout):: anode
  integer(i_kind) , intent(in   ):: iunit
  integer(i_kind) , intent(  out):: istat
  type(obs_diags) , intent(in   ):: diaglookup
  logical,optional, intent(in   ):: skip

  character(len=*),parameter:: myname_=myname//'.obsnode_xread_'
  integer(i_kind):: k,mloz,mlev,mlevp,meff,msig
  integer(i_kind),allocatable,dimension(:):: ich
  logical:: skip_
_ENTRY_(myname_)

  skip_=.false.
  if(present(skip)) skip_=skip

  istat=0
  if(skip_) then
     read(iunit,iostat=istat)
     if (istat/=0) then
        call perr(myname_,'skipping read(%nloz), istat =',istat)
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
     read(iunit,iostat=istat)    anode%nloz,mlev,mlevp,meff,msig
     if (istat/=0) then
        call perr(myname_,'read(%nloz), istat =',istat)
        _EXIT_(myname_)
        return
     end if

     ASSERT(mlev==anode%nloz+1)
     ASSERT(msig==nsig)

     if(associated(anode%diags  )) deallocate(anode%diags  )
     if(associated(anode%res    )) deallocate(anode%res    )
     if(associated(anode%err2   )) deallocate(anode%err2   )
     if(associated(anode%raterr2)) deallocate(anode%raterr2)
     if(associated(anode%prs    )) deallocate(anode%prs    )
     if(associated(anode%ipos   )) deallocate(anode%ipos   )
     if(associated(anode%apriori)) deallocate(anode%apriori)
     if(associated(anode%efficiency)) deallocate(anode%efficiency)
     if(associated(anode%dprsi  )) deallocate(anode%dprsi  )
     if(associated(anode%wij    )) deallocate(anode%wij    )

     allocate( anode%diags  (mlev ), &
               anode%res    (mlev ), &
               anode%err2   (mlev ), &
               anode%raterr2(mlev ), &
               anode%prs    (mlevp), &
               anode%ipos   (mlev ), &
               anode%apriori(meff ), &
               anode%efficiency(meff), &
               anode%wij  (4,msig ), &
               anode%dprsi(  msig )  )

     allocate(ich(mlev))

     read(iunit,iostat=istat)          ich    , &
                                 anode%res    , &
                                 anode%err2   , &
                                 anode%raterr2, &
                                 anode%prs    , &
                                 anode%ipos   , &
                                 anode%apriori, &
                                 anode%efficiency, &
                                 anode%rozcon , &
                                 anode%dprsi  , &
                                 anode%wij    , &
                                 anode%ij
     if (istat/=0) then
        call perr(myname_,'read(%(res,err2,...)), istat =',istat)
        call perr(myname_,'                        meff =',meff)
        call perr(myname_,'                        mloz =',mloz)
        call perr(myname_,'                        mlev =',mlev)
        call perr(myname_,'                       mlevp =',mlevp)
        call perr(myname_,'                        msig =',msig)
        _EXIT_(myname_)
        return
     end if

     do k=1,mlev
        anode%diags(k)%ptr => obsdiaglookup_locate(diaglookup,anode%idv,anode%iob,ich(k))
        if(.not.associated(anode%diags(k)%ptr)) then
           call perr(myname_,'obsdiaglookup_locate(k), k =',k)
           call perr(myname_,'                      %idv =',anode%idv)
           call perr(myname_,'                      %iob =',anode%iob)
           call perr(myname_,'                       ich =',ich(k))
           call  die(myname_)
        endif
     enddo
     deallocate(ich)
  endif
_EXIT_(myname_)
  return
end subroutine obsnode_xread_

subroutine obsnode_xwrite_(anode,junit,jstat)
  implicit none
  class(oznode),intent(in):: anode
  integer(i_kind),intent(in   ):: junit
  integer(i_kind),intent(  out):: jstat

  character(len=*),parameter:: myname_=myname//'.obsnode_xwrite_'
  integer(i_kind):: k,mlev,mlevp,meff,msig
_ENTRY_(myname_)

  mlev =size(anode%diags)
  mlevp=size(anode%prs)
  if(mlev/=anode%nloz+1) then
     call perr(myname_,'mlev/=anode%nloz+1, mlev =',mlev)
     call perr(myname_,'                   mlevp =',mlevp)
     call perr(myname_,'                   %nloz =',anode%nloz)
     call  die(myname_)
  endif
  ASSERT(mlev==anode%nloz+1)

  meff =size(anode%apriori)
  msig =size(anode%dprsi)

  write(junit,iostat=jstat) anode%nloz,mlev,mlevp,meff,msig
  if(jstat/=0) then
     call perr(myname_,'write(%nloz,...)), jstat =',jstat)
     _EXIT_(myname_)
     return
  endif

  write(junit,iostat=jstat) (/ (anode%diags(k)%ptr%ich,k=1,mlev)/), &
                                anode%res    , & ! nlev
                                anode%err2   , & ! nlev
                                anode%raterr2, & ! nlev
                                anode%prs    , & ! nlevp
                                anode%ipos   , & !
                                anode%apriori, &
                                anode%efficiency, &
                                anode%rozcon , &
                                anode%dprsi  , &
                                anode%wij    , &
                                anode%ij
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
  use gridmod, only: nsig
  implicit none
  class(oznode),intent(inout):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_sethop_'
  real(r_kind),dimension(4):: twij
  integer(i_kind):: k
_ENTRY_(myname_)
  ASSERT(size(anode%wij,2)==nsig)
  ASSERT(nsig>0)

  call cvgridlookup_getiw(anode%elat,anode%elon,anode%ij,twij(:))
  do k=1,nsig
     anode%wij(:,k)=twij(:)*anode%rozcon*anode%dprsi(k)
  enddo
_EXIT_(myname_)
  return
end subroutine obsnode_sethop_

function obsnode_isvalid_(anode) result(isvalid_)
  implicit none
  logical:: isvalid_
  class(oznode),intent(in):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_isvalid_'
  integer(i_kind):: k
_ENTRY_(myname_)
  isvalid_= all ( (/ (associated(anode%diags(k)%ptr), k=1,anode%nloz) /) )
_EXIT_(myname_)
  return
end function obsnode_isvalid_

pure subroutine gettlddp_(anode,jiter,tlddp,nob)
  use kinds, only: r_kind
  implicit none
  class(oznode), intent(in):: anode
  integer(kind=i_kind),intent(in):: jiter
  real(kind=r_kind),intent(inout):: tlddp
  integer(kind=i_kind),optional,intent(inout):: nob

  integer(kind=i_kind):: k
  do k=1,anode%nloz
     tlddp = tlddp + anode%diags(k)%ptr%tldepart(jiter)*anode%diags(k)%ptr%tldepart(jiter)
  enddo
  if(present(nob)) nob=nob+anode%nloz
  return
end subroutine gettlddp_

end module m_oznode
