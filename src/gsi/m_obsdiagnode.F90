module m_obsdiagnode
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module m_obsdiagnode
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2016-05-18
!
! abstract: module of node type obs_diag and linked-list type obs_diags.
!
! program history log:
!   2016-05-18  j guo   - added this document block for the initial implementation.
!   2016-06-24  j.guo   - Added support of using m_latlonrange to find a cluster
!                         latlonrange from (elat,elon) values of observations.
!                       . cleaned out some components from obsdiagnode, which
!                         were put in for debugging purposes. (%dlat,%dlon).
!                       . removed some earlier routines for debuggings and
!                         testings.  e.g. lmock_() and obsnode_mock_().
!                       . use a fixed miter size for both write_() and read_(),
!                         for a simpler control in the future.
!                       . renamed lsize_() to lcount_().  Then reimplemented a
!                         new lsize_() to separate different functionalities.
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
  use kinds , only: i_kind,r_kind
  use mpeu_util, only: assert_,tell,warn,perr,die
  implicit none

  private

  public:: obs_diag
  public:: obs_diags
  public:: fptr_obsdiagnode

        ! Primery behaviors:
  public:: obsdiagllist_reset   ! destructor + initializer
  public:: obsdiagllist_appendnode
  public:: obsdiagllist_rewind  ! rewind an obsdiagllist
  public:: obsdiagllist_nextnode

  public:: obsdiagllist_headnode
  public:: obsdiagllist_tailnode

  public:: obsdiagllist_read    ! reader, for input
  public:: obsdiagllist_write   ! writer, for otuput
  public:: obsdiagllist_lsize   ! size inquiry
  public:: obsdiagllist_lcount  ! size inquiry with recount
  public:: obsdiagllist_lsort   ! sort nodes according to their keys
  public:: obsdiagllist_checksum! size consistency checking
  public:: obsdiagllist_summary ! status report

  interface obsdiagllist_reset ; module procedure  lreset_; end interface
  interface obsdiagllist_rewind; module procedure lrewind_; end interface
  interface obsdiagllist_read  ; module procedure   lread_; end interface
  interface obsdiagllist_checksum; module procedure &
     lchecksum_  , &
     lchecksum1_ , &
     lchecksum2_ ; end interface
  interface obsdiagllist_lsize  ; module procedure lsize_   ; end interface
  interface obsdiagllist_lcount ; module procedure lcount_  ; end interface
  interface obsdiagllist_lsort  ; module procedure lsort_   ; end interface
  interface obsdiagllist_write  ; module procedure lwrite_  ; end interface
  interface obsdiagllist_summary; module procedure lsummary_; end interface

  interface obsdiagllist_appendnode; module procedure obsnode_append_; end interface
  interface obsdiagllist_nextnode  ; module procedure &
     obsnode_next_, &
     make_or_next_; end interface

  interface obsdiagllist_headnode  ; module procedure lheadnode_  ; end interface
  interface obsdiagllist_tailnode  ; module procedure ltailnode_  ; end interface

        ! node lookup, secondary function with its searching component
  public:: obsdiaglookup_build  ! setup, its searching component
  public:: obsdiaglookup_locate ! node lookup, with the searching component
  public:: obsdiaglookup_clean  ! clean, its searching component

  interface obsdiaglookup_build ; module procedure lbuild_; end interface
  interface obsdiaglookup_locate; module procedure locate_; end interface
  interface obsdiaglookup_clean ; module procedure lclean_; end interface

  public:: obsdiagllist_dump
  interface obsdiagllist_dump; module procedure ldump_; end interface

  !public:: obsdiagnode_append
  !interface obsdiagnode_append; module procedure obsnode_append_; end interface
  !public:: obsdiagnode_first
  !interface obsdiagnode_first ; module procedure  obsnode_first_; end interface
  !public:: obsdiagnode_next
  !interface obsdiagnode_next  ; module procedure   obsnode_next_; end interface
  public:: obsdiagnode_init
  public:: obsdiagnode_assert
  public:: obsdiagnode_set
  public:: obsdiagnode_get
  interface obsdiagnode_init  ; module procedure obsnode_init_; end interface
  interface obsdiagnode_assert; module procedure anode_assert_; end interface
  interface obsdiagnode_set   ; module procedure obsnode_set_ ; end interface
  interface obsdiagnode_get   ; module procedure obsnode_get_ ; end interface

  type obs_diag
     type(obs_diag), pointer :: next => null()
     real(r_kind), pointer :: nldepart(:) => null()    ! (miter+1)
     real(r_kind), pointer :: tldepart(:) => null()    ! (miter)
     real(r_kind), pointer :: obssen(:)   => null()    ! (miter)
     real(r_kind) :: wgtjo
     real(r_kind) :: elat, elon         ! earth lat-lon for redistribution
     integer(i_kind) :: idv,iob,ich     ! device, obs., and channel indices
     logical, pointer :: muse(:)          => null()    ! (miter+1), according the setup()s
     logical :: luse
  end type obs_diag

  type fptr_obsdiagnode         ! Fortran array element of a type(obs_diag) pointer
     type(obs_diag),pointer:: ptr => null()
  end type fptr_obsdiagnode

  type:: obs_diags
     integer(i_kind):: n_alloc=0
     type(obs_diag), pointer :: head => null()
     type(obs_diag), pointer :: tail => null()
     type(fptr_obsdiagnode), allocatable, dimension(:):: lookup
  end type obs_diags

#include "myassert.H"
#include "mytrace.H"

  character(len=*),parameter:: myname="m_obsdiagnode"

#define _obsnode_   obs_diag
#define _obsllist_  obs_diags

contains
subroutine lgotonode_(headll,thisnode)
! Move the tail pointer to thisnode. 
! It is assumed that given thisnode is one of nodes in the list.  Otherwise
! this function would break the list.
  implicit none
  type(_obsllist_),target,intent(inout):: headll
  type(_obsnode_ ),target,intent(in   ):: thisnode
  headll%tail => thisnode
end subroutine lgotonode_

function lheadnode_(headll) result(here_)
! Return the head node
  implicit none
  type(_obsnode_),pointer:: here_
  type(_obsllist_),target,intent(in):: headll
  here_ => headll%head
end function lheadnode_

function ltailnode_(headll) result(here_)
! Return the current tail node
  implicit none
  type(_obsnode_ ),pointer:: here_
  type(_obsllist_),target,intent(in):: headll
  here_ => headll%tail
end function ltailnode_

subroutine lwrite_(diagll,iunit,luseonly,jiter,miter,jj_type,ii_bin,luserange)
  use m_latlonrange  , only: latlonrange
  use m_latlonrange  , only: latlonrange_enclose
  use mpeu_util, only: stdout
  use mpeu_util, only: stdout_lead
  implicit none
  type(_obsllist_)    ,intent(inout):: diagll   ! the linked list of data
  integer(kind=i_kind),intent(in   ):: iunit    ! the output unit
  logical             ,intent(in   ):: luseonly ! write only if(luse)
  integer(kind=i_kind),intent(in   ):: jiter    ! diag width for the IO (or this iter)
  integer(kind=i_kind),intent(in   ):: miter    ! diag width of the memory
  integer(kind=i_kind),intent(in   ):: jj_type, ii_bin
  type(latlonrange),optional,intent(inout):: luserange

  character(len=*),parameter:: myname_=myname//"::lwrite_"
  integer(kind=i_kind):: iobs,kobs,lobs,mobs
  integer(kind=i_kind):: istat
  type(_obsnode_), pointer:: inode
  logical:: isluse_
_ENTRY_(myname_)
!_timer_on_(myname_)

  lobs=obsdiagllist_lcount(diagll,luseonly=luseonly)
  mobs=lobs
  if(.not.luseonly) mobs=obsdiagllist_lsize(diagll)

  call obsheader_write_(iunit,ii_bin,jj_type,lobs,jiter,miter,istat)
  if(istat/=0) then
     call perr(myname_,'obsheader_write_(), istat =',istat)
     call perr(myname_,'                    iunit =',iunit)
     call perr(myname_,'                   ii_bin =',ii_bin)
     call perr(myname_,'                    jtype =',jj_type)
     call perr(myname_,'                    jiter =',jiter)
     call perr(myname_,'                    miter =',miter)
     call perr(myname_,'    total-luse-node, lobs =',lobs)
     call perr(myname_,'     total-all-node, mobs =',mobs)
     call perr(myname_,'                 luseonly =',luseonly)
     call  die(myname_)
  endif

  _TRACE_(myname_,'looping through obshead pointers')

  if(lobs<=0) then
     !_timer_off_(myname_)
     _EXIT_(myname_)
     return
  endif

  iobs=0
  kobs=0
  inode => obsnode_first_(diagll)
  do while(associated(inode))
     iobs=iobs+1
     isluse_=obsnode_isluse_(inode)
     if(isluse_ .or. .not.luseonly) then

                ! Update luserange with a luse observation, for the lat-lon-
                ! range on the current pe.

        if(isluse_ .and. present(luserange)) &
           call latlonrange_enclose(luserange,inode%elat,inode%elon)

                ! Count it, then write the node out.  Use of miter suggests a
                ! fixed output size.
        kobs=kobs+1
        call obsnode_write_(inode,iunit,miter,istat)
        if(istat/=0) then
           call perr(myname_,'obsnode_write_(), istat =',istat)
           call perr(myname_,'                  iunit =',iunit)
           call perr(myname_,'                  jiter =',jiter)
           call perr(myname_,'                  miter =',miter)
           call perr(myname_,'                 ii_bin =',ii_bin)
           call perr(myname_,'                  jtype =',jj_type)
           call perr(myname_,'current-luse-node, kobs =',kobs)
           call perr(myname_,' current-all-node, iobs =',iobs)
           call perr(myname_,'  total-luse-node, lobs =',lobs)
           call perr(myname_,'   total-all-node, mobs =',mobs)
           call perr(myname_,'               luseonly =',luseonly)
           call  die(myname_)
        endif
     endif
     inode => obsnode_next_(diagll)
  enddo

  ASSERT(kobs==lobs)
  ASSERT(iobs==mobs)

!_timer_off_(myname_)
_EXIT_(myname_)
  return
end subroutine lwrite_

subroutine ldump_(diagll,jiter)
  use mpeu_util, only: stdout
  implicit none
  type(_obsllist_),         intent(inout):: diagll  ! the list to dump
  integer(i_kind ),optional,intent(in   ):: jiter   ! jiter of diagll

  character(len=*),parameter:: myname_=myname//"::ldump_"
  integer(kind=i_kind):: iobs,lobs,mobs
  integer(kind=i_kind):: jiter_
  type(_obsnode_), pointer:: inode
  logical:: isluse_,ismuse_
_ENTRY_(myname_)
!_timer_on_(myname_)
  jiter_=0
  if(present(jiter)) jiter_=jiter

  call lbuild_(diagll)        ! create a pointer array %lookup, sorted by (idv,iob,ich)

  lobs=0
  mobs=0
  do iobs=1,size(diagll%lookup(:))
     inode => diagll%lookup(iobs)%ptr

     isluse_=obsnode_isluse_(inode)
     if(isluse_) lobs=lobs+1

     ismuse_=jiter_>=1.and.jiter_<=size(inode%muse)
     if(ismuse_) ismuse_=inode%muse(jiter_)
     if(ismuse_) mobs=mobs+1

     write(stdout,'(2x,2l1,3i8,2x,2f12.4)') isluse_,ismuse_, &
        inode%idv,inode%iob,inode%ich, inode%elat,inode%elon
  enddo
  write(stdout,'(2x,a,4i8)') '***',jiter_,size(diagll%lookup(:)),lobs,mobs
  call lclean_(diagll)        ! destroy the pointer array %lookup.

!_timer_off_(myname_)
_EXIT_(myname_)
  return
end subroutine ldump_

subroutine lread_(diagll,iunit,redistr,jiter,miter,jj_type,ii_bin,jread,leadnode,jiter_expected)
!_timer_use_
  implicit none
  type(_obsllist_),intent(inout):: diagll
  integer(kind=i_kind),intent(in   ):: iunit
  logical        ,intent(in   ):: redistr
  integer(kind=i_kind),intent(in   ):: jiter
  integer(kind=i_kind),intent(in   ):: miter
  integer(kind=i_kind),intent(in   ):: jj_type, ii_bin
  integer(kind=i_kind),intent(  out):: jread
  type(_obsnode_), pointer, intent(out):: leadnode

  integer(kind=i_kind),intent(in),optional:: jiter_expected

  character(len=*),parameter:: myname_=myname//"::lread_"
  integer(kind=i_kind):: ki,kj,kobs
  integer(kind=i_kind):: kiter,miter_read
        ! jiter : current iter count
        ! miter : maximum iter size
        ! kiter(read): current iter count as it was written
        ! miter_read : maximum iter size as it was written
  integer(kind=i_kind):: kk,istat
  type(_obsnode_), pointer:: anode
_ENTRY_(myname_)
!_timer_on_(myname_)
!call timer_ini(myname_)

  call obsheader_read_(iunit,ki,kj,kobs,kiter,miter_read,istat)
  if(istat/=0) then
     call perr(myname_,'obsheader_read_(), istat =',istat)
     call perr(myname_,'                   iunit =',iunit)
     call  die(myname_)
  endif

  if(ki/=ii_bin .or. kj/=jj_type .or. miter/=miter_read) then
     call perr(myname_,'obsheader_read_(), unexpected header values (ii,jj,miter)')
     call perr(myname_,'  expecting miter =',miter)
     call perr(myname_,'     actual miter =',miter_read)
     call perr(myname_,'     expecting ii =',ii_bin)
     call perr(myname_,'        actual ii =',ki)
     call perr(myname_,'     expecting jj =',jj_type)
     call perr(myname_,'        actual jj =',kj)
     call  die(myname_)
  endif

  if(present(jiter_expected)) then
     if(jiter_expected>=0) then
        if(kiter/=jiter_expected) then
           call perr(myname_,'obsheader_read_(), unexpected input jiter =',kiter)
           call perr(myname_,'                         with input miter =',miter_read)
           call perr(myname_,'                    expecting input jiter =',jiter_expected)
           call perr(myname_,'                                    miter =',miter)
           call perr(myname_,'                                    jiter =',jiter)
           call  die(myname_)
        endif
     endif
  endif
  jread=kiter

                !-- construct an an_obsnode
  leadnode => null()
  anode => obsnode_alloc_(miter)
  do kk=1,kobs
                !-- initialize an_obsnode from a file (iunit).  Use of miter suggests a
                !-- fixed input size.
     call obsnode_read_(anode,iunit,miter,istat,redistr=redistr)
     if(istat<0) then
        call perr(myname_,'obsnode_read_(), istat =',istat)
        call perr(myname_,'               redistr =',redistr)
        call  die(myname_)
     endif

                ! istat <0:     a failed read(anode)
                !      ==0:     passed, thus an incomplete anode
                !       >0:     a good anode to keep
     if(istat==0) cycle
     if(redistr) call obsnode_setluse_(anode)

                ! keep this obsnode in its linked-list, diagll := obsdiags(jj,ii)
     call obsnode_append_(diagll,anode)
                !-- mark the beginning of this linked-list segment
     if(.not.associated(leadnode)) leadnode => anode

                !-- drop this anode, to construct a new.  This _alloc_
                !   ensures an anode is not in anyway referencible to
                !   the one that has been appended to the linked-list.
                !   Then, a deep-deallocation of anode is alwasy safe.
     anode => obsnode_alloc_(miter)
  enddo  ! < kobs >
  call obsnode_dealloc_(anode,deep=.true.)  ! Clean up the malloc of anode

! ----------------------------------------------------------
!call timer_fnl(myname_)
!_timer_off_(myname_)
_EXIT_(myname_)
  return
end subroutine lread_

subroutine lreset_(diagll)
  implicit none
  type(_obsllist_), intent(inout):: diagll

  character(len=*),parameter:: myname_=myname//"::lreset_"
  type(_obsnode_),pointer:: l_obsnode
  type(_obsnode_),pointer:: n_obsnode
  integer(kind=i_kind):: ip
_ENTRY_(myname_)

  l_obsnode => obsnode_first_(diagll)
  ip=0
  do while(associated(l_obsnode))
     ip=ip+1
     !_TRACEV_(myname_,'deallocating at ip =',ip)
     !call obsnode_check_(myname_,l_obsnode)
        ! Steps of forward resetting,
        !   (1) hold the %next node,
        !   (2) clean (leaving the %next node untouched,
        !   (3) deallocate the current node,
        !   (4) point the starting point to the %next node.
     n_obsnode => obsnode_next_(diagll)
     call obsnode_dealloc_(l_obsnode,deep=.true.)
     l_obsnode => n_obsnode
  enddo
  !n_obsnode   => null()
  !l_obsnode   => null()

  diagll%n_alloc = 0
  diagll%head => null()
  diagll%tail => null()
  if(allocated(diagll%lookup)) deallocate(diagll%lookup)

_EXIT_(myname_)
  return
end subroutine lreset_
subroutine lrewind_(diagll)
  implicit none
  type(_obsllist_),target,intent(inout):: diagll
  diagll%tail => null()
  return
end subroutine lrewind_

subroutine lchecksum_(diagll,leadnode,itype,ibin,sorted)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    lchecksum_
!   prgmmr:      J. Guo
!
! abstract: check the size values against a known counts.
!
! program history log:
!   2015-06-26  guo     - 
!
!   input argument list: (see Fortran declarations below)
!
!   output argument list: (see Fortran declarations below)
!
! attributes:
!   language: f90/f95/f2003/f2008
!   machine:
!
!$$$ end documentation block
  use mpeu_util, only: stdout
  use mpeu_util, only: stdout_lead
  implicit none
  type(_obsllist_), intent(in):: diagll
  type(_obsnode_ ), pointer, optional, intent(in):: leadnode
  integer(kind=i_kind),optional,intent(in ):: itype
  integer(kind=i_kind),optional,intent(in ):: ibin
  logical             ,optional,intent(out):: sorted

  character(len=*),parameter:: myname_=myname//"::lchecksum_"
  integer(kind=i_kind):: jtype,jbin
  integer(kind=i_kind):: mcount
  integer(kind=i_kind):: nuse,nooo,ndup
  integer(kind=i_kind),dimension(3):: ksum
!jtest
!  logical:: lasso,lhead

_ENTRY_(myname_)
!jtest
!  ASSERT(present(leadnode))
!  lasso=associated(leadnode)
!  lhead=associated(diagll%head,leadnode)

  mcount=lcount_(diagll,recount=.true.,nuse=nuse,nooo=nooo,ndup=ndup,ksum=ksum,leadnode=leadnode)
  if(present(sorted)) sorted = nooo==0.and.ndup==0

!jtest
!  if(mcount/=diagll%n_alloc) then
!     call perr(myname_,'checksum failed, mcount =',mcount)
!     call perr(myname_,'      diagllist%n_alloc =',diagll%n_alloc)
!     if(present(itype)) &
!     call perr(myname_,'                  itype =',itype)
!     if(present(ibin)) &
!     call perr(myname_,'                   ibin =',ibin)
!     call  die(myname_)
!  endif

  if(present(itype)) jtype=itype
  if(present(ibin))  jbin =ibin
_EXIT_(myname_)
  return
end subroutine lchecksum_
subroutine lchecksum1_(diagll,leadnode,itype)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    lchecksum_
!   prgmmr:      J. Guo
!
! abstract: check the size values against a known counts.
!
! program history log:
!   2015-06-26  guo     - 
!
!   input argument list: (see Fortran declarations below)
!
!   output argument list: (see Fortran declarations below)
!
! attributes:
!   language: f90/f95/f2003/f2008
!   machine:
!
!$$$ end documentation block
  implicit none
  type(_obsllist_), dimension(:),intent(in):: diagll
  integer(kind=i_kind),optional,intent(in):: itype
  type(fptr_obsdiagnode),optional,dimension(:),intent(in):: leadnode

  character(len=*),parameter:: myname_=myname//"::lchecksum1_"
  integer(kind=i_kind):: i
_ENTRY_(myname_)
  if(present(leadnode)) then
     ASSERT(size(diagll)==size(leadnode))
     do i=1,size(diagll)
        call lchecksum_(diagll(i),itype=itype,ibin=i,leadnode=leadnode(i)%ptr)
     enddo
  else
     do i=1,size(diagll)
        call lchecksum_(diagll(i),itype=itype,ibin=i)
     enddo
  endif
_EXIT_(myname_)
  return
end subroutine lchecksum1_
subroutine lchecksum2_(diagll)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    lchecksum_
!   prgmmr:      J. Guo
!
! abstract: check the size values against a known counts.
!
! program history log:
!   2015-06-26  guo     - 
!
!   input argument list: (see Fortran declarations below)
!
!   output argument list: (see Fortran declarations below)
!
! attributes:
!   language: f90/f95/f2003/f2008
!   machine:
!
!$$$ end documentation block
  implicit none
  type(_obsllist_), dimension(:,:),intent(in):: diagll

  character(len=*),parameter:: myname_=myname//"::lchecksum2_"
  integer(kind=i_kind):: it,ib
_ENTRY_(myname_)
  do it=1,size(diagll,1)
     do ib=1,size(diagll,2)
        call lchecksum_(diagll(it,ib),itype=it,ibin=ib)
     enddo
  enddo
_EXIT_(myname_)
  return
end subroutine lchecksum2_

subroutine lsummary_(diagll,verbose)
  implicit none
  type(_obsllist_), intent(in):: diagll
  logical,optional, intent(in):: verbose

  character(len=*),parameter:: myname_=myname//"::lsummary_"
  type(_obsnode_ ), pointer:: inode
  type(_obsllist_), target :: templl
  integer(kind=i_kind):: iobs_
  logical:: verbose_
  verbose_=.false.
  if(present(verbose)) verbose_=verbose
_ENTRY_(myname_)

  if(verbose_) then
     templl = diagll
     iobs_ = 0
     inode => obsnode_first_(templl)
     do while(associated(inode))
        iobs_=iobs_+1
        call obsnode_show_(inode,iobs_)
        inode => obsnode_next_(templl)
     enddo
  endif
_EXIT_(myname_)
  return
end subroutine lsummary_

function lsize_(diagll) result(lobs_)
  implicit none
  integer(kind=i_kind):: lobs_
  type(_obsllist_),   target, intent(in):: diagll
  lobs_=diagll%n_alloc
end function lsize_

function lcount_(diagll,luseonly,recount,nuse,nooo,ndup,ksum,leadnode) result(lobs_)
  use mpeu_util, only: assert_
  implicit none
  integer(kind=i_kind):: lobs_
  type(_obsllist_),   target, intent(in):: diagll
  logical         , optional, intent(in):: luseonly
  logical         , optional, intent(in):: recount
  integer(kind=i_kind),optional,intent(out):: nuse      ! no. of luse
  integer(kind=i_kind),optional,intent(out):: nooo      ! no. out-of-orders
  integer(kind=i_kind),optional,intent(out):: ndup      ! no. duplicates
  integer(kind=i_kind),optional,dimension(:),intent(out):: ksum ! key value sum
  type(_obsnode_ ), pointer, optional, intent(in):: leadnode

  character(len=*),parameter:: myname_=myname//"::lcount_"
  type(_obsnode_ ), pointer:: inode
  type(_obsllist_), target :: templl
  integer(kind=i_kind):: nuse_
  integer(kind=i_kind):: k
  integer(kind=i_kind),dimension(3) :: kprev
  logical:: luseonly_,recount_,checksum_
_ENTRY_(myname_)

  luseonly_=.false.
  if(present(luseonly)) luseonly_=luseonly
  recount_ =.false.
  if(present(recount )) recount_ =recount
  if(present(leadnode)) recount_ =.true.

  checksum_= present(nuse).or.present(nooo).or.present(ndup).or.present(ksum)
  recount_ = recount_ .or. checksum_
  !if(.not.recount_) recount_ = checksum_

  if(present(ksum)) then
     ALWAYS_ASSERT( size(ksum)==size(kprev) )
  endif

  if(.not.(luseonly_.or.recount_)) then
     lobs_=diagll%n_alloc

  else  ! recount through the list
     templl = diagll    ! A copy of diagll, such that diagll can remain intent(in)

     lobs_ = 0
     nuse_ = 0

     if(checksum_) call checksum_init_(kprev,nooo=nooo,ndup=ndup,ksum=ksum)

     inode => obsnode_first_(templl,atnode=leadnode)
     do while(associated(inode))
        if(obsnode_isluse_(inode)) nuse_=nuse_+1
        if(.not.luseonly_ .or. obsnode_isluse_(inode)) lobs_=lobs_+1

        if(checksum_) call checksum_add_(kprev, &
           (/inode%idv,inode%iob,inode%ich/),nooo=nooo,ndup=ndup,ksum=ksum)

        inode => obsnode_next_(templl)
     enddo
     if(present(nuse)) nuse=nuse_
  endif

_EXIT_(myname_)
  return
contains
subroutine checksum_init_(kprev,nooo,ndup,ksum)
  implicit none
  integer(kind=i_kind),dimension(:),intent(out):: kprev
  integer(kind=i_kind),optional,intent(out):: nooo
  integer(kind=i_kind),optional,intent(out):: ndup
  integer(kind=i_kind),optional,dimension(:),intent(out):: ksum

  kprev(:)= 0
  if(present(nooo)) nooo=0
  if(present(ndup)) ndup=0
  if(present(ksum)) ksum(:)=0
end subroutine checksum_init_
subroutine checksum_add_(kprev,knext,nooo,ndup,ksum)
  implicit none
  integer(kind=i_kind),dimension(:),intent(inout):: kprev
  integer(kind=i_kind),dimension(:),intent(in   ):: knext
  integer(kind=i_kind),optional,intent(inout):: nooo
  integer(kind=i_kind),optional,intent(inout):: ndup
  integer(kind=i_kind),optional,dimension(:),intent(inout):: ksum

  k=compare_(kprev,knext)
  if(present(nooo).and.k> 0) nooo=nooo+1
  if(present(ndup).and.k==0) ndup=ndup+1
  if(present(ksum)) ksum(:)=ksum(:)+knext(:)
  kprev(:)=knext(:)
end subroutine checksum_add_
end function lcount_

function obsnode_first_(diagll,atnode) result(here_)
  implicit none
  type(_obsnode_ ), pointer     :: here_
  type(_obsllist_), target, intent(inout):: diagll
  type(_obsnode_ ), optional, pointer,intent(in):: atnode

  character(len=*),parameter:: myname_=myname//"::obsnode_first_"
_ENTRY_(myname_)
  !_TRACEV_(myname_,'%n_alloc =',diagll%n_alloc)
  !_TRACEV_(myname_,'associated(%head) =',associated(diagll%head))
  here_ => diagll%head
  if(present(atnode)) here_=>atnode
  diagll%tail => here_  ! update the tail-node

  if(associated(here_)) call obsnode_check_(myname_,here_)
_EXIT_(myname_)
  return
end function obsnode_first_

function obsnode_next_(diagll) result(next_)
  implicit none
  type(_obsnode_ ), pointer      :: next_
  type(_obsllist_), target, intent(inout):: diagll

  character(len=*),parameter:: myname_=myname//"::obsnode_next_"
_ENTRY_(myname_)
  next_ => diagll%head
  if(associated(diagll%tail)) next_ => diagll%tail%next
  diagll%tail => next_  ! update the tail-node
_EXIT_(myname_)
  return
end function obsnode_next_

function make_or_next_(diagll,create,idv,iob,ich,elat,elon,luse,miter) result(next_)
  implicit none
  type(_obsnode_ ), pointer      :: next_
  type(_obsllist_), target, intent(inout):: diagll

  logical             , intent(in):: create     ! make or next
  integer(kind=i_kind), intent(in):: idv,iob,ich
  real   (kind=r_kind), intent(in):: elat,elon
  logical             , intent(in):: luse
  integer(kind=i_kind), intent(in):: miter

  character(len=*),parameter:: myname_=myname//"::make_or_next_"
  logical:: matched
_ENTRY_(myname_)

  if(create) then
     allocate(next_)
     call obsnode_append_(diagll,next_)
     call obsnode_init_(next_,idv,iob,ich,elat,elon,luse,miter)

  else
     next_ => diagll%head
     if(associated(diagll%tail)) next_ => diagll%tail%next
     diagll%tail => next_  ! update the tail-node

     ! Check the next node against (idv,iob,ich)
     matched = associated(next_)
     if(matched) matched = next_%idv==idv .and. &
                           next_%iob==iob .and. &
                           next_%ich==ich

     if(.not.matched) then
        call   perr(myname_,"unexpected node, associated(next) =", associated(next_))
        call   perr(myname_,"          expecting (idv,iob,ich) =", (/idv,iob,ich/))
        call   perr(myname_,"                             elat =", elat)
        call   perr(myname_,"                             elon =", elon)
        if(associated(next_)) then
           call perr(myname_,"               next%(idv,iob,ich) =", (/next_%idv,next_%iob,next_%ich/))
           call perr(myname_,"                        next%elat =", next_%elat)
           call perr(myname_,"                        next%elon =", next_%elon)
           call perr(myname_,"                        next%luse =", next_%luse)
           call perr(myname_,"                  size(next%muse) =", size(next_%muse))
        endif
        call die(myname_)
     endif
  endif ! (create)
_EXIT_(myname_)
  return
end function make_or_next_

subroutine obsnode_append_(diagll,targetnode)
        ! Link the next node of the list to the given targetnode.  The return
        ! result is a pointer associated to the same targetnode.
!--  use jfunc, only: miter
  implicit none
  type(_obsllist_), intent(inout):: diagll
  type(_obsnode_ ), pointer, intent(in):: targetnode

  character(len=*),parameter:: myname_=myname//"::obsnode_append_"
!-- type(_obsnode_ ),pointer:: anode
_ENTRY_(myname_)
  if(.not.associated(diagll%head)) then
                ! this is a fresh starting -node- for this linked-list ...
     diagll%n_alloc = 1
     diagll%head => targetnode
     diagll%tail => diagll%head

  else
                ! this is for a new next -node- from here ...
     diagll%n_alloc = diagll%n_alloc +1
     diagll%tail%next => targetnode
     diagll%tail      => diagll%tail%next

     !diagll%tail%append(next_)
     !    append(t,next_)
     !           t%next => next_
     !           t      => t%next
  endif
  if(associated(diagll%tail)) diagll%tail%next => null()

!--  anode => diagll%tail
!--  ASSERT(lbound(anode%muse    ,1)==1.and.ubound(anode%muse    ,1)==miter+1)
!--  ASSERT(lbound(anode%nldepart,1)==1.and.ubound(anode%nldepart,1)==miter+1)
!--  ASSERT(lbound(anode%tldepart,1)==1.and.ubound(anode%tldepart,1)==miter  )
!--  ASSERT(lbound(anode%obssen  ,1)==1.and.ubound(anode%obssen  ,1)==miter  )
!--  anode => null()

_EXIT_(myname_)
  return
end subroutine obsnode_append_

subroutine obsnode_insert_(diagll,targetnode)
        ! Insert targetnode to diagll's current location, mostly %tail.  At the
        ! return, diagll%tail is associated to targetnode.
!--  use jfunc, only: miter
  implicit none
  type(_obsllist_), intent(inout):: diagll
  type(_obsnode_ ), pointer, intent(in):: targetnode

  character(len=*),parameter:: myname_=myname//"::obsnode_insert_"
  type(_obsnode_),pointer:: next_
_ENTRY_(myname_)
  if(.not.associated(diagll%head)) then
                ! This is a fresh start case: insert a node as append
     diagll%n_alloc = 1
     diagll%head => targetnode
     diagll%tail => diagll%head            ! now the current node
     diagll%tail%next => null()            ! set %next to nothing there before

  elseif(.not.associated(diagll%tail)) then
                ! This is a rewound case: insert a node as the new %head
     next_ => diagll%head
     diagll%n_alloc = diagll%n_alloc +1
     diagll%head      => targetnode
     diagll%tail      => diagll%head       ! now the current node
     diagll%tail%next => next_             ! set %next to the original %head

  else
                ! This is a normal case: insert a node in between %tail and
                ! %tail%next.
     next_ => diagll%tail%next
     diagll%n_alloc = diagll%n_alloc +1
     diagll%tail%next => targetnode
     diagll%tail      => diagll%tail%next  ! now the current node.
     diagll%tail%next => next_             ! set %next to the original %tail%next
        ! Note in the last stateument, targetnode%next has been implicitly modifed.
  endif

!--  associate(anode => diagll%tail)
!--    ASSERT(lbound(anode%muse    ,1)==1.and.ubound(anode%muse    ,1)==miter+1)
!--    ASSERT(lbound(anode%nldepart,1)==1.and.ubound(anode%nldepart,1)==miter+1)
!--    ASSERT(lbound(anode%tldepart,1)==1.and.ubound(anode%tldepart,1)==miter  )
!--    ASSERT(lbound(anode%obssen  ,1)==1.and.ubound(anode%obssen  ,1)==miter  )
!--  end associate ! (anode => diagll%tail)

_EXIT_(myname_)
  return
end subroutine obsnode_insert_

subroutine lsort_(diagll,itype,ibin)
!       lsort_: node-sort diagll, to line-up nodes according to their keys
!_timer_use_
!  use timermod , only: timer_ini,timer_fnl
  !use mpeu_util, only: indexSet
  !use mpeu_util, only: indexSort
  !use mpeu_util, only: die
  implicit none
  type(_obsllist_) , intent(inout):: diagll
  integer(kind=i_kind),optional,intent(in):: itype,ibin

  character(len=*),parameter:: myname_=myname//'::lsort_'
  integer(kind=i_kind):: i,nobs,mobs
  logical:: sorted
_ENTRY_(myname_)
!_timer_on_(myname_)
!  call timer_ini(myname_)

  call lchecksum_(diagll,itype=itype,ibin=ibin,sorted=sorted)
  if(sorted) then
     _EXIT_(myname_)
     return
  endif

        ! created a sorted table
  call lbuild_(diagll)

  nobs = diagll%n_alloc
  mobs = size(diagll%lookup(:))
  ASSERT(nobs==mobs)

        ! rebuild the linked-list
  diagll%n_alloc=0
  diagll%head => null()
  diagll%tail => null()

        ! rebuild the list according to the sorted table
  do i=1,mobs
     call obsnode_append_(diagll,diagll%lookup(i)%ptr)
  enddo
  ASSERT(nobs==diagll%n_alloc)
  if(associated(diagll%tail)) then
     ASSERT(.not.associated(diagll%tail%next))
  endif

        ! discard the sorted table
  call lclean_(diagll)

  call lchecksum_(diagll,itype=itype,ibin=ibin,sorted=sorted)
  if(.not.sorted) then
     call perr(myname_,'failed post-sorting lchecksum_(diagll), sorted =',sorted)
     if(present(itype)) &
     call perr(myname_,'                                         itype =',itype)
     if(present(ibin )) &
     call perr(myname_,'                                          ibin =',ibin )
     call  die(myname_)
  endif

!  call timer_fnl(myname_)
!_timer_off_(myname_)
_EXIT_(myname_)
  return
end subroutine lsort_

subroutine lbuild_(diagll,leadnode,jiter)
!_timer_use_
!  use timermod , only: timer_ini,timer_fnl
  use mpeu_util, only: indexset
  use mpeu_util, only: indexsort
  !use mpeu_util, only: die
  implicit none
  type(_obsllist_), intent(inout):: diagll
  type(_obsnode_ ), pointer, optional, intent(in):: leadnode
  integer(i_kind) , optional, intent(in):: jiter

  character(len=*),parameter:: myname_=myname//'::lbuild_'
  type(_obsnode_),pointer:: inode,pnode
  integer(kind=i_kind),allocatable,dimension(:):: indx,idv_,iob_,ich_
  integer(kind=i_kind):: i,m,n
  integer(kind=i_kind):: idum
  logical:: good
_ENTRY_(myname_)
!_timer_on_(myname_)
!  call timer_ini(myname_)
  if(present(jiter)) idum=jiter

        ! Mark the leading node
  inode => null()
  if(present(leadnode)) inode => leadnode
  if(.not.associated(inode)) inode => diagll%head

  m=diagll%n_alloc
  if(m<0) call die(myname_,'unexpected diagll, %n_alloc =',m)

        ! Count, starting from the leading node
  n=0
  pnode => inode
  do while(associated(pnode))
     n=n+1
     pnode => pnode%next
  enddo

  if(n>diagll%n_alloc) then
     call perr(myname_,'unexpected diagll, %n_alloc =',m)
     call  die(myname_,'               actual count =',n)
  endif

  allocate(diagll%lookup(n))
  allocate(indx(n),idv_(n),iob_(n),ich_(n))

  associate(lookup => diagll%lookup(:))
        ! loop over the linked-list, to get keys.
     i=0
     pnode => inode
     do while(associated(pnode))
        i=i+1
        if(i<=n) then
           lookup(i)%ptr => pnode
             idv_(i)     =  pnode%idv
             iob_(i)     =  pnode%iob
             ich_(i)     =  pnode%ich
           !call obsnode_get(idv=idv_(i),iob=iob_(i),ich=ich_(i))
        endif
        pnode => pnode%next
     enddo
  end associate

        ! sort %lookup(1:n), by its (idv,iob,ich) values
  call indexset (indx)
  call indexsort(indx,ich_)
  call indexsort(indx,iob_)
  call indexsort(indx,idv_)

  associate(lookup => diagll%lookup(:))
     lookup(1:n) = lookup(indx(1:n))
  end associate

  idv_(1:n) = idv_(indx(1:n))
  iob_(1:n) = iob_(indx(1:n))
  ich_(1:n) = ich_(indx(1:n))

  associate(lookup => diagll%lookup(:))
     good = .true.
     do i=1,n
        good = lookup(i)%ptr%idv==idv_(i) .and. &
               lookup(i)%ptr%iob==iob_(i) .and. &
               lookup(i)%ptr%ich==ich_(i)
        if(.not.good) exit
     enddo

     if(.not.good) then
        call perr(myname_,'verification failed at %lookup(i)%ptr,  i =',i)
        call perr(myname_,'                                 %ptr%idv =',lookup(i)%ptr%idv)
        call perr(myname_,'                                      idv_=',idv_(i))
        call perr(myname_,'                                 %ptr%iob =',lookup(i)%ptr%iob)
        call perr(myname_,'                                      iob_=',iob_(i))
        call perr(myname_,'                                 %ptr%ich =',lookup(i)%ptr%ich)
        call perr(myname_,'                                      ich_=',ich_(i))
        call die(myname_)
     endif
  end associate

  deallocate(indx,idv_,iob_,ich_)

!  call timer_fnl(myname_)
!_timer_off_(myname_)
_EXIT_(myname_)
  return
end subroutine lbuild_

subroutine lclean_(diagll)
  implicit none
  type(_obsllist_), intent(inout):: diagll

  character(len=*),parameter:: myname_=myname//'::lclean_'
  integer(kind=i_kind):: ier,i
_ENTRY_(myname_)
  associate(lookup => diagll%lookup(:))
     do i=1,size(lookup)
        lookup(i)%ptr => null()
     end do
  end associate
  deallocate(diagll%lookup,stat=ier)
  if(ier/=0) call die(myname_,'deallocate(diagll%lookup), stat =',ier)
_EXIT_(myname_)
  return
end subroutine lclean_

function locate_(diagll,idv,iob,ich) result(here_)
  use timermod , only: timer_ini,timer_fnl
  implicit none
  type(_obsnode_ ), pointer:: here_
  type(_obsllist_), intent(in):: diagll
  integer(kind=i_kind), intent(in):: idv,iob,ich

  character(len=*),parameter:: myname_=myname//"::locate_"
  type(_obsnode_ ),pointer:: idiag
  integer(kind=i_kind):: m,i,lb,ub
  logical:: done
_ENTRY_(myname_)
  call timer_ini(myname_)

  here_ => null()     ! return null() if the key is not located.

  associate(lookup => diagll%lookup(:))
     lb=lbound(lookup,1)
     ub=ubound(lookup,1)
     done=.false.
     do while(.not.done)
        i=(lb+ub)/2
        idiag => lookup(i)%ptr

        m=compare_((/idiag%idv,idiag%iob,idiag%ich/),(/idv,iob,ich/))
        done = m==0
        if(done) exit

        ! We are searching for equal, so skip the i-th point if not equal.
        if(m<0) then
           ! if idiag%(idv,iob,ich) < (/idv,iob,ich/), move the lower range (lb) up
           ! to continue the search above i
           lb=i+1
        else
           ! if idiag%(idv,iob,ich) > (/idv,iob,ich/), move the upper range (ub) down
           ! to continue the search below i.
           ub=i-1
        endif

        if(ub<lb) exit      ! termionate the search
     enddo
  end associate

  if(done) then
     here_ => idiag
  endif

  call timer_fnl(myname_)
_EXIT_(myname_)
  return
end function locate_

function compare_(key1,key2) result (m)
  implicit none
  integer(kind=i_kind):: m
  integer(kind=i_kind),dimension(:),intent(in):: key1,key2

  integer(kind=i_kind):: n,i
  m=0
  n=min(size(key1),size(key2))
  do i=1,n
     if    (key1(i)<key2(i)) then
        m=-1; exit
     elseif(key1(i)>key2(i)) then
        m=+1; exit
     endif
  enddo
end function compare_

!-------------------
function obsnode_islocal_(anode) result(islocal_)
  use mpimod, only: mype
  use m_cvgridlookup, only: cvgridlookup_islocal
  implicit none
  logical:: islocal_
  type(_obsnode_),intent(in):: anode

  character(len=*),parameter:: myname_=myname//"::obsnode_islocal_"
_ENTRY_(myname_)
  islocal_=cvgridlookup_islocal(anode%elat,anode%elon,mype)
_EXIT_(myname_)
  return
end function obsnode_islocal_

function obsnode_isluse_(anode) result(isluse_)
  implicit none
  logical:: isluse_
  type(_obsnode_),intent(in):: anode

  character(len=*),parameter:: myname_=myname//"::obsnode_isluse_"
_ENTRY_(myname_)
  isluse_=anode%luse
_EXIT_(myname_)
  return
end function obsnode_isluse_

subroutine obsnode_setluse_(anode)
  use mpimod, only: mype
  use m_cvgridlookup, only: cvgridlookup_isluse
  implicit none
  type(_obsnode_),intent(inout):: anode

  character(len=*),parameter:: myname_=myname//"::obsnode_setluse_"
_ENTRY_(myname_)
  anode%luse=cvgridlookup_isluse(anode%elat, anode%elon, mype)
  !    call obstype_setluse(anode%luse, anode%elat, anode%elon, mype)
_EXIT_(myname_)
  return
end subroutine obsnode_setluse_

subroutine obsheader_read_(iunit,ii_bin,jj_type,lobs,jiter,miter,istat)
  implicit none
  integer(kind=i_kind),intent(in ):: iunit
  integer(kind=i_kind),intent(out):: ii_bin,jj_type,lobs,jiter,miter
  integer(kind=i_kind),intent(out):: istat

  character(len=*),parameter:: myname_=myname//"::obsheader_read_"
_ENTRY_(myname_)
  read(iunit,iostat=istat) ii_bin,jj_type,lobs,jiter,miter
_EXIT_(myname_)
  return
end subroutine obsheader_read_

subroutine obsheader_write_(iunit,ii_bin,jj_type,lobs,jiter,miter,istat)
  implicit none
  integer(kind=i_kind),intent(in ):: iunit
  integer(kind=i_kind),intent(in ):: ii_bin,jj_type,lobs,jiter,miter
  integer(kind=i_kind),intent(out):: istat

  character(len=*),parameter:: myname_=myname//"::obsheader_write_"
_ENTRY_(myname_)
  write(iunit,iostat=istat) ii_bin,jj_type,lobs,jiter,miter
_EXIT_(myname_)
  return
end subroutine obsheader_write_

subroutine obsnode_check_(who,anode)
!--  use jfunc, only: miter        ! for debugging
  implicit none
  character(len=*),intent(in):: who
  type(_obsnode_),intent(in):: anode

  logical:: equival
  character(len=256)::mywho

  mywho=who
    !_TRACEV_(who,'associated(anode%muse    ) =',associated(anode%muse    ))
    !_TRACEV_(who,'associated(anode%nldepart) =',associated(anode%nldepart))
    !_TRACEV_(who,'associated(anode%tldepart) =',associated(anode%tldepart))
    !_TRACEV_(who,'associated(anode%obssen  ) =',associated(anode%obssen  ))

  equival = associated(anode%nldepart) .eqv. associated(anode%muse    )
  if(equival) equival = associated(anode%tldepart) .eqv. associated(anode%nldepart)
  if(equival) equival = associated(anode%obssen  ) .eqv. associated(anode%tldepart)
  if(equival) equival = associated(anode%muse)

  ASSERT(equival)

!--  ASSERT(lbound(anode%muse    ,1)==1.and.ubound(anode%muse    ,1)==miter+1)
!--  ASSERT(lbound(anode%nldepart,1)==1.and.ubound(anode%nldepart,1)==miter+1)
!--  ASSERT(lbound(anode%tldepart,1)==1.and.ubound(anode%tldepart,1)==miter  )
!--  ASSERT(lbound(anode%obssen  ,1)==1.and.ubound(anode%obssen  ,1)==miter  )

  return
end subroutine obsnode_check_

function obsnode_alloc_(miter) result(anode_)
  implicit none
  type(_obsnode_), pointer   :: anode_
  integer(kind=i_kind), intent(in):: miter

  character(len=*),parameter:: myname_=myname//"::obsnode_alloc_"
_ENTRY_(myname_)
  allocate(anode_)
  anode_%next => null()

  allocate(anode_%muse    (miter+1), &
           anode_%nldepart(miter+1), &
           anode_%tldepart(miter  ), &
           anode_%obssen  (miter  )  )

  anode_%luse = .false.
  anode_%elat = 0._r_kind
  anode_%elon = 0._r_kind
  anode_%idv  =-1
  anode_%iob  =-1
  anode_%ich  =-1

  anode_%muse    (:)= .false.
  anode_%nldepart(:)=-huge(0._r_kind)
  anode_%tldepart(:)= 0._r_kind
  anode_%wgtjo      =-huge(0._r_kind)
  anode_%obssen  (:)= 0._r_kind

  call obsnode_check_(myname_,anode_)
_EXIT_(myname_)
  return
end function obsnode_alloc_

subroutine obsnode_init_(anode,idv,iob,ich,elat,elon,luse,miter)
  implicit none
  type(_obsnode_),intent(inout):: anode
  integer(kind=i_kind), intent(in):: idv,iob,ich
  real   (kind=r_kind), intent(in):: elat,elon
  logical, intent(in):: luse
  integer(kind=i_kind), intent(in):: miter

  character(len=*),parameter:: myname_=myname//"::obsnode_init_"
_ENTRY_(myname_)

  anode%next => null()
  anode%idv   = idv
  anode%iob   = iob
  anode%ich   = ich
  anode%elat  = elat
  anode%elon  = elon
  anode%luse  = luse



  anode%wgtjo      =-huge(0._r_kind)

  allocate(anode%muse    (miter+1), &
           anode%nldepart(miter+1), &
           anode%tldepart(miter  ), &
           anode%obssen  (miter  )  )

  anode%muse    (:)= .false.
  anode%nldepart(:)=-huge(0._r_kind)
  anode%tldepart(:)= 0._r_kind
  anode%obssen  (:)= 0._r_kind

  call obsnode_check_(myname_,anode)
_EXIT_(myname_)
  return
end subroutine obsnode_init_

subroutine anode_assert_(anode,idv,iob,ich,who,what)
  implicit none
  type(_obsnode_),intent(in):: anode
  integer(kind=i_kind), intent(in):: idv,iob,ich
  character(len=*),intent(in):: who
  character(len=*),intent(in):: what

  character(len=*),parameter:: myname_=myname//"::anode_assert_"
  logical:: valid
  character(len=:),allocatable:: what_
_ENTRY_(myname_)
  valid = &
        anode%idv == idv .and. &
        anode%iob == iob .and. &
        anode%ich == ich

  if(.not.valid) then
     what_=repeat(" ",len(trim(what)))
     call perr(who,trim(what)//", %(idv,iob,ich) =",(/anode%idv,anode%iob,anode%ich/))
     call perr(who,     what_//"   (idv,iob,ich) =",(/      idv,      iob,      ich/))
     call  die(who)
  endif

_EXIT_(myname_)
  return
end subroutine anode_assert_

subroutine obsnode_set_(anode, &
        idv,iob,ich,elat,elon,luse,wgtjo, &
        jiter,muse,nldepart,tldepart,obssen)
  implicit none
  type(_obsnode_),intent(inout):: anode
  integer(kind=i_kind),optional,intent(in):: idv,iob,ich
  real   (kind=r_kind),optional,intent(in):: elat,elon
  logical             ,optional,intent(in):: luse
  real   (kind=r_kind),optional,intent(in):: wgtjo

  integer(kind=i_kind),optional,intent(in):: jiter
  logical             ,optional,intent(in):: muse
  real   (kind=r_kind),optional,intent(in):: nldepart
  real   (kind=r_kind),optional,intent(in):: tldepart
  real   (kind=r_kind),optional,intent(in):: obssen

  character(len=*),parameter:: myname_=myname//"::obsnode_set_"
_ENTRY_(myname_)

  if(present(idv )) anode%idv =idv
  if(present(iob )) anode%iob =iob
  if(present(ich )) anode%ich =ich
  if(present(elat)) anode%elat=elat
  if(present(elon)) anode%elon=elon
  if(present(luse)) anode%luse=luse

  if(present(wgtjo )) anode%wgtjo =wgtjo


  if(present(jiter)) then
     if(present(muse  ).or.present(nldepart)) then
        ASSERT(jiter>=lbound(anode%muse    ,1))
        ASSERT(jiter<=ubound(anode%muse    ,1))
        ASSERT(jiter>=lbound(anode%nldepart,1))
        ASSERT(jiter<=ubound(anode%nldepart,1))
     endif
     if(present(obssen).or.present(tldepart)) then
        ASSERT(jiter>=lbound(anode%obssen  ,1))
        ASSERT(jiter<=ubound(anode%obssen  ,1))
        ASSERT(jiter>=lbound(anode%tldepart,1))
        ASSERT(jiter<=ubound(anode%tldepart,1))
     endif

     if(present(muse    )) anode%muse    (jiter) = muse
     if(present(nldepart)) anode%nldepart(jiter) = nldepart
     if(present(tldepart)) anode%tldepart(jiter) = tldepart
     if(present(obssen  )) anode%obssen  (jiter) = obssen
  endif

  !call obsnode_check_(myname_,anode_)
_EXIT_(myname_)
  return
end subroutine obsnode_set_

subroutine obsnode_get_(anode, &
        idv,iob,ich,elat,elon,luse,wgtjo, &
        jiter,muse,nldepart,tldepart,obssen)
  implicit none
  type(_obsnode_),intent(inout):: anode
  integer(kind=i_kind),optional,intent(out):: idv,iob,ich
  real   (kind=r_kind),optional,intent(out):: elat,elon
  logical             ,optional,intent(out):: luse
  real   (kind=r_kind),optional,intent(out):: wgtjo

  integer(kind=i_kind),optional,intent(in ):: jiter
  logical             ,optional,intent(out):: muse
  real(kind=r_kind)   ,optional,intent(out):: nldepart
  real(kind=r_kind)   ,optional,intent(out):: tldepart
  real(kind=r_kind)   ,optional,intent(out):: obssen

  character(len=*),parameter:: myname_=myname//"::obsnode_get_"
_ENTRY_(myname_)

  if(present(idv )) idv  = anode%idv
  if(present(iob )) iob  = anode%iob
  if(present(ich )) ich  = anode%ich
  if(present(elat)) elat = anode%elat
  if(present(elon)) elon = anode%elon
  if(present(luse)) luse = anode%luse

  if(present(wgtjo )) wgtjo  = anode%wgtjo

  if(present(jiter)) then
     if(present(muse  ).or.present(nldepart)) then
        ASSERT(jiter>=lbound(anode%muse    ,1))
        ASSERT(jiter<=ubound(anode%muse    ,1))
        ASSERT(jiter>=lbound(anode%nldepart,1))
        ASSERT(jiter<=ubound(anode%nldepart,1))
     endif
     if(present(obssen).or.present(tldepart)) then
        ASSERT(jiter>=lbound(anode%obssen  ,1))
        ASSERT(jiter<=ubound(anode%obssen  ,1))
        ASSERT(jiter>=lbound(anode%tldepart,1))
        ASSERT(jiter<=ubound(anode%tldepart,1))
     endif

     if(present(muse    )) muse     = anode%muse    (jiter)
     if(present(nldepart)) nldepart = anode%nldepart(jiter)
     if(present(tldepart)) tldepart = anode%tldepart(jiter)
     if(present(obssen  )) obssen   = anode%obssen  (jiter)
  endif

  !call obsnode_check_(myname_,anode_)
_EXIT_(myname_)
  return
end subroutine obsnode_get_

subroutine obsnode_read_(anode,iunit,kiter,istat,redistr)
  implicit none
  type(_obsnode_), intent(inout):: anode
  integer(kind=i_kind), intent(in   ):: iunit
  integer(kind=i_kind), intent(in   ):: kiter   ! input size
  integer(kind=i_kind), intent(out  ):: istat
  logical        , intent(in   ):: redistr

  character(len=*),parameter:: myname_=myname//'::obsnode_read_'
  integer(kind=i_kind):: ier
  !real(kind=r_kind),dimension(1:kiter):: zobssen
_ENTRY_(myname_)

  istat=0
  read(iunit,iostat=ier) anode%luse,anode%elat,anode%elon, &
                         anode%idv ,anode%iob ,anode%ich
  if(ier/=0) then
     call perr(myname_,'read(%luse,%elat,%elon,...), iostat =',ier)
     istat=-1
     _EXIT_(myname_)
     return
  endif

  istat=1
  if(redistr) then
     istat=0
     if(anode%luse) then
        if(obsnode_islocal_(anode)) istat=1
     endif
  endif

  if(istat==0) then
     read(iunit,iostat=ier)
     if(ier/=0) then
        call perr(myname_,'skipping read(%nchanl,%muse,...), iostat =',ier)
        istat=-2
        _EXIT_(myname_)
        return
     endif

  else
     read(iunit,iostat=ier)       &
         anode%muse    (1:kiter+1), &    ! = lmuse(1:kiter)
         anode%nldepart(1:kiter+1), &    ! = znldepart(1:kiter)
         anode%tldepart(1:kiter), &    ! = ztldepart(1:kiter)
         anode%wgtjo,             &    ! = zwgtjo
         anode%obssen  (1:kiter)       ! = zobssen(1:kiter)
     if(ier/=0) then
        call perr(myname_,'read(%nchanl,%muse,...), iostat =',ier)
        istat=-3
        _EXIT_(myname_)
        return
     endif

!     if    (lobsensfc.and..not.lsensrecompute) then
!        anode%obssen(jiter+1:miter  )=zobssen(jiter+1:miter  )
!     elseif(lobserver) then
!        anode%obssen(      1:jiter-1)=zobssen(      1:jiter-1)
!     else
!        anode%obssen(      1:miter  )=zobssen(      1:miter  )
!     endif
  endif

  call obsnode_check_(myname_,anode)
_EXIT_(myname_)
  return
end subroutine obsnode_read_

subroutine obsnode_write_(anode,iunit,jiter,istat)
  implicit none
  type(_obsnode_), intent(in   ):: anode
  integer(kind=i_kind), intent(in   ):: iunit
  integer(kind=i_kind), intent(in   ):: jiter   ! the output size
  integer(kind=i_kind), intent(inout):: istat

  character(len=*),parameter:: myname_=myname//'::obsnode_write_'
_ENTRY_(myname_)

  write(iunit,iostat=istat) anode%luse,anode%elat,anode%elon, &
                            anode%idv,anode%iob,anode%ich
  if(istat/=0) then
     call perr(myname_,'write(%luse,%elat,%elon,...), iostat =',istat)
     _EXIT_(myname_)
    return
  endif

  write(iunit,iostat=istat)     &
        anode%muse    (1:jiter+1),&
        anode%nldepart(1:jiter+1),&
        anode%tldepart(1:jiter),&
        anode%wgtjo,            &
        anode%obssen(1:jiter)

        if(istat/=0) then
           call perr(myname_,'write(%nchanl,%muse,...), iostat =',istat)
           _EXIT_(myname_)
           return
        endif
  call obsnode_check_(myname_,anode)
_EXIT_(myname_)
  return
end subroutine obsnode_write_

subroutine obsnode_dealloc_(anode,deep)
  implicit none
  type(_obsnode_),pointer,intent(inout):: anode
  logical,optional,intent(in):: deep

  character(len=*),parameter:: myname_=myname//'::obsnode_dealloc_'
  logical:: deep_
_ENTRY_(myname_)
  call obsnode_check_(myname_,anode)

  deep_=.false.
  if(present(deep)) deep_=deep
  ASSERT(associated(anode))

!  _TRACEV_(myname_,'if(deep_), deep_ =',deep_)
  if(deep_) then
!     _TRACEV_(myname_,'associated(anode%nldepart) =',associated(anode%nldepart))
     if(associated(anode%nldepart)) deallocate(anode%nldepart)
!     _TRACEV_(myname_,'associated(anode%tldepart) =',associated(anode%tldepart))
     if(associated(anode%tldepart)) deallocate(anode%tldepart)
!     _TRACEV_(myname_,'associated(anode%obssen  ) =',associated(anode%obssen  ))
     if(associated(anode%obssen  )) deallocate(anode%obssen  )
!     _TRACEV_(myname_,'associated(anode%muse    ) =',associated(anode%muse    ))
     if(associated(anode%muse    )) deallocate(anode%muse    )
  endif
    ! This is not a recursive dealloc_().  Therefore, the associated target of
    ! %next is not deallocated, but only %next itself is nullified.
!  _TRACEV_(myname_,'associated(%next) =',associated(anode%next))
  anode%next => null()
!  _TRACEV_(myname_,'associated(%next) =',associated(anode%next))
  deallocate(anode)
!  _TRACEV_(myname_,'associated(anode) =',associated(anode))
_EXIT_(myname_)
  return
end subroutine obsnode_dealloc_

subroutine obsnode_show_(anode,iob)
  use mpeu_util, only: stdout
  implicit none
  type(_obsnode_),intent(in):: anode
  integer(kind=i_kind),intent(in):: iob

  character(len=*),parameter:: myname_=myname//'::obsnode_show_'
_ENTRY_(myname_)
  write(stdout,'(2a,5i4,l4,2f8.2)') myname,":: iob,ity,%(idv,iob,ich,luse,elat,elon) =", &
        iob,0,anode%idv,anode%iob,anode%ich,anode%luse,anode%elat,anode%elon
  call obsnode_check_(myname_,anode)
_EXIT_(myname_)
  return
end subroutine obsnode_show_

end module m_obsdiagnode
