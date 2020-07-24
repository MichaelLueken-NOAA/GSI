module m_obsllist
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module m_obsllist
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2016-05-18
!
! abstract: class-module of linked-list of polymorphic obsnode.
!
! program history log:
!   2016-05-18  j guo   - added this document block for the initial polymorphic
!                         implementation.
!   2016-06-24  j.guo   - added support of using m_latlonrange to find a cluster
!                         latlonrange from (elat,elon) values of observations.
!   2016-07-25  j.guo   - added gettlddotprod, to accumulate obsnode tld-dot_produst
!   2016-09-19  j.guo   - added function lincr_() to extend []_lsize().
!   2017-08-26  G.Ge    - change allocate(headll%mold,mold=mold)
!                             to allocate(headll%mold,source=mold)
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
  use kinds , only: i_kind
  use mpeu_util, only: assert_,die,perr,warn,tell
  use m_obsnode, only: obsnode
  implicit none
  private

  public:: obsllist

  type obsllist
     private
     integer(i_kind):: n_alloc    =0

     integer(i_kind):: my_obstype =0
     class(obsnode),pointer:: mold => null()    ! a mold for the nodes

     class(obsnode),pointer:: head => null()    ! 
     class(obsnode),pointer:: tail => null()

     integer(i_kind):: l_alloc    =0            ! previous n_alloc, see showincr
  end type obsllist

  public:: obsllist_mold        ! get the mold of the obsllist
  interface obsllist_mold; module procedure lmold_; end interface

  public:: obsllist_reset       ! reset or finalize obsllist to its empty state.
  public:: obsllist_appendnode  ! append a node to obsllist

  interface obsllist_reset     ; module procedure lreset_     ; end interface
  interface obsllist_appendnode; module procedure lappendnode_; end interface

  public:: obsllist_rewind      ! rewind obsllist
  public:: obsllist_nextnode    ! move obsllist to its next node

  interface obsllist_rewind    ; module procedure lrewind_    ; end interface
  interface obsllist_nextnode  ; module procedure lnextnode_  ; end interface

  public:: obsllist_headnode    ! locate the head node of obsllist
  public:: obsllist_tailnode    ! locate the tail node of obsllist

  interface obsllist_headnode  ; module procedure lheadnode_  ; end interface
  interface obsllist_tailnode  ; module procedure ltailnode_  ; end interface

  public:: obsllist_lsize        ! get the size of a llist
  public:: obsllist_lcount       ! get the size of a llist
  public:: obsllist_lsort        ! sort nodes according to their keys
  public:: obsllist_write        ! output a llist to a file unit
  public:: obsllist_read         ! input from a file created by _write()
  public:: obsllist_checksum     ! size consistency checking
  public:: obsllist_summary      ! show some information about the llist

  interface obsllist_lsize  ; module procedure lsize_, &
                                               lincr_   ; end interface
  interface obsllist_lcount ; module procedure lcount_  ; end interface
  interface obsllist_lsort  ; module procedure lsort_   ; end interface
  interface obsllist_write  ; module procedure lwrite_  ; end interface
  interface obsllist_read   ; module procedure lread_   ; end interface
  interface obsllist_checksum; module procedure &
     lchecksum_, &
     lchecksum1_ ; end interface
  interface obsllist_summary; module procedure lsummary_; end interface

  public:: obsllist_gettlddotprod     ! get "LHS" (dot-product of (:)%diags%tldepar, plus count)
  interface obsllist_gettlddotprod ; module procedure ltlddotprod_ ; end interface

  character(len=*),parameter:: myname="m_obsllist"

#include "myassert.H"
#include "mytrace.H"
contains

subroutine ltlddotprod_(headll,jiter,tlddp,nnode,nob)
!-- get "lhs" of the given linked-list
  use kinds, only: i_kind,r_kind
  use m_obsnode, only: obsnode_next, obsnode_isluse
  implicit none
  type(obsllist),target, intent(in):: headll    ! a linked-list 
  integer(kind=i_kind) , intent(in):: jiter     ! for this iteration
  real   (kind=r_kind) , intent(inout):: tlddp  ! dot_product((:)%tld)
  integer(kind=i_kind) , optional, intent(inout):: nnode ! node count
  integer(kind=i_kind) , optional, intent(inout):: nob   ! obs. count

  class(obsnode),pointer:: inode
  inode => lheadnode_(headll)
  do while(associated(inode))
     if(obsnode_isluse(inode)) then
        call inode%gettlddp(jiter,tlddp,nob=nob)
        if(present(nnode)) nnode=nnode+1
     endif
     inode => obsnode_next(inode)
  enddo
end subroutine ltlddotprod_

function lmold_(headll) result(ptr_)
  implicit none
  class(obsnode),pointer:: ptr_
  type(obsllist),target,intent(in):: headll
  ptr_ => null()
  if(associated(headll%mold)) ptr_ => headll%mold
end function lmold_

!--------------------------- will go to m_obsllist ----------------------
subroutine lrewind_(headll)
  implicit none
  type(obsllist),target,intent(inout):: headll
  headll%tail => null()
end subroutine lrewind_

function lnextnode_(headll) result(here_)
  use m_obsnode, only: obsnode_next
  implicit none
  class(obsnode),pointer:: here_
  type(obsllist),target,intent(inout):: headll

  if(associated(headll%tail)) then
        ! when not the first time lnextnode_(), after call lrewind_()
     headll%tail => obsnode_next(headll%tail)
  else
        ! When the first time lnextnode_(), after call lrewind_()
     headll%tail => lheadnode_(headll)
  endif
  here_ => headll%tail
end function lnextnode_

function lheadnode_(headll) result(here_)
  implicit none
  class(obsnode),pointer:: here_
  type(obsllist),target,intent(in):: headll
  here_ => headll%head
end function lheadnode_

function ltailnode_(headll) result(here_)
  implicit none
  class(obsnode),pointer:: here_
  type(obsllist),target,intent(in):: headll
  here_ => headll%tail
end function ltailnode_

function lsize_(headll)
  implicit none
  integer(i_kind):: lsize_
  type(obsllist),intent(in):: headll
  lsize_=headll%n_alloc
end function lsize_
function lincr_(headll,incr)
  implicit none
  integer(i_kind):: lincr_
  type(obsllist),intent(inout):: headll
  logical,intent(in):: incr
  lincr_=headll%n_alloc
  if(incr) then
     lincr_=lincr_-headll%l_alloc
     headll%l_alloc=headll%n_alloc
  endif
end function lincr_

subroutine lreset_(headll,mold,stat)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    lreset_
!   prgmmr:      J. Guo
!
! abstract: reset a linked-list to empty.
!
! program history log:
!   2015-01-12  guo     - reset headll for a generic obsnode
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
  use m_obsnode, only: obsnode_next
  use m_obsnode, only: obsnode_clean
  use m_obsnode, only: obsnode_type => obsnode_mytype
  implicit none
  type(obsllist), intent(inout):: headll
  class(obsnode), intent(in   ):: mold
  integer(i_kind),optional,intent(out):: stat

  character(len=*),parameter:: myname_=myname//"::lreset_"
  character(len=:),allocatable:: mymold_
  integer(i_kind):: n
  integer(i_kind):: ier
_ENTRY_(myname_)

  if(present(stat)) stat=0

  call obsnode_clean(headll%head,deep=.true.,depth=n,stat=ier)
  if(ier/=0.or.n/=0) then
     call perr(myname_,'obsnode_clean(.deep.), stat =',ier)
     call perr(myname_,'                      depth =',n)
     call perr(myname_,'              lsize(headll) =',lsize_(headll))
     call perr(myname_,'       headll%head%mytype() =',obsnode_type(headll%head))
     call perr(myname_,'       headll%mold%mytype() =',obsnode_type(headll%mold))
     if(.not.present(stat)) call die(myname_)
     stat=ier
     _EXIT_(myname_)
     return
  endif

  call nodedestroy_(headll%head)

  headll%n_alloc = 0
  headll%l_alloc = 0
  headll%head    => null()
  headll%tail    => null()

  if(associated(headll%mold)) then
     mymold_ = obsnode_type(headll%mold)
     deallocate(headll%mold,stat=ier)
     if(ier/=0) then
        call perr(myname_,'deallocate(headll%mold), stat =',ier)
        call perr(myname_,'    obsnode_type(headll%mold) =',mymold_)
        if(.not.present(stat)) call die(myname_)
        stat=ier
        _EXIT_(myname_)
        return
     endif
  endif

  allocate(headll%mold, mold=mold)
_EXIT_(myname_)
  return
end subroutine lreset_

subroutine lappendnode_(headll,targetnode)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    lappendnode_
!   prgmmr:      J. Guo
!
! abstract: append a node to the given linked-list
!
! program history log:
!   2015-01-12  guo     - constructed for generic _obsnode_
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
        ! Link the next node of the list to the given targetnode.  The return
        ! result is a pointer associated to the same targetnode.
  use m_obsnode, only: obsnode_append
  implicit none
  type(obsllist), intent(inout):: headll
  !class(obsnode), target, intent(in):: targetnode
  class(obsnode), pointer, intent(in):: targetnode

  character(len=*),parameter:: myname_=myname//'::lappendnode_'
!_ENTRY_(myname_)
  ASSERT(associated(targetnode))

  if(.not.associated(headll%head)) then
                ! this is a fresh starting -node- for this linked-list ...
     call obsnode_append(headll%head,targetnode)
     headll%tail => headll%head
     headll%n_alloc = 1

  else
     ASSERT(associated(headll%tail))
     ASSERT(.not.associated(headll%tail,targetnode))

     call obsnode_append(headll%tail,targetnode)
     headll%n_alloc = headll%n_alloc + 1

  endif

!_EXIT_(myname_)
  return
end subroutine lappendnode_

!--------------------------- will go to m_obsllistio ----------------------
subroutine lread_(headll,iunit,redistr,diaglookup,jtype)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    lread_
!   prgmmr:      todling
!   prgmmr:      J. Guo
!
! abstract: read obs-specific data structure from file.
!
! program history log:
!   2007-10-03  todling - (original read_obsdiags::read_${obstype}head_()
!   2008-12-08  todling - update to May08 version
!   2015-01-12  guo     - restructured for generic _obsnode_, with redistributions
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

  !use obsmod, only: obs_diags
  use m_obsdiagnode, only: obs_diags
  use m_obsnode, only: obsnode_read
  use m_obsnode, only: obsnode_setluse
  implicit none
  type(obsllist), intent(inout):: headll
  integer(i_kind), intent(in  ):: iunit
  logical        , intent(in  ):: redistr
  type(obs_diags), intent(in  ):: diaglookup
  integer(i_kind), intent(in  ):: jtype

  character(len=*),parameter:: myname_=myname//"::lread_"
  class(obsnode),pointer :: anode => null()
  integer(i_kind) :: kk,istat,mobs,jread
_ENTRY_(myname_)

! read in an obs-specific header of the next block
! >>>>>----------------------------
! obsheader is the information about an obs-block, where an obs-block
! a collection of nodes of the same _obsnode_ type,
! !-- not about the corresponding linked-list.

  ASSERT(associated(headll%mold))

  call obsheader_read_(headll%mold,iunit,mobs,jread,istat)

  if(istat/=0) then
     call perr(myname_,'%obsheader_read_(mobs,jread), istat =',istat)
     call perr(myname_,'                              iunit =',iunit)
     call  die(myname_)
  endif

  if(jtype/=jread) then
     call perr(myname_,'unexpected record type, jread =',jread)
     call perr(myname_,'              expecting jtype =',jtype)
     call perr(myname_,'                         mobs =',mobs)
     call perr(myname_,'                        iunit =',iunit)
     call  die(myname_)
  end if
!   ----------------------------<<<<<

  if(mobs==0) then
     ! No more record to read
     _EXIT_(myname_)
     return
  endif

  !-- construct an anode
  anode => alloc_nodecreate_(mold=headll%mold)
  do kk=1,mobs
     !-- initialize anode from a file (iunit)
     call obsnode_read(anode,iunit,istat,redistr=redistr,diaglookup=diaglookup)
     if(istat<0) then
        call perr(myname_,'obsnode_read(), istat =',istat)
        call perr(myname_,'                iunit =',iunit)
        call perr(myname_,'                   kk =',kk)
        call perr(myname_,'                 mobs =',mobs)
        call perr(myname_,'              redistr =',redistr)
        call perr(myname_,'                jtype =',jtype)
        call  die(myname_)
     endif

     if(istat==0) cycle

     !-- If this anode is to be kept ...
     if(redistr) then
        ! recompute its %luse and %hop for the redistributed grid partition,

        call obsnode_setluse(anode)     ! reset %luse for subdomain ownership
        call anode%sethop()             ! recompute %hop for the new grid
     endif

                !-- keep this obsnode in its linked-list, obsllist := obsdiags(jtype,ibin)
     call lappendnode_(headll,targetnode=anode)

                !-- Drop the earlier object, contruct a new anode.
                !-- No deep deallocation is needed for anode, since its
                !-- associated target has been passed to headll
     anode => null()
     anode => alloc_nodecreate_(mold=headll%mold)

  enddo  ! < mobs >

  call nodedestroy_(anode)  ! Clean up the working-space an_onsnode
    
_EXIT_(myname_)
  return
end subroutine lread_

subroutine lwrite_(headll,iunit,luseonly,jtype,luserange)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    lwrite_
!   prgmmr:      todling
!   prgmmr:      J. Guo
!
! abstract: write obs-specific data structure from file.
!
! program history log:
!   2007-10-03  todling - (original write_obsdiags::write_${obstype}head_()
!   2008-12-08  todling - update to May08 version
!   2015-01-12  guo     - restructured for generic _obsnode_, with redistributions
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

  use m_obsnode, only: obsnode_next
  use m_obsnode, only: obsnode_isluse
  use m_obsnode, only: obsnode_write
  use m_latlonrange, only: latlonrange
  use m_latlonrange, only: latlonrange_enclose
  implicit none
  type(obsllist), intent(in):: headll
  integer(i_kind ), intent(in):: iunit       ! unit for output
  logical         , intent(in):: luseonly
  integer(i_kind ), intent(in):: jtype
  type(latlonrange),optional,intent(inout):: luserange

  character(len=*),parameter:: myname_=myname//"::lwrite_"
  class(obsnode), pointer :: inode => null()
  integer(i_kind) :: istat
  integer(i_kind) :: mobs,lobs,iobs,kobs
  logical:: isluse_
_ENTRY_(myname_)

! if(jtype/=iobstype) then
!    call perr(myname_,'unexpected record type, jtype =',jtype)
!    call perr(myname_,'           expecting iobstype =',iobstype)
!    call perr(myname_,'                        iunit =',iunit)
!    call  die(myname_)
! end if

! read in an obs-specific header of the next block
! >>>>>----------------------------
! !-- A header is about a collection of nodes of the same obsnode type,
! !-- not about the corresponding linked-list.

  ASSERT(associated(headll%mold))

  lobs = lcount_(headll,luseonly=luseonly)      ! actual count of write
  mobs = lobs
  if(.not.luseonly) mobs = lsize_(headll)       ! actual count of nodes

  call obsheader_write_(headll%mold,iunit,lobs,jtype,istat)

  if(istat/=0) then
     call perr(myname_,'obsheader_write_(), istat =',istat)
     call perr(myname_,'                    iunit =',iunit)
     call perr(myname_,'                    jtype =',jtype)
     call perr(myname_,'        no. node of write =',lobs)
     call perr(myname_,'         no. node of data =',mobs)
     call perr(myname_,'                 luseonly =',luseonly)
     call  die(myname_)
  endif
! ----------------------------<<<<<

  if(lobs==0) then
     ! No more record to write
     _EXIT_(myname_)
     return
  endif

!-- looping over the linked-list for every obsnode,

  inode => lheadnode_(headll)
  iobs=0
  kobs=0
  do while(associated(inode))
     iobs=iobs+1
     isluse_=obsnode_isluse(inode)
     if(isluse_ .or. .not.luseonly) then
        if(isluse_.and.present(luserange)) &
                call latlonrange_enclose(luserange,inode%elat,inode%elon)
        kobs=kobs+1
        call obsnode_write(inode,iunit,istat)
        if(istat/=0) then
           call perr(myname_,' obsnode_write(), istat =',istat)
           call perr(myname_,'                  iunit =',iunit)
           call perr(myname_,'                  jtype =',jtype)
           call perr(myname_,'current-luse-node, kobs =',kobs)
           call perr(myname_,' current-all-node, iobs =',iobs)
           call perr(myname_,'  total-luse-node-count =',lobs)
           call perr(myname_,'   total-all-node-count =',mobs)
           call perr(myname_,'               luseonly =',luseonly)
           call  die(myname_)
        endif
     endif
     inode => obsnode_next(inode)
  enddo

  ASSERT(iobs==mobs)
  ASSERT(kobs==lobs)
_EXIT_(myname_)
  return
end subroutine lwrite_

subroutine lchecksum_(headll,itype,ibin,leadnode,sorted)
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
  use mpeu_util, only: stdout,stdout_lead
  implicit none
  type(obsllist), intent(in):: headll
  integer(kind=i_kind),optional,intent(in ):: itype,ibin
  class(obsnode),optional,pointer, intent(in):: leadnode
  logical             ,optional,intent(out):: sorted

  character(len=*),parameter:: myname_=myname//"::lchecksum_"
  integer(kind=i_kind):: lrecount
  integer(kind=i_kind):: jtype,jbin
  integer(kind=i_kind):: nuse,nooo,ndup,ksum(2)
_ENTRY_(myname_)
  lrecount=lcount_(headll,recount=.true.,nuse=nuse,nooo=nooo,ndup=ndup,ksum=ksum,leadnode=leadnode)
  if(present(sorted)) sorted = nooo==0.and.ndup==0

  jtype=itype
  jbin =ibin
_EXIT_(myname_)
  return
end subroutine lchecksum_
subroutine lchecksum1_(headll,itype)
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
  type(obsllist), dimension(:),intent(in):: headll
  integer(kind=i_kind),optional ,intent(in):: itype

  character(len=*),parameter:: myname_=myname//"::lchecksum1_"
  integer(kind=i_kind):: i
_ENTRY_(myname_)
  do i=1,size(headll)
     call lchecksum_(headll(i),itype=itype,ibin=i)
  enddo
_EXIT_(myname_)
  return
end subroutine lchecksum1_

subroutine lsummary_(headll,verbose)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    lsummary_
!   prgmmr:      J. Guo
!
! abstract: summarize for the contents of a linked-list.
!
! program history log:
!   2015-01-12  guo     - constructed for generic _obsnode_
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
  use m_obsnode, only: obsnode_next
  use m_obsnode, only: obsnode_show
  implicit none
  type(obsllist), intent(in):: headll
  logical,optional, intent(in):: verbose

  character(len=*),parameter:: myname_=myname//"::lsummary_"
  class(obsnode), pointer:: inode
  integer(i_kind):: iobs_
  logical:: verbose_
  verbose_=.false.
  if(present(verbose)) verbose_=verbose
_ENTRY_(myname_)
  !call tell(myname_,' headllist%n_alloc =',headll%n_alloc)

  if(verbose_) then
     iobs_ = 0
     inode => lheadnode_(headll)
     do while(associated(inode))
        iobs_=iobs_+1
        call obsnode_show(inode,iobs_)
        inode => obsnode_next(inode)
     enddo
  endif
_EXIT_(myname_)
  return
end subroutine lsummary_

function lcount_(headll,luseonly,recount,nuse,nooo,ndup,ksum,leadnode) result(lobs_)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    lcount_
!   prgmmr:      J. Guo
!
! abstract: inquire for the size information about the linked-list
!
! program history log:
!   2015-01-12  guo     - constructed for generic _obsnode_
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
  use m_obsnode, only: obsnode_next
  use m_obsnode, only: obsnode_isluse
  implicit none
  integer(kind=i_kind):: lobs_
  type(obsllist), target, intent(in):: headll
  logical,optional,intent(in):: luseonly        ! count only luse data
  logical,optional,intent(in):: recount
  integer(kind=i_kind),optional,intent(out):: nuse      ! no. luse
  integer(kind=i_kind),optional,intent(out):: nooo      ! no. out-of-orders
  integer(kind=i_kind),optional,intent(out):: ndup      ! no. duplicates
  integer(kind=i_kind),optional,dimension(:),intent(out):: ksum ! key value sum
  class(obsnode), pointer, optional, intent(in):: leadnode

  character(len=*),parameter:: myname_=myname//"::lcount_"
  class(obsnode), pointer:: inode
  integer(i_kind):: nuse_
  integer(kind=i_kind),dimension(2) :: kprev
  logical:: luseonly_,recount_,checksum_
_ENTRY_(myname_)

  luseonly_=.false.
  if(present(luseonly)) luseonly_=luseonly
  recount_ =.false.
  if(present(recount )) recount_ =recount
  if(present(leadnode)) recount_ =.true.

  checksum_= present(nuse).or.present(nooo).or.present(ndup).or.present(ksum)
  if(.not.recount_) recount_ = checksum_

  if(present(ksum)) then
     ALWAYS_ASSERT(size(ksum)==size(kprev))
  endif

  if(.not.(luseonly_.or.recount_)) then
     lobs_=headll%n_alloc

  else
     lobs_ = 0
     nuse_ = 0

     if(checksum_) call checksum_init_(kprev,nooo=nooo,ndup=ndup,ksum=ksum)

     inode => lheadnode_(headll)
     do while(associated(inode))
        if(obsnode_isluse(inode)) nuse_=nuse_+1
        if(.not.luseonly_ .or. obsnode_isluse(inode)) lobs_=lobs_+1

        if(checksum_) call checksum_add_(kprev, &
           knext=(/inode%idv,inode%iob/),nooo=nooo,ndup=ndup,ksum=ksum)

        inode => obsnode_next(inode)
     enddo
     if(present(nuse)) nuse=nuse_
  endif

_EXIT_(myname_)
  return
end function lcount_

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

  integer(kind=i_kind):: k
  k=compare_(kprev,knext)
  if(present(nooo).and.k> 0) nooo=nooo+1
  if(present(ndup).and.k==0) ndup=ndup+1
  if(present(ksum)) ksum(:)=ksum(:)+knext(:)
  kprev(:)=knext(:)
end subroutine checksum_add_

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

subroutine lsort_(headll,itype,ibin)
!       lsort_: node-sort diagll, to line-up nodes according to their keys
!_timer_use_
!  use timermod , only: timer_ini,timer_fnl
  use mpeu_util, only: indexset
  use mpeu_util, only: indexsort
  use m_obsnode, only: obsnode_next
  !use mpeu_util, only: die
  implicit none
  type(obsllist), intent(inout):: headll
  integer(kind=i_kind),optional,intent(in):: itype,ibin

  character(len=*),parameter:: myname_=myname//'::lsort_'
  class(obsnode),pointer:: pnode
  integer(kind=i_kind),allocatable,dimension(:):: indx,idv_,iob_
  integer(kind=i_kind):: i,n
  logical:: sorted

  type fptr_of_obsnode
     class(obsnode),pointer:: ptr
  end type fptr_of_obsnode
  type(fptr_of_obsnode),allocatable,dimension(:):: lookup
_ENTRY_(myname_)
!_timer_on_(myname_)
!  call timer_ini(myname_)

  call lchecksum_(headll,itype=itype,ibin=ibin,sorted=sorted)
  if(sorted) then
     _EXIT_(myname_)
     return
  endif

  n=lsize_(headll)

  allocate(lookup(n))
  allocate(indx(n),idv_(n),iob_(n))

        ! Loop over the linked-list, to get keys.
  i=0
  pnode => lheadnode_(headll)
  do while(associated(pnode))
     i=i+1
     if(i<=n) then
        lookup(i)%ptr => pnode
          idv_(i)     =  pnode%idv
          iob_(i)     =  pnode%iob
     endif
     pnode => obsnode_next(pnode)
  enddo
 
  ASSERT(i==n)

        ! sort %lookup(1:n), by its (idv,iob) values
  call indexset (indx)
  call indexsort(indx,iob_)
  call indexsort(indx,idv_)
  lookup(1:n) = lookup(indx(1:n))

  deallocate(indx,idv_,iob_)

        ! Rebuild the linked-list from lookup(1:n)%ptr
  headll%n_alloc = 0
  headll%head    => null()
  headll%tail    => null()

        ! rebuild the list according to the sorted table
  do i=1,n
     call lappendnode_(headll,lookup(i)%ptr)
  enddo
  ASSERT(n==headll%n_alloc)
  if(associated(headll%tail)) then
     ASSERT(.not.associated(headll%tail%llpoint))
  endif

        ! discard the table
  deallocate(lookup)

  call lchecksum_(headll,itype=itype,ibin=ibin,sorted=sorted)
  if(.not.sorted) then
     call perr(myname_,'failed post-sorting lchecksum_(), sorted =',sorted)
     if(present(itype)) &
     call perr(myname_,'                                   itype =',itype)
     if(present(ibin )) &
     call perr(myname_,'                                    ibin =',ibin)
     call  die(myname_)
  endif

!  call timer_fnl(myname_)
!_timer_off_(myname_)
_EXIT_(myname_)
  return
end subroutine lsort_

function alloc_nodecreate_(mold) result(ptr_)
!-- allocate() + init()
  implicit none
  class(obsnode),pointer:: ptr_
  class(obsnode),target,intent(in):: mold
  allocate(ptr_,mold=mold)
  call ptr_%init()
  return
end function alloc_nodecreate_

subroutine nodedestroy_(node)
!-- clean() + deallocate()
  use m_obsnode, only: obsnode_type => obsnode_mytype
  implicit none
  class(obsnode),pointer,intent(inout):: node
  character(len=*),parameter:: myname_=myname//'::nodedestroy_'
  integer(i_kind):: ier
  if(associated(node)) then
     call node%clean()
     deallocate(node,stat=ier)
     if(ier/=0) then
        call perr(myname_,'can not deallocate(node), stat =',ier)
        call perr(myname_,'            obsnode_type(node) =',obsnode_type(node))
        call die(myname_)
     endif
  endif
  return
end subroutine nodedestroy_

subroutine obsheader_read_(anode,iunit,iobs,itype,istat)
!-- read header of some type
  use m_obsnode, only: obsnode
  implicit none
  class(obsnode) ,intent(in ):: anode
  integer(i_kind),intent(in ):: iunit
  integer(i_kind),intent(out):: iobs,itype
  integer(i_kind),intent(out):: istat
  call anode%headerread(iunit,iobs,itype,istat)
end subroutine obsheader_read_

subroutine obsheader_write_(anode,junit,mobs,mtype,istat)
!-- write header of some type
  use m_obsnode, only: obsnode
  implicit none
  class(obsnode) ,intent(in ):: anode
  integer(i_kind),intent(in ):: junit
  integer(i_kind),intent(in ):: mobs,mtype
  integer(i_kind),intent(out):: istat
  call anode%headerwrite(junit,mobs,mtype,istat)
end subroutine obsheader_write_
end module m_obsllist
