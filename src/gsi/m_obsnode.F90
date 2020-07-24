module m_obsnode
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module m_obsnode
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2015-01-12
!
! abstract: basic obsnode functionalities interfacing the distributed grid
!
! program history log:
!   2015-01-12  j guo   - added this document block.
!   2016-05-18  j guo   - finished its initial polymorphic implementation.
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

  use kinds, only: i_kind,r_kind
  use mpeu_util, only: tell,perr,die
  use mpeu_util, only: assert_
  use m_obsdiagnode, only: obs_diag
  use m_obsdiagnode, only: obs_diags
  implicit none
  private       ! except
  public:: obsnode              ! data structure

  type, abstract:: obsnode
     ! private

     ! - Being not "private", type(obsnode) allowes its type extentions
     !   to access its components without additional interfaces.
     ! - On the other hand, by turning private on, one can use the
     !   compiler to report where the components of this type have been
     !   used.

     class(obsnode),pointer :: llpoint => null()

     logical         :: luse =.false.         ! flag indicating if ob is used in pen.
     real(r_kind)    :: time = 0._r_kind      ! observation time in sec, relative to the time window
     real(r_kind)    :: elat = 0._r_kind      ! earth lat-lon for redistribution
     real(r_kind)    :: elon = 0._r_kind      ! earth lat-lon for redistribution

     integer(i_kind) :: idv  =-1              ! device id
     integer(i_kind) :: iob  =-1              ! initial obs sequential id

#ifdef _TO_DO_
     integer(i_kind):: nprof   ! count of corresponding profile locations
     integer(i_kind):: idspl   ! cross referencing index to profile locations
                                ! given i-th observation, corresponding profile
                                ! is block ([]%idspl+1 : []%idspl+[]%nprof)
#endif
  contains

        !----------- overrideable procedures -----------------------------------
     procedure, nopass:: headerread  => obsheader_read_         ! read a header
     procedure, nopass:: headerwrite => obsheader_write_        ! write a header

     procedure:: init  => init_                                  ! initialize a node
     procedure:: clean => clean_                                 ! clean a node

        !----------- procedures must be defined by extensions ------------------
     procedure(intrfc_mytype_ ),nopass,deferred:: mytype     ! return my type name
     procedure(intrfc_sethop_ ), deferred:: sethop     ! re-construct H
     procedure(intrfc_xread_  ), deferred:: xread      ! read extensions
     procedure(intrfc_xwrite_ ), deferred:: xwrite     ! write extensions
     procedure(intrfc_isvalid_), deferred:: isvalid    ! validate extensions
 
     procedure(intrfc_gettlddp_), deferred:: gettlddp  ! (tlddp,nob)=(sum(%tld*%tld),sum(1)
        !--------- non_overrideable procedures are implemented statically ------
  end type obsnode

!-- module procedures, such as base-specific operations

  public:: obsnode_clean
  interface obsnode_clean; module procedure deepclean_; end interface

        ! nodes operations
  public:: obsnode_next         ! nextnode => obsnode_next (thisnode)
  public:: obsnode_append       ! call obsnode_append(thisnode,targetnode)

  interface obsnode_next  ; module procedure next_  ; end interface
  interface obsnode_append; module procedure append_; end interface

        ! Getters-and-setters
  public:: obsnode_islocal      ! is anode local? -- obsnode_islocal(anode)
  public:: obsnode_isluse       ! is anode luse?  -- obsnode_isluse(anode)
  public:: obsnode_setluse      ! set anode%luse. -- call obsnode_setluse(anode)

  interface obsnode_islocal; module procedure islocal_ ; end interface
  interface obsnode_isluse ; module procedure isluse_  ; end interface
  interface obsnode_setluse; module procedure setluse_ ; end interface

!-- module procedures, requiring base-specific operations

        ! reader-and-writer
  public:: obsnode_read         ! call obsnode_read(anode, ...)
  public:: obsnode_write        ! call obsnode_write(anode, ...)

  interface obsnode_read   ; module procedure read_   ; end interface
  interface obsnode_write  ; module procedure write_  ; end interface

  public:: obsnode_show         ! call obsnode_show(anode)
  interface obsnode_show   ; module procedure show_   ; end interface

  public:: obsnode_mytype       ! call obsnode_type(anode)
  interface obsnode_mytype ; module procedure nodetype_   ; end interface

  abstract interface
     subroutine intrfc_xread_(anode,iunit,istat,diaglookup,skip)
        use kinds,only: i_kind
        use m_obsdiagnode, only: obs_diags
        import:: obsnode
        implicit none
        class(obsnode), intent(inout):: anode
        integer(kind=i_kind), intent(in ):: iunit
        integer(kind=i_kind), intent(out):: istat
        type(obs_diags)     , intent(in ):: diaglookup
        logical,optional    , intent(in ):: skip
     end subroutine intrfc_xread_
  end interface

  abstract interface
     subroutine intrfc_xwrite_(anode,junit,jstat)
        use kinds,only: i_kind
        import:: obsnode
        implicit none
        class(obsnode), intent(in):: anode
        integer(kind=i_kind), intent(in ):: junit
        integer(kind=i_kind), intent(out):: jstat
     end subroutine intrfc_xwrite_
  end interface

  abstract interface
     function intrfc_isvalid_(anode) result(isvalid_)
        import:: obsnode
        implicit none
        logical:: isvalid_
        class(obsnode), intent(in):: anode
     end function intrfc_isvalid_
  end interface

  abstract interface
     subroutine intrfc_sethop_(anode)
        use kinds, only: r_kind
        import:: obsnode
        implicit none
        class(obsnode), intent(inout):: anode
     end subroutine intrfc_sethop_
  end interface

  abstract interface
     function intrfc_mytype_()
        import:: obsnode
        implicit none
        character(len=:),allocatable:: intrfc_mytype_
     end function intrfc_mytype_
  end interface

  abstract interface
     pure subroutine intrfc_gettlddp_(anode,jiter,tlddp,nob)
        use kinds, only: i_kind,r_kind
        import:: obsnode
        implicit none
        class(obsnode),intent(in):: anode
        integer(kind=i_kind),intent(in):: jiter
        real(kind=r_kind),intent(inout):: tlddp
        integer(kind=i_kind),optional,intent(inout):: nob
     end subroutine intrfc_gettlddp_
  end interface

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname='m_obsnode'

#include "mytrace.H"
#include "myassert.H"

contains
function next_(anode) result(here_)
!-- associate to thisnode%llpoint.
  implicit none
  class(obsnode),pointer:: here_
  class(obsnode),target,intent(in):: anode

  character(len=*),parameter :: myname_=myname//'::next_'
_ENTRY_(myname_)
        !!! trying to go next on a null reference is a serious logical error.
  here_ => anode%llpoint
_EXIT_(myname_)
  return
end function next_

subroutine append_(thisnode,targetnode,follow)
!-- append targetnode to thisnode%llpoint, or thisnode if .not.associated(thisnode)
  implicit none
  class(obsnode),pointer ,intent(inout):: thisnode
  class(obsnode),pointer ,intent(in   ):: targetnode
  logical       ,optional,intent(in):: follow  ! Follow targetnode%llpoint to its last node.
                                               ! The default is to nullify(thisnode%llpoint)

  character(len=*),parameter:: myname_=myname//"::append_"
  logical:: follow_
_ENTRY_(myname_)
  ASSERT(associated(targetnode))     ! verify for any exception.

  follow_=.false.
  if(present(follow)) follow_=follow

  if(.not.associated(thisnode)) then
     thisnode => targetnode              ! as the first node

  else
     thisnode%llpoint => targetnode      ! as an additional node
     thisnode => thisnode%llpoint

  endif

  if(follow_) then
     ! Follow thisnode to thisnode%llpoint, till its end, as targetnode is a
     ! valid linked-list.  The risk is the possibility of some circular
     ! association, evenif both linked-lists, thisnode and targetnode are given
     ! clean.

     do while(associated(thisnode%llpoint))
        ASSERT(.not.associated(thisnode%llpoint,targetnode))
        ! This assertion tries to identify possible circular association between
        ! linked-list::thisnode and linked-list::targetnode.

        thisnode => thisnode%llpoint
     enddo

  else
     ! Nullify(thisnode%llpoint) to avoid any possibility of circular
     ! association.  Note this action will touch the input target argument
     ! (targetnode) indirectly.

     thisnode%llpoint => null()
  endif
_EXIT_(myname_)
  return
end subroutine append_

function islocal_(anode)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    islocal_
!   prgmmr:      J. Guo
!
! abstract: check if this node is for the local grid partition.
!
! program history log:
!   2015-01-12  guo     - constructed for generic obsnode
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
  use mpimod, only: mype
  use m_cvgridlookup, only: cvgridlookup_islocal
  implicit none
  logical:: islocal_
  class(obsnode),intent(in):: anode
  character(len=*),parameter:: myname_=myname//'::islocal_'
_ENTRY_(myname_)
  islocal_=cvgridlookup_islocal(anode%elat,anode%elon,mype)
_EXIT_(myname_)
  return
end function islocal_

function isluse_(anode)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    isluse_
!   prgmmr:      J. Guo
!
! abstract: check the %luse value of this node
!
! program history log:
!   2015-01-12  guo     - constructed for generic obsnode
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
  logical:: isluse_
  class(obsnode),intent(in):: anode
  character(len=*),parameter:: myname_=myname//'::isluse_'
_ENTRY_(myname_)
  isluse_=anode%luse
_EXIT_(myname_)
  return
end function isluse_

subroutine setluse_(anode)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    lsummary_
!   prgmmr:      J. Guo
!
! abstract: set %luse value for locally-owned node.
!
! program history log:
!   2015-01-12  guo     - constructed for generic obsnode
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
  use mpimod, only: mype
  use m_cvgridlookup, only: cvgridlookup_isluse
  implicit none
  class(obsnode),intent(inout):: anode
  character(len=*),parameter:: myname_=myname//'::setluse_'
_ENTRY_(myname_)
  anode%luse = cvgridlookup_isluse(anode%elat, anode%elon, mype)
_EXIT_(myname_)
  return
end subroutine setluse_

!===================================================================
! Routines below are default code to be used, if they are not override
! by the code invoked this include-file.
subroutine obsheader_read_(iunit,mobs,jread,istat)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    obsheader_read_
!   prgmmr:      J. Guo
!
! abstract: read the jtype-block header record.
!
! program history log:
!   2015-01-12  guo     - constructed for generic obsnode
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
  integer(i_kind),intent(in ):: iunit
  integer(i_kind),intent(out):: mobs
  integer(i_kind),intent(out):: jread
  integer(i_kind),intent(out):: istat
  
  character(len=*),parameter:: myname_=myname//'::obsheader_read_'
_ENTRY_(myname_)
  read(iunit,iostat=istat) mobs,jread
_EXIT_(myname_)
  return
end subroutine obsheader_read_

subroutine obsheader_write_(junit,mobs,jwrite,jstat)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    obsheader_write_
!   prgmmr:      J. Guo
!
! abstract: write the jtype-block header record.
!
! program history log:
!   2015-01-12  guo     - constructed for generic obsnode
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
  integer(i_kind),intent(in ):: junit
  integer(i_kind),intent(in ):: mobs
  integer(i_kind),intent(in ):: jwrite
  integer(i_kind),intent(out):: jstat
  
  character(len=*),parameter:: myname_=myname//'::obsheader_write_'
_ENTRY_(myname_)
  write(junit,iostat=jstat) mobs,jwrite
_EXIT_(myname_)
  return
end subroutine obsheader_write_

subroutine init_(anode)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    init_
!   prgmmr:      J. Guo
!
! abstract: allocate a node.
!
! program history log:
!   2015-01-12  guo     - constructed for generic obsnode
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
  class(obsnode),intent(out):: anode

  character(len=*),parameter:: myname_=myname//'::init_'
_ENTRY_(myname_)
!_TRACEV_(myname_,'%mytype() =',anode%mytype())
  anode%llpoint => null()
  anode%luse = .false.
  anode%time = 0._r_kind
  anode%elat = 0._r_kind
  anode%elon = 0._r_kind
  anode%idv  =-1
  anode%iob  =-1
_EXIT_(myname_)
  return
end subroutine init_

subroutine clean_(anode)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    clean_
!   prgmmr:      J. Guo
!
! abstract: a shallow node clean
!
! program history log:
!   2015-01-12  guo     - constructed for generic obsnode
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
  class(obsnode),intent(inout):: anode

  character(len=*),parameter:: myname_=myname//'::clean_'
_ENTRY_(myname_)
!_TRACEV_(myname_,'%mytype() =',anode%mytype())
  call anode%init()
_EXIT_(myname_)
  return
end subroutine clean_

subroutine deepclean_(anode,deep,depth,stat)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 subroutine deepclean_
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2018-04-11
!
! abstract: a deep node clean
!
! program history log:
!   2018-04-11  j guo   - added this document block
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

  implicit none
  class(obsnode ),pointer ,intent(inout):: anode
  logical        ,optional,intent(in ):: deep   ! with deep=.true., the full
                                                ! linked-list headed by anode
                                                ! will be "deep" cleaned.
  integer(i_kind),optional,intent(out):: depth  ! depth of deep-cleaned nodes at
                                                ! the return.  zero is expected
                                                ! unless in an error.
  integer(i_kind),optional,intent(out):: stat   ! status return.

  character(len=*),parameter:: myname_=myname//'::deepclean_'
  integer(i_kind):: ier,depth_
  logical:: deep_

  if(present(depth)) depth=0
  if(present(stat )) stat=0

  if(.not.associated(anode)) return

  deep_=.false.
  if(present(deep )) deep_=deep

  if(deep_) then
     depth_=0
     call recurs_nodeclean_(anode,depth_,ier)
     if(present(depth)) depth=depth_

     if(ier/=0) then
        call perr(myname_,'recurs_nodeclean_(), stat =',ier)
        call perr(myname_,'                    depth =',depth_)
        call perr(myname_,'           anode%mytype() =',nodetype_(anode))
        if(.not.present(stat)) call die(myname_)
        stat=ier
        return
     endif

  else
     ! Full-clean anode itself, but not %llpoint.  This includes any dynamic
     ! component of anode defined in its type/endtype block.
     call anode%clean()
  endif

  return
end subroutine deepclean_

recursive subroutine recurs_nodeclean_(anode,depth,stat)
  implicit none
  class(obsnode),pointer,intent(inout):: anode
        ! This routine intends to fully erase the contents of argument anode,
        ! but not the storage of it.  A target attribute is used to prevent any
        ! attempt to deallocate.  Also see step (2) and (4) below.
  integer(i_kind),intent(inout):: depth
  integer(i_kind),intent(  out):: stat

  character(len=*),parameter:: myname_=myname//"::recurs_nodeclean_"

  stat=0
  if(associated(anode)) then

     if(associated(anode%llpoint)) then
        depth=depth+1
 
     ! (1) deep-clean the target of %llpoint, a level deeper than anode.

        call recurs_nodeclean_(anode%llpoint,depth,stat)
        if(stat/=0) return

     ! (2) deallocate %llpoint to release the memory associated with it.  This is
     !     in concert with step (4) below.

        deallocate(anode%llpoint,stat=stat)
        if(stat/=0) then
           call perr(myname_,"deallocate(anode%llpoint), stat =",stat)
           call perr(myname_,'                          depth =',depth)
           return
        endif

        depth=depth-1
     endif

     ! (3) full-clean anode itself other than %llpoint, including any its dynamic
     !     component defined in its type/endtype block.

     call anode%clean()

     ! (4) memory storage of anode itself is NOT expected to be deallocated.
     !     This is in concert with step (2) above.
  endif
  return
end subroutine recurs_nodeclean_

subroutine read_(anode,iunit,istat,redistr,diaglookup)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    read_
!   prgmmr:      J. Guo
!
! abstract: read the input for a node.
!
! program history log:
!   2015-01-12  guo     - constructed for generic obsnode
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
  use m_obsdiagnode, only: obsdiaglookup_locate
  use m_obsdiagnode, only: obs_diag
  use m_obsdiagnode, only: obs_diags
  implicit none
  class(obsnode),intent(inout):: anode
  integer(i_kind),intent(in   ):: iunit
  integer(i_kind),intent(  out):: istat
  logical        ,intent(in   ):: redistr
  type(obs_diags),intent(in   ):: diaglookup

  character(len=*),parameter:: myname_=myname//'::read_'
  integer(i_kind):: ier
_ENTRY_(myname_)

  istat=0
  read(iunit,iostat=ier) anode%luse,anode%time,anode%elat,anode%elon, &
                         !anode%dlat,anode%dlon, &
                         anode%idv ,anode%iob
  if(ier/=0) then
     call perr(myname_,'read(%(luse,time,elat,elon,...)), iostat =',ier)
     istat=-1
     _EXIT_(myname_)
     return
  endif

  istat=1                ! Now a complete xread(anode) is expected.
  if(redistr) then       ! Or additional conditions must be considered.
     istat=0             ! A complete xread(anode) is not expected, unless
     if(anode%luse) then ! ... .and. ...
        if(islocal_(anode)) istat=1
     endif
  endif

  call anode%xread(iunit,ier,diaglookup,skip=istat==0)
  if(ier/=0) then
     call perr(myname_,'anode%xread(), iostat =',ier)
     call perr(myname_,'                 skip =',istat==0)
     call perr(myname_,'                istat =',istat)
     istat=-2
     _EXIT_(myname_)
     return
  endif

_EXIT_(myname_)
  return
end subroutine read_

subroutine write_(anode,junit,jstat)
  implicit none
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    write_
!   prgmmr:      J. Guo
!
! abstract: write a node for output.
!
! program history log:
!   2015-01-12  guo     - constructed for generic obsnode
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
  class(obsnode),intent(in):: anode
  integer(i_kind),intent(in   ):: junit
  integer(i_kind),intent(  out):: jstat

  character(len=*),parameter:: myname_=myname//'::write_'
_ENTRY_(myname_)

  jstat=0
  write(junit,iostat=jstat) anode%luse,anode%time,anode%elat,anode%elon, &
                            anode%idv,anode%iob
  if(jstat/=0) then
     call perr(myname_,'write(%(luse,elat,elon,...)), jstat =',jstat)
     _EXIT_(myname_)
     return
  endif

  call anode%xwrite(junit,jstat)
  if (jstat/=0) then
     call perr(myname_,'anode%xwrite(), jstat =',jstat)
     _EXIT_(myname_)
     return
  end if
_EXIT_(myname_)
  return
end subroutine write_

subroutine show_(anode,iob)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    show_
!   prgmmr:      J. Guo
!
! abstract: show selected obsnode data.
!
! program history log:
!   2015-01-12  guo     - constructed for generic obsnode
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
  implicit none
  class(obsnode),intent(inout):: anode
  integer(i_kind),intent(in   ):: iob

  character(len=*),parameter:: myname_=myname//'::show_'
  logical:: isvalid_
_ENTRY_(myname_)
  isvalid_=anode%isvalid()
  write(stdout,"(2a,3i4,2x,2l1,3f8.2)") myname,":: iob,%(idv,iob,luse,vald,time,elat,elon) =", &
        iob,anode%idv,anode%iob,anode%luse,isvalid_,anode%time,anode%elat,anode%elon
_EXIT_(myname_)
  return
end subroutine show_

function nodetype_(anode)
!-- Return its type information, even when the argument is a null.
  implicit none
  character(len=:),allocatable:: nodetype_
  class(obsnode),pointer,intent(in):: anode
  nodetype_=".null.[obsnode]"
  if(associated(anode)) nodetype_=anode%mytype()
end function nodetype_

end module m_obsnode
