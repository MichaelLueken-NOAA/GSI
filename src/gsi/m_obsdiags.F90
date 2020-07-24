module m_obsdiags
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    module m_obsdiags
!   prgmmr:      j guo <jguo@nasa.gov>
!      org:      NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:      2015-02-04
!
! abstract: a bundle of GSI multiple obstypes and the obsdiags linked-lists
!
! program history log:
!   2015-02-04  j guo   - Re-implemented read_obsdiags() and write_obsdiags(),
!                         to support reconfigurable observation operators.  This
!                         implemenstation uses an obsllist template to support,
!                         in ceterian degree, static polymoprhism for different
!                         observation types.
!   2015-10-09  j guo   - By using Fortran 2003 dynamic polymorphism, this
!                         version has removed many ugly type dispatching select
!                         case constructs, by using an obsllist, a linked-list
!                         of dynamic polymorphic observation type (obsnode),
!                         which replaced the earlier obsllist template.
!   2016-06-22  j guo   - Added latlonrange for selected file readings, to let
!                         []_mread() to skip unnecessary read() of some files
!                         containing no relevant observations.
!                       . Added obsdiags_alwayslocal, as a user controlable
!                         switch to allow or to bypass selected file readings.
!                       . Added check_sizes_ outputs to allow size checkings.
!                       . Added #define show_llrange, for text-dumping of latlonranges.
!                       . Added #define debug_obsdiags, for text-dumping
!                         specific sections of obsdiags(:,:).
!                       . locally renamed mpi_comm_world to gsi_comm_world.
!   2018-01-23  k apodaca - Add a new observation type i.e. lightning (light) 
!                           suitable for the goes/glm instrument
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

  use kinds, only: i_kind, r_kind
  use mpeu_util, only: tell,warn,perr,die
  use mpeu_util, only: assert_
  use mpeu_util, only: stdout_open,stdout_close,stdout
  use mpeu_mpif, only: gsi_comm_world => mpi_comm_world

  use gsi_oboper, only: oboper

  use m_obsllist, only: obsllist
  use m_obsdiagnode, only: obs_diag
  use m_obsdiagnode, only: obs_diags

  use gsi_obopertypemanager, only: nobs_type => oboper_count
  use gsi_4dvar            , only: nobs_bins

  !use obsmod, only: obsdiags     ! (nobs_type,nobs_bins)
  implicit none
  private       ! except

  public:: obopers_config
  interface obopers_config; module procedure config_; end interface

        ! oboper instance creater with initialization to objects corresponding
        ! linked-list data instances.
  public:: oboper_create
  public:: oboper_headnode
  public:: oboper_destroy
  interface oboper_create; module procedure &
     createbydtype_, &
     createbyindex_, &
     createbyvmold_
  end interface
  interface oboper_headnode; module procedure headnode_; end interface
  interface oboper_destroy ; module procedure destroy_ ; end interface

  public:: obsdiags_reset
  public:: obsdiags_write
  public:: obsdiags_read
  public:: obsdiags_sort

  interface obsdiags_reset; module procedure reset_; end interface
  interface obsdiags_write; module procedure write_; end interface
  interface obsdiags_read ; module procedure mread_; end interface
  interface obsdiags_sort ; module procedure lsort_; end interface

  public:: obsdiags_create
  public:: obsdiags_destroy
  public:: obsdiags_inquire
  interface obsdiags_create ; module procedure  create_obsmod_vars; end interface
  interface obsdiags_destroy; module procedure destroy_obsmod_vars; end interface
  interface obsdiags_inquire; module procedure inquire_obsdiags   ; end interface

  public:: obsdiags_summary

  interface obsdiags_summary ; module procedure summary_ ; end interface

  public:: obsdiags_alwayslocal
  logical,save:: obsdiags_alwayslocal = .false.

! Note: User configurables
!
!   (1) obsdiags_mread(..,mpes,..) via /setup/:mpes_observer:
!
!       mpes==0, for reading "my own data";
!       mpes=>0, reading "all data", from pe 0 to mpes-1, but only up to the
!                highest count of the actually accessible files.
!
!       This value is configured through gsimod namelist/setup/:mpes_observer,
!       with the default value set to 0, to behave as it was ("my own data").
!       Otherwise, a simple usage is to let mpes_observer=1024, or other large
!       enough value, such that the solver-mode will try to determine how many
!       files created by the observer-mode are actually there to read.
!
!   (2) obsdiags_alwayslocal via /setup/:alwayslocal:
!       
!       obsdiags_alwayslocal sets an alternative default value of the optional
!       argument of obsdiags_mread(..,alwayslocal).
!
!       obsdiags_alwayslocal==.false., its default value.
!               It let obsdiags_mread() to check the locality of a file first,
!               using latlonrange_islocal(ipe), to avoid unnecessary opening+
!               reading some files.
!       obsdiags_alwayslocal==.true., override latlonrange_islocal(ipe).
!               It let obsdiags_mread() to always open+read all file.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname='m_obsdiags'
  logical,save:: lobsdiags_allocated_ = .false.
  logical,save:: lobstypes_allocated_ = .false.

  logical,parameter:: all_pes =.false.  ! report status on all pes or root only
  !logical,parameter:: all_pes =.true.  ! report status on all pes or root only
  logical,parameter:: do_summary =.false.  ! report status on all pes or root only
  !logical,parameter:: do_summary =.true.  ! report status on all pes or root only

        ! synch_messages is a flag to invoke mpi_barrier() calls before some
        ! status messages.  These calls are otherwise not necessary for the
        ! functionalities, but used here to ensure those messages mean what they
        ! intent to mean, in case that only the root pe is used to report some
        ! all pe status.

  !logical,parameter:: synch_messages = .true.          ! turn synch on
  !logical,parameter:: synch_messages = .not. all_pes   ! conditionally turn on
  logical,parameter:: synch_messages = .false.          ! turn synch off

  public :: obsdiags
  public :: obsllists

  type(obsllist ),save,dimension(:,:),pointer :: obsllists => null()
  type(obs_diags),save,dimension(:,:),pointer :: obsdiags  => null()  ! (nobs_type,nobs_bins)
  
  integer(i_kind),save:: jfunc__jiter = -1
  integer(i_kind),save:: jfunc__miter = -1
  integer(i_kind),save:: jfunc__jiterstart = -1

  integer(i_kind),save:: gsi_4dvar__nobs_bins  = -1
  integer(i_kind),save:: gsi_4dvar__min_offset = -1
  real   (r_kind),save:: gsi_4dvar__hr_obsbin  = -999._r_kind

!#define debug_trace
!#define debug_verbose
#include "mytrace.H"
#include "myassert.H"

#define _timer_on_
#ifdef  _timer_on_
#undef  _timer_on_
#undef  _timer_off_
#undef  _timer_use_
#define _timer_on_(id)  call timer_ini(id)
#define _timer_off_(id) call timer_fnl(id)
#define _timer_use_     use timermod, only: timer_ini,timer_fnl
#else
#define _timer_on_(id)
#define _timer_off_(id)
#define _timer_use_
#endif

  logical,parameter:: check_sizes_=.false.
  !logical,parameter:: check_sizes_=.true.

  !-- if(check_sizes_) then 
  !--   these size counters,

  integer(i_kind),allocatable,dimension(:),save:: lsize_type    !  luse counts of ob_type
  integer(i_kind),allocatable,dimension(:),save:: nsize_type    ! total counts of ob_type
  integer(i_kind),allocatable,dimension(:),save:: lsize_diag    !  luse counts of obs_diags
  integer(i_kind),allocatable,dimension(:),save:: msize_diag    !  muse counts of obs_diags
  integer(i_kind),allocatable,dimension(:),save:: nsize_diag    ! total counts of obs_diags

  !--   will be used to generate extra log-information, reporting different
  !--   size-counts of linked-lists, of all j-type, i-bin, on all pes.  Search
  !--   "check_sizes_" here for details.
  !-- endif

contains

subroutine config_()
!> Coupling external configurations (through external modules) to obopers' own
!> module configurations
  implicit none

!> For all obopers, import external configurations
  call     jfunc__import_()
  call gsi_4dvar__import_()

!> For specific obopers, import specific configurations
  call  lwcpoper__config_()

  return
contains
subroutine jfunc__import_()
!> jfunc parameters imported
  use jfunc, only: jiter
  use jfunc, only: miter
  use jfunc, only: jiterstart
  implicit none
  jfunc__jiter = jiter
  jfunc__miter = miter
  jfunc__jiterstart = jiterstart
  return
end subroutine jfunc__import_
subroutine gsi_4dvar__import_()
!> gsi4dvar parameters imported
  use gsi_4dvar, only: nobs_bins
  use gsi_4dvar, only: min_offset
  use gsi_4dvar, only: hr_obsbin
  implicit none
  gsi_4dvar__nobs_bins  = nobs_bins
  gsi_4dvar__min_offset = min_offset
  gsi_4dvar__hr_obsbin  = hr_obsbin
  return
end subroutine gsi_4dvar__import_
subroutine lwcpoper__config_()
!> gsi_lwcpoper parameters for configuration
!> gfs_stratosphere imports
  use gfs_stratosphere, only: use_gfs_stratosphere
  use gfs_stratosphere, only: nsig_save
!> lwcpoper
  use gsi_lwcpoper    , only: lwcpoper_config
  implicit none
  
  call lwcpoper_config()   ! reset to default
!> From gfs_stratosphere to gsi_lwcpoper, and expected to be refactored into an attribute of profile-vectors objects)
  if(use_gfs_stratosphere) call lwcpoper_config(nsig_save=nsig_save)
  return
end subroutine lwcpoper__config_
end subroutine config_

function createbydtype_(dtype) result(self)
!>> create an oboper to its components instanciated in this data module, with
!>> a given oboper registered dtype
  use gsi_obopertypemanager, only: oboper_typemold     ! (dtype)
  implicit none
  class(oboper),pointer:: self
  character(len=*),intent(in):: dtype
  character(len=*),parameter:: myname_=myname//"::createbydtype_"

  self => createbyvmold_(oboper_typemold(dtype))

#ifdef debug_verbose
! show status of the object for debugging
  call tell(myname_,'--- argument dtype =',trim(dtype))
  call tell(myname_,'associated(return) =',associated(self))
  !if(associated(self)) call oboper_show_(myname_,self)
#endif
end function createbydtype_

function createbyindex_(ioper) result(self)
!>> create an oboper to its components instanciated in this data module, with
!>> a given oboper registered index.
  use gsi_obopertypemanager, only: oboper_typemold     ! (ioper)
  use gsi_obopertypemanager, only: oboper_lbound
  use gsi_obopertypemanager, only: oboper_ubound
  implicit none
  class(oboper),pointer:: self
  integer(kind=i_kind),intent(in):: ioper

  character(len=*),parameter:: myname_=myname//"::createbyindex_"
  class(oboper),pointer:: mold_

  mold_ => oboper_typemold(ioper)
  if(associated(mold_)) then
     allocate(self,mold=mold_)

     if(ioper<lbound(obsllists,1) .or. ioper>ubound(obsllists,1)) then
        call perr(myname_,'unexpected value, ioper =',ioper)
        call perr(myname_,'    lbound(obsllists,1) =',lbound(obsllists,1))
        call perr(myname_,'    ubound(obsllists,1) =',ubound(obsllists,1))
        call perr(myname_,'              %mytype() =',self%mytype())
        call perr(myname_,'      %mytype(nodetype) =',self%mytype(nodetype=.true.))
        call  die(myname_)
     endif
     if(ioper<lbound( obsdiags,1) .or. ioper>ubound( obsdiags,1)) then
        call perr(myname_,'unexpected value, ioper =',ioper)
        call perr(myname_,'    lbound( obsdiags,1) =',lbound( obsdiags,1))
        call perr(myname_,'    ubound( obsdiags,1) =',ubound( obsdiags,1))
        call perr(myname_,'              %mytype() =',self%mytype())
        call perr(myname_,'      %mytype(nodetype) =',self%mytype(nodetype=.true.))
        call  die(myname_)
     endif

     call self%init(obsllists(ioper,:), &
                      obsdiags(ioper,:)  )
     mold_ => null()

  else
     call perr(myname_,'.not.associated, ioper =',ioper)
     call  die(myname_)
  endif

#ifdef debug_verbose
!>> show status of the object for debugging
  call tell(myname_,'--- argument ioper =',ioper)
  call tell(myname_,'associated(return) =',associated(self))
  !if(associated(self)) call oboper_show_(myname_,self)
#endif
end function createbyindex_

function createbyvmold_(mold) result(self)
!>> initialize an oboper to its components (linked-lists)
  use gsi_obopertypemanager, only: oboper_typeindex      ! to type-index
  use gsi_obopertypemanager, only: oboper_typeindex      ! to type-index
  implicit none
  class(oboper),pointer:: self
  class(oboper),target,intent(in):: mold

  character(len=*),parameter:: myname_=myname//"::createbyvmold_"
  integer(kind=i_kind):: itype  ! for a registered obsnode type index

  self => mold
  if(associated(self)) then
     allocate(self,mold=mold)

     ! Get its corresponding obsnode type name, then convert to its type-index
     itype=oboper_typeindex(self)

     if(itype<lbound(obsllists,1) .or. itype>ubound(obsllists,1)) then
        call perr(myname_,'unexpected value, itype =',itype)
        call perr(myname_,'    lbound(obsllists,1) =',lbound(obsllists,1))
        call perr(myname_,'    ubound(obsllists,1) =',ubound(obsllists,1))
        call perr(myname_,'              %mytype() =',self%mytype())
        call perr(myname_,'      %mytype(nodetype) =',self%mytype(nodetype=.true.))
        call  die(myname_)
     endif
     if(itype<lbound( obsdiags,1) .or. itype>ubound( obsdiags,1)) then
        call perr(myname_,'unexpected value, itype =',itype)
        call perr(myname_,'    lbound( obsdiags,1) =',lbound( obsdiags,1))
        call perr(myname_,'    ubound( obsdiags,1) =',ubound( obsdiags,1))
        call perr(myname_,'              %mytype() =',self%mytype())
        call perr(myname_,'      %mytype(nodetype) =',self%mytype(nodetype=.true.))
        call  die(myname_)
     endif

     call self%init(obsllists(itype,:), &
                     obsdiags(itype,:)  )
  endif

#ifdef debug_verbose
! show status of the object for debugging
  call tell(myname_,'--- argument mold%mytype() =',mold%mytype())
  call tell(myname_,'     mold%mytype(nodetype) =',mold%mytype(nodetype=.true.))
  call tell(myname_,'        associated(return) =',associated(self))
  if(associated(self)) call oboper_show_(myname_,self)
#endif
end function createbyvmold_

subroutine oboper_show_(mname,self)
  use gsi_oboper, only: oboper
  use gsi_obopertypemanager, only: oboper_typeindex
  use gsi_obopertypemanager, only: oboper_typeinfo
  use m_obsnodetypemanager , only: obsnode_typeindex     ! to type-index
  use mpeu_util, only: tell
  implicit none
  character(len=*),intent(in):: mname
  class(oboper),target,intent(in):: self

  call tell(mname,' oboper_typeindex(%) =',oboper_typeindex(self))
  call tell(mname,'  oboper_typeinfo(%) =',oboper_typeinfo(self))
  call tell(mname,'  associated(%obsll) =',associated(self%obsll))
  call tell(mname,'associated(%odiagll) =',associated(self%odiagll))
  call tell(mname,'     self%nodetype() =', self%mytype(nodetype=.true.))
end subroutine oboper_show_

subroutine destroy_(self)
  implicit none
  class(oboper),pointer,intent(inout):: self
  if(associated(self)) then
     call self%clean()
     deallocate(self)
  endif
end subroutine destroy_

function headnode_(ioboper,ibin) result(anode)
!>> Example: -- get the head node of an oboper%obsll(ibin)
!>>   psptr => psnode_typecast(headnode(ioboper_ps))
  use m_obsnode , only: obsnode
  use m_obsllist, only: obsllist_headnode
  use gsi_oboper, only: oboper
  implicit none
  integer(kind=i_kind),intent(in):: ioboper
  integer(kind=i_kind),intent(in):: ibin
  class(obsnode),pointer:: anode

  character(len=*),parameter:: myname_=myname//"::headnode_"
  class(oboper),pointer:: oboper_

  oboper_ => createbyindex_(ioboper)
  if(.not.associated(oboper_)) then
     call perr(myname_,'createbuindex_(), associated(oboper_) =',associated(oboper_))
     call perr(myname_,'                                ioper =',ioboper)
     call perr(myname_,'                                 ibin =',ibin)
     call  die(myname_)
  endif

  anode => obsllist_headnode(oboper_%obsll(ibin))
  call destroy_(oboper_)
end function headnode_

subroutine lobsdiags_statuscheck_(who,allocated)
!-- check the allocation status of basic obsdiags components.
  use obsmod, only: luse_obsdiag
  implicit none
  character(len=*),intent(in):: who
  logical,intent(in):: allocated

  if(.not.luse_obsdiag) return
  if(allocated) then
     if( .not.lobsdiags_allocated_ .or. &
         .not.lobstypes_allocated_ ) then
        if(.not.lobsdiags_allocated_) call perr(who,'.not.lobsdiags_allocated_')
        if(.not.lobstypes_allocated_) call perr(who,'.not.lobstypes_allocated_')
        call die(who)
     endif

  else
     if( lobsdiags_allocated_ .or. &
         lobstypes_allocated_ ) then
        if(lobsdiags_allocated_) call perr(who,'lobsdiags_allocated_ already')
        if(lobstypes_allocated_) call perr(who,'lobstypes_allocated_ already')
        call die(who)
     endif
  endif
end subroutine lobsdiags_statuscheck_

subroutine mread_(cdfile,mpes,force,jiter_expected,alwayslocal)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    m_obdiags::mread_
!   prgmmr:      tremolet
!
! abstract: Read obsdiags data structure from file.
!
! program history log:
!   2007-07-05  tremolet
!   2007-08-04  todling  - using get_lun to determine file unit number
!   2007-10-03  todling  - expanded to account for full observer 
!   2009-01-08  todling  - remove reference to ozohead
!   2009-01-23  todling  - add read_gpshead
!   2009-04-02  meunier  - add read_laghead
!   2010-04-27  tangborn - addded read_colvkhead
!   2010-05-26  treadon  - add read_tcphead
!   2011-05-18  todling  - aero, aerol, and pm2_5
!   2011-09-20  hclin    - 1d wij for aero
!   2015-02-04  j guo   - Re-implemented to support re-configurable observation
!                         operators.  read_() is split to read_() for a single
!                         file, and mread_() for one file only or all files for
!                         redistribution
!   2015-10-09  j guo   - Now it uses Fortran 2003 dynamic polymorphism.
!
!   input argument list:
!     cdfile - filename to read data from
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

  use mpeu_util, only: tell,perr,die,stdout_open,stdout_close,stdout
  _timer_use_
  use kinds, only: r_kind,i_kind

  use obsmod, only: lobserver
  use mpimod, only: mype
  use m_latlonrange, only: latlonrange
  use m_latlonrange, only: latlonrange_reset
  use m_latlonrange, only: latlonrange_islocal
  use m_latlonrange, only: latlonrange_readbcast
  use m_latlonrange, only: latlonrange_alldump

  use m_obsdiagnode, only: obsdiagllist_dump
  implicit none
  character(len=*), intent(in) :: cdfile          ! prefix, "obsdiags.<miter>"
  integer(i_kind),optional,intent(in):: mpes      ! number of files, from 0 to mpes-1
  logical        ,optional,intent(in):: force     ! force to read ob_types, regardless l4dvar etc.
  integer(i_kind),optional,intent(in):: jiter_expected  ! expected input jiter
  logical        ,optional,intent(in):: alwayslocal ! read all files

! ----------------------------------------------------------
  character(len=*),parameter:: myname_=myname//"::mread_"
  logical:: redistr,exist_
  integer(i_kind):: lpe,upe,ipe,ier
  integer(i_kind):: jtyp,jread
  logical:: force_read
  logical:: alwayslocal_
  logical:: fileislocal
  type(latlonrange),allocatable,dimension(:):: allranges
_ENTRY_(myname_)
_timer_on_(myname_)
!call stdout_open("obsdiags_mread")
  force_read=.false.
  if(present(force)) force_read=force
  alwayslocal_=obsdiags_alwayslocal
  if(present(alwayslocal)) alwayslocal_=alwayslocal

  call lobsdiags_statuscheck_(myname_,allocated=.true.)

        ! Determine the configuration, either read-my-own-data-only, or
        ! try-to-read-all-data-available.

  lpe=mype
  upe=lpe
  redistr=.false.
  if(present(mpes)) then
     if(mpes>0) then
        redistr=.true.
        lpe=0
        upe=-1
        do ipe=lpe,mpes-1
           inquire(file=trim(filename_(cdfile,ipe)), exist=exist_)
           if(exist_) upe=ipe
        enddo
     endif
  endif

        ! Reset components of obsdiags, for their re-construction from files
  call reset_()

  if(check_sizes_) then
     allocate(lsize_type(nobs_type))
     allocate(nsize_type(nobs_type))
     allocate(lsize_diag(nobs_type))
     allocate(nsize_diag(nobs_type))
     allocate(msize_diag(nobs_type))

     lsize_type(:)=0
     nsize_type(:)=0
     lsize_diag(:)=0
     nsize_diag(:)=0
     msize_diag(:)=0
  endif

        ! mpi_barrier() calls are not necessary.  They are used here to ensure
        ! the log-messages mean what they really mean, if only the root is used to
        ! report the all-pe status.

  if(synch_messages) call mpi_barrier(gsi_comm_world,ier)

  if(redistr) then
     if(mype==0) then
        call tell(myname_,'Reading obsdiags files for redistribution, npes =',upe-lpe+1)
        call tell(myname_,'                    prefix of the files, cdfile =',trim(cdfile))
        call tell(myname_,'                                            lpe =',lpe)
        call tell(myname_,'                                            upe =',upe)
        call tell(myname_,'                                    alwayslocal =',alwayslocal_)
     endif

     allocate(allranges(0:upe))
     call latlonrange_reset(allranges)
     call latlonrange_readbcast(hdfilename_(cdfile),allranges,root=0,comm=gsi_comm_world)

!#define show_llrange
#ifdef show_llrange
     call latlonrange_alldump(allranges,"obsllrange")
#endif


     jread=-1    ! checker of the input jiter values
     do ipe=lpe,upe
        fileislocal=latlonrange_islocal(allranges(ipe))
        if(alwayslocal_.or.fileislocal) then
           call read_(cdfile,ipe,redistr,fileislocal=fileislocal, &
                   force=force, &
                   jiter_expected=jiter_expected, &
                   verbose=.not.alwayslocal_.or.mype==0, &
                   jread=jread)
        endif
     enddo

!#define debug_obsdiags
#ifdef debug_obsdiags
        ! This is an example of dumping information for debugging, on selected
        ! pes, for specific jtyp and ibin.
        !
        ! This example is on pe #1, for (jtype==3 .and. ibin==3).

     if(mype==1) then
        call tell(myname_)
        call tell(myname_,'dumping obsdiags(), jtyp =',3)
        call tell(myname_,'                    ibin =',3)
        call tell(myname_,'                   jread =',jread)
        call obsdiagllist_dump(obsdiags(3,3),jiter=jread)
     endif
#endif

        ! Sort to ensure the ordering is unique.
     call lsort_()

     call latlonrange_reset(allranges)
     deallocate(allranges)

  else  ! of if(redistr)
     call read_(cdfile,mype,redistr,fileislocal=.true., &
         force=force, &
         jiter_expected=jiter_expected, &
         verbose=.true.)

  endif ! of if(redistr)

  if(mype==0) then
     call tell(myname_,'Finished reading of all obsdiags files, npes =',upe-lpe+1)
  endif
  
  if(check_sizes_) then
     do jtyp=lbound(lsize_type,1),ubound(lsize_type,1)
        if( msize_diag(jtyp)>0.or.lsize_diag(jtyp)>0.or.nsize_diag(jtyp)>0 .or. &
                                  lsize_type(jtyp)>0.or.nsize_type(jtyp)>0 ) then
           write(stdout,'(i5.3,i5,7x,5i8,2x,l1)')   mype,jtyp ,lsize_type(jtyp),nsize_type(jtyp), &
                                              msize_diag(jtyp),lsize_diag(jtyp),nsize_diag(jtyp)
        endif
     enddo

     call impi_reducesum_(lsize_type,root=0,comm=gsi_comm_world)
     call impi_reducesum_(nsize_type,root=0,comm=gsi_comm_world)
     call impi_reducesum_(lsize_diag,root=0,comm=gsi_comm_world)
     call impi_reducesum_(nsize_diag,root=0,comm=gsi_comm_world)
     call impi_reducesum_(msize_diag,root=0,comm=gsi_comm_world)

     if(mype==0) then
        do jtyp=lbound(lsize_type,1),ubound(lsize_type,1)
           if( msize_diag(jtyp)>0.or.lsize_diag(jtyp)>0.or.nsize_diag(jtyp)>0 .or. &
                                     lsize_type(jtyp)>0.or.nsize_type(jtyp)>0 ) then
              write(stdout,'(2x,a,i5,7x,5i8,2x,l1)')    '***',jtyp ,lsize_type(jtyp),nsize_type(jtyp), &
                                                   msize_diag(jtyp),lsize_diag(jtyp),nsize_diag(jtyp)
           endif
        enddo
     endif

     deallocate(lsize_type)
     deallocate(nsize_type)
     deallocate(lsize_diag)
     deallocate(nsize_diag)
     deallocate(msize_diag)
  endif

  if(synch_messages) call mpi_barrier(gsi_comm_world,ier)
  if(do_summary) call summary_(myname_)

  if(lobserver) then
     if(.not.force_read) then
        !call destroyobs(   skipit=.true.)
        call reset_(obsdiags_keep=.true.)
     endif
  endif
!call stdout_close()
_timer_off_(myname_)
_EXIT_(myname_)
  return
end subroutine mread_

subroutine reset_(obsdiags_keep)
  !use obsmod, only: obsdiags
  use obsmod, only: luse_obsdiag
  use obsmod, only: lobsdiag_allocated

  use m_obsdiagnode, only: obsdiagllist_reset
  use m_obsdiagnode, only: obsdiagllist_rewind
  use m_obsllist, only: obsllist_reset
  use m_obsnode , only: obsnode
  use gsi_obopertypemanager, only: oboper_typemold
  use gsi_obopertypemanager, only: oboper_lbound
  use gsi_obopertypemanager, only: oboper_ubound
  use m_obsnodetypemanager , only: obsnode_typemold

  _timer_use_
  implicit none
  logical,optional,intent(in):: obsdiags_keep
  character(len=*),parameter:: myname_=myname//'::reset_'
  integer(i_kind):: ii,jj
  logical:: obsdiags_keep_
  integer(i_kind):: ier
  class(obsnode),pointer:: mnode_
  class(oboper ),pointer:: moper_
_ENTRY_(myname_)
_timer_on_(myname_)

_TRACEV_(myname_,'lobsdiag_allocated   =',lobsdiag_allocated)
_TRACEV_(myname_,'lobsdiags_allocated_ =',lobsdiags_allocated_)

  ASSERT(nobs_type>0)
  ASSERT(nobs_bins>0)

  ! Both objects, obsdiags and obsllists are checked for their associated sizes
  ! and allocated shapes, regardless luse_obsdiag or not.  This is to simplify
  ! the algorithm logic.  The enforcements of (luse_obsdiag) are done on lower
  ! levels only.

  if(.not.lobstypes_allocated_) then
     lobstypes_allocated_=.true.
     if(.not.associated(obsllists)) call die(myname_,'unexpectedly, .not.associated(obsllists)')
  endif

  if(.not.lobsdiags_allocated_) then
     lobsdiags_allocated_=.true.
     if(.not.associated(obsdiags )) call die(myname_,'unexpectedly, .not.associated(obsdiags)')
  endif

  ASSERT(all(shape(obsdiags  )==shape(obsllists  )))
  ASSERT(     size(obsdiags,1)== size(obsllists,1) )
  ASSERT(     size(obsdiags,2)== size(obsllists,2) )

  obsdiags_keep_=.false.
  if(present(obsdiags_keep)) obsdiags_keep_=obsdiags_keep

  do ii=1,size(obsllists,2)      ! nobs_bins
     do jj=1,size(obsllists,1)    ! nobs_type
        if(luse_obsdiag) then
           if(.not.obsdiags_keep_) then
              call obsdiagllist_reset(obsdiags(jj,ii))
              lobsdiag_allocated=.false.

           else
              call obsdiagllist_rewind(obsdiags(jj,ii))

              ! In cases of rewinding without resetting, an obsdiagllist can
              ! be either initialized (lobsdiag_allocated), or not initialized
              ! (.not.lobsdiag_allocated).  So the code here should not try
              ! to alter the value of lobsdiag_allocated.
           endif
        endif

!++++
        moper_ => oboper_typemold(jj)
        if(.not.associated(moper_)) then
           call perr(myname_,'oboper_typemold(j) not associated, j =',jj)
           call perr(myname_,'                       oboper_lbound =',oboper_lbound)
           call perr(myname_,'                       oboper_ubound =',oboper_ubound)
           call  die(myname_)
        endif

        mnode_ => moper_%nodemold()
        if(.not.associated(mnode_)) then
           call perr(myname_,'moper_%nodemold() not associated, j =',jj)
           call perr(myname_,'                    moper_%mytype() =',moper_%mytype())
           call  die(myname_)
        endif
      
        moper_ => null()

!++++

        call obsllist_reset(obsllists(jj,ii),mold=mnode_, stat=ier)
        if(ier/=0) then
           call perr(myname_,'call obsllist_reset(), stat =',ier)
           call perr(myname_,'                       ibin =',ii)
           call perr(myname_,'                      jtype =',jj)
           call perr(myname_,'              mold%mytype() =',mnode_%mytype())
           call  die(myname_)
        endif

        mnode_ => null()
     enddo
  enddo
_timer_off_(myname_)
_EXIT_(myname_)
  return
end subroutine reset_

subroutine lsort_()
!$$$  subprogram documentation block
!
! abstract: sort entries of obsdiags(:,:) and obsllists(:,:)
!
! program history log:
!
!   input argument list:
!
!$$$

  use gsi_unformatted, only: unformatted_open
  use obsmod, only: luse_obsdiag

  use m_obsllist, only: obsllist_lsort
  use m_obsdiagnode, only: obsdiagllist_lsize
  use m_obsdiagnode, only: obsdiagllist_lsort

  _timer_use_
  implicit none

  character(len=*), parameter :: myname_=myname//"::lsort_"

  integer(i_kind) :: ii,jj !,iobs,lobs,ierr
_ENTRY_(myname_)
_timer_on_(myname_)
! ----------------------------------------------------------
  call lobsdiags_statuscheck_(myname_,allocated=.true.)

  if (luse_obsdiag) then

     ASSERT(all(shape(obsdiags)==shape(obsllists)))
     ASSERT(size(obsdiags,1)==size(obsllists,1))
     ASSERT(size(obsdiags,2)==size(obsllists,2))

  endif

  do jj=1,size(obsdiags,1)
     do ii=1,size(obsdiags,2)
        call obsdiagllist_lsort(obsdiags(jj,ii),itype=jj,ibin=ii)
     enddo
  enddo

  do jj=1,size(obsllists,1)
     do ii=1,size(obsllists,2)
        call obsllist_lsort(obsllists(jj,ii),itype=jj,ibin=ii)
     enddo
  enddo

! ----------------------------------------------------------
_timer_off_(myname_)
_EXIT_(myname_)
  return
end subroutine lsort_

subroutine write_(cdfile,luseonly,force)
!$$$  subprogram documentation block
!
! abstract: write obsdiags data structure to file.
!
! program history log:
!   2007-07-05  tremolet
!   2007-10-03  todling - expanded to account for full observer
!   2007-10-24  todling - add parameter nchnperobs to obsdiag
!   2009-01-08  todling - remove reference to ozohead
!   2009-01-27  todling - add gps write
!   2010-05-26  treadon - add write_tcphead
!   2010-06-03  todling - add write_colvkhead
!   2011-05-18  todling - aero, aerol, and pm2_5
!   2015-02-04  j guo   - Re-implemented to support re-configurable observation
!                         operators.
!   2015-10-09  j guo   - Now it uses Fortran 2003 dynamic polymorphism.
!
!   input argument list:
!     cdfile - filename to write data
!
!$$$

  use mpeu_util, only: tell,die,perr,stdout_open,stdout_close
_timer_use_

  use gsi_unformatted, only: unformatted_open
  use mpimod, only: mype
  use gsi_4dvar, only: l4dvar
  use  jfunc, only: jiter, miter

  use m_obsllist, only: obsllist_write
  use m_obsdiagnode, only: obsdiagllist_lsize
  use m_obsdiagnode, only: obsdiagllist_write

  use m_latlonrange, only: latlonrange
  use m_latlonrange, only: latlonrange_reset
  use m_latlonrange, only: latlonrange_gatherwrite
  use m_latlonrange, only: latlonrange_gatherdump

  implicit none
  character(len=*), intent(in) :: cdfile        ! := "obsdiags.<miter>"
  logical,optional, intent(in) :: luseonly      ! output only if(%luse)
  logical,optional, intent(in) :: force         ! write all out regardlessly

  character(len=*), parameter :: myname_=myname//"::write_"

  integer(i_kind) :: iunit,istat
  integer(i_kind) :: ii,jj,ier
  logical :: luseonly_
  logical :: force_write
  type(latlonrange):: luserange
! ----------------------------------------------------------
_ENTRY_(myname_)
_timer_on_(myname_)
!call stdout_open("obsdiags_write")
  force_write=.false.
  if(present(force)) force_write=force
  call lobsdiags_statuscheck_(myname_,allocated=.true.)

  ASSERT(all(shape(obsdiags)==shape(obsllists)))
  ASSERT(size(obsdiags,1)==size(obsllists,1))
  ASSERT(size(obsdiags,2)==size(obsllists,2))

  luseonly_=.false.
  if(present(luseonly)) luseonly_=luseonly
 
  call unformatted_open( unit=iunit, &
        file=trim(filename_(cdfile,mype)), &
        class='.obsdiags.', &
        action='write', &
        status='unknown', &
        newunit=.true., &       ! with newunit=.true., unit returns a value assigned by Fortran.
        iostat=istat,silent=.true.)
  if(istat/=0) then
     call perr(myname_,'unformatted_open(), file =',filename_(cdfile,mype))
     call perr(myname_,'                 newunit =',iunit)
     call perr(myname_,'                  iostat =',istat)
     call  die(myname_)
  endif

  if(do_summary) call summary_(myname_)

  do ii=1,size(obsdiags,2)
     do jj=1,size(obsdiags,1)
        call obsdiagllist_write(obsdiags(jj,ii),iunit,luseonly_,jiter,miter,jj,ii,luserange=luserange)

        if (force_write .or. l4dvar) then
           call obsllist_write(obsllists(jj,ii),iunit,luseonly_,jj,luserange=luserange)
        endif

        write(iunit)ii,jj   ! a jj_obstype-block trailer
     enddo
  enddo

  close(iunit)

        ! latlonrange_gatherwrite() implies a mpi_barrier() action.
  call latlonrange_gatherwrite(luserange,hdfilename_(cdfile),root=0,comm=gsi_comm_world)

#ifdef show_llrange
        ! Text-dump to diagnose the values
  call latlonrange_gatherdump(          "cvgllrange",root=0,comm=gsi_comm_world)
  call latlonrange_gatherdump(luserange,"obsllrange",root=0,comm=gsi_comm_world)
#endif

  call latlonrange_reset(luserange)

  if(synch_messages) call mpi_barrier(gsi_comm_world,ier)
  if (mype==0) call tell(myname_,'Finish writing obsdiags to file ',filename_(cdfile,mype))

! ----------------------------------------------------------
!call stdout_close()
_timer_off_(myname_)
_EXIT_(myname_)
  return
end subroutine write_

subroutine read_(cdfile,ipe,redistr,fileislocal,force,jiter_expected,verbose,jread)
  use mpeu_util, only: tell,perr,die
  use mpeu_util, only: stdout
  use mpimod, only: mype
  use gsi_4dvar, only: l4dvar
  use gsi_unformatted, only: unformatted_open
  use  jfunc, only: jiter,miter
  _timer_use_

  use obsmod, only: lobserver

  use m_obsllist, only: obsllist_read
  use m_obsllist, only: obsllist_lsize
  use m_obsllist, only: obsllist_lcount

  use m_obsdiagnode, only: obs_diag
  use m_obsdiagnode, only: obsdiagllist_read
  use m_obsdiagnode, only: obsdiagllist_lsize
  use m_obsdiagnode, only: obsdiagllist_lcount
  use m_obsdiagnode, only: obsdiaglookup_build
  use m_obsdiagnode, only: obsdiaglookup_clean

  implicit none
  character(len=*), intent(in ):: cdfile        ! prefix of the input file
  integer(i_kind ), intent(in ):: ipe           ! ipe of the input file
  logical         , intent(in ):: redistr       ! data redistribution is expected
  logical         , intent(in ):: fileislocal   ! the file to read, is known local

  logical,optional, intent(in ):: force         ! (force to read ob_type data
  integer(i_kind ), optional, intent(in   ):: jiter_expected    ! expecte input jiter
  logical,optional, intent(in ):: verbose       ! report each reading
  integer(i_kind ), optional, intent(inout):: jread     ! jiter read from the input

  character(len=*),parameter:: myname_=myname//'::read_'
  character(len=*),parameter:: diag_timer_=myname_//'.obsdiagllist_read'
  character(len=*),parameter:: list_timer_=myname_//'.obsllist_read'
  integer(i_kind):: ii,jj
  integer(i_kind):: ki,kj
  integer(i_kind):: iunit,istat
  integer(i_kind):: jread_
  integer(i_kind):: lsize_type_,nsize_type_
  integer(i_kind):: lsize_diag_,nsize_diag_,msize_diag_
  type(obs_diag),pointer:: leadnode => null()
  logical:: force_read
  logical:: verbose_
_ENTRY_(myname_)
_timer_on_(myname_)

  call lobsdiags_statuscheck_(myname_,allocated=.true.)
  force_read=.false.
  if(present(force)) force_read=force

  verbose_=.false.
  if(present(verbose)) verbose_=verbose

  ASSERT(all(shape(obsdiags)==shape(obsllists)))
  ASSERT(size(obsdiags,1)==size(obsllists,1))
  ASSERT(size(obsdiags,2)==size(obsllists,2))
  if(check_sizes_) then
     ASSERT(size(obsdiags,1)==size(lsize_type ))
     ASSERT(size(obsdiags,1)==size(nsize_type ))
     ASSERT(size(obsdiags,1)==size(lsize_diag ))
     ASSERT(size(obsdiags,1)==size(nsize_diag ))
  endif

  call unformatted_open( unit=iunit, &
        file=trim(filename_(cdfile,ipe)), &
        class='.obsdiags.', &
        action='read', &
        status='old', &
        newunit=.true., &       ! with newunit=.true., unit returns a value assigned by Fortran.
        iostat=istat,silent=.true.)
  if(istat/=0) then
     call perr(myname_,'unformatted_open(), file =',trim(filename_(cdfile,ipe)))
     call perr(myname_,'                    mype =',mype)
     call perr(myname_,'                     ipe =',ipe)
     call perr(myname_,'                   miter =',miter)
     call perr(myname_,'                 redistr =',redistr)
     call perr(myname_,'                 newunit =',iunit)
     call perr(myname_,'                  iostat =',istat)
     call  die(myname_)
  endif

  if(verbose_) call tell(myname_,'Reading obsdiags, file =',trim(filename_(cdfile,ipe)))

  leadnode => null()
  do ii=1,size(obsdiags,2)
     do jj=1,size(obsdiags,1)
        if(check_sizes_) then
           lsize_type_= obsllist_lcount(obsllists(jj,ii),luseonly=.true.,recount=.true.)
           nsize_type_= obsllist_lsize (obsllists(jj,ii)                )

           lsize_diag_= obsdiagllist_lcount(obsdiags(jj,ii),luseonly=.true.,recount=.true.)
           !msize_diag_= obsdiagllist_lcount(obsdiags(jj,ii),museonly=.true.)
           nsize_diag_= obsdiagllist_lsize (obsdiags(jj,ii)                )
        endif

        call obsdiagllist_read(obsdiags(jj,ii),iunit,redistr,jiter,miter,jj,ii,jread_,leadnode=leadnode, &
           jiter_expected=jiter_expected)
        if(present(jread)) then
           if(jread/=jread_) then
              if(jread>0) then
                 call perr(myname_,'not the same iteration, jiter =',jiter)
                 call perr(myname_,'                  saved jread =',jread)
                 call perr(myname_,'                current jread =',jread_)
                 call  die(myname_)
              endif
              jread=jread_
           endif
        endif

        call obsdiaglookup_build(obsdiags(jj,ii),leadnode=leadnode,jiter=jread)
        leadnode => null()     ! nullified after its use, to avoid leadnode dangling arround.

        if (force_read .or. l4dvar.and..not.(lobserver.and.jiter==1)) then
           call obsllist_read(obsllists(jj,ii),iunit,redistr,obsdiags(jj,ii),jj)
        endif

        call obsdiaglookup_clean(obsdiags(jj,ii))

        read(iunit)ki,kj
        if(ki/=ii .or. kj/=jj) then
           call perr(myname_,'mismatched block id, file =',filename_(cdfile,ipe))
           if(kj/=jj) then
              call perr(myname_,'               reading kj =',kj)
              call perr(myname_,'             expecting jj =',jj)
           endif
           if(ki/=ii) then
              call perr(myname_,'               reading ki =',ki)
              call perr(myname_,'             expecting ii =',ii)
           endif
           call perr(myname_,'                     file =',filename_(cdfile,ipe))
           call perr(myname_,'                   cdfile =',cdfile)
           call perr(myname_,'                     mype =',mype)
           call perr(myname_,'                      ipe =',ipe)
           call perr(myname_,'                    miter =',miter)
           call perr(myname_,'                  redistr =',redistr)
           call perr(myname_,'                  newunit =',iunit)
           call perr(myname_,'                   iostat =',istat)
           call  die(myname_)
        endif

        ASSERT(1<=jj.and.jj<=nobs_type)

        if(check_sizes_) then
           lsize_type_= obsllist_lcount(obsllists(jj,ii),luseonly=.true.)-lsize_type_
           nsize_type_= obsllist_lsize (obsllists(jj,ii)                )-nsize_type_

           lsize_diag_= obsdiagllist_lcount(obsdiags(jj,ii),luseonly=.true.)-lsize_diag_
           !msize_diag_= obsdiagllist_lcount(obsdiags(jj,ii),museonly=.true.)-msize_diag_
           nsize_diag_= obsdiagllist_lsize (obsdiags(jj,ii)                )-nsize_diag_

           if( fileislocal  .or. lsize_type_>0.or.nsize_type_>0 .or. &
               msize_diag_>0.or. lsize_diag_>0.or.nsize_diag_>0 ) then
              write(stdout,'(i5.3,2i5,2x,5i6,2x,l1)') ipe,jj,ii,lsize_type_,nsize_type_, &
                                                    msize_diag_,lsize_diag_,nsize_diag_,fileislocal
           endif
      
           lsize_type(jj)= lsize_type(jj) +lsize_type_
           nsize_type(jj)= nsize_type(jj) +nsize_type_

           lsize_diag(jj)= lsize_diag(jj) +lsize_diag_
           !msize_diag(jj)= msize_diag(jj) +msize_diag_
           nsize_diag(jj)= nsize_diag(jj) +nsize_diag_
        endif

     enddo       ! jj=1,size(obsdiags,1)
  enddo         ! ii=1,size(obsdiags,2)

  close(iunit)
! ----------------------------------------------------------
_timer_off_(myname_)
_EXIT_(myname_)
  return
end subroutine read_

function filename_(prefix,ipe)
!>> name of partitioned (obsdiags,obsllists) files
  implicit none
  character(len=:),allocatable:: filename_
  character(len=*)    , intent(in ):: prefix
  integer(kind=i_kind), intent(in ):: ipe

  character(len=4):: chpe
  write(chpe,'(i4.4)') ipe
  filename_=trim(adjustl(prefix))//'.'//trim(chpe)
end function filename_

function hdfilename_(prefix)
!>> name of the header file
  use kinds, only: i_kind
  implicit none
  character(len=:),allocatable:: hdfilename_
  character(len=*)    , intent(in ):: prefix
  hdfilename_=trim(adjustl(prefix))//'.headers'
end function hdfilename_

subroutine summary_(title)
!-- get a summary of obsdiags(:,:) and obsllists(:,:)
  use obsmod, only: luse_obsdiag
  use mpeu_util, only: tell,die,perr,stdout_open,stdout_close
_timer_use_

  use gsi_unformatted, only: unformatted_open
  use gsi_4dvar, only: nobs_bins

  use m_obsllist, only: obsllist_lsize => obsllist_lcount
  use m_obsdiagnode, only: obsdiagllist_lsize => obsdiagllist_lcount

  implicit none
  character(len=*), intent(in) :: title

  character(len=*), parameter :: myname_=myname//"::summary_"

  integer(i_kind) :: ii,jj
  integer(i_kind),dimension(nobs_type,nobs_bins):: ldiag,ndiag
  integer(i_kind),dimension(nobs_type,nobs_bins):: lobss,nobss
_ENTRY_(myname_)
_timer_on_(myname_)
! ----------------------------------------------------------

  call lobsdiags_statuscheck_(myname_,allocated=.true.)

  if (luse_obsdiag) then
     ASSERT(all(shape(obsdiags)==shape(obsllists)))
     ASSERT(size(obsdiags,1)==size(obsllists,1))
     ASSERT(size(obsdiags,2)==size(obsllists,2))
  endif

  do ii=1,size(obsdiags,2)
     do jj=1,size(obsdiags,1)
        ldiag(jj,ii) = obsdiagllist_lsize(obsdiags(jj,ii),luseonly=.true. ,recount=.true.)
        ndiag(jj,ii) = obsdiagllist_lsize(obsdiags(jj,ii),luseonly=.false.,recount=.true.)
     enddo
  enddo

  do ii=1,size(obsllists,2)
     do jj=1,size(obsllists,1)
        lobss(jj,ii) = obsllist_lsize(obsllists(jj,ii),luseonly=.true. ,recount=.true.)
        nobss(jj,ii) = obsllist_lsize(obsllists(jj,ii),luseonly=.false.,recount=.true.)
     enddo
  enddo

  call gather_write_(title,lobss,ldiag,nobss,ndiag,root=0,comm=gsi_comm_world)

! ----------------------------------------------------------
_timer_off_(myname_)
_EXIT_(myname_)
  return
end subroutine summary_

subroutine gather_write_(title,lobss,ldiag,nobss,ndiag,root,comm)
  use mpimod   , only: mype,npe
  use kinds    , only: i_kind
  use mpeu_mpif, only: mpi_ikind
  _timer_use_
  implicit none
  character(len=*),intent(in):: title
  integer(kind=i_kind),dimension(:,:),intent(in):: lobss,ldiag
  integer(kind=i_kind),dimension(:,:),intent(in):: nobss,ndiag
  integer(kind=mpi_ikind),intent(in):: root
  integer(kind=mpi_ikind),intent(in):: comm

  character(len=*),parameter:: myname_=myname//'::gather_write_'
  integer(kind=i_kind):: jj,ii,ipe
  integer(kind=i_kind) :: mtyp,mbin,mpes
  integer(kind=i_kind),allocatable,dimension(:,:,:):: ldiagm,ndiagm
  integer(kind=i_kind),allocatable,dimension(:,:,:):: lobssm,nobssm

_ENTRY_(myname_)
_timer_on_(myname_)
  mtyp=size(lobss,1)
  mbin=size(lobss,2)
  ASSERT(mtyp==size(nobss,1))
  ASSERT(mbin==size(nobss,2))
  ASSERT(mtyp==size(ldiag,1))
  ASSERT(mbin==size(ldiag,2))
  ASSERT(mtyp==size(ndiag,1))
  ASSERT(mbin==size(ndiag,2))

  mpes=0        ! its value is significant only on root
  if(mype==root) mpes=npe

  allocate(lobssm(mtyp,mbin,0:mpes-1))
  allocate(ldiagm(mtyp,mbin,0:mpes-1))
  allocate(nobssm(mtyp,mbin,0:mpes-1))
  allocate(ndiagm(mtyp,mbin,0:mpes-1))

  call impi_gather_(lobss,lobssm,root,comm)
  call impi_gather_(nobss,nobssm,root,comm)
  call impi_gather_(ldiag,ldiagm,root,comm)
  call impi_gather_(ndiag,ndiagm,root,comm)

  if(mype==root) then
     do ipe=0,npe-1
        write(stdout,'(2a,i6)'     ) title,'(): local obs/diag counts, ipe =',ipe
        write(stdout,'(2a,9(1x,a))') title,'(): typ', ('|  -----lo -----ld  -----no -----nd',ii=1,mbin)
        do jj=1,mtyp
           write(stdout,'(2a,i3,9(1x,a,2(1x,2i8)))') &
                                      title,'(): ',jj , &
             ("|",lobssm(jj,ii,ipe),ldiagm(jj,ii,ipe), &
                  nobssm(jj,ii,ipe),ndiagm(jj,ii,ipe), ii=1,mbin)
        enddo
     enddo
  endif

  deallocate(lobssm)
  deallocate(ldiagm)
  deallocate(nobssm)
  deallocate(ndiagm)
_timer_off_(myname_)
_EXIT_(myname_)
end subroutine gather_write_

subroutine impi_barrier_(comm)
  use mpeu_mpif, only: mpi_ikind
  use mpeu_util, only: die
  implicit none
  integer(kind=mpi_ikind),intent(in):: comm

  character(len=*),parameter:: myname_=myname//"::impi_barrier_"
  integer(kind=mpi_ikind):: ier

  call mpi_barrier(comm,ier)
  if(ier/=0) call die(myname_,'mpi_barrier() error, ierror =',ier)
end subroutine impi_barrier_

subroutine impi_gather_(isend,irecv,root,comm)
  use mpeu_mpif,only: mpi_ikind,mpi_type
  use mpeu_util, only: die
  use kinds, only: i_kind
  implicit none
  integer(kind=i_kind),dimension(:,:  ),intent(in ):: isend
  integer(kind=i_kind),dimension(:,:,:),intent(out):: irecv
  integer(kind=mpi_ikind),intent(in):: root
  integer(kind=mpi_ikind),intent(in):: comm

  character(len=*),parameter:: myname_=myname//"::impi_gather_"
  integer(kind=mpi_ikind):: itype,isize,ierr

  isize=size(isend)
  itype=mpi_type(isend)
  call mpi_gather(isend,isize,itype, &
                  irecv,isize,itype, root,comm,ierr)
  if(ierr/=0) call die(myname_,'mpi_gather() error, ierror =',ierr)
end subroutine impi_gather_

subroutine impi_reducesum_(iredu,root,comm)
  use mpeu_mpif,only: mpi_ikind,mpi_type,mpi_sum
  use mpeu_util, only: die
  use kinds, only: i_kind
  implicit none
  integer(kind=i_kind),dimension(:),intent(inout):: iredu
  integer(kind=mpi_ikind),intent(in):: root
  integer(kind=mpi_ikind),intent(in):: comm

  character(len=*),parameter:: myname_=myname//"::impi_reducesum_"
  integer(kind=mpi_ikind):: itype,isize,ierr
  !integer(kind=kind(iredu)),dimension(size(iredu)):: irecv

  isize=size(iredu)
  itype=mpi_type(iredu)
  call mpi_reduce((iredu),iredu,isize,itype, mpi_sum, root,comm,ierr)
  if(ierr/=0) call die(myname_,'mpi_reduce(mpi_sum) error, ierror =',ierr)
  !iredu(:)=irecv(:)
end subroutine impi_reducesum_

subroutine create_obsmod_vars()
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    create_obsmod_vars
!     prgmmr:    derber            org: np23           date: 2003-09-25
!
! abstract:  allocate arrays to hold observation related information
!
! program history log:
!   2003-09-25  derber
!   2004-05-13  treadon, documentation
!   2015-10-09  j guo   - moved here from module obsmod with modifcations
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!$$$ end documentation block
  use gsi_4dvar, only: nobs_bins
  implicit none
  lobstypes_allocated_=.true.
  lobsdiags_allocated_=.true.
  allocate(obsllists(nobs_type,nobs_bins))
  allocate(obsdiags (nobs_type,nobs_bins))
  return
end subroutine create_obsmod_vars

subroutine destroy_obsmod_vars()
!-- Created to pair with create_obsmod_vars().
  implicit none
  deallocate(obsllists)
  deallocate(obsdiags )
  lobstypes_allocated_=.false.
  lobsdiags_allocated_=.false.
  return
end subroutine destroy_obsmod_vars

! ----------------------------------------------------------------------
subroutine inquire_obsdiags(kiter)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    inquire_obsdiags
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-08-07  lueken - added  subprogram doc block
!
!   input argument list:
!    kiter
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

  use constants, only:  one,two,three,four,five
  use mpimod, only: mpi_max,mpi_comm_world,ierror,mype
  use mpeu_mpif, only: mpi_type, mpi_ikind
  implicit none

  integer(i_kind), intent(in   ) :: kiter

  real(r_kind) :: sizei, sizer, sizel, sizep, ziter, zsize, ztot
  integer(i_kind) :: ii,jj,iobsa(2),iobsb(2)
  type(obs_diag), pointer :: obsptr => null()

! Any better way to determine size or i_kind, r_kind, etc... ?
  sizei=four
  sizer=8.0_r_kind
  sizel=one
  sizep=four

  iobsa(:)=0
  do ii=1,size(obsdiags,2)
     do jj=1,size(obsdiags,1)
        obsptr => obsdiags(jj,ii)%head
        do while (associated(obsptr))
           iobsa(1)=iobsa(1)+1
           if (any(obsptr%muse(:))) iobsa(2)=iobsa(2)+1
           obsptr => obsptr%next
        enddo
     enddo
  enddo

  call mpi_reduce(iobsa,iobsb,2_mpi_ikind,mpi_type(iobsa),mpi_max,0_mpi_ikind,mpi_comm_world,ierror)

  if (mype==0) then
     ziter=real(kiter,r_kind)
     zsize = sizer*(three*ziter+two) + sizei + sizel*(ziter+one) + sizep*five
     ztot=real(iobsb(1),r_kind)*zsize
     ztot=ztot/(1024.0_r_kind*1024.0_r_kind)
 
     write(6,*)'obsdiags: Bytes per element=',nint(zsize)
     write(6,*)'obsdiags: length total, used=',iobsb(1),iobsb(2)
     write(6,'(A,F8.1,A)')'obsdiags: Estimated memory usage= ',ztot,' Mb'
  endif

end subroutine inquire_obsdiags

end module m_obsdiags
