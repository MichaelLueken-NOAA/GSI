module gsi_unformatted
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    module gsi_unformatted
!   prgmmr:      j guo <jguo@nasa.gov>
!      org:      NASA/GSFC, Global Modeling and Assimilation Office, 900.3
!     date:      2013-05-06
!
! abstract: open() for named sequential unformatted files with "convert" specifier.
!
! program history log:
!   2013-05-06  j guo   - initial implementation.
!                       - added this document block
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

  use kinds    , only: i_kind
  use mpeu_util, only: mpeu_luavail => luavail
  use mpeu_util, only: mpeu_getarec => getarec
  use mpeu_util, only: die,perr,warn,tell
  use mpeu_util, only: indexSet, indexSort
  implicit none
  private       ! except

  public :: unformatted_open    ! (unit,file[,class][,action][,position][,status][,iostat])
                                !  -- an interface to open a Fortran sequential unformatted file
  public :: fileinfo_lookup     ! (class,convert)
                                !  -- look up class in a fileinfo table, for a convert definition
  public :: fileinfo_reset      ! ([fileinfo])
                                !  -- deallocate(fileinfo_xx) or set an alternate fileinfo filename.
  public :: fileinfo_len        ! the internal lenth of fields class and convert.

  interface unformatted_open; module procedure   open_; end interface
  interface fileinfo_lookup ; module procedure lookup_; end interface
  interface fileinfo_reset  ; module procedure  reset_; end interface

!!! For implementor: For compilers do not support "convert", this is
!!! the place to #ifdef them out.
#ifdef __X_OR_Y_OR_Z_FORTRAN_COMPILERS__
#define _DO_NOT_SUPPORT_OPEN_WITH_CONVERT_
#endif

!!! Usage:
!!!
!!!   ! lookup convert value for a user defined class '.bufr.'
!!!     use gsi_unformatted, only: fileinfo_lookup, fileinfo_len
!!!     character(len=fileinfo_len):: convert
!!!     convert='native'
!!!     call fileinfo_lookup('.bufr.',convert)
!!! or
!!!   ! lookup convert value for a given filename
!!!     use gsi_unformatted, only: fileinfo_lookup, fileinfo_len
!!!     character(len=fileinfo_len):: convert
!!!     convert=''
!!!     call fileinfo_lookup(filename,convert)
!!! or
!!!   ! open an existed bufr file for input
!!!     use gsi_unformatted, only: unformatted_open
!!!     call unformatted_open(unit,file,class='.bufr.',status='old',iostat=ier)
!!!
!!# This is an example of file 'unformatted_fileinfo', a user defined database
!!# file supporting module gsi_unformatted, for user specific local GSI system
!!# configuration.
!!#
!!# Note: convert="" and convert="native" may be different, when this code is
!!#     compiled with a compiler "--convert <keyword>" flag.  "native" refers
!!#     to the platform default, while "" refers to <keyword> of flag "--convert".
!!#
!!# Note: Just like a filename, class is case sensitive, at least for this
!!#     implementation.
!!#
!!# Note: Reserved values in this implementation of this module.
!!#     class   == ".default."          -- for all files
!!#     convert == "_not_supported_"    -- flag compilers not supporting "convert"
!!#     convert == "_not_found_"        -- flag a failed lookup() call.
!!#
!!#     class/file      convert
!!#--------------------------------------
!!      .default.       ""              # a class name reserved for all files
!!      .bufr.          little_endian   # for all bufr files
!!      prepbufr        native          # an exception from .bufr.
!!      .berror.        native          # for files grouped under .berror.
!!      berror_stats    big_endian      # an exception from .berror.
!!      .diag.          big_endian      # for files grouped under .diag.


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname='gsi_unformatted'

  integer(i_kind),parameter:: fileinfo_len=64
  integer(i_kind),parameter:: fileinfo_inc=32

  integer(i_kind),parameter:: fileinfo_rec=256
  integer(i_kind),parameter:: fileinfo_fnl=512  ! in case it has a very long pathname

#ifdef _DO_NOT_SUPPORT_OPEN_WITH_CONVERT_
  logical,parameter:: convert_supported_ = .false.
#else
  logical,parameter:: convert_supported_ = .true.
#endif

! Declare a "_fileinfo_" data structure, defining a class-vs-convert table.

  character(len=*),parameter:: default_fileinfo_name ='unformatted_fileinfo'

  logical,save:: fileinfo_initialized_ = .false.

  character(len=fileinfo_fnl),save:: fileinfo_name_=default_fileinfo_name
  integer(i_kind),save:: fileinfo_msize_=fileinfo_inc   ! allocated size
  integer(i_kind),save:: fileinfo_lsize_=-1             ! actual size
  character(len=fileinfo_len),dimension(:),pointer,save:: fileinfo_class_
  character(len=fileinfo_len),dimension(:),pointer,save:: fileinfo_cnvrt_
  integer(i_kind)            ,dimension(:),pointer,save:: fileinfo_index_

! name/class    convert

contains
subroutine open_(unit,file,class,newunit,action,position,status,iostat,silent)
  implicit none
  integer(i_kind) ,intent(inout):: unit                 ! logical unit
  character(len=*),intent(in   ):: file                 ! filename
  character(len=*),optional,intent(in ):: class         ! file class of convert specifier
  logical         ,optional,intent(in ):: newunit       ! file class of convert specifier
  character(len=*),optional,intent(in ):: action        ! 'read', 'write', or 'readwrite'
  character(len=*),optional,intent(in ):: position      ! 'rewind', 'append', or 'asis'
  character(len=*),optional,intent(in ):: status        ! 'old', 'new', 'unknown', or 'scratch'
  integer(i_kind) ,optional,intent(out):: iostat        ! the return status
  logical         ,optional,intent(in ):: silent        ! not complain on missing fileinfo

  integer(i_kind):: iostat_
  logical        :: newunit_
  character(len=fileinfo_len):: class_
  character(len=fileinfo_len):: action_
  character(len=fileinfo_len):: position_
  character(len=fileinfo_len):: status_
  character(len=fileinfo_len):: convert_
  character(len=*),parameter:: myname_=myname//'::open_'

  class_   ='.default.'; if(present(class   )) class_   =class
  action_  ='readwrite'; if(present(action  )) action_  =action
  position_='rewind'   ; if(present(position)) position_=position
  newunit_ =.false.    ; if(present(newunit )) newunit_ =newunit
  status_  ='unknown'  ; if(present(status  )) status_  =status

  if(present(iostat)) iostat=0

#ifdef _DO_NOT_SUPPORT_OPEN_WITH_CONVERT_
  convert_="_not_supported_"    ! open(file with the compiler default convert
  if(newunit_) then
     open(newunit=unit,file=file,access='sequential',form='unformatted', &
        action=action_,position=position_,status=status_,iostat=iostat_)
  else
     open(   unit=unit,file=file,access='sequential',form='unformatted', &
        action=action_,position=position_,status=status_,iostat=iostat_)
  endif

#else
  convert_="_not_found_"        ! set a difault value
  call lookup_(class_,convert_,silent=silent) ! may override convert value, if an entry of class_ is found.
  call lookup_(file  ,convert_,silent=silent) ! may override convert value, if an entry of file is found.

  select case(convert_)
     case("","_not_found_")        ! open(file) with the compiler default convert
        if(newunit_) then
           open(newunit=unit,file=file,access='sequential',form='unformatted', &
              action=action_,position=position_,status=status_,iostat=iostat_)
        else
           open(   unit=unit,file=file,access='sequential',form='unformatted', &
              action=action_,position=position_,status=status_,iostat=iostat_)
        endif

     case default                  ! open(file) with user specified convert
        if(newunit_) then
           open(newunit=unit,file=file,access='sequential',form='unformatted', &
              action=action_,position=position_,status=status_,iostat=iostat_, &
              convert=convert_)
        else
           open(   unit=unit,file=file,access='sequential',form='unformatted', &
              action=action_,position=position_,status=status_,iostat=iostat_, &
              convert=convert_)
        endif
  end select
#endif

  if(iostat_/=0) then
     call perr(myname_,'open() error, iostat =',iostat_)
     call perr(myname_,'                unit =',unit)
     call perr(myname_,'                file =',trim(file))
     call perr(myname_,'              status =',trim(status_))
     call perr(myname_,'              action =',trim(action_))
     call perr(myname_,'            position =',trim(position_))
     call perr(myname_,'       FILEINFO_NAME =',trim(fileinfo_name_))
     call perr(myname_,'      fileinfo.class =',trim(class_))
     call perr(myname_,'    fileinfo.convert =',trim(convert_))

     if(.not.present(iostat)) call die(myname_)
     iostat=iostat_
     return
  endif
end subroutine open_

subroutine lookup_(class,convert,silent)
  implicit none
  character(len=*),intent(in   ):: class        ! class or filename itself
  character(len=*),intent(inout):: convert      ! may be override if an entry of class is found
  logical,optional,intent(in   ):: silent       ! not complain on missing fileinfo file

  character(len=*),parameter:: myname_=myname//'::lookup_'
  integer(i_kind):: l
  logical:: silent_
  silent_=.false.
  if(present(silent)) silent_=silent
  if(.not.fileinfo_initialized_) call init_(.not.silent_)

  if(fileinfo_lsize_<=0) return

  l=fileinfo_lsize_
  call lookitup_(l, fileinfo_index_(1:l), &
                    fileinfo_class_(1:l), &
                    fileinfo_cnvrt_(1:l), &
                    class, convert        )
end subroutine lookup_
!!!!!!!
subroutine reset_(fileinfo)
  implicit none
  character(len=*),optional,intent(in):: fileinfo       ! an alternate fileinfo name

        ! Reset fileinfo_name_, even if the fileinfo part has not been not
        ! initialized_.  So one can lookup() from a different fileinfo.
  fileinfo_name_= default_fileinfo_name
  if(present(fileinfo)) fileinfo_name_= fileinfo

        ! Initialization (init_()) is defered to the time an actual lookup().
  if(.not.fileinfo_initialized_) return

        ! Reset to the pre-init_() state, except fileinfo_name_
  fileinfo_initialized_ = .false.
  fileinfo_msize_       = fileinfo_inc
  fileinfo_lsize_       = -1
  deallocate( fileinfo_class_, &
              fileinfo_cnvrt_, &
              fileinfo_index_  )
end subroutine reset_
subroutine init_(verbose)
  implicit none
  logical,intent(in):: verbose

! local variables
  integer(i_kind):: lu,ier,i,n

  character(len=fileinfo_len):: classi
  character(len=fileinfo_len):: cnvrti
  character(len=fileinfo_rec):: arec

  character(len=fileinfo_len),pointer,dimension(:):: p_class
  character(len=fileinfo_len),pointer,dimension(:):: p_cnvrt
  character(len=*),parameter:: myname_=myname//'::init_'

  fileinfo_initialized_=.true.

  if(.not.convert_supported_.and.verbose) call warn(myname_,'Not supported, open(convert=..)')

        ! read in the fileinfo table anyway
  lu=mpeu_luavail()
  open(lu,file=fileinfo_name_,status='old',form='formatted',iostat=ier)
  if(ier/=0) then
#ifndef _DO_NOT_SUPPORT_OPEN_WITH_CONVERT_
     if(verbose) then
        call warn(myname_,'Can not open, file =',trim(fileinfo_name_))
        call warn(myname_,'Will use default convert values in code')
     endif
#endif

     fileinfo_lsize_=0
     allocate( fileinfo_class_(0), &
               fileinfo_cnvrt_(0), &
               fileinfo_index_(0)  )
     return
  endif

  n=fileinfo_msize_
  allocate( fileinfo_class_(n), &
            fileinfo_cnvrt_(n))

  i=0
  call mpeu_getarec(lu,arec,ier,commchar='#!')
  do while(ier==0)
     read(arec,*,iostat=ier) classi,cnvrti
     if(ier/=0) cnvrti=""

     i=i+1
     if(i>fileinfo_msize_) then  ! realloc()
        p_class => fileinfo_class_
        p_cnvrt => fileinfo_cnvrt_
        n=fileinfo_msize_+fileinfo_inc
        allocate( fileinfo_class_(n), &
                  fileinfo_cnvrt_(n))

        fileinfo_class_(1:fileinfo_msize_)=p_class(:)
        fileinfo_cnvrt_(1:fileinfo_msize_)=p_cnvrt(:)

        fileinfo_msize_=n
        deallocate(p_class,p_cnvrt)
     endif

     fileinfo_class_(i)=classi
     fileinfo_cnvrt_(i)=cnvrti

     call mpeu_getarec(lu,arec,ier,commchar='#!')
  enddo

  fileinfo_lsize_=i
  close(lu)

  allocate(fileinfo_index_(fileinfo_lsize_))
  call indexset (fileinfo_index_(1:fileinfo_lsize_))
  call indexsort(fileinfo_index_(1:fileinfo_lsize_), &
                 fileinfo_class_(1:fileinfo_lsize_), descend=.false.)

end subroutine init_

subroutine lookitup_(lsize_,index_,class_,cnvrt_,classi,convert)
  implicit none
  integer(i_kind) ,intent(in):: lsize_
  integer(i_kind) ,dimension(:),intent(in):: index_
  character(len=*),dimension(:),intent(in):: class_
  character(len=*),dimension(:),intent(in):: cnvrt_

  character(len=*),intent(in   ):: classi
  character(len=*),intent(inout):: convert

  logical:: done
  integer(i_kind):: lb,ub,i,l
  character(len=*),parameter:: myname_=myname//'::init_'

  done=.false.
  lb=1; ub=lsize_
  do while( .not. (done .or. ub<lb) )
     i=(lb+ub)/2
     l=index_(i)
     if    (classi<class_(l)) then
        ub=i-1
     elseif(classi>class_(l)) then
        lb=i+1
     else
        done=.true.
        convert=cnvrt_(l)
        exit
     endif
  enddo
  return
end subroutine lookitup_

end module gsi_unformatted
