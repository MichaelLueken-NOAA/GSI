module m_latlonrange
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    module m_latlonrange
!   prgmmr:      j guo <jguo@nasa.gov>
!      org:      NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:      2016-02-09
!
! abstract: lat-lon range information for subdomains and observation clusters
!
! program history log:
!   2016-02-09  j guo   - initial code
!   2016-02-09  j guo   - added this document block
!   2016-06-22  j guo   - refined the definition for local subdomain ranges
!                       . added text-dump routines for diagnosis
!                       . added "pure" qualifier to some routines.
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
  use mpeu_mpif, only: mpi_ikind
  use mpeu_util, only: assert_
  use mpeu_util, only: tell,perr,die
  use timermod, only: timer_ini,timer_fnl
  implicit none
  private       ! except
  public :: latlonrange          ! data structure

  public :: latlonrange_reset    ! the initializer or the finalizer
  public :: latlonrange_enclose  ! updating the range by including a (elat,elon)
  public :: latlonrange_islocal  ! falls (partially) in the "local" subdomain.

  public :: latlonrange_set      ! set the range for a known lat-lon domain
  public :: latlonrange_overlaps ! []_overlaps(ar, br) is "(br) overlaps (ar)"

  public :: latlonrange_gatherwrite    ! gather then root-write, for all-ranges
  public :: latlonrange_readbcast      ! root-read then bcast, for all-ranges

        ! Interfaces for debugging purposes

  public :: latlonrange_alldump    ! all-text-dump a local latlonrange array.
     interface latlonrange_alldump; module procedure alldump_; end interface

  public :: latlonrange_gatherdump ! gather local scalar latlonranges, then root-text-dump all.
     interface latlonrange_gatherdump; module procedure &
             gatherdump_local_, &            ! gather-dump local cvgrid lat-lon-range
             gatherdump_; end interface      ! gather-dump user provided lat-lon-range

  type latlonrange
     private
     integer(kind=i_kind):: nuse=0       !       > 0; count of luse==.true., or
                                         !         0; for no data (so no range defined), or
                                         !        -1; for a subdomain grid range.

        ! latitude range value defined within an absolute domain,
        !       [-90.,+90.]
     real(kind=r_kind):: alat_min=0._r_kind     ! the lower bound of lat
     real(kind=r_kind):: alat_max=0._r_kind     ! the upper bound of lat

        ! longitude range value is defined within a relative domain,
        !       [alon_ref-180.,alon_ref+180.]
     real(kind=r_kind):: alon_ref=0._r_kind     ! a reference lon used for all data
     real(kind=r_kind):: alon_min=0._r_kind     ! the lower bound of lon
     real(kind=r_kind):: alon_max=0._r_kind     ! the upper bound of lon
  end type latlonrange

  interface latlonrange_reset; module procedure &
     reset_,  &                  ! rank-0
     reset_r1_; end interface    ! rank-1

  interface latlonrange_enclose    ; module procedure enclose_ ; end interface
  interface latlonrange_islocal    ; module procedure islocal_ ; end interface

  interface latlonrange_set        ; module procedure set_     ; end interface
  interface latlonrange_overlaps   ; module procedure overlaps_; end interface

  interface latlonrange_gatherwrite; module procedure gatherwrite_; end interface
  interface latlonrange_readbcast  ; module procedure readbcast_  ; end interface

! Usecase I: Create a file of all-ranges, for meta info of the all-obsdiags outputs.
!
!       use m_latlonrange, only: latlonrange
!       use m_latlonrange, only: latlonrange_reset
!       use m_latlonrange, only: latlonrange_enclose
!       use m_latlonrange, only: latlonrange_gatherwrite()
!       type(latlonrange):: orange
!
!       call latlonrange_reset(orange)
!       do .. in an obstypes node entries
!          if(obsnode%isluse()) &
!            call latlonrange_enclose(orange, obsnode%elat,obsnode%elon)
!       enddo
!  
!       call latlonrange_gatherwrite(orange,'obsdiags.01',root,comm)
!       call latlonrange_reset(orange)
!
! Usecase II: Read in a file of all-ranges, to identify the locality before the
!       read of any obsdiags file.
!
!       type(latlonrange):: myrange
!       type(latlonrange),allocatable,dimension(:):: allranges
!
!       allocatge(allranges(0:mpes-1)
!       call latlonrange_reset(allranges(:))
!       call latlonrange_readbcast('obsdiags.01.headers',allranges(:),comm,root)
!       call latlonrange_set(myrange, ...)
!
!       do ipe=lpe,upe
!          if(latlonrange_overlaps(myrange,allranges(ipe)) then
!             call read_(...)
!          endif
!       enddo
!
!       call latlonrange_reset(allranges)
!       call latlonrange_reset( myrange )
!       deallocate(allranges)

  real(r_kind),parameter:: alat_np=+90._r_kind
  real(r_kind),parameter:: alat_sp=-90._r_kind

  real(r_kind),parameter:: deg_in_rad=4._r_kind*atan(1._r_kind)/180._r_kind
  real(r_kind),parameter:: rad_in_deg=1._r_kind/deg_in_rad
  real(r_kind),parameter:: r360deg   =360._r_kind

  integer(i_kind),parameter:: root=0

  integer(kind=mpi_ikind),save:: mpi_type_=0
  logical,save:: mpi_type_defined_=.false.

  type(latlonrange),save:: localrange_ = latlonrange()
  logical,save:: localrange_defined_   =.false.

  logical,parameter:: alwayslocal_   =.true.

#include "myassert.H"

#define _TIMER_ON_
#ifdef  _TIMER_ON_
#undef  _TIMER_ON_
#undef  _TIMER_OFF_
#define _TIMER_ON_(id)  call timer_ini(id)
#define _TIMER_OFF_(id) call timer_fnl(id)
#else
#define _TIMER_ON_(id)
#define _TIMER_OFF_(id)
#endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname='m_latlonrange'

contains
pure function londiff_(clon,glon) result(dlon)
        ! londiff_ = dlon := (glon.minus.clon +180.) .modulo. 360. -180.
        ! Note: it is not clon.minus.glon, but glon.minus.clon.
  implicit none
  real(kind=r_kind):: dlon      ! is defined in range [-180.,180.].
                ! The actual result for an argument at an end point could
                ! be either -180 or +180, depending on the value of glon
                ! relative to clon.

  real(kind=r_kind),intent(in):: clon,glon

  real(kind=kind(dlon)):: rad_clon,cos_clon,sin_clon
  real(kind=kind(dlon)):: rad_glon,cos_glon,sin_glon

  !print'(a,2f10.3)','deg_in_rad=',deg_in_rad,3.1415926/180.
  !print'(a,2f10.3)','rad_in_deg=',rad_in_deg,180.00/3.1415926

  rad_clon = clon *deg_in_rad
  cos_clon = cos(rad_clon)
  sin_clon = sin(rad_clon)

  rad_glon = glon *deg_in_rad
  cos_glon = cos(rad_glon)
  sin_glon = sin(rad_glon)

  ! dlon = glon .minus. clon, in a trignomitry based algorithm.

  dlon = atan2( cos_clon*sin_glon - sin_clon*cos_glon,  &   ! sin(g-c) = sin(g)*cos(c)-cos(g)*sin(c)
                cos_clon*cos_glon + sin_clon*sin_glon)      ! cos(g-c) = cos(g)*cos(c)+sin(g)*sin(c)
  dlon = dlon *rad_in_deg                                   ! d := "g-c" = atan2(sin(g-c),cos(g-c))
end function londiff_

pure subroutine reset_(llrange)
        ! set llrange to null.
  implicit none
  type(latlonrange),intent(out):: llrange
  llrange=latlonrange()
end subroutine reset_
pure subroutine reset_r1_(llranges)
  implicit none
  type(latlonrange),dimension(:),intent(out):: llranges
  llranges(:)=latlonrange()
end subroutine reset_r1_

pure subroutine set_(llrange,alat1,alon1,alat2,alon2)
        ! set llrange for the local lat-lon subdomain
  implicit none
  type(latlonrange),intent(out):: llrange
  real(kind=r_kind),intent( in):: alat1,alon1   ! corner #1
  real(kind=r_kind),intent( in):: alat2,alon2   ! corner #2

  llrange=latlonrange() ! reinitialize to zero
  call enclose_(llrange,alat1,alon1)
  call enclose_(llrange,alat2,alon2)
  llrange%nuse=-1       ! it suggests this range is for a grid.
end subroutine set_

pure subroutine enclose_(llrange,alat,alon)
  implicit none
  type(latlonrange),intent(inout):: llrange
  real(kind=r_kind),intent(in   ):: alat,alon

  real(kind=kind(alat)):: alat_,alon_
  real(kind=kind(alat)),parameter:: alon_halfinterval=360._r_kind*epsilon(1._r_kind)

  if(llrange%nuse==0) then
     llrange%alat_min=alat
     llrange%alat_max=alat

     llrange%alon_ref=alon
     llrange%alon_min=alon -alon_halfinterval
     llrange%alon_max=alon +alon_halfinterval

  else
     alat_=max(alat_sp,min(alat,alat_np))
     llrange%alat_min=min(llrange%alat_min,alat_)
     llrange%alat_max=max(llrange%alat_max,alat_)
    
     alon_ = llrange%alon_ref + londiff_(llrange%alon_ref,alon)
     llrange%alon_min=min(llrange%alon_min,alon_-alon_halfinterval)
     llrange%alon_max=max(llrange%alon_max,alon_+alon_halfinterval)
  endif

  llrange%nuse=llrange%nuse+1
end subroutine enclose_

subroutine localrange_config_()
!-- Config localrange_ for its lat-lon range.
  use m_cvgridlookup, only: cvgridlookup_sdget
  implicit none

  real(r_kind):: sdlat_ref ,sdlon_ref    ! a reference lat-lon grid point
  real(r_kind):: sdlat_lbnd,sdlon_lbnd   ! lower-left corner
  real(r_kind):: sdlat_ubnd,sdlon_ubnd   ! upper-right corner

  call cvgridlookup_sdget( sdlat_ref =sdlat_ref , sdlon_ref =sdlon_ref , &
                           sdlat_lbnd=sdlat_lbnd, sdlon_lbnd=sdlon_lbnd, &
                           sdlat_ubnd=sdlat_ubnd, sdlon_ubnd=sdlon_ubnd  )

  localrange_=latlonrange() ! reinitialize
  call enclose_(localrange_,sdlat_ref ,sdlon_ref )      ! include a reference grid point
  call enclose_(localrange_,sdlat_lbnd,sdlon_lbnd)      ! include the lower-left corner
  call enclose_(localrange_,sdlat_ubnd,sdlon_ubnd)      ! include the upper-right corner

                ! This algorithm does not assume sdlat_lbnd < sdlat_ubnd or
                ! sdlon_lbnd < sdlon_ubnd.

  localrange_defined_=.true.
end subroutine localrange_config_

function islocal_(r)
!-- is the given latlonrange::r local, against *localrange_*?
  implicit none
  logical:: islocal_
  type(latlonrange),intent(in):: r

  if(.not.localrange_defined_) call localrange_config_()

  islocal_=overlaps_(localrange_,r)

    !-- Or for debugging purposes (see module parameter alwayslocal_ defined on
    !-- the top), let
    !   islocal_ = alwayslocal_.or.overlaps_(localrange_,r)
end function islocal_

function overlaps_(arange,brange)
        ! overlaps_(ar,br) := br .overlaps. ar
  implicit none
  logical:: overlaps_
  type(latlonrange),intent(in):: arange
  type(latlonrange),intent(in):: brange

        ! An overlaping status of a pair of ranges, requires neither to be null
  overlaps_ = arange%nuse /= 0 .and. brange%nuse /= 0 

        ! In some cases, both boundaries of the larger range may be outside
        ! the smaler range.  So we need to swap the two ranges to find the
        ! right result.
  if(overlaps_) overlaps_ = lat_overlaps_(arange,brange) .or. &
                            lat_overlaps_(brange,arange)

        ! an overlaping status is both overlaping-in-lat and overlaping-in-lon.
  if(overlaps_) overlaps_ = lon_overlaps_(arange,brange) .or. &
                            lon_overlaps_(brange,arange)

contains
pure function lat_overlaps_(ar,br) result(overlaps_)
        ! at least one of boundary points of br, falls in the range of ar
  implicit none
  logical:: overlaps_
  type(latlonrange),intent(in):: ar,br

  real(kind=r_kind):: alat_

  alat_       = br%alat_min
  overlaps_   = alat_>=ar%alat_min .and. &
                alat_<=ar%alat_max
  if(.not.overlaps_) then
     alat_     = br%alat_max
     overlaps_ = alat_>=ar%alat_min .and. &
                 alat_<=ar%alat_max
  endif
end function lat_overlaps_

pure function lon_overlaps_(ar,br) result(overlaps_)
        ! at least one of boundary points of br, falls in the range of ar
  implicit none
  logical:: overlaps_
  type(latlonrange),intent(in):: ar,br

  real(kind=r_kind):: aref_,alon_

  aref_ = ar%alon_ref
  alon_       = aref_ + londiff_(aref_,br%alon_min)
  overlaps_   = alon_ >= ar%alon_min .and. &
                alon_ <= ar%alon_max
  if(.not.overlaps_) then
     alon_     = aref_ + londiff_(aref_,br%alon_max)
     overlaps_ = alon_ >= ar%alon_min .and. &
                 alon_ <= ar%alon_max
  endif
end function lon_overlaps_
end function overlaps_

subroutine gatherwrite_(llrange,hdfile,root,comm)
        ! gather all llranges from all pes of this comm, then write them all to
        ! a "header" file (hdfile="obsdiags.##.headers") by the root pe.

  use mpeu_mpif, only: mpi_type
  use gsi_unformatted, only: unformatted_open
  use mpimod, only: npes => npe
  use mpimod, only: mype
  implicit none
  type(latlonrange   ),intent(in):: llrange     ! a laglonrange to write
  character(len=*    ),intent(in):: hdfile      ! a filename to write to
  integer(kind=i_kind),intent(in):: root
  integer(kind=i_kind),intent(in):: comm
  
  character(len=*),parameter:: myname_=myname//"::gatherwrite_"
  integer(kind=i_kind):: ier,lu
  integer(kind=i_kind):: lsize,isize,rsize,mtype
  integer(kind=i_kind):: irec,nrec
  integer(kind=i_kind),            dimension(1  ):: isend
  integer(kind=i_kind),allocatable,dimension(:,:):: irecv
  real   (kind=r_kind),            dimension(5  ):: rsend
  real   (kind=r_kind),allocatable,dimension(:,:):: rrecv
_TIMER_ON_(myname_)

  lsize=0
  if(mype==root) lsize=npes
  isize   =size(isend)
  rsize   =size(rsend)
  allocate( irecv(isize,0:lsize-1), &
            rrecv(rsize,0:lsize-1)  )

  isend(1)=llrange%nuse
  mtype   =mpi_type(isend)
  call mpi_gather(isend,isize,mtype, irecv,isize,mtype, root,comm,ier)
  ASSERT(ier==0)

  rsend(:) = (/ llrange%alat_min, &
                llrange%alat_max, &
                llrange%alon_ref, & 
                llrange%alon_min, &
                llrange%alon_max  /)
  mtype   =mpi_type(rsend)
  call mpi_gather(rsend,rsize,mtype, rrecv,rsize,mtype, root,comm,ier)
  ASSERT(ier==0)

  if(mype==root) then
     call unformatted_open(unit=lu, &
                           file=trim(hdfile), &
                           class='.obsdiags.', &
                           action='write', &
                           status='unknown', &
                           newunit=.true., &       ! with newunit=.true., unit returns a value assigned by Fortran.
                           iostat=ier,silent=.true.)
     if(ier/=0) then
        call perr(myname_,'unformatted_open(), file =',trim(hdfile))
        call perr(myname_,'                 newunit =',lu)
        call perr(myname_,'                  iostat =',ier)
        call  die(myname_)
     endif

     nrec=size(irecv,2)
     write(lu) nrec
     do irec=0,nrec-1
        write(lu) irec, irecv(1:isize,irec),rrecv(1:rsize,irec)
     enddo
     close(lu)
  endif

  deallocate(irecv,rrecv)
_TIMER_OFF_(myname_)
end subroutine gatherwrite_

subroutine readbcast_(hdfile,allranges,root,comm)
  use gsi_unformatted, only: unformatted_open
  use mpeu_mpif, only: mpi_type
  use mpimod, only: mype
  implicit none
  character(len=*    ),intent(in):: hdfile      ! input file
  type(latlonrange),dimension(0:),intent(out):: allranges ! data received by all pes.
  integer(kind=i_kind),intent(in):: root        ! the root pe does the read().
  integer(kind=i_kind),intent(in):: comm        ! the communicator of my "world".

  character(len=*),parameter:: myname_=myname//"::readbcast_"
  integer(kind=i_kind):: ier,lu,irec,jrec,nrec
  integer(kind=i_kind),allocatable,dimension(:,:):: ibufr
  real   (kind=r_kind),allocatable,dimension(:,:):: rbufr
_TIMER_ON_(myname_)
  
  nrec=0
  if(mype==root) then
     call unformatted_open(unit=lu, &
                           file=trim(hdfile), &
                           class='.obsdiags.', &
                           action='read', &
                           status='old', &
                           newunit=.true., &       ! with newunit=.true., unit returns a value assigned by Fortran.
                           iostat=ier,silent=.true.)
     if(ier/=0) then
        call perr(myname_,'unformatted_open(), file =',trim(hdfile))
        call perr(myname_,'                 newunit =',lu)
        call perr(myname_,'                  iostat =',ier)
        call  die(myname_)
     endif

     read(lu) nrec
  endif

  call mpi_bcast(nrec,1,mpi_type(nrec),root,comm,ier)
  ASSERT(ier==0)

  allocate( ibufr(1,0:nrec-1), &
            rbufr(5,0:nrec-1)  )

  if(mype==root) then
     do irec=0,nrec-1
        read(lu) jrec, ibufr(1,irec),rbufr(1:5,irec)
        ASSERT(irec==jrec)
     enddo
     close(lu)
  endif

  call mpi_bcast(ibufr,size(ibufr),mpi_type(ibufr),root,comm,ier)
  ASSERT(ier==0)

  call mpi_bcast(rbufr,size(rbufr),mpi_type(rbufr),root,comm,ier)
  ASSERT(ier==0)

  ASSERT(nrec==size(allranges))
  do irec=0,nrec-1
     allranges(irec)=latlonrange( nuse=ibufr(1,irec), &
                              alat_min=rbufr(1,irec), &
                              alat_max=rbufr(2,irec), &
                              alon_ref=rbufr(3,irec), &
                              alon_min=rbufr(4,irec), &
                              alon_max=rbufr(5,irec)  )
  enddo
  deallocate(ibufr,rbufr)
_TIMER_OFF_(myname_)
end subroutine readbcast_

subroutine alldump_(allranges,varname)
!-- simply text-dump a local allranges to stdout
  use mpeu_util, only: stdout
  use mpeu_util, only: stdout_lead
  use mpeu_util, only: stdout_open
  use mpeu_util, only: stdout_close
  implicit none
  type(latlonrange), dimension(0:), intent(in):: allranges
  character(len=* ), intent(in):: varname

  character(len=*),parameter:: myname_=myname//"::alldump_"
  integer(i_kind):: irec
  character(len=:), allocatable:: varlead_
_TIMER_ON_(myname_)

  varlead_=stdout_lead(varname)

  write  (stdout,'(a,i4,l2,i8,5f12.4)') varlead_,   -1,islocal_(localrange_      ),localrange_
  write  (stdout,'(a,1x,a,i8)')         varlead_,'NPES =',size(allranges)
  do irec=0,size(allranges)-1
     write(stdout,'(a,i4,l2,i8,5f12.4)') varlead_, irec,islocal_(  allranges(irec)),allranges(irec)
  enddo
  write  (stdout,'(a,1x,a)')            varlead_,'========'
_TIMER_OFF_(myname_)
end subroutine alldump_

subroutine gatherdump_local_(varname,root,comm)
!-- gather-dump the internal lat-lon-range of the cvgrid, localrange_
  implicit none
  character(len=*    ),intent(in):: varname     ! identity of the range
  integer(kind=i_kind),intent(in):: root
  integer(kind=i_kind),intent(in):: comm

  character(len=*),parameter:: myname_=myname//"::gatherdump_local_"
_TIMER_ON_(myname_)
  if(.not.localrange_defined_) call localrange_config_()
  call gatherdump_(localrange_,varname,root,comm)
_TIMER_OFF_(myname_)
end subroutine gatherdump_local_

subroutine gatherdump_(llrange,varname,root,comm)
!-- gather all llranges from all pes of this comm, then text dump the gathered
!-- array of llranges to stdout.

  use mpeu_mpif, only: mpi_type
  use mpeu_util, only: stdout
  use mpeu_util, only: stdout_lead
  use mpimod, only: npes => npe
  use mpimod, only: mype
  implicit none
  type(latlonrange   ),intent(in):: llrange     ! a laglonrange to write
  character(len=*    ),intent(in):: varname       ! identity of the range
  integer(kind=i_kind),intent(in):: root
  integer(kind=i_kind),intent(in):: comm
  
  character(len=*),parameter:: myname_=myname//"::gatherdump_"
  integer(kind=i_kind):: ier
  integer(kind=i_kind):: lsize,isize,rsize,mtype
  integer(kind=i_kind):: irec,nrec
  integer(kind=i_kind),            dimension(1  ):: isend
  integer(kind=i_kind),allocatable,dimension(:,:):: irecv
  real   (kind=r_kind),            dimension(5  ):: rsend
  real   (kind=r_kind),allocatable,dimension(:,:):: rrecv
  character(len=:),allocatable:: varlead_
_TIMER_ON_(myname_)

  lsize=0
  if(mype==root) lsize=npes
  isize   =size(isend)
  rsize   =size(rsend)
  allocate( irecv(isize,0:lsize-1), &
            rrecv(rsize,0:lsize-1)  )

  isend(1)=llrange%nuse
  mtype   =mpi_type(isend)
  call mpi_gather(isend,isize,mtype, irecv,isize,mtype, root,comm,ier)
  ASSERT(ier==0)

  rsend(:) = (/ llrange%alat_min, &
                llrange%alat_max, &
                llrange%alon_ref, & 
                llrange%alon_min, &
                llrange%alon_max  /)
  mtype   =mpi_type(rsend)
  call mpi_gather(rsend,rsize,mtype, rrecv,rsize,mtype, root,comm,ier)
  ASSERT(ier==0)

  varlead_=stdout_lead(myname,varname)
  if(mype==root) then
     nrec=size(irecv,2)
     write  (stdout,'(a,1x,a,i8)')     varlead_,'nrec =',nrec
     do irec=0,nrec-1
        write(stdout,'(a,2i8,5f12.4)')  varlead_,irec, irecv(1:isize,irec),rrecv(1:rsize,irec)
     enddo
     write  (stdout,'(a,1x,a)')        varlead_,'=========='
  endif

  deallocate(irecv,rrecv)
_TIMER_OFF_(myname_)
end subroutine gatherdump_

end module m_latlonrange
!.
