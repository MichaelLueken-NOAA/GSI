module m_extozone
!#define NO_EXTRA_
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    module m_extozone
!   prgmmr:      j guo <jguo@nasa.gov>
!      org:      NASA/GSFC, Global Modeling and Assimilation Office, 900.3
!     date:      2013-09-27
!
! abstract: a module for reading extra ozone observation data
!
! program history log:
!   2012-09-10  wargan  - add omi text option for omi with efficiency factors
!   2012-12-06  jin     - read mls o3lev data from bufr files.
!                         keep's kris's ascii reading and rename the type to 'o3levtext'.
!   2012-12-14  jin/Wargan - fixed mflg error under "o3levtext" handling block.
!   2013-01-18  Wargan  - read omieffnc from a netcdf file, and removed omieff text.
!   2013-02-05  guo     - stop statements were replaced with die() to signal an _abort_.
!                       - dfile_format() is used to handle formats of "omieff" and "o3lev".
!                       - code is slightly cleaned up for possible restructuring
!                         with alternative format handling.
!   2013-09-27  guo     - added this document block.
!   2014-02-03  guo     - restructured from gmao extended read_ozone() routine
!                         as a separated module for extended ozone obstypes, in
!                         particular, if not in bufr format.
!                       - removed "text" option of "o3lev", for it is not been
!                         used anymore.
!                       - moved history log messages here from read_ozone.
!   2015-09-17  Thomas  - add l4densvar and thin4d to data selection procedure
!   2015-10-01  guo     - consolidate use of ob location (in deg)
!   2016-09-19  guo     - moved function dfile_format() here from obsmod.F90.
!   2018-05-21  j.jin   - added time-thinning. Moved the checking of thin4d into satthin.F90.  
!   2018-05-25  wargan  - added omps and hooks for lims, uars mls, mipas
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

  use kinds, only: i_kind,r_kind,r_double
  use obsmod, only: time_window_max
  implicit none
  private     ! except

  public:: is_extozone
  public:: extozone_read

!!  public:: extozone_setupoz
!!  public:: extozone_setupozlev
!!  public:: extozone_setupoztot

!!  public:: extozone_intoz
!!  public:: extozone_intozlev
!!  public:: extozone_intoztot

  interface is_extozone; module procedure is_extozone_; end interface
  interface extozone_read; module procedure read_; end interface


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*) , parameter:: myname='m_extozone'

  integer(kind=i_kind), parameter:: maxobs_ = 1000000

  real   (kind=r_kind), parameter:: rmiss = -9999.9_r_kind
  real   (kind=r_kind), parameter:: badoz = 10000.0_r_kind
  real   (kind=r_kind), parameter:: r6    =     6.0_r_kind
  real   (kind=r_kind), parameter:: r360  =   360.0_r_kind
  real   (kind=r_kind) :: ptime,timeinflat,crit0
  integer(kind=i_kind) :: ithin_time,n_tbin,it_mesh

contains
function is_extozone_(dfile,dtype,dplat,class)

  use mpeu_util, only: die,perr
  implicit none
  logical:: is_extozone_              ! this is the function return variable
  character(len=*),intent(in):: dfile   ! observation input filename (see gsiparm.anl:&obsinput/)
  character(len=*),intent(in):: dtype   ! observation type (see gsiparm.anl:&obsinput/)
  character(len=*),intent(in):: dplat   ! platform (see gsiparm.anl:&obsinput/)
  character(len=*),optional,intent(in):: class  ! specify either "level" or "layer" for sub-class

! Use case 1: (in read_obs())
!
!   elseif(is_extozone(dfile,obstype,dplat)) then
!     call read_ozone(obstype,dplat,dsis,dfile,...)
!
!
! Use case 2: (in setuprhsall())
!
!   elseif(is_extozone(dfile,obstype,dplat,class='level')) then
!     call setupozlev(obstype,dplat,dsis,...)
!   elseif(is_extozone(dfile,obstype,dplat,class='layer')) then
!     call setupozlay(obstype,dplat,dsis,...)
!   elseif(is_extozone(dfile,obstype,dplat,class='total')) then
!     call setupozlay(obstype,dplat,dsis,...)
!
! Use case 3: (where intoz() or stpoz() are called)
!
!   elseif(is_extozone(dfile,obstype,dplat)) then
!     call intoz(obstype,dplat,dsis,...)

  character(len=*), parameter:: myname_=myname//'::is_extozone_'

  integer(kind=i_kind),parameter:: iany    = 0
  integer(kind=i_kind),parameter:: iunknown=-1

  integer(kind=i_kind),parameter:: ilevel  = 1
  integer(kind=i_kind),parameter:: ilayer  = 2
  integer(kind=i_kind),parameter:: itotal  = 3

  integer(kind=i_kind),parameter:: itext   = 1
  integer(kind=i_kind),parameter:: ibufr   = 2
  integer(kind=i_kind),parameter:: inc     = 3

  integer(kind=i_kind):: class_,ifile_

  ifile_=iunknown
  select case(dfile_format(dfile))
     case("text")
        ifile_=itext
     case("bufr")
        ifile_=ibufr
     case("nc")
        ifile_=inc
  end select

  class_=iany
  if(present(class)) then
     select case(class)
        case('level')
           class_=ilevel
        case('layer')
           class_=ilayer
        case('total')
           class_=itotal
        case default
           class_=iunknown
           call perr(myname_,'unknown ozone class, class =',class)
           call perr(myname_,'                     dfile =',dfile)
           call perr(myname_,'                     dtype =',dtype)
           call perr(myname_,'                     dplat =',dplat)
           call  die(myname_)
     end select
  endif
    
  is_extozone_= .false.
#ifndef NO_EXTRA_
  select case(class_)
     case(iany)
        is_extozone_= &
          ifile_==ibufr .and. dtype == 'o3lev'     .or. &
          ifile_==inc   .and. dtype == 'mls55'     .or. &
          ifile_==inc   .and. dtype == 'ompslpvis' .or. &
          ifile_==inc   .and. dtype == 'ompslpuv'  .or. &
          ifile_==inc   .and. dtype == 'ompslp'    .or. &
          ifile_==inc   .and. dtype == 'lims'      .or. &
          ifile_==inc   .and. dtype == 'uarsmls'   .or. &
          ifile_==inc   .and. dtype == 'mipas'     .or. &
          ifile_==inc   .and. dtype == 'omieff'    .or. &
          ifile_==inc   .and. dtype == 'tomseff'

     case(ilevel)
        is_extozone_= &
          ifile_==ibufr .and. dtype == 'o3lev'     .or. &
          ifile_==inc   .and. dtype == 'mls55'     .or. &
          ifile_==inc   .and. dtype == 'ompslpvis' .or. &
          ifile_==inc   .and. dtype == 'ompslpuv'  .or. &
          ifile_==inc   .and. dtype == 'ompslp'    .or. &
          ifile_==inc   .and. dtype == 'lims'      .or. &
          ifile_==inc   .and. dtype == 'uarsmls'   .or. &
          ifile_==inc   .and. dtype == 'mipas'

     case(ilayer)
        is_extozone_= .false.

     case(itotal)
        is_extozone_= &
          ifile_==inc   .and. dtype == 'omieff'  .or. &
          ifile_==inc   .and. dtype == 'tomseff'

     case default
        is_extozone_= .false.
  end select
#endif

end function is_extozone_

! ----------------------------------------------------------------------
function dfile_format(dfile) result(dform)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    function dfile_format
!   prgmmr:      j guo <jguo@nasa.gov>
!      org:      NASA/GSFC, Global Modeling and Assimilation Office, 610.1
!     date:      2013-02-04
!
! abstract: - check filename suffix to guess its format
!
! program history log:
!   2013-02-04  j guo   - added this document block
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

! function interface:

  implicit none

  character(len=len('unknown')):: dform ! a 2-4 byte code for a format guess,
  ! from a list of filename suffixes, 'bufr', 'text' (also for 'txt',
  ! 'tcvitle', or 'vitl'), 'nc', or return a default value 'unknown'.  One
  ! can extend the list to other suffixes, such as 'hdf', 'hdf4', 'hdf5',
  ! etc., if they are needed in the future.
  character(len=*),intent(in):: dfile  ! a given filename

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname_=myname//'::dfile_format'
  integer(i_kind):: i,l

  dform='unknown'
  l=len_trim(dfile)

  i=max(0,l-6)+1      ! 6 byte code?
  select case(dfile(i:l))
     case ('tcvitl')
        dform='text'
  end select
  if(dform/='unknown') return

  i=max(0,l-4)+1! 4 byte code?
  select case(dfile(i:l))
     case ('bufr','avhr')
        dform='bufr'
     case ('text','vitl')
        dform='text'
  end select
  if(dform/='unknown') return

  i=max(0,l-3)+1   ! 3 byte code?
  select case(dfile(i:l))
     case ('txt')    ! a short
        dform='text'
  end select
  if(dform/='unknown') return

  i=max(0,l-2)+1    ! 2 byte code?
  select case(dfile(i:l))
     case ('nc')
        dform='nc'
  end select
  if(dform/='unknown') return

  dform='bufr'
  return
end function dfile_format

subroutine read_(dfile,dtype,dplat,dsis, &      ! intent(in), keys for type managing
  nread,npuse,nouse, &                          ! intent(out), beside other implicit output variables
  jsatid,gstime,lunout,twind,ithin,rmesh,nobs)       ! intent(in), beside other implicit input variables

  use constants, only: zero
  use mpeu_util, only: die,perr,tell
  use mpimod, only: npe
  use satthin, only: satthin_time_info => radthin_time_info
!  use mpeu_util, only: mprefix,stdout
!     nobs     - array of observations on each subdomain for each processor

  implicit none
  character(len=*), intent(in):: dfile   ! obs_input filename
  character(len=*), intent(in):: dtype   ! observation type (see gsiparm.anl:&obsinput/)
  character(len=*), intent(in):: dplat   ! platform (see gsiparm.anl:&obsinput/)
  character(len=*), intent(in):: dsis    ! sensor/instrument/satellite tag (see gsiparm.anl:&obsinput/ and ozinfo)

! Use case: (in read_ozone())
!
!   elseif(is_extozone(dfile,dtype,dplat)) then
!     call extozone_read(dtype,dplat,dsis,dfile,...)
!

  character(len=*), parameter:: myname_=myname//'::read_'

  integer(kind=i_kind), intent(out):: nread     ! number of obs record reads in this call
  integer(kind=i_kind),dimension(npe), intent(inout):: nobs     ! number of obs record reads in this call
  integer(kind=i_kind), intent(out):: npuse     ! nnmber of "preofiles" retained in this call
  integer(kind=i_kind), intent(out):: nouse     ! nnmber of obs retained in this call

  ! jsatid,gstime,lunout,twind,ithin,rmesh
  character(len=*)    , intent(in ):: jsatid    ! platform id (verification)
  real   (kind=r_kind), intent(in ):: gstime    ! analysis time (minute) from reference date
  integer(kind=i_kind), intent(in ):: lunout    ! logical unit to send obs output.
  real   (kind=r_kind), intent(in ):: twind     ! input group time window (hour)

  integer(kind=i_kind), intent(in ):: ithin     ! flag to thin data
  real   (kind=r_kind), intent(in ):: rmesh     ! thining mesh size (km)

  real(kind=r_kind),pointer,dimension(:,:):: p_out

  integer(kind=i_kind):: nreal,nchan,ilat,ilon
  integer(kind=i_kind):: maxobs
  integer(kind=i_kind):: i,k

  nread=-1
  npuse=-1
  nouse=-1

  if(.not.is_extozone(dfile,dtype,dplat)) then
     call perr(myname_,'unexpected use, dfile =',dfile)
     call perr(myname_,'                dtype =',dtype)
     call perr(myname_,'                dplat =',dplat)
     call die (myname_)
  endif
  call satthin_time_info(dtype, dplat, dsis, ptime, ithin_time)
  if( ptime > 0.0_r_kind) then
     n_tbin=nint(2.0_r_kind*time_window_max/ptime)
  else
     n_tbin=1
  endif

  select case(dtype)
     case('omieff','tomseff')       ! layer-ozone or total-ozone types
        select case(dfile_format(dfile))
           case('nc')
              call oztot_ncinquire_(nreal,nchan,ilat,ilon, &
                                    ithin,rmesh,maxobs)

              allocate(p_out(nreal+nchan,maxobs))
              p_out(:,:)=rmiss

              call oztot_ncread_(dfile,dtype,dsis, p_out,nread,npuse,nouse, &
                                 gstime, twind, ithin)

                ! Skip all "thinned" data records, and reset values of
                ! nread, npuse, and nouse, as they are required by
                ! upper level read_obs().

              k=0
              do i=1,maxobs
                 if(p_out(1,i)>zero) then
                    k=k+1
                    if(i>k) p_out(:,k)=p_out(:,i)
                 endif
              enddo
              nouse=k
              npuse=k
        end select

     case('o3lev')         ! level-ozone types
        select case(dfile_format(dfile))
!           case('text')
!              call ozlev_textinquire(dfile,dtype,dplat,  &
!                                     nreal,nchan,ilat,ilon, maxobs)
!
!              allocate(p_out(nreal+nchan,maxobs))
!              p_out(:,:)=rmiss
!
!              call ozlev_textread_(dfile,dtype,dplat,dsis, p_out,nread,npuse,nouse, &
!                                   gstime,twind, nreal,nchan,ilat,ilon)
!
           case('bufr')
              call ozlev_bufrinquire_(nreal,nchan,ilat,ilon,maxobs)

              allocate(p_out(nreal+nchan,maxobs))
              p_out(:,:)=rmiss

              call ozlev_bufrread_(dfile,dtype,dsis, p_out,nread,npuse,nouse, &
                                   jsatid, gstime,twind)
        end select

     case('mls55','ompslpvis','ompslpuv','ompslp','lims','uarsmls','mipas')
        select case(dfile_format(dfile))
           case('nc')
              call ozlev_ncinquire_( nreal,nchan,ilat,ilon,maxobs)

              allocate(p_out(nreal+nchan,maxobs))
              p_out(:,:)=rmiss
 
              call ozlev_ncread_(dfile,dtype, p_out,nread,npuse,nouse, gstime,twind)
 
        end select

  end select

  if(nouse<0 .or. .not.associated(p_out)) then
     call perr(myname_,'can not process, dtype =',trim(dtype))
     call perr(myname_,'                 dplat =',trim(dplat))
     call perr(myname_,'                  dsis =',trim(dsis))
     call perr(myname_,'                 dfile =',trim(dfile))
     call perr(myname_,'     associated(p_out) =',associated(p_out))
     if(associated(p_out)) then
        call perr(myname_,'  actual size(p_out,1) =',size(p_out,1))
        call perr(myname_,'  actual size(p_out,2) =',size(p_out,2))
        call perr(myname_,'                 nread =',nread)
        call perr(myname_,'                 npuse =',npuse)
        call perr(myname_,'                 nouse =',nouse)
     endif
     call  die(myname_)
  endif

  if(nreal+nchan/=size(p_out,1)) then
     call perr(myname_,'unexpected size(p_out,1), nreal+nchan =',nreal+nchan)
     call perr(myname_,'                                nreal =',nreal)
     call perr(myname_,'                                nchan =',nchan)
     call perr(myname_,'                 actual size(p_out,1) =',size(p_out,1))
     call perr(myname_,'                 actual size(p_out,2) =',size(p_out,2))
     call  die(myname_)
  endif

! Output candidate observations from _dfile_ to _obsfile_ (pre-opened in lunout).

!  write(stdout,'(3a,3i8,f8.2)') mprefix('read_ozone'), &
!     ' obstype,nread,npuse,nouse,no/npuse = ',dtype,nread,npuse,nouse,real(nouse)/npuse

        ! While nreal+nchan defines the leading dimension of p_out(:,:), it is
        ! obviously a bad idea that n_out has been missing from the header for
        ! its second dimension.

  call count_obs(npuse,nreal+nchan,ilat,ilon,p_out,nobs)
  write(lunout) dtype, dsis, nreal, nchan, ilat, ilon
  write(lunout) p_out(:,1:npuse)

  deallocate(p_out)

end subroutine read_

subroutine oztot_ncinquire_( nreal,nchan,ilat,ilon, ithin,rmesh,maxrec)
  use satthin, only: satthin_itxmax    => itxmax
  use satthin, only: satthin_makegrids => makegrids
  implicit none

  integer(kind=i_kind), intent(out):: nreal  ! number of real parameters per record
  integer(kind=i_kind), intent(out):: nchan  ! number of channels or levels per record
  integer(kind=i_kind), intent(out):: ilat   ! index to latitude in nreal parameters.
  integer(kind=i_kind), intent(out):: ilon   ! index to longitude in nreal parameters.

  integer(kind=i_kind), intent(in ):: ithin     ! flag to thin data
  real   (kind=r_kind), intent(in ):: rmesh     ! size (km) of the thinning mesh

  integer(kind=i_kind), intent(out):: maxrec    ! extimated input record count

  character(len=*), parameter:: myname_=myname//'::oztot_ncinquire_'

! Configure the record buffer for this obs class
  nreal=36
  nchan=1
  ilat=4
  ilon=3

! Make thinning grids, and to define the total record size.
  call satthin_makegrids(rmesh,ithin,n_tbin=n_tbin)
  maxrec=satthin_itxmax

end subroutine oztot_ncInquire_

!..................................................................................
subroutine oztot_ncread_(dfile,dtype,dsis, ozout,nmrecs,ndata,nodata, &
                         gstime,twind, ithin)
!..................................................................................

  use netcdf, only: nf90_open
  use netcdf, only: nf90_nowrite
  use netcdf, only: nf90_noerr
  use netcdf, only: nf90_inq_dimid
  use netcdf, only: nf90_inquire_dimension
  use netcdf, only: nf90_inq_varid
  use netcdf, only: nf90_get_var
  use netcdf, only: nf90_close

  use satthin, only: satthin_makegrids    => makegrids
  use satthin, only: satthin_tdiff2crit   => tdiff2crit
  use satthin, only: satthin_map2tgrid    => map2tgrid
  use satthin, only: satthin_finalcheck   => finalcheck
  use satthin, only: satthin_destroygrids => destroygrids

  use gridmod, only: nlat,nlon,regional,tll2xy,rlats,rlons
  use gsi_4dvar, only: l4dvar,iwinbgn,winlen,l4densvar
  use obsmod, only: nloz_omi

  use constants, only: deg2rad,zero,r60inv
! use mpeu_util, only: mprefix,stdout

  implicit none
  character(len=*), intent(in):: dfile   ! obs_input filename
  character(len=*), intent(in):: dtype   ! observation type (see gsiparm.anl:&obsinput/)
  character(len=*), intent(in):: dsis    ! sensor/instrument/satellite tag (see gsiparm.anl:&obsinput/ and ozinfo)

  real   (kind=r_kind), dimension(:,:), intent(out):: ozout
  integer(kind=i_kind), intent(out):: nmrecs ! read count in ozout
  integer(kind=i_kind), intent(out):: ndata  ! "good" profile count in ozout
  integer(kind=i_kind), intent(out):: nodata ! "good" data count in ozout

  real   (kind=r_kind), intent(in):: gstime ! analysis time (minute) from reference date
  real   (kind=r_kind), intent(in):: twind  ! input group time window (hour)


  integer(kind=i_kind), intent(in ):: ithin     ! flag to thin data

  character(len=*), parameter:: myname_=myname//'::oztot_ncread_'

! parameters for output bookkeeping
  integer(kind=i_kind):: i, irec
  integer(kind=i_kind):: itx,itt

! parameters for NetCDF arrays
  integer(kind=i_kind), allocatable :: iya(:),ima(:),idda(:),ihha(:),imina(:),fovna(:)
  integer(kind=i_kind), allocatable :: toqfa(:),alqfa(:)
  real   (kind=r_kind), allocatable :: rseca(:),slatsa(:),slonsa(:),totoza(:),szaa(:)
  real   (kind=r_kind), allocatable :: aprioria(:,:),efficiencya(:,:)
  integer(kind=i_kind)  nrecdimid,nomilevsdimid,lonvarid,latvarid,yyvarid,mmvarid
  integer(kind=i_kind)  ddvarid,hhvarid,minvarid,ssvarid,fovnvarid,toqfvarid,alqfvarid
  integer(kind=i_kind)  szavarid,totozvarid,apriorivarid,efficiencyvarid
  integer(kind=i_kind)  ier, ncid, nomilevs
  integer(kind=i_kind):: iy,im,idd,ihh,imin,nmind
  integer(kind=i_kind):: toqf,alqf,fovn
  real   (kind=r_kind):: rsec,slats,slons
  real   (kind=r_kind),dimension(nloz_omi):: apriori, efficiency
  integer(kind=i_kind):: binary(17)

  real   (kind=r_kind):: dlon,dlon_earth,dlon_earth_deg
  real   (kind=r_kind):: dlat,dlat_earth,dlat_earth_deg
  real   (kind=r_kind):: tdiff,sstime,t4dv,crit1,dist1,rsat
  integer(kind=i_kind):: idate5(5)

  real   (kind=r_double):: totoz, sza

  logical:: outside, removescans, iuse
  integer(kind=i_kind):: maxobs


  removescans = .true. ! remove omi scan numbers >= 25?
  maxobs=size(ozout,2)
  rsat=999._r_kind

! Using omi/toms with efficience factors in netcdf format
  nodata = 0
  ndata  = 0

  call check(nf90_open(trim(dfile),nf90_nowrite,ncid),stat=ier)

  ! ignore if the file is not actually present.
  if(ier/=nf90_noerr) then
     call satthin_destroygrids()
     return
  end if

  ! Get dimensions from omi input file
  call check(nf90_inq_dimid(ncid, "nrec", nrecdimid),stat=ier)

  ! ignore if the file header is empty
  if(ier/=nf90_noerr) then
     call check(nf90_close(ncid),stat=ier)
     call satthin_destroygrids()
     return
  endif

  ! Get dimensions from omi/toms input file
!  nmrecs=0
  call check(nf90_inquire_dimension(ncid, nrecdimid, len = nmrecs),stat=ier)

  ! ignore if the file header is empty
  if(ier/=nf90_noerr .or. nmrecs==0) then
     call check(nf90_close(ncid),stat=ier)
     call satthin_destroygrids()
     return
  endif

  ! Continue the input
  call check(nf90_inq_dimid(ncid, "nlevs", nomilevsdimid))
  call check(nf90_inquire_dimension(ncid, nomilevsdimid, len = nomilevs))
     
  ! We have dimensions so we can allocate arrays
  allocate(iya(nmrecs),ima(nmrecs),idda(nmrecs),ihha(nmrecs),imina(nmrecs), &
       rseca(nmrecs),fovna(nmrecs),slatsa(nmrecs),slonsa(nmrecs),totoza(nmrecs), &
       toqfa(nmrecs),alqfa(nmrecs),szaa(nmrecs))
  allocate(aprioria(nomilevs,nmrecs),efficiencya(nomilevs,nmrecs))

  ! Read variables and store them in these arrays
  call check(nf90_inq_varid(ncid, "lon", lonvarid))
  call check(nf90_get_var(ncid, lonvarid, slonsa))

  call check(nf90_inq_varid(ncid, "lat", latvarid))
  call check(nf90_get_var(ncid, latvarid, slatsa))

  call check(nf90_inq_varid(ncid, "yy", yyvarid))
  call check(nf90_get_var(ncid, yyvarid, iya))

  call check(nf90_inq_varid(ncid, "mm", mmvarid))
  call check(nf90_get_var(ncid, mmvarid, ima))

  call check(nf90_inq_varid(ncid, "dd", ddvarid))
  call check(nf90_get_var(ncid, ddvarid, idda))

  call check(nf90_inq_varid(ncid, "hh", hhvarid))
  call check(nf90_get_var(ncid, hhvarid, ihha))

  call check(nf90_inq_varid(ncid, "min", minvarid))
  call check(nf90_get_var(ncid, minvarid, imina))

  call check(nf90_inq_varid(ncid, "ss", ssvarid))
  call check(nf90_get_var(ncid, ssvarid, rseca))

  call check(nf90_inq_varid(ncid, "fov", fovnvarid))
  call check(nf90_get_var(ncid, fovnvarid, fovna))

  call check(nf90_inq_varid(ncid, "sza", szavarid))
  call check(nf90_get_var(ncid, szavarid, szaa))

  call check(nf90_inq_varid(ncid, "toqf", toqfvarid))
  call check(nf90_get_var(ncid, toqfvarid, toqfa))

  call check(nf90_inq_varid(ncid, "alqf", alqfvarid))
  call check(nf90_get_var(ncid, alqfvarid, alqfa))

  call check(nf90_inq_varid(ncid, "toz", totozvarid))
  call check(nf90_get_var(ncid, totozvarid, totoza))

  call check(nf90_inq_varid(ncid, "apriori", apriorivarid))
  call check(nf90_get_var(ncid, apriorivarid, aprioria))

  call check(nf90_inq_varid(ncid, "efficiency", efficiencyvarid))
  call check(nf90_get_var(ncid, efficiencyvarid, efficiencya))

  ! close the data file
  call check(nf90_close(ncid))
     
  ! now screen the data and put them into the right places
  recloop: do irec = 1, nmrecs
     iy = iya(irec)
     im = ima(irec)
     idd = idda(irec)
     ihh = ihha(irec)
     imin = imina(irec)
     rsec = rseca(irec)
     fovn = fovna(irec)
     slats = slatsa(irec)
     slons = slonsa(irec)
     totoz = totoza(irec)
     toqf = toqfa(irec)
     alqf = alqfa(irec)
     sza = szaa(irec)
     do i = 1, nomilevs
        apriori(i) = aprioria(i, irec)
        ! Reported efficiencies from layer 4 up (counting from the surface
        ! so from 1 - 8 here) are incorrect (too large) for high sza
        ! Setting them all to 1.0 per pk bhartia''s recommendation
        if (i <= 8) then
           efficiency(i) = 1._r_kind
        else
           efficiency(i) = efficiencya(i, irec) 
        endif
     end do
        
!     nmrecs=nmrecs+1
     !      Convert observation location to radians
     if(abs(slats)>90._r_kind .or. abs(slons)>r360) cycle recloop
     if(slons< zero) slons=slons+r360
     if(slons==r360) slons=zero
     dlat_earth_deg = slats
     dlon_earth_deg = slons
     dlat_earth = slats * deg2rad
     dlon_earth = slons * deg2rad

     if(regional)then
        call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)
        if(outside) cycle recloop
     else
        dlat = dlat_earth
        dlon = dlon_earth
        call grdcrd1(dlat,rlats,nlat,1)
        call grdcrd1(dlon,rlons,nlon,1)
     endif

!    convert observation time to relative time
     idate5(1) = iy    !year 
     idate5(2) = im    !month
     idate5(3) = idd   !day
     idate5(4) = ihh   !hour
     idate5(5) = imin  !minute
     call w3fs21(idate5,nmind)

     t4dv=real((nmind-iwinbgn),r_kind)*r60inv
     sstime=real(nmind,r_kind)
     tdiff=(sstime-gstime)*r60inv
     if (l4dvar.or.l4densvar) then
        if (t4dv<zero .OR. t4dv>winlen) cycle recloop
     else
        if(abs(tdiff) > twind) cycle recloop
     end if

     if (totoz > badoz ) cycle recloop
        
! Apply data screening based on quality flags
! Bit 10 (from the left) in TOQF represents row anomaly.  All 17 bits in toqf is
! supposed to converted into array elements of binary(:), either for "tomseff" or
! "omieff".
     binary(:)=0
     select case(dtype)
        case('omieff')

           if (toqf /= 0 .and. toqf /= 1) cycle recloop
 
!       Remove obs at high solar zenith angles
           if (sza > 84.0_r_kind) cycle recloop

!       remove the bad scan position data: fovn beyond 25
           if (removescans) then
              if (fovn >=25) cycle recloop
           endif
           if (fovn <=2 .or. fovn >=58) cycle recloop

!       remove the data in which the c-pair algorithm ((331 and 360 nm) is used. 
           if (alqf == 3 .or. alqf == 13) cycle recloop

        case('tomseff')
        ! The meaning of quality flags for toms version 8 is similar to that
        ! for sbuv:
        ! 0 - good data, 1 - good data with sza > 84 deg
           if (toqf /= 0) cycle recloop

        case default
     end select

!    thin omi/toms data

     crit0 = 0.01_r_kind 
     timeinflat=r6
     call satthin_tdiff2crit(tdiff,ptime,ithin_time,timeinflat,crit0,crit1,it_mesh)
     if (ithin /= 0) then
        call satthin_map2tgrid(dlat_earth,dlon_earth,dist1,crit1,itx,ithin,itt,iuse,dsis,it_mesh=it_mesh)
        if(.not. iuse)cycle recloop
        call satthin_finalcheck(dist1,crit1,itx,iuse)
        if(.not. iuse)cycle recloop
        ndata=ndata+1
        nodata=ndata
     else
        ndata=ndata+1
        nodata=ndata
        itx= ndata
     end if

!     ASSERT_(size(ozout,2)>=itx)

     if(itx <= maxobs) then
        ozout(1,itx)=rsat
        ozout(2,itx)=t4dv
        ozout(3,itx)=dlon               ! grid relative longitude
        ozout(4,itx)=dlat               ! grid relative latitude
        ozout(5,itx)=dlon_earth_deg     ! earth relative longitude (degrees)
        ozout(6,itx)=dlat_earth_deg     ! earth relative latitude (degrees)
        ozout(7,itx)=real(toqf)         ! - total ozone quality code (not used)
        ozout(8,itx)=real(sza)          ! solar zenith angle
        ozout(9,itx)=binary(10)         ! row anomaly flag, is actually fixed to 0
        ozout(10,itx)=0._r_kind         ! - cloud amount (not used)
        ozout(11,itx)=0._r_kind         ! - vzan (not used)
        ozout(12,itx)=0._r_kind         ! - aerosol index (not used)
        ozout(13,itx)=0._r_kind         ! - ascending/descending (not used)
        ozout(14,itx)=real(fovn)        ! scan position
                                        ! "(not used)" flags above imply that they
                                        ! are not used in setupozlay().

! Added apriori and efficiency profiles 
        ozout(15:25,itx)=apriori        
        ozout(26:36,itx)=efficiency
        ozout(37,itx)=totoz 
     endif

!     ASSERT_(size(ozout,1)==37)

  end do recloop

  deallocate(iya,ima,idda,ihha,imina, &
       rseca,fovna,slatsa,slonsa,totoza, &
       toqfa,alqfa,szaa,aprioria,efficiencya)
! end
!    End of loop over observations
! End of omi block with efficiency factors in netcdf format

  call satthin_destroygrids()

  nodata=min(ndata,maxobs)
!  write(stdout,'(3a,3i8,f8.2)') mprefix('read_ozone'), &
!     ' obstype,nmrecs,ndata,nodata,no/ndata = ',dtype,nmrecs,ndata,nodata,real(nodata)/ndata
  return

end subroutine oztot_ncread_
!..................................................................................


subroutine ozlev_ncinquire_( nreal,nchan,ilat,ilon, maxrec)
  implicit none

  integer(kind=i_kind), intent(out):: nreal  ! number of real parameters per record
  integer(kind=i_kind), intent(out):: nchan  ! number of channels or levels per record
  integer(kind=i_kind), intent(out):: ilat   ! index to latitude in nreal parameters.
  integer(kind=i_kind), intent(out):: ilon   ! index to longitude in nreal parameters.

  integer(kind=i_kind), intent(out):: maxrec    ! extimated input record count

  character(len=*), parameter:: myname_=myname//'::ozlev_ncinquire_'

! Configure the record, they are not (dfile,dtype,dplat) dependent in this case.
  nreal = 12
  nchan =  1   ! There are 'levs' levels but each is treated 
               ! as a separate observation so that nchanl = 1
  ilat=4
  ilon=3

  maxrec = maxobs_
end subroutine ozlev_ncinquire_

!..................................................................................
subroutine ozlev_ncread_(dfile,dtype,ozout,nmrecs,ndata,nodata, gstime,twind)
!..................................................................................
  use netcdf, only: nf90_open
  use netcdf, only: nf90_nowrite
  use netcdf, only: nf90_noerr
  use netcdf, only: nf90_inq_dimid
  use netcdf, only: nf90_inquire_dimension
  use netcdf, only: nf90_inq_varid
  use netcdf, only: nf90_get_var
  use netcdf, only: nf90_close

  use gridmod, only: nlat,nlon,regional,tll2xy,rlats,rlons
  use gsi_4dvar, only: l4dvar,iwinbgn,winlen,l4densvar

  use constants, only: deg2rad,zero,one_tenth,r60inv
  use ozinfo, only: jpch_oz,nusis_oz,iuse_oz
  use mpeu_util, only: perr,die
!  use mpeu_util, only: mprefix,stdout

  implicit none
  character(len=*), intent(in):: dfile   ! obs_input filename
  character(len=*), intent(in):: dtype   ! obs_input dtype

  real   (kind=r_kind), dimension(:,:), intent(out):: ozout
  integer(kind=i_kind), intent(out):: nmrecs ! count of records read
  integer(kind=i_kind), intent(out):: ndata  ! count of processed
  integer(kind=i_kind), intent(out):: nodata ! count of retained

  real   (kind=r_kind), intent(in):: gstime ! analysis time (minute) from reference date
  real   (kind=r_kind), intent(in):: twind  ! input group time window (hour)


  character(len=*), parameter:: myname_=myname//'::ozlev_ncread_'

  integer(kind=i_kind):: ier, iprof, nprofs, maxobs
  integer(kind=i_kind):: i, ilev, ikx, ncid, k0
  integer(kind=i_kind),allocatable,dimension(:):: ipos

  integer(kind=i_kind):: nrecdimid,lonvarid,latvarid,yyvarid,mmvarid
  integer(kind=i_kind):: ddvarid,hhvarid,minvarid,ssvarid
  integer(kind=i_kind):: pressvarid
  integer(kind=i_kind):: convvarid, errvarid, ozonevarid
  integer(kind=i_kind):: levsdimid,levs

  integer(kind=i_kind), allocatable :: iya(:),ima(:),idda(:),ihha(:),imina(:),iseca(:)
  real   (kind=r_kind), allocatable :: slatsa(:),slonsa(:)
  real   (kind=r_kind), allocatable :: press(:), ozone(:,:), qual(:), press2d(:,:)
  real   (kind=r_kind), allocatable :: conv(:), err(:,:)

  integer(kind=i_kind):: nmind

  real   (kind=r_kind):: slons0,slats0
  real   (kind=r_kind):: ppmv, pres, pob, obserr, usage

  real   (kind=r_kind):: dlon,dlon_earth,dlon_earth_deg
  real   (kind=r_kind):: dlat,dlat_earth,dlat_earth_deg
  real   (kind=r_kind):: tdiff,sstime,t4dv,rsat
  integer(kind=i_kind):: idate5(5)

  logical:: outside
  logical:: first

  maxobs=size(ozout,2)
  rsat=999._r_kind
  ndata = 0
  nmrecs=0
  nodata=-1
!..................................................................................
  ! ---------------mls and other limb ozone data in netcdf format ----------
  ! Open file and read dimensions
  call check(nf90_open(trim(dfile),nf90_nowrite,ncid),stat=ier)

  ! ignore if the file is not actually present.
  if(ier/=nf90_noerr) return                 

  ! Get dimensions from omi input file
  call check(nf90_inq_dimid(ncid, "nprofiles", nrecdimid),stat=ier)

  ! ignore if the file header is empty
  if(ier/=nf90_noerr) then
     call check(nf90_close(ncid),stat=ier)
     return               
  endif

  ! Get dimensions from the input file: # of profiles and # of levels
  nprofs=0
  call check(nf90_inquire_dimension(ncid, nrecdimid, len = nprofs),stat=ier)

  ! ignore if the file header is empty
  if(ier/=nf90_noerr) then
     call check(nf90_close(ncid),stat=ier)
     return                   
  endif

  if(nprofs==0) then
     nodata=0
     call check(nf90_close(ncid),stat=ier)
     return               
  endif

  ! Continue the input
  call check(nf90_inq_dimid(ncid, "nlevs", levsdimid))
  call check(nf90_inquire_dimension(ncid, levsdimid, len = levs))

  !  Note: Make sure that 'ozinfo' has the same number of levels
  ! for nrt it is 55
  allocate(ipos(levs))
  ipos=999
     
  ! Process limb data in netcdf format
  ikx = 0 
  first=.false.
  do i=1,jpch_oz
     if( (.not. first) .and. index(nusis_oz(i), trim(dtype))/=0) then
        k0=i
        first=.true.
     end if
     if(first.and.index(nusis_oz(i),trim(dtype))/=0) then 
        ikx=ikx+1
        ipos(ikx)=k0+ikx-1
     end if
  end do
    
  if (ikx/=levs) call die(myname_//': inconsistent levs for '//dtype)
     
  nmrecs=0
  ! Allocate space and read data

  allocate(iya(nprofs),ima(nprofs),idda(nprofs),ihha(nprofs),imina(nprofs), &
       iseca(nprofs),slatsa(nprofs),slonsa(nprofs), ozone(levs,nprofs), &
       err(levs,nprofs),qual(nprofs), conv(nprofs))
  if (index(dtype, 'ompslp') == 0) then 
     allocate(press(levs))
  else
     allocate(press2d(levs,nprofs))
  endif

  ! Read variables and store them in these arrays
  call check(nf90_inq_varid(ncid, "lon", lonvarid))
  call check(nf90_get_var(ncid, lonvarid, slonsa))

  call check(nf90_inq_varid(ncid, "lat", latvarid))
  call check(nf90_get_var(ncid, latvarid, slatsa))

  call check(nf90_inq_varid(ncid, "yy", yyvarid))
  call check(nf90_get_var(ncid, yyvarid, iya))

  call check(nf90_inq_varid(ncid, "mm", mmvarid))
  call check(nf90_get_var(ncid, mmvarid, ima))

  call check(nf90_inq_varid(ncid, "dd", ddvarid))
  call check(nf90_get_var(ncid, ddvarid, idda))

  call check(nf90_inq_varid(ncid, "hh", hhvarid))
  call check(nf90_get_var(ncid, hhvarid, ihha))

  call check(nf90_inq_varid(ncid, "min", minvarid))
  call check(nf90_get_var(ncid, minvarid, imina))

  call check(nf90_inq_varid(ncid, "ss", ssvarid))
  call check(nf90_get_var(ncid, ssvarid, iseca))

  if (index(dtype, 'ompslp') == 0) then
     call check(nf90_inq_varid(ncid, "press", pressvarid))
     call check(nf90_get_var(ncid, pressvarid, press))
  else
     call check(nf90_inq_varid(ncid, "press", pressvarid))
     call check(nf90_get_var(ncid, pressvarid, press2d))
  endif

  if (index(dtype, 'ompslp') == 0) then
     call check(nf90_inq_varid(ncid, "conv", convvarid))
     call check(nf90_get_var(ncid, convvarid, conv))
  end if

  !call check(nf90_inq_varid(ncid, "qual", qualvarid))
  !call check(nf90_get_var(ncid, qualvarid, qual))
     
  call check(nf90_inq_varid(ncid, "oberr", errvarid))
  call check(nf90_get_var(ncid, errvarid, err))

  call check(nf90_inq_varid(ncid, "ozone", ozonevarid))
  call check(nf90_get_var(ncid, ozonevarid, ozone))

  ! close the data file
  call check(nf90_close(ncid))
     
  ! 'Unpack' the data
  nmrecs = 0
  nodata = 0
  do iprof = 1, nprofs
     do ilev = 1, levs
        ! note that most of the quality control is done at the 
        ! pre-processing stage
!        if (press(ilev) > 262.0_r_kind .or. press(ilev) < 0.1_r_kind ) cycle ! undefined
        if (ozone(ilev, iprof) < -900.0_r_kind) cycle            ! undefined
        if (err(ilev, iprof) < -900.0_r_kind) cycle              ! undefined
        if (iuse_oz(ipos(ilev)) < 0) then
           usage = 100._r_kind
        else
           usage = zero
        endif
        nmrecs=nmrecs+1

        !       convert observation location to radians
        slons0=slonsa(iprof)
        slats0=slatsa(iprof)
        if(abs(slats0)>90._r_kind .or. abs(slons0)>r360) cycle
        if(slons0< zero) slons0=slons0+r360
        if(slons0==r360) slons0=zero
        dlat_earth_deg = slats0
        dlon_earth_deg = slons0
        dlat_earth = slats0 * deg2rad
        dlon_earth = slons0 * deg2rad

        if(regional)then
           call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)
           if(outside) cycle    
        else
           dlat = dlat_earth
           dlon = dlon_earth
           call grdcrd1(dlat,rlats,nlat,1)
           call grdcrd1(dlon,rlons,nlon,1)
        endif

        idate5(1) = iya(iprof) !year
        idate5(2) = ima(iprof) !month
        idate5(3) = idda(iprof) !day
        idate5(4) = ihha(iprof) !hour
        idate5(5) = imina(iprof) !minute
        call w3fs21(idate5,nmind)
        t4dv=real((nmind-iwinbgn),r_kind)*r60inv
        if (l4dvar.or.l4densvar) then
           if (t4dv<zero .OR. t4dv>winlen) then
              write(6,*)'read_ozone: ', dtype,' obs time idate5=',idate5,', t4dv=',&
                   t4dv,' is outside time window, sstime=',sstime*r60inv
              cycle
           end if
        else
           sstime=real(nmind,r_kind)
           tdiff=(sstime-gstime)*r60inv
           if(abs(tdiff) > twind)then
              write(6,*)'read_ozone: ',dtype,' obs time idate5=',idate5,', tdiff=',&
                   tdiff,' is outside time window=',twind
              cycle
           end if
        end if

        obserr = err(ilev, iprof)
        ppmv = ozone(ilev, iprof)
        if (index(dtype, 'ompslp') == 0) then
           pres = press(ilev) 
        else
           pres = press2d(ilev, iprof)
        end if
        pob = log(pres * one_tenth)
        ndata  = ndata+1
        if(ndata<=maxobs) then
           nodata = nodata + 1
           ozout(1,ndata)=rsat
           ozout(2,ndata)=t4dv
           ozout(3,ndata)=dlon                 ! grid relative longitude
           ozout(4,ndata)=dlat                 ! grid relative latitude
           ozout(5,ndata)=dlon_earth_deg       ! earth relative longitude (degrees)
           ozout(6,ndata)=dlat_earth_deg       ! earth relative latitude (degrees)

           ozout(7,ndata)=rmiss                ! used to be solar zenith angle
           ozout(8,ndata)=usage
           ozout(9,ndata)=pob                  ! pressure 
           ozout(10,ndata)=obserr              ! ozone mixing ratio precision in ppmv
           ozout(11,ndata)=real(ipos(ilev))   ! pointer of obs level index in ozinfo.txt
           ozout(12,ndata)=levs                ! # of  vertical levels
           ozout(13,ndata)=ppmv                ! ozone mixing ratio in ppmv
        endif
           
     end do
  end do


  deallocate(iya,ima,idda,ihha,imina,iseca,slatsa,slonsa, ozone, &
       err,qual, conv)
  if (index(dtype, 'ompslp') == 0) then
     deallocate(press)
  else
     deallocate(press2d)
  end if
  deallocate(ipos)
     
  !     write(stdout,'(3a,3i8,f8.2)') mprefix('read_ozone'), &
  !        ' obstype,nmrecs,ndata,nodata,no/ndata = ',dtype,nmrecs,ndata,nodata,real(nodata)/ndata
     
  if (ndata > maxobs) then 
     call perr(myname_,'Number of limb obs reached maxobs = ', maxobs)
     call perr(myname_,'                            ndata = ', ndata)
     call perr(myname_,'                           nodata = ', nodata)
     call die(myname_)
  endif

  !---------------End mls nrt netcdf---------------------------
end subroutine ozlev_ncread_

subroutine ozlev_bufrinquire_(nreal,nchan,ilat,ilon,maxrec)
  implicit none

  integer(kind=i_kind), intent(out):: nreal  ! number of real parameters per record
  integer(kind=i_kind), intent(out):: nchan  ! number of channels or levels per record
  integer(kind=i_kind), intent(out):: ilat   ! index to latitude in nreal parameters.
  integer(kind=i_kind), intent(out):: ilon   ! index to longitude in nreal parameters.

  integer(kind=i_kind), intent(out):: maxrec    ! extimated input record count

  character(len=*), parameter:: myname_=myname//'::ozlev_bufrinquire_'

! Configure the record, they are not (dfile,dtype,dplat) dependent in this case.
  nreal = 12
  nchan =  1
  ilat=4
  ilon=3

  maxrec = maxobs_
end subroutine ozlev_bufrinquire_

subroutine ozlev_bufrread_(dfile,dtype,dsis, ozout,nmrecs,ndata,nodata, &
                           jsatid, gstime,twind)

  use gridmod, only: nlat,nlon,regional,tll2xy,rlats,rlons
  use gsi_4dvar, only: l4dvar,iwinbgn,winlen,l4densvar

  use constants, only: deg2rad,zero,r60inv
  use ozinfo, only: jpch_oz,nusis_oz,iuse_oz
  use radinfo, only: dec2bin

  use mpeu_util, only: warn,tell
!  use mpeu_util, only: mprefix,stdout

  implicit none
  character(len=*), intent(in):: dfile   ! obs_input filename
  character(len=*), intent(in):: dtype   ! observation type (see gsiparm.anl:&obsinput/)
  character(len=*), intent(in):: dsis    ! sensor/instrument/satellite tag (see gsiparm.anl:&obsinput/ and ozinfo)

  real   (kind=r_kind), dimension(:,:), intent(out):: ozout
  integer(kind=i_kind), intent(out):: nmrecs ! count of actual records read
  integer(kind=i_kind), intent(out):: ndata  ! count of processed data
  integer(kind=i_kind), intent(out):: nodata ! count of retained data

  character(len=*)    , intent(in):: jsatid ! platform id (verification)
  real   (kind=r_kind), intent(in):: gstime ! analysis time (minute) from reference date
  real   (kind=r_kind), intent(in):: twind  ! input group time window (hour)


  character(len=*),parameter:: myname_=myname//'::ozlev_bufrread_'

  integer(kind=i_kind),parameter:: nloz =37
  integer(kind=i_kind),parameter:: lunin=10

  character(len=*), parameter:: mlstr  ='SAID CLAT CLON YEAR MNTH DAYS HOUR MINU SECO SOZA CONV MLST PCCF'
  character(len=*), parameter:: mlstrl2='PRLC OZMX OZMP OZME'

  integer(kind=i_kind):: idate,jdate,ksatid,iy,iret,im,ihh,idd
  integer(kind=i_kind):: mpos
  character(len= 8):: subset

  integer(kind=i_kind):: maxobs
  integer(kind=i_kind):: nmind
  integer(kind=i_kind):: k
  integer(kind=i_kind):: kidsat
  integer(kind=i_kind):: idate5(5)
  integer(kind=i_kind):: ikx
  integer(kind=i_kind):: decimal,binary_mls(18)

  real(kind=r_kind):: tdiff,sstime,dlon,dlat,t4dv
  real(kind=r_kind):: slons0,slats0,rsat,dlat_earth,dlon_earth
  real(kind=r_kind):: dlat_earth_deg,dlon_earth_deg
  real(kind=r_double),dimension(13):: hdrmls
  real(kind=r_double),dimension(4,37):: hdrmlsl2
  real(kind=r_double):: hdrmls13
  real(kind=r_kind),allocatable,dimension(:):: mlspres,mlsoz,mlsozpc,usage1
  integer(kind=i_kind),allocatable,dimension(:):: ipos
  integer(kind=i_kind):: iprofs,nprofs
  logical:: outside

  maxobs=size(ozout,2)
  rsat=999._r_kind

  open(lunin,file=dfile,form='unformatted')

  call openbf(lunin,'IN',lunin)
  call datelen(10)
  call readmg(lunin,subset,idate,iret)

  ! If it has failed at the first read(), ...
  if (iret/=0 .or. subset /= 'GM008015') then
     call closbf(lunin)
     close(lunin)

     call warn(myname_,'Failed at reading BUFR file, dfile =',trim(dfile))
     call warn(myname_,'                             dtype =',trim(dtype))
     call warn(myname_,'                            jsatid =',trim(jsatid))
     call warn(myname_,'                             lunin =',lunin)
     call warn(myname_,'                              iret =',iret)
     call warn(myname_,'                            subset =',trim(subset))

     nprofs = 0
     nodata = 0

     return
  endif

  call tell(myname_,'MLS o3lev data type, subset =',trim(subset))
  write(6,*)'read_ozone:  GMAO MLS o3lev data type, subset=',subset

  !    Q: o3lev data has 37 levels. Then, why is size(ipos) of 44?
  !    A: Because there are 44 entries in ozinfo.txt at the time

  mpos=max(nloz,jpch_oz)
  allocate (ipos(mpos))      ! 44? 37?
  allocate (usage1(nloz))

  ipos=999
  nmrecs=0
  nprofs=0
  nodata=0

! Set dependent variables and allocate arrays

  allocate (mlspres(nloz))
  allocate (mlsoz(nloz))
  allocate (mlsozpc(nloz))

  ikx=0
  do k=1,jpch_oz
     if(index(nusis_oz(k),'o3lev')/=0) then  ! all "o3lev" in ozinfo.txt
        ikx=ikx+1
        ipos(ikx)=k
     end if
  end do

  iy=0
  im=0
  idd=0
  ihh=0

! This is the top of the profile loop
  ndata=0
  iprofs=0
  obsloop: do 
     call readsb(lunin,iret)
     if (iret/=0) then           !JJJ, end of the subset
        call readmg(lunin,subset,jdate,iret)  !JJJ  open a new mg
        if (iret/=0) exit obsloop    !JJJ, no more  mg,  EOF      
        cycle obsloop
     endif

     do k=1,nloz
        if (iuse_oz(ipos(k)) < 0) then
           usage1(k) = 100._r_kind
        else
           usage1(k) = zero
        endif
     end do

!    extract header information
     call ufbint(lunin,hdrmls,13,1,iret,mlstr)
     rsat = hdrmls(1); ksatid=rsat

     if(jsatid == 'aura')kidsat = 785
     if (ksatid /= kidsat) cycle obsloop

     nmrecs=nmrecs+nloz

!    Convert observation location to radians
     slats0= hdrmls(2)
     slons0= hdrmls(3)
     if(abs(slats0)>90._r_kind .or. abs(slons0)>r360) cycle obsloop
     if(slons0< zero) slons0=slons0+r360
     if(slons0==r360) slons0=zero
     dlat_earth_deg = slats0
     dlon_earth_deg = slons0
     dlat_earth = slats0 * deg2rad
     dlon_earth = slons0 * deg2rad

     if(regional)then
        call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)
        if(outside) cycle obsloop
     else
        dlat = dlat_earth
        dlon = dlon_earth
        call grdcrd1(dlat,rlats,nlat,1)
        call grdcrd1(dlon,rlons,nlon,1)
     endif

! convert observation time to relative time
     idate5(1) = hdrmls(4)  !year
     idate5(2) = hdrmls(5)  !month
     idate5(3) = hdrmls(6)  !day
     idate5(4) = hdrmls(7)  !hour
     idate5(5) = hdrmls(8)  !minute
     call w3fs21(idate5,nmind)

     t4dv=real((nmind-iwinbgn),r_kind)*r60inv
     if (l4dvar.or.l4densvar) then
        if (t4dv<zero .or. t4dv>winlen) then
           write(6,*)'read_ozone: mls obs time idate5=',idate5,', t4dv=',&
              t4dv,' is outside time window, sstime=',sstime*r60inv
           cycle obsloop
        endif
     else
        sstime=real(nmind,r_kind)
        tdiff=(sstime-gstime)*r60inv
        if (abs(tdiff) > twind) then
           write(6,*)'read_ozone: mls obs time idate5=',idate5,', tdiff=',&
              tdiff,' is outside time window=',twind
           cycle obsloop
        endif
     end if

!    v2.2 data screening, only accept:
!    Pressure range:       215-0.02mb
!    Precision:            positive ozmp;    
!    Status flag:          only use even number
!    Quality(pccf):        use >1.2 for data at 215-100mb & low latitude, 
!                          use >0.4 for data elsewhere
!    Convergence:          use <1.8

!    Bit 1 in mlst represents data should not be used
!    Note: in bufr bits are defined from left to right as: 123456789...
!    whereas in hdf5 (and the nasa document) bits are defined from right to left as: ...876543210

     decimal=int(hdrmls(12))
     call dec2bin(decimal,binary_mls,18)
     if (binary_mls(1) == 1 ) cycle obsloop

     if(hdrmls(11) >= 1.8_r_kind) cycle obsloop

!    extract pressure, ozone mixing ratio and precision
     call ufbrep(lunin,hdrmlsl2,4,nloz,iret,mlstrl2)

     iprofs=iprofs+1    ! counting the profiles
     do k=1,nloz
        mlspres(k)=log(hdrmlsl2(1,k)*0.001_r_kind)    ! mls pressure in Pa, coverted to log(cb)
        mlsoz(k)=hdrmlsl2(2,k)                     ! ozone mixing ratio in ppmv
        mlsozpc(k)=hdrmlsl2(3,k)                   ! ozone mixing ratio precision in ppmv
        if (dsis /= 'mls_aura_ozpc') mlsozpc(k)=hdrmlsl2(4,k)   ! use obserr if (dsis /= 'mls_aura_ozpc')
     end do

     do k=1,nloz
        if(hdrmlsl2(1,k)>21600._r_kind .or. hdrmlsl2(1,k)<2._r_kind) usage1(k)=1000._r_kind
        if(hdrmlsl2(3,k)<=0._r_kind) usage1(k)=1000._r_kind
     end do

     hdrmls13=hdrmls(13)*0.1_r_kind
     if (abs(slats0)<30._r_kind) then
        do k=1,nloz
           if(hdrmlsl2(1,k)>10000._r_kind .and. hdrmlsl2(1,k)<21600._r_kind) then
              if(hdrmls13 <= 1.2_r_kind) usage1(k)=1000._r_kind
           else
              if(hdrmls13 <= 0.4_r_kind) usage1(k)=1000._r_kind
           endif
        end do
     else
        if(hdrmls13 <= 0.4_r_kind) then
           do k=1,nloz
              usage1(k)=1000._r_kind
           end do
        end if
     end if

     do k=1,nloz

        ndata=ndata+1

        if(ndata <= maxobs) then
           ozout( 1,ndata)=rsat
           ozout( 2,ndata)=t4dv
           ozout( 3,ndata)=dlon               ! grid relative longitude
           ozout( 4,ndata)=dlat               ! grid relative latitude
           ozout( 5,ndata)=dlon_earth_deg     ! earth relative longitude (degrees)
           ozout( 6,ndata)=dlat_earth_deg     ! earth relative latitude (degrees)
           ozout( 7,ndata)=hdrmls(10)         ! solar zenith angle

           ozout( 8,ndata)=usage1(k)          ! 
           ozout( 9,ndata)=mlspres(k)          ! mls pressure in log(cb)
           ozout(10,ndata)=mlsozpc(k)   ! ozone mixing ratio precision in ppmv
           ozout(11,ndata)=real(ipos(k))       ! pointer of obs level index in ozinfo.txt
           ozout(12,ndata)=nloz         ! # of mls vertical levels
           ozout(13,ndata)=mlsoz(k)     ! ozone mixing ratio in ppmv
        endif
     end do

  end do obsloop
! End of o3lev bufr loop

  nodata=min(ndata,maxobs)   ! count of retained data
  if(nodata < ndata) then    ! warning if all are not properly retained
     call warn(myname_,'insufficient buffer space, expecting =',ndata)
     call warn(myname_,'                     actual retained =',nodata)
     call warn(myname_,'                       size(ozout,2) =',maxobs)
  endif

!  write(stdout,'(3a,3i8,f8.2)') mprefix('read_ozone'), &
!     ' obstype,nmrecs,ndata,nodata,no/ndata = ',dtype,nmrecs,ndata,nodata,real(nodata)/ndata

end subroutine ozlev_bufrread_

!..................................................................................
subroutine check(status,stat)
!..................................................................................
  use netcdf, only: nf90_noerr
  use netcdf, only: nf90_strerror
  use mpeu_util, only: die,perr,warn
  implicit none
  integer(i_kind), intent (in) :: status
  integer(i_kind), optional, intent(out):: stat
  character(len=*),parameter:: myname_=myname//'::check'

  if(present(stat)) stat=status

  if(status /= nf90_noerr) then 
     if(present(stat)) then
        call warn(myname_,'ignored, nf90_strerror =',trim(nf90_strerror(status)))
     else
        call perr(myname_,'nf90_strerror =',trim(nf90_strerror(status)))
        call die(myname_)
     endif
  endif
end subroutine check
end module m_extozone
