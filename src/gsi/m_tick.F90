module m_tick
!$$$ module documentation block
!           .      .    .                                       .
! module:   m_tick
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-01-08 todling - encapsulated in this module
!   2009-08-06 lueken  - added module doc block
!
! subroutines included:
!   sub tick
!   sub incymd
!   sub leap_year
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block
  use kinds, only: i_kind

  implicit none

! set default to private
  private
! set subroutines to public
  public :: tick

contains

#ifdef ibm_sp

subroutine tick (nymd, nhms, ndt)
!$$$ subprogram documentation block
!               .      .    .                                       .
! subprogram:   tick
!   prgmmr:
!
! abstract: Dummy routine for ibm_sp.
!
! program history log:
!   2009-08-06 lueken  - added subprogram doc block
!
!   input argument list:
!    ndt
!    nymd
!    nhms
!
!   output argument list:
!    nymd
!    nhms
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block
  implicit none

  integer(i_kind),intent(in   ) :: ndt
  integer(i_kind),intent(inout) :: nymd
  integer(i_kind),intent(inout) :: nhms

  end subroutine tick

#else

subroutine tick (nymd, nhms, ndt)
!$$$ subprogram documentation block
!               .      .    .                                       .
! subprogram:   tick
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-08-06 lueken  - added subprogram doc block
!
!   input argument list:
!    ndt
!    nymd
!    nhms
!
!   output argument list:
!    nymd
!    nhms
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block
  implicit none

! Input:
  integer(i_kind),intent(in   ) :: ndt                  ! Time-step
! Input/Output:
  integer(i_kind),intent(inout) :: nymd                 ! Current yyyymmdd
  integer(i_kind),intent(inout) :: nhms                 ! Current hhmmss
! Local:
  integer(i_kind) :: nsecf, nhmsf, nsec, n

! Origin:     L.L. Takacs
! Revision:   S.-J. Lin Mar 2000

  nsecf(n)   = n/10000*3600 + mod(n,10000)/100* 60 + mod(n,100)
  nhmsf(n)   = n/3600*10000 + mod(n,3600 )/ 60*100 + mod(n, 60)

  nsec = nsecf(nhms) + ndt

  if (nsec>86400)  then
     do while (nsec>86400)
        nsec = nsec - 86400
        nymd = invymd (nymd,1)
     enddo
  endif

  if (nsec==86400)  then
     nsec = 0
     nymd = invymd (nymd,1)
  endif

  if (nsec < 0)  then
     do while (nsec < 0)
        nsec = 86400 + nsec
        nymd = invymd (nymd,-1)
     enddo
  endif

  nhms = nhmsf (nsec)
  return
end subroutine tick

integer(i_kind) function invymd (nymd,m)
!$$$ subprogram documentation block
!               .      .    .                                       .
! subprogram:   invymd
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-08-06 lueken  - added subprogram doc block
!
!   input argument list:
!    nymd
!    M
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block
  implicit none

!  Purpose
!     invymd:  nymd changed by one day
!     modymd:  nymd converted to julian date
!  Description of parameters
!     nymd     current date in yymmdd format
!     m        +/- 1 (day adjustment)

  integer(i_kind) ndpm(12)
  data    ndpm /31, 28, 31, 30, 31, 30, &
                31, 31, 30, 31, 30, 31/
  integer(i_kind),intent(in):: nymd, m
  integer(i_kind) ny, nm, nd

  ny = nymd / 10000
  nm = mod(nymd,10000) / 100
  nd = mod(nymd,100) + m

  if (nd==0) then
     nm = nm - 1
     if (nm==0) then
        nm = 12
        ny = ny - 1
     endif
     nd = ndpm(nm)
     if (nm==2 .and. leap_year(ny))  nd = 29
  endif

  if (.not. (nd==29 .and. nm==2 .and. leap_year(ny)))  then

     if (nd>ndpm(nm)) then
        nd = 1
        nm = nm + 1
        if (nm>12) then
           nm = 1
           ny = ny + 1
        endif
     endif

  end if
  invymd = ny*10000 + nm*100 + nd
  return
end function invymd

logical function leap_year(ny)
!$$$ subprogram documentation block
!               .      .    .                                       .
! subprogram:   leap_year
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-08-06 lueken  - added subprogram doc block
!
!   input argument list:
!    ny
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block
!
! Determine if year ny is a leap year
!
! Author: S.-J. Lin
  implicit none

  integer(i_kind),intent(in   ) :: ny

  integer(i_kind) ny00

!
! No leap years prior to 1900
!
  parameter ( ny00 = 1900 )   ! The threshold for starting leap-year 

  if( mod(ny,4) == 0 .and. ny >= ny00 ) then
     leap_year = .true.
  else
     leap_year = .false.
  endif

  return 
end function leap_year
#endif
end module m_tick
