module mpeu_mpif
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    module mpeu_mpif
!   prgmmr:      j guo <jguo@nasa.gov>
!      org:      NASA/GSFC, Global Modeling and Assimilation Office, 900.3
!     date:      2010-03-22
!
! abstract: a portable interface to include "mpif.h" for mpi.
!
! program history log:
!   2010-03-22  j guo   - added this document block
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

   implicit none
   private          ! except

   public :: mpi_type
   interface mpi_type; module procedure &
      type_r0_of_integer1, &
      type_r0_of_integer2, &
      type_r0_of_integer4, &
      type_r0_of_integer8, &
      type_r0_of_real4   , &
      type_r0_of_real8   , &
      type_r0_of_real16  , &
      type_r1_of_integer1, &
      type_r1_of_integer2, &
      type_r1_of_integer4, &
      type_r1_of_integer8, &
      type_r1_of_real4   , &
      type_r1_of_real8   , &
      type_r1_of_real16  , &
      type_r2_of_integer1, &
      type_r2_of_integer2, &
      type_r2_of_integer4, &
      type_r2_of_integer8, &
      type_r2_of_real4   , &
      type_r2_of_real8   , &
      type_r2_of_real16  , &
      type_r3_of_integer1, &
      type_r3_of_integer2, &
      type_r3_of_integer4, &
      type_r3_of_integer8, &
      type_r3_of_real4   , &
      type_r3_of_real8   , &
      type_r3_of_real16
   end interface

   public :: mpi_ikind
   public :: mpi_address_kind

   public :: mpi_integer1
   public :: mpi_integer2
   public :: mpi_integer4
   public :: mpi_integer8

   public :: mpi_integer
   public :: mpi_real
   public :: mpi_double_precision
   public :: mpi_logical
   public :: mpi_character

   public :: mpi_2integer
   public :: mpi_2real
   public :: mpi_2double_precision

   public :: mpi_real4
   public :: mpi_real8
   public :: mpi_real16

   public :: mpi_comm_world
   public :: mpi_comm_null
   public :: mpi_request_null

   public :: mpi_sum
   public :: mpi_prod
   public :: mpi_min
   public :: mpi_max
   public :: mpi_minloc
   public :: mpi_maxloc

   public :: mpi_max_error_string
   public :: mpi_status_size
   public :: mpi_error

!#if !defined(sysLinux)
   public :: mpi_offset_kind
   public :: mpi_info_null
   public :: mpi_mode_rdonly
   public :: mpi_mode_rdwr
   public :: mpi_mode_wronly
   public :: mpi_mode_create
   public :: mpi_seek_set
!#endif
   public :: mpi_byte

#ifdef MPICH_
   public :: mpipriv     ! the common block name
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mpeu_mpif - a portable interface to the mpi "mpif.h" COMMONs.
!
! !DESCRIPTION:
!
!   The purpose of \verb"mpeu_mpif" module is to provide a portable
!   interface of \verb"mpif.h" with different mpi implementation.
!   By combining module \verb"mpeu_mpif" and \verb"mpeu_mpif90", it may be
!   possible to build a fortran 90 mpi binding module graduately.
!
!   Although it is possible to use \verb'include "mpif.h"' directly
!   in individual modules, it has several problems:
!   \begin{itemize}
!   \item It may conflict with either the source code of a {\sl fixed}
!      format or the code of a {\sl free} format;
!   \item It does not provide the protection and the safety of using
!     these variables as what a \verb"module" would provide.
!   \end{itemize}
!
!   More information may be found in the module \verb"mpeu_mpif90".
!
! !INTERFACE:

include "mpif.h"

! !REVISION HISTORY:
!    01Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!    16Feb05 - Todling - added a number of vars needed by gsi
!    18Mar05 - Todling - added a few more mpi vars used by gsi
!    20Mar05 - Todling - see no reason for sysLinux ifdef (commented out)
!    18Feb09 - Jing Guo <Jing.Guo@nasa.gov>
!     - Copied from gmao_mpeu/m_mpif.F here to avoid the
!       dependency of gsi_gridcomp/ on gmao_mpeu.
!     - Renamed to mpeu_mpif to avoid possible name conflict
!       with gmao_mpeu/.
!    23Feb09 - Jing Guo <Jing.Guo@nasa.gov>
!     - Renamed to *.F90, by assuming the system mpif.h is
!       compatible with free-format Fortran.
!    22Feb17 - Todling - add dummy vars to satisfy emc requirement
!                        of all vars being used - this is counter
!                        advanced fortran interface practices.
!EOP
!_______________________________________________________________________

    integer, parameter:: mpi_ikind=kind(mpi_comm_world)

    character(len=*),parameter :: myname='mpeu_mpif'

contains
function type_r0_of_integer1(mold) result(mtype)
  use kinds, only: i_byte
  implicit none
  integer(kind=mpi_ikind):: mtype
  integer(i_byte),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_integer1
end function type_r0_of_integer1
function type_r0_of_integer2(mold) result(mtype)
  use kinds, only: i_short
  implicit none
  integer(kind=mpi_ikind):: mtype
  integer(i_short),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_integer2
end function type_r0_of_integer2
function type_r0_of_integer4(mold) result(mtype)
  use kinds, only: i_long
  implicit none
  integer(kind=mpi_ikind):: mtype
  integer(i_long),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_integer4
end function type_r0_of_integer4
function type_r0_of_integer8(mold) result(mtype)
  use kinds, only: i_llong
  implicit none
  integer(kind=mpi_ikind):: mtype
  integer(i_llong),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_integer8
end function type_r0_of_integer8

function type_r0_of_real4(mold) result(mtype)
  use kinds, only: r_single
  implicit none
  integer(kind=mpi_ikind):: mtype
  real(r_single),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_real4
end function type_r0_of_real4
function type_r0_of_real8(mold) result(mtype)
  use kinds, only: r_double
  implicit none
  integer(kind=mpi_ikind):: mtype
  real(r_double),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_real8
end function type_r0_of_real8
function type_r0_of_real16(mold) result(mtype)
  use kinds, only: r_quad
  implicit none
  integer(kind=mpi_ikind):: mtype
  real(r_quad),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_real16
end function type_r0_of_real16

function type_r1_of_integer1(mold) result(mtype)
  use kinds, only: i_byte
  implicit none
  integer(kind=mpi_ikind):: mtype
  integer(i_byte),dimension(:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_integer1
end function type_r1_of_integer1
function type_r1_of_integer2(mold) result(mtype)
  use kinds, only: i_short
  implicit none
  integer(kind=mpi_ikind):: mtype
  integer(i_short),dimension(:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_integer2
end function type_r1_of_integer2
function type_r1_of_integer4(mold) result(mtype)
  use kinds, only: i_long
  implicit none
  integer(kind=mpi_ikind):: mtype
  integer(i_long),dimension(:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_integer4
end function type_r1_of_integer4
function type_r1_of_integer8(mold) result(mtype)
  use kinds, only: i_llong
  implicit none
  integer(kind=mpi_ikind):: mtype
  integer(i_llong),dimension(:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_integer8
end function type_r1_of_integer8

function type_r1_of_real4(mold) result(mtype)
  use kinds, only: r_single
  implicit none
  integer(kind=mpi_ikind):: mtype
  real(r_single),dimension(:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_real4
end function type_r1_of_real4
function type_r1_of_real8(mold) result(mtype)
  use kinds, only: r_double
  implicit none
  integer(kind=mpi_ikind):: mtype
  real(r_double),dimension(:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_real8
end function type_r1_of_real8
function type_r1_of_real16(mold) result(mtype)
  use kinds, only: r_quad
  implicit none
  integer(kind=mpi_ikind):: mtype
  real(r_quad),dimension(:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_real16
end function type_r1_of_real16

function type_r2_of_integer1(mold) result(mtype)
  use kinds, only: i_byte
  implicit none
  integer(kind=mpi_ikind):: mtype
  integer(i_byte),dimension(:,:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_integer1
end function type_r2_of_integer1
function type_r2_of_integer2(mold) result(mtype)
  use kinds, only: i_short
  implicit none
  integer(kind=mpi_ikind):: mtype
  integer(i_short),dimension(:,:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_integer2
end function type_r2_of_integer2
function type_r2_of_integer4(mold) result(mtype)
  use kinds, only: i_long
  implicit none
  integer(kind=mpi_ikind):: mtype
  integer(i_long),dimension(:,:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_integer4
end function type_r2_of_integer4
function type_r2_of_integer8(mold) result(mtype)
  use kinds, only: i_llong
  implicit none
  integer(kind=mpi_ikind):: mtype
  integer(i_llong),dimension(:,:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_integer8
end function type_r2_of_integer8

function type_r2_of_real4(mold) result(mtype)
  use kinds, only: r_single
  implicit none
  integer(kind=mpi_ikind):: mtype
  real(r_single),dimension(:,:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_real4
end function type_r2_of_real4
function type_r2_of_real8(mold) result(mtype)
  use kinds, only: r_double
  implicit none
  integer(kind=mpi_ikind):: mtype
  real(r_double),dimension(:,:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_real8
end function type_r2_of_real8
function type_r2_of_real16(mold) result(mtype)
  use kinds, only: r_quad
  implicit none
  integer(kind=mpi_ikind):: mtype
  real(r_quad),dimension(:,:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_real16
end function type_r2_of_real16

function type_r3_of_integer1(mold) result(mtype)
  use kinds, only: i_byte
  implicit none
  integer(kind=mpi_ikind):: mtype
  integer(i_byte),dimension(:,:,:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_integer1
end function type_r3_of_integer1
function type_r3_of_integer2(mold) result(mtype)
  use kinds, only: i_short
  implicit none
  integer(kind=mpi_ikind):: mtype
  integer(i_short),dimension(:,:,:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_integer2
end function type_r3_of_integer2
function type_r3_of_integer4(mold) result(mtype)
  use kinds, only: i_long
  implicit none
  integer(kind=mpi_ikind):: mtype
  integer(i_long),dimension(:,:,:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_integer4
end function type_r3_of_integer4
function type_r3_of_integer8(mold) result(mtype)
  use kinds, only: i_llong
  implicit none
  integer(kind=mpi_ikind):: mtype
  integer(i_llong),dimension(:,:,:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_integer8
end function type_r3_of_integer8

function type_r3_of_real4(mold) result(mtype)
  use kinds, only: r_single
  implicit none
  integer(kind=mpi_ikind):: mtype
  real(r_single),dimension(:,:,:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_real4
end function type_r3_of_real4
function type_r3_of_real8(mold) result(mtype)
  use kinds, only: r_double
  implicit none
  integer(kind=mpi_ikind):: mtype
  real(r_double),dimension(:,:,:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_real8
end function type_r3_of_real8
function type_r3_of_real16(mold) result(mtype)
  use kinds, only: r_quad
  implicit none
  integer(kind=mpi_ikind):: mtype
  real(r_quad),dimension(:,:,:),intent(in):: mold
  integer(mpi_ikind):: dummymold
  dummymold=kind(mold)
  mtype=mpi_real16
end function type_r3_of_real16

end module mpeu_mpif
!.
