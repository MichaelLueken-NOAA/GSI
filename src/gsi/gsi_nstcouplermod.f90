!----------------------------------------------------------------------------
!BOP
!  
! !MODULE:  gsi_nstcouplermod ---
!
! !INTERFACE:
!
! !DESCRIPTION: This module provides general interface for
!               nst- sea surface skin temp analysis
!
! !REVISION HISTORY:
!
!  2011-10-20 RT/ Akella- Initial code
!  2012-03-05 SA-         _full fields: tref, dt_cool, dt_warm, z_c, z_w, ... are declared here instead of satthin     
!  2015-05-01 Li-         Change the nst fields to be single precision
!  2017-09-14 LI-         Change the default value to be 1 for fac_dtl & fac_tsl
!
!EOP
!-------------------------------------------------------------------------
!  def tref_full      - sea surface reference temperature-- foundation sst
!  def dt_cool_full   - sea cooling amount across sub-layer (or, cool-layer)
!  def z_c_full       - sub-layer thickness
!  def dt_warm_full   - sea diurnal warming amount
!  def z_w_full       - diurnal warming layer thickness
!  *********************************************************
!   Following 4 fields are for gfs
!  *********************************************************
!   def c_0_full       - coefficient 1 to calculate d(tz)/d(tr)
!   def c_d_full       - coefficient 2 to calculate d(tz)/d(tr)
!   def w_0_full       - coefficient 3 to calculate d(tz)/d(tr)
!   def w_d_full       - coefficient 4 to calculate d(tz)/d(tr)
!  *********************************************************
!-------------------------------------------------------------------------

module gsi_nstcouplermod

! !USES:
use kinds,         only: r_single, r_kind, i_kind

implicit none
private

!
! !PUBLIC MEMBER FUNCTIONS:
!
public gsi_nstcoupler_init_nml
public gsi_nstcoupler_init
public gsi_nstcoupler_read
public gsi_nstcoupler_skindepth
public gsi_nstcoupler_deter
public gsi_nstcoupler_final

public :: nst_gsi,nstinfo,zsea1,zsea2,fac_dtl,fac_tsl
public :: tref_full,dt_cool_full,z_c_full,dt_warm_full,z_w_full
public :: c_0_full,c_d_full,w_0_full,w_d_full

integer(i_kind) :: nst_gsi   ! indicator of tr analysis
integer(i_kind) :: nstinfo   ! number of nst variables
integer(i_kind) :: zsea1     ! upper depth (in mm) to do the mean
integer(i_kind) :: zsea2     ! lower depth (in mm) to do the mean
integer(i_kind) :: fac_dtl   ! indicator of dtl
integer(i_kind) :: fac_tsl   ! indicator of tsl

real(r_single),allocatable,dimension(:,:,:):: tref_full,dt_cool_full,z_c_full,dt_warm_full,z_w_full
real(r_single),allocatable,dimension(:,:,:):: c_0_full,c_d_full,w_0_full,w_d_full

!-------------------
interface gsi_nstcoupler_init
  subroutine nst_init_()
     implicit none
  end subroutine nst_init_
end interface
!-------------------

interface gsi_nstcoupler_read
  subroutine nst_read_(mype_io)
     use kinds,         only: i_kind
     implicit none

     integer(i_kind), intent(in   ) :: mype_io
     
  end subroutine nst_read_
end interface
!-------------------

interface gsi_nstcoupler_final
  subroutine nst_final_()
     implicit none
  end subroutine nst_final_
end interface
!-------------------

interface gsi_nstcoupler_skindepth
  subroutine skindepth_(obstype, zob)
     use kinds,   only: r_kind
     implicit none

     character(10), intent(in)  :: obstype
     real(r_kind),  intent(out) :: zob
  end subroutine skindepth_
end interface
!-------------------

interface gsi_nstcoupler_deter
  subroutine deter_nst_(dlat_earth,dlon_earth,obstime,zob,tref,dtw,dtc,tz_tr)
     use kinds,   only: r_kind
     implicit none

     real(r_kind), intent(in ) :: dlat_earth,dlon_earth,obstime,zob
     real(r_kind), intent(out) :: tref,dtw,dtc,tz_tr
  end subroutine deter_nst_
end interface
!-------------------

contains

!-------------------
subroutine gsi_nstcoupler_init_nml

  use mpimod, only: mype

  implicit none
  
  if ( mype == 0 ) &
    write(6,*)'NST_INIT_NML_: Initializing default NST namelist variables'

  nst_gsi   = 0          ! 0 = no nst info at all in gsi
                         ! 1 = read nst info but not applied
                         ! 2 = read nst info, applied to tb simulation but no tr analysis
                         ! 3 = read nst info, applied to tb simulation and do tr Analysis
  nstinfo   = 0          ! number of nst fields used in tr analysis
  zsea1     = 0          ! upper depth to do the mean
  zsea2     = 0          ! lower depth to do the mean
  fac_dtl   = 1          ! indicator to apply dtl model
  fac_tsl   = 1          ! indicator to apply tsl model

  return

end subroutine gsi_nstcoupler_init_nml
!-------------------

end module gsi_nstcouplermod

