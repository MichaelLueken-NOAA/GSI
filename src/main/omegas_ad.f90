subroutine omegas_ad( im, ix, km, dphi_i, dlam_i, u_i, v_i, div_i, &
     ps_i, rcl, del, sl, vvel_o, u_i_ad, v_i_ad, div_i_ad, vvel_o_ad, &
     adjoint)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    omegas_ad    forward & adjoint model for GFS vertical velocity calculation
!     prgmmr:    treadon     org: np23                date: 2003-12-18
!
! abstract:  This subroutine contains the forward and ajoint models for the
!            GFS vertical velocity calculation
!
! program history log:
!   2003-12-18  treadon - initial routine
!   2004-06-15  treadon - reformat documenation
!   2006-04-12  treadon - change del and sl from 1d to 2d arrays
!   2006-08-10  treadon - change dphi_i,dlam_i from ln(ps) to ps
!
!   input argument list:
!     im      - actual number of profiles to be processed
!     ix      - maximum number of profiles to process (array dimension)
!     km      - number of levels in vertical profile
!     dphi_i  - partial derivative of ps with respect to latitude
!     dlam_i  - partial derivative of ps with respect to longitude
!     u_i     - zonal wind
!     v_i     - meridional wind
!     div_i   - horizontal divergence
!     ps_i    - surface pressure
!     rcl     - 1/sin(lat)**2
!     del     - "sigma" thickness of layers
!     sl      - "sigma" value at layer midpoints
!     vvel_o_ad- vertical velocity perturbation
!     adjoint - logical flag (.false.=forward model only, .true.=forward and ajoint)
!
!   output argument list:
!     vvel_o  - vertical velocity
!     u_i_ad   - partial derivative of vertical velocity with respect to zonal wind
!     v_i_ad   - partial derivative of vertical velocity with respect to meridional wind
!     div_i_ad - partial derivative of vertical velocity with respect to horizontal divergence

!
! remarks:
!    The adjoint portion of this routine was generated by the
!    Tangent linear and Adjoint Model Compiler,  TAMC 5.3.0
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
!==============================================
! all entries are defined explicitly
!==============================================
  use kinds, only: r_kind,i_kind
  use constants, only: zero, half
  implicit none

!==============================================
! define arguments
!==============================================
  logical adjoint
  integer(i_kind) ix
  integer(i_kind) km
  real(r_kind) div_i_ad(km,ix)
  real(r_kind) u_i_ad(km,ix)
  real(r_kind) v_i_ad(km,ix)
  real(r_kind) vvel_o_ad(km,ix)
  real(r_kind) del(km,ix)
  real(r_kind) div_i(km,ix)
  real(r_kind) dlam_i(ix)
  real(r_kind) dphi_i(ix)
  integer(i_kind) im
  real(r_kind) ps_i(ix)
  real(r_kind) rcl(ix)
  real(r_kind) sl(km,ix)
  real(r_kind) u_i(km,ix)
  real(r_kind) v_i(km,ix)
  real(r_kind) vvel_o(km,ix)

!==============================================
! define local variables
!==============================================
  real(r_kind) cb_ad(km,ix)
  real(r_kind) cg_ad(km,ix)
  real(r_kind) db_ad(km,ix)
  real(r_kind) dot_ad(km+1,ix)
  real(r_kind) cb(km,ix)
  real(r_kind) cg(km,ix)
  real(r_kind) db(km,ix)
  real(r_kind) dlam(ix)
  real(r_kind) dot(km+1,ix)
  real(r_kind) dphi(ix)
  integer(i_kind) i
  integer(i_kind) ip1
  integer(i_kind) ip2
  integer(i_kind) k

!----------------------------------------------
! RESET LOCAL ADJOINT VARIABLES
!----------------------------------------------
  do ip2 = 1, ix
     do ip1 = 1, km
        cb_ad(ip1,ip2) = zero
     end do
  end do
  do ip2 = 1, ix
     do ip1 = 1, km
        cg_ad(ip1,ip2) = zero
     end do
  end do
  do ip2 = 1, ix
     do ip1 = 1, km
        db_ad(ip1,ip2) = zero
     end do
  end do
  do ip2 = 1, ix
     do ip1 = 1, km+1
        dot_ad(ip1,ip2) = zero
     end do
  end do

!----------------------------------------------
! ROUTINE BODY
!----------------------------------------------
!----------------------------------------------
! FUNCTION AND TAPE COMPUTATIONS
!----------------------------------------------
  do i = 1, im
     do k = 1, km+1
        dot(k,i) = zero
     end do
  end do
  do i = 1, im
     dphi(i) = dphi_i(i)*rcl(i) / ps_i(i)
     dlam(i) = dlam_i(i)*rcl(i) / ps_i(i)
  end do
  do i = 1, im
     do k = 1, km
        cg(k,i) = u_i(k,i)*dlam(i)+v_i(k,i)*dphi(i)
     end do
  end do
  do i = 1, im
     db(1,i) = del(1,i)*div_i(1,i)
     cb(1,i) = del(1,i)*cg(1,i)
  end do
  do i = 1, im
     do k = 1, km-1
        db(k+1,i) = db(k,i)+del(k+1,i)*div_i(k+1,i)
        cb(k+1,i) = cb(k,i)+del(k+1,i)*cg(k+1,i)
     end do
  end do
  do i = 1, im
     do k = 1, km-1
        dot(k+1,i) = dot(k,i)+del(k,i)*(db(km,i)+cb(km,i)-div_i(k,i)- &
             cg(k,i))
     end do
  end do
  do i = 1, im
     do k = 1, km
        vvel_o(k,i) = sl(k,i)*(cg(k,i)-cb(km,i)-db(km,i))-half*(dot(k+1, &
             i)+dot(k,i))
        vvel_o(k,i) = vvel_o(k,i)*ps_i(i)
     end do
  end do
  
  if (.not.adjoint) return

!----------------------------------------------
! ADJOINT COMPUTATIONS
!----------------------------------------------
  do i = 1, im
     do k = 1, km
        vvel_o_ad(k,i) = vvel_o_ad(k,i)*ps_i(i)
        cb_ad(km,i) = cb_ad(km,i)-vvel_o_ad(k,i)*sl(k,i)
        cg_ad(k,i) = cg_ad(k,i)+vvel_o_ad(k,i)*sl(k,i)
        db_ad(km,i) = db_ad(km,i)-vvel_o_ad(k,i)*sl(k,i)
        dot_ad(k+1,i) = dot_ad(k+1,i)-half*vvel_o_ad(k,i)
        dot_ad(k,i) = dot_ad(k,i)-half*vvel_o_ad(k,i)
        vvel_o_ad(k,i) = zero
     end do
  end do
  do i = 1, im
     do k = km-1, 1, -1
        cb_ad(km,i) = cb_ad(km,i)+dot_ad(k+1,i)*del(k,i)
        cg_ad(k,i) = cg_ad(k,i)-dot_ad(k+1,i)*del(k,i)
        db_ad(km,i) = db_ad(km,i)+dot_ad(k+1,i)*del(k,i)
        div_i_ad(k,i) = div_i_ad(k,i)-dot_ad(k+1,i)*del(k,i)
        dot_ad(k,i) = dot_ad(k,i)+dot_ad(k+1,i)
        dot_ad(k+1,i) = zero
     end do
  end do
  do i = 1, im
     do k = km-1, 1, -1
        cg_ad(k+1,i) = cg_ad(k+1,i)+cb_ad(k+1,i)*del(k+1,i)
        cb_ad(k,i) = cb_ad(k,i)+cb_ad(k+1,i)
        cb_ad(k+1,i) = zero
        div_i_ad(k+1,i) = div_i_ad(k+1,i)+db_ad(k+1,i)*del(k+1,i)
        db_ad(k,i) = db_ad(k,i)+db_ad(k+1,i)
        db_ad(k+1,i) = zero
     end do
  end do
  do i = 1, im
     cg_ad(1,i) = cg_ad(1,i)+cb_ad(1,i)*del(1,i)
     cb_ad(1,i) = zero
     div_i_ad(1,i) = div_i_ad(1,i)+db_ad(1,i)*del(1,i)
     db_ad(1,i) = zero
  end do
  do i = 1, im
     do k = 1, km
        u_i_ad(k,i) = u_i_ad(k,i)+cg_ad(k,i)*dlam(i)
        v_i_ad(k,i) = v_i_ad(k,i)+cg_ad(k,i)*dphi(i)
        cg_ad(k,i) = zero
     end do
  end do
  
  return
end subroutine omegas_ad


