module lagmod
!$$$   module documentation block
!                .      .    .                                       .
! module:  lagmod
! prgmmr:  purser/cucurull           org: np23            date: 2005-12-01
!
! abstract: This module contains routines (including tangent linear and adjoint code)
!           to perform basic operations related to Lagrange polynomial interpolation.
!
! program history log:
!   2003        purser   - routines for the forward code
!   2005-12-01  cucurull - implement tangent linear and adjoint code 
!                        - adapt the code to GSI standards
!   2008-05-09  safford  - complete documentation block
! 
! subroutines included:
!   setq
!   setq_tl
!   setq_ad
!   lagdw
!   lagdw_tl
!   lagdw_ad
!   slagdw
!   slagdw_tl
!   slagdw_ad
!
! variables defined:
!
! attributes:
!   language:  f90
!   machine:   ibm RS/6000 SP
!
!$$$ end documentation block

!============================================================================
use kinds, only: r_kind,i_kind
use constants, only: zero,one
implicit none

! set default to private
private
! set subroutines to public
public :: setq
public :: setq_tl
public :: setq_ad
public :: lagdw
public :: lagdw_tl
public :: lagdw_ad
public :: slagdw
public :: slagdw_tl
public :: slagdw_ad

contains


subroutine setq(q,x,n)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    setq
!
!   prgrmmr:
!
! abstract:      Precompute the n constant denominator factors of the n-point 
!                Lagrange polynomial interpolation formula.
!
! program history log:
!   2008-05-09  safford - add subprogram doc block, rm unused uses
!
!   input argument list:
!     x -    The n abscissae.
!     n -    The number of points involved.
!
!   output argument list:
!     q -    The n denominator constants.
!
! attributes:
!   language:  f90
!   machine:   ibm RS/6000 SP
!
!$$$ end documentation block

implicit none

integer(i_kind)          ,intent(in   ) :: n
real(r_kind),dimension(n),intent(  out) :: q
real(r_kind),dimension(n),intent(in   ) :: x
!-----------------------------------------------------------------------------
integer(i_kind)                      :: i,j
!=============================================================================
do i=1,n
   q(i)=one
   do j=1,n
      if(j /= i)q(i)=q(i)/(x(i)-x(j))
   enddo
enddo
end subroutine setq


!============================================================================
subroutine setq_tl(q,q_tl,x,x_tl,n)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    setq_tl
!
!   prgrmmr:
!
! abstract:     
!
! program history log:
!   2008-05-09  safford - add subprogram doc block, rm unused uses
!
!   input argument list:
!     x     - 
!     x_tl  - 
!     n     -    The number of points involved.
!
!   output argument list:
!     q    -   
!     q_tl
!
! attributes:
!   language:  f90
!   machine:   ibm RS/6000 SP
!
!$$$ end documentation block

!============================================================================
implicit none

integer(i_kind)          ,intent(in   ) :: n
real(r_kind),dimension(n),intent(  out) :: q,q_tl
real(r_kind),dimension(n),intent(in   ) :: x,x_tl
!-----------------------------------------------------------------------------
integer(i_kind)                      :: i,j
real(r_kind)                         :: rat
!=============================================================================
do i=1,n
   q(i)=one
   q_tl(i)=zero
   do j=1,n
      if(j /= i) then
         rat=one/(x(i)-x(j))
         q_tl(i)=(q_tl(i)-q(i)*(x_tl(i)-x_tl(j))*rat)*rat
         q(i)=q(i)*rat
      endif
   enddo
enddo
end subroutine setq_tl


!============================================================================
subroutine setq_ad(q_ad,x,x_ad,n)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    setq_ad
!
!   prgrmmr:
!
! abstract:     
!
! program history log:
!   2008-05-09  safford - add subprogram doc block, rm unused uses
!
!   input argument list:
!     q_ad - 
!     x    -
!     x_ad -
!     n    -    The number of points involved.
!
!   output argument list:
!     x_ad -
!     q_ad - 
!
! attributes:
!   language:  f90
!   machine:   ibm RS/6000 SP
!
!$$$ end documentation block

!============================================================================
implicit none

integer(i_kind)          ,intent(in   ) :: n
real(r_kind),dimension(n),intent(in   ) :: x
real(r_kind),dimension(n),intent(inout) :: x_ad
real(r_kind),dimension(n),intent(inout) :: q_ad
!-----------------------------------------------------------------------------
integer(i_kind)                      :: i,j
real(r_kind),dimension(n)            :: q
real(r_kind),dimension(n,n)          :: jac
!=============================================================================
call setq(q,x,n)
jac=zero
do i=1,n
   do j=1,n
      if(j /= i) then
         jac(j,i)=q(i)/(x(i)-x(j))
         jac(i,i)=jac(i,i)-jac(j,i)
      endif
   enddo
enddo
x_ad=x_ad+matmul(jac,q_ad)
end subroutine setq_ad


!============================================================================
subroutine lagdw(x,xt,q,w,dw,n)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    lagdw
!
!   prgrmmr:
!
! abstract:      Construct the lagrange weights and their derivatives when 
!                target abscissa is known and denominators q have already 
!                been precomputed
!
! program history log:
!   2008-05-09  safford - add subprogram doc block, rm unused uses
!
!   input argument list:
!     x   - Grid abscissae
!     xt  - Target abscissa
!     q   - Q factors (denominators of the lagrange weight formula)
!     n   - Number of grid points involved in the interpolation
!
!   output argument list:
!     w   - Lagrange weights
!     dw  - Derivatives, dw/dx, of Lagrange weights w
!
! attributes:
!   language:  f90
!   machine:   ibm RS/6000 SP
!
!$$$ end documentation block

!============================================================================
implicit none

integer(i_kind)          ,intent(in   ) :: n
real(r_kind)             ,intent(in   ) :: xt
real(r_kind),dimension(n),intent(in   ) :: x,q
real(r_kind),dimension(n),intent(  out) :: w,dw
!-----------------------------------------------------------------------------
real(r_kind),dimension(n)            :: d,pa,pb,dpa,dpb
integer(i_kind)                      :: j
!============================================================================
pa(1)=one
dpa(1)=zero
do j=1,n-1
   d(j)=xt-x(j)
   pa (j+1)=pa (j)*d(j)
   dpa(j+1)=dpa(j)*d(j)+pa(j)
enddo
d(n)=xt-x(n)

pb(n)=one
dpb(n)=zero
do j=n,2,-1
   pb (j-1)=pb (j)*d(j)
   dpb(j-1)=dpb(j)*d(j)+pb(j)
enddo
do j=1,n
   w (j)= pa(j)*pb (j)*q(j)
   dw(j)=(pa(j)*dpb(j)+dpa(j)*pb(j))*q(j)
enddo
end subroutine lagdw


!============================================================================
subroutine lagdw_tl(x,x_tl,xt,q,q_tl,w,w_tl,dw,dw_tl,n)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    lagdw_tl
!
!   prgrmmr:
!
! abstract:     
!
! program history log:
!   2008-05-09  safford - add subprogram doc block, rm unused uses
!
!   input argument list:
!     x     -
!     xt    - 
!     q     -
!     q_tl  -
!     n     - The number of points involved.
!
!   output argument list:
!     w     - 
!     dw    -
!     w_tl  -
!     dw_tl -
!
! attributes:
!   language:  f90
!   machine:   ibm RS/6000 SP
!
!$$$ end documentation block

!============================================================================
implicit none

integer(i_kind)          ,intent(in   ) :: n
real(r_kind)             ,intent(in   ) :: xt
real(r_kind),dimension(n),intent(in   ) :: x,q,x_tl,q_tl
real(r_kind),dimension(n),intent(  out) :: w,dw,w_tl,dw_tl
!-----------------------------------------------------------------------------
real(r_kind),dimension(n)            :: d,pa,pb,dpa,dpb
real(r_kind),dimension(n)            :: d_tl,pa_tl,pb_tl,dpa_tl,dpb_tl
integer(i_kind)                      :: j
!============================================================================
pa(1)=one
dpa(1)=zero
pa_tl(1)=zero
dpa_tl(1)=zero

do j=1,n-1
   d(j)=xt-x(j)
   d_tl(j)=-x_tl(j)
   pa    (j+1)=pa    (j)*d(j)
   pa_tl (j+1)=pa_tl (j)*d(j)+pa(j)*d_tl(j)
   dpa   (j+1)=dpa   (j)*d(j)+pa(j)
   dpa_tl(j+1)=dpa_tl(j)*d(j)+dpa(j)*d_tl(j)+pa_tl(j)
enddo
d(n)=xt-x(n)
d_tl(n)=-x_tl(n)

pb(n)=one
dpb(n)=zero
pb_tl(n)=zero
dpb_tl(n)=zero
do j=n,2,-1
   pb    (j-1)=pb    (j)*d(j)
   pb_tl (j-1)=pb_tl (j)*d(j)+pb (j)*d_tl(j)
   dpb   (j-1)=dpb   (j)*d(j)+pb (j)
   dpb_tl(j-1)=dpb_tl(j)*d(j)+dpb(j)*d_tl(j)+pb_tl(j)
enddo
do j=1,n
   w    (j)= pa   (j)*pb (j)*q(j)
   dw   (j)=(pa   (j)*dpb(j)+dpa(j)*pb    (j))*q(j)
   w_tl (j)=(pa_tl(j)*pb (j)+pa (j)*pb_tl (j))*q(j)+pa(j)*pb(j)*q_tl(j)
   dw_tl(j)=(pa_tl(j)*dpb(j)+pa (j)*dpb_tl(j)+dpa_tl(j)*pb(j)+dpa(j)*pb_tl(j))*q(j)+ &
            (pa   (j)*dpb(j)+dpa(j)*pb    (j))*q_tl(j)
enddo
end subroutine lagdw_tl


!============================================================================
subroutine lagdw_ad(x,x_ad,xt,q,q_ad,w_ad,dw_ad,n)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    lagdw_ad
!
!   prgrmmr:
!
! abstract:     
!
! program history log:
!   2008-05-09  safford - add subprogram doc block, rm unused uses
!
!   input argument list:
!     n     -
!     xt    -
!     x     -
!     q     -
!     q_ad  -
!     x_ad  -
!     w_ad  -
!     dw_ad -
!
!   output argument list:
!     q_ad  -
!     x_ad  -
!     w_ad  -
!     dw_ad -
!
! attributes:
!   language:  f90
!   machine:   ibm RS/6000 SP
!
!$$$ end documentation block

!============================================================================
implicit none

integer(i_kind)          ,intent(in   ) :: n
real(r_kind)             ,intent(in   ) :: xt
real(r_kind),dimension(n),intent(in   ) :: x,q
real(r_kind),dimension(n),intent(inout) :: q_ad,x_ad,w_ad,dw_ad
!-----------------------------------------------------------------------------
real(r_kind),dimension(n)            :: d,pa,pb,dpa,dpb
real(r_kind),dimension(n)            :: d_ad,pa_ad,pb_ad,dpa_ad,dpb_ad
integer(i_kind)                      :: j
!============================================================================

pa_ad=zero; dpb_ad=zero; dpa_ad=zero; pb_ad=zero
d_ad=zero

! passive variables
pa(1)=one
dpa(1)=zero
do j=1,n-1
   d(j)=xt-x(j)
   pa (j+1)=pa (j)*d(j)
   dpa(j+1)=dpa(j)*d(j)+pa(j)
enddo
d(n)=xt-x(n)
pb(n)=one
dpb(n)=zero
do j=n,2,-1
   pb (j-1)=pb (j)*d(j)
   dpb(j-1)=dpb(j)*d(j)+pb(j)
enddo
!
do j=n,1,-1
   pa_ad (j)=pa_ad (j)+ dpb(j)*q  (j)*dw_ad(j)
   dpb_ad(j)=dpb_ad(j)+ pa (j)*q  (j)*dw_ad(j)
   dpa_ad(j)=dpa_ad(j)+ pb (j)*q  (j)*dw_ad(j)
   pb_ad (j)=pb_ad (j)+ dpa(j)*q  (j)*dw_ad(j)
   q_ad  (j)=q_ad  (j)+(pa (j)*dpb(j)+dpa  (j)*pb(j))*dw_ad(j)
   pa_ad (j)=pa_ad (j)+ pb (j)*q  (j)*w_ad (j)
   pb_ad (j)=pb_ad (j)+ pa (j)*q  (j)*w_ad (j)
   q_ad  (j)=q_ad  (j)+ pa (j)*pb (j)*w_ad (j)
end do
do j=2,n
   dpb_ad(j)=dpb_ad(j)+d  (j)*dpb_ad(j-1)
   d_ad  (j)=d_ad  (j)+dpb(j)*dpb_ad(j-1)
   pb_ad (j)=pb_ad (j)       +dpb_ad(j-1)
   pb_ad (j)=pb_ad (j)+d  (j)*pb_ad (j-1)
   d_ad  (j)=d_ad  (j)+pb (j)*pb_ad (j-1)
enddo

dpb_ad(n)=zero ! not sure about it
pb_ad(n) =zero ! not sure about it

x_ad(n)=x_ad(n)-d_ad(n)
do j=n-1,1,-1
   dpa_ad(j)=dpa_ad(j)+d  (j)*dpa_ad(j+1)
   d_ad  (j)=d_ad  (j)+dpa(j)*dpa_ad(j+1)
   pa_ad (j)=pa_ad (j)       +dpa_ad(j+1)
   pa_ad (j)=pa_ad (j)+d  (j)*pa_ad (j+1)
   d_ad  (j)=d_ad  (j)+pa (j)*pa_ad (j+1)
   x_ad  (j)=x_ad  (j)       -d_ad  (j  )
enddo

end subroutine lagdw_ad


!============================================================================
subroutine slagdw(x,xt,aq,bq,w,dw,n)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    slagdw
!
!   prgrmmr:
!
! abstract:      Construct weights and their derivatives for interpolation 
!                to a given target based on a linearly weighted mixture of 
!                the pair of lagrange interpolators which omit the respective 
!                end points of the source nodes provided. The number of source 
!                points provided must be even and the nominal target interval 
!                is the unique central one. The linear weighting pair is 
!                determined by the relative location of the target within 
!                this central interval, or else the extreme values, 0 and 1, 
!                when target lies outside this interval.  The objective is to 
!                provide an interpolator whose derivative is continuous.
!
! program history log:
!   2008-05-09  safford - add subprogram doc block, rm unused uses
!
!   input argument list:
!     x     - Grid abscissae (n points)
!     xt    - Target abscissa
!     aq,bq - Q factors (denominators of the lagrange weight formula for
!             the two respective consecutive subsets of n-1 grid points)
!     n     - Even number of grid points involved in the interpolation
!
!   output argument list:
!     w     - Final weights (n values)
!     dw    - Final derivatives, dw/dx, of weights w (n values)
!
! attributes:
!   language:  f90
!   machine:   ibm RS/6000 SP
!
!$$$ end documentation block

!============================================================================
implicit none

integer(i_kind)            ,intent(in   ) :: n
real(r_kind)               ,intent(in   ) :: xt
real(r_kind),dimension(n)  ,intent(in   ) :: x
real(r_kind),dimension(n-1),intent(in   ) :: aq,bq
real(r_kind),dimension(n)  ,intent(  out) :: w,dw
!-----------------------------------------------------------------------------
real(r_kind),dimension(n)               :: aw,bw,daw,dbw
real(r_kind)                            :: xa,xb,dwb,wb
integer(i_kind)                         :: na
!============================================================================
call lagdw(x(1:n-1),xt,aq,aw(1:n-1),daw(1:n-1),n-1)
call lagdw(x(2:n  ),xt,bq,bw(2:n  ),dbw(2:n  ),n-1)
aw(n)=zero
daw(n)=zero
bw(1)=zero
dbw(1)=zero
na=n/2
if(na*2 /= n)stop 'In slagdw; n must be even'
xa =x(na     )
xb =x(na+1)
dwb=one/(xb-xa)
wb =(xt-xa)*dwb
if    (wb>one )then
   wb =one
   dwb=zero
elseif(wb<zero)then
   wb =zero
   dwb=zero
endif
bw =bw -aw
dbw=dbw-daw

w =aw +wb*bw
dw=daw+wb*dbw+dwb*bw
end subroutine slagdw


!============================================================================
subroutine slagdw_tl(x,x_tl,xt,aq,aq_tl,bq,bq_tl,dw,dw_tl,n)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    slagdw_tl
!
!   prgrmmr:
!
! abstract:      
!
! program history log:
!   2008-05-09  safford - add subprogram doc block, rm unused uses
!
!   input argument list:
!     n
!     xt
!     x,x_tl
!     aq,bq,aq_tl,bq_tl
!     dw_tl,dw
!
!   output argument list:
!     dw_tl,dw
!
! attributes:
!   language:  f90
!   machine:   ibm RS/6000 SP
!
!$$$ end documentation block

!============================================================================
implicit none

integer(i_kind)            ,intent(in   ) :: n
real(r_kind)               ,intent(in   ) :: xt
real(r_kind),dimension(n)  ,intent(in   ) :: x,x_tl
real(r_kind),dimension(n-1),intent(in   ) :: aq,bq,aq_tl,bq_tl
real(r_kind),dimension(n)  ,intent(  out) :: dw_tl,dw
!-----------------------------------------------------------------------------
real(r_kind),dimension(n)               :: aw,bw,daw,dbw
real(r_kind),dimension(n)               :: aw_tl,bw_tl,daw_tl,dbw_tl
real(r_kind)                            :: xa,xb,dwb,wb
real(r_kind)                            :: xa_tl,xb_tl,dwb_tl,wb_tl
integer(i_kind)                         :: na
!============================================================================
call lagdw_tl(x(1:n-1),x_tl(1:n-1),xt,aq,aq_tl,aw(1:n-1),aw_tl(1:n-1),&
             daw(1:n-1),daw_tl(1:n-1),n-1)
call lagdw_tl(x(2:n  ),x_tl(2:n  ),xt,bq,bq_tl,bw(2:n  ),bw_tl(2:n  ),&
             dbw(2:n  ),dbw_tl(2:n  ),n-1)
aw(n) =zero
daw(n)=zero
bw(1) =zero
dbw(1)=zero
!
aw_tl(n) =zero
daw_tl(n)=zero
bw_tl(1) =zero
dbw_tl(1)=zero
na=n/2
if(na*2 /= n)stop 'In slagdw; n must be even'
xa   =x   (na)
xa_tl=x_tl(na)
xb   =x   (na+1)
xb_tl=x_tl(na+1)
dwb   = one          /(xb-xa)
dwb_tl=-(xb_tl-xa_tl)/(xb-xa)**2
wb   =             (xt-xa)*dwb
wb_tl=(-xa_tl)*dwb+(xt-xa)*dwb_tl
if    (wb > one)then
   wb    =one
   dwb   =zero
   wb_tl =zero
   dwb_tl=zero
elseif(wb < zero)then
   wb    =zero
   dwb   =zero
   wb_tl =zero
   dwb_tl=zero

endif

bw    =bw    -aw
bw_tl =bw_tl -aw_tl
dbw   =dbw   -daw
dbw_tl=dbw_tl-daw_tl

!w=aw+wb*bw
dw   =daw   + wb   *dbw           + dwb   *bw
dw_tl=daw_tl+(wb_tl*dbw+wb*dbw_tl)+(dwb_tl*bw+dwb*bw_tl)
end subroutine slagdw_tl


!============================================================================
subroutine slagdw_ad(x,x_ad,xt,aq,aq_ad,bq,bq_ad,w_ad,dw,dw_ad,n)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    slagdw_ad
!
!   prgrmmr:
!
! abstract:      
!
! program history log:
!   2008-05-09  safford - add subprogram doc block, rm unused uses
!
!   input argument list:
!     n
!     xt
!     x
!     aq,bq
!     aq_ad,bq_ad
!     x_ad,dw_ad,w_ad
!
!   output argument list:
!     aq_ad,bq_ad
!     dw
!     x_ad,dw_ad,w_ad
!
! attributes:
!   language:  f90
!   machine:   ibm RS/6000 SP
!
!$$$ end documentation block

!============================================================================
implicit none
integer(i_kind)            ,intent(in   ) :: n
real(r_kind)               ,intent(in   ) :: xt
real(r_kind),dimension(n)  ,intent(in   ) :: x
real(r_kind),dimension(n-1),intent(in   ) :: aq,bq
real(r_kind),dimension(n-1),intent(inout) :: aq_ad,bq_ad
real(r_kind),dimension(n)  ,intent(  out) :: dw
real(r_kind),dimension(n)  ,intent(inout) :: x_ad,dw_ad,w_ad
!-----------------------------------------------------------------------------
real(r_kind),dimension(n)                 :: aw,bw,daw,dbw
real(r_kind),dimension(n)                 :: aw_ad,bw_ad,daw_ad,dbw_ad
real(r_kind)                              :: xa,xb,dwb,wb
real(r_kind)                              :: xa_ad,xb_ad,wb_ad,dwb_ad
integer(i_kind)                           :: na
!============================================================================
daw_ad=zero;wb_ad=zero;dbw_ad=zero;dwb_ad=zero;bw_ad=zero
aw_ad =zero;xb_ad=zero;xa_ad =zero

! passive variables needed (from forward code)
call lagdw(x(1:n-1),xt,aq,aw(1:n-1),daw(1:n-1),n-1)
call lagdw(x(2:n  ),xt,bq,bw(2:n  ),dbw(2:n  ),n-1)

aw(n) =zero
daw(n)=zero
bw(1) =zero
dbw(1)=zero
na=n/2
if(na*2 /= n)stop 'In slagdw; n must be even'
xa =x(na     )
xb =x(na+1)
dwb=one/(xb-xa)
wb =(xt-xa)*dwb
if    (wb>one)then
   wb =one
   dwb=zero
elseif(wb<zero)then
   wb =zero
   dwb=zero
endif
bw =bw -aw
dbw=dbw-daw
!w=aw+wb*bw ! not needed
dw=daw+wb*dbw+dwb*bw ! not needed
!
!
daw_ad=daw_ad+dw_ad
wb_ad =wb_ad +dot_product(dbw,dw_ad)
dbw_ad=dbw_ad+wb*dw_ad
dwb_ad=dwb_ad+dot_product(bw,dw_ad)
bw_ad =bw_ad +dwb*dw_ad

aw_ad=aw_ad+w_ad
wb_ad=wb_ad+dot_product(bw,w_ad)
bw_ad=bw_ad+wb*w_ad

daw_ad=daw_ad-dbw_ad
aw_ad =aw_ad -bw_ad

if(wb>one)then
   wb_ad =zero
   dwb_ad=zero
elseif(wb<zero)then
   wb_ad =zero
   dwb_ad=zero
endif

xa_ad     =xa_ad-wb_ad*dwb
dwb_ad    =dwb_ad+(xt-xa)*wb_ad
xb_ad     =xb_ad-(dwb_ad/(xb-xa)**2)
xa_ad     =xa_ad+(dwb_ad/(xb-xa)**2)
x_ad(na+1)=x_ad(na+1)+xb_aD
x_ad(na  )=x_ad(na  )+xa_ad

dbw_ad(1)=zero;bw_ad(1)=zero
daw_ad(n)=zero;aw_ad(n)=zero

call lagdw_ad(x(2:n  ),x_ad(2:n  ),xt,bq,bq_ad,bw_ad(2:n  ),&
             dbw_ad(2:n  ),n-1)

call lagdw_ad(x(1:n-1),x_ad(1:n-1),xt,aq,aq_ad,aw_ad(1:n-1),&
             daw_ad(1:n-1),n-1)

end subroutine slagdw_ad
!============================================================================
end module lagmod

