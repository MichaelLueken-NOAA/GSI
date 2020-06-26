module intswcpmod

!$$$ module documentation block
!           .      .    .                                       .
! module:   intswcpmod    module for intswcp and its tangent linear intswcp_tl
!   prgmmr:
!
! abstract: module for intswcp and its tangent linear intswcp_tl
!
! program history log:
!
! subroutines included:
!   sub intswcp_
!
! variable definitions:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

use m_obsnode, only: obsnode
use m_swcpnode, only: swcpnode
use m_swcpnode, only: swcpnode_typecast
use m_swcpnode, only: swcpnode_nextcast
use m_obsdiagnode, only: obsdiagnode_set
implicit none

private
public intswcp

interface intswcp; module procedure &
   intswcp_
end interface

contains

subroutine intswcp_(swcphead,rval,sval)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    intswcp       apply nonlin qc obs operator for swcp
!   prgmmr: Ting-Chi Wu        org: CIRA/CSU                date: 2017-06-28 
!
! abstract: apply observation operator and adjoint for solid-water content path
!           with addition of nonlinear qc.
!
! program history log:
!   2017-06-28  Ting-Chi Wu - mimic the structure in intpw.f90 and intgps.f90 
!                           - intswcp.f90 includes 2 tl/adj options
!                             1) when l_wcp_cwm = .false.: 
!                                operator = f(t,p,q)
!                             2) when l_wcp_cwm = .true. and cwm partition6: 
!                                 operator = f(qi,qs,qg,qh) partition6
!
!   input argument list:
!     swcphead   - obs type pointer to obs structure
!     st       - t increment in grid space
!     sp       - p increment in grid space
!     sq       - q increment in grid space
!     sqi      - qi increment in grid space
!     sqs      - qs increment in grid space
!     sqg      - qg increment in grid space
!     sqh      - qh increment in grid space
!
!   output argument list:
!     rt       - results of t from swcp observation operator 
!     rp       - results of p from swcp observation operator 
!     rq       - results of q from swcp observation operator 
!     rqi      - results of qi from swcp observation operator 
!     rqs      - results of qs from swcp observation operator 
!     rqg      - results of qg from swcp observation operator 
!     rqh      - results of qh from swcp observation operator 
!
!   comments:
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind
  use obsmod, only: lsaveobsens,l_do_adjoint,luse_obsdiag
  use obsmod, only: l_wcp_cwm
  use gridmod, only: nsig
  use qcmod, only: nlnqc_iter,varqc_iter
  use constants, only: zero,half,one,tiny_r_kind,cg_term,r3600
  use jfunc, only: jiter
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundlegetpointer
  use gsi_4dvar, only: ladtest_obs
  implicit none

! Declare passed variables
  class(obsnode), pointer, intent(in   ) :: swcphead
  type(gsi_bundle)        ,intent(in   ) :: sval
  type(gsi_bundle)        ,intent(inout) :: rval

! Declare local variables
  integer(i_kind) k,ier,istatus
  integer(i_kind),dimension(nsig):: i1,i2,i3,i4
! real(r_kind) penalty
  real(r_kind) :: t_tl,p_tl,q_tl
  real(r_kind) :: t_ad,p_ad,q_ad
  real(r_kind) :: qi_tl,qs_tl,qg_tl,qh_tl
  real(r_kind) :: qi_ad,qs_ad,qg_ad,qh_ad
  real(r_kind) val,w1,w2,w3,w4
  real(r_kind) cg_swcp,grad,p0,wnotgross,wgross,pg_swcp
  real(r_kind),pointer,dimension(:) :: st, sp, sq
  real(r_kind),pointer,dimension(:) :: sqi, sqs, sqg, sqh
  real(r_kind),pointer,dimension(:) :: rt, rp, rq
  real(r_kind),pointer,dimension(:) :: rqi, rqs, rqg, rqh
  type(swcpnode), pointer :: swcpptr

! If no swcp data return
  if(.not. associated(swcphead))return
! Retrieve pointers
! Simply return if any pointer not found
  ier=0

  if (.not.l_wcp_cwm) then

     call gsi_bundlegetpointer(sval,'tsen',st,istatus);ier=istatus+ier
     call gsi_bundlegetpointer(sval,'prse',sp,istatus);ier=istatus+ier
     call gsi_bundlegetpointer(sval,'q'   ,sq,istatus);ier=istatus+ier
     call gsi_bundlegetpointer(rval,'tsen',rt,istatus);ier=istatus+ier
     call gsi_bundlegetpointer(rval,'prse',rp,istatus);ier=istatus+ier
     call gsi_bundlegetpointer(rval,'q'   ,rq,istatus);ier=istatus+ier

  else

     call gsi_bundlegetpointer(sval,'qi',sqi,istatus);ier=istatus+ier
     call gsi_bundlegetpointer(sval,'qs',sqs,istatus);ier=istatus+ier
     call gsi_bundlegetpointer(sval,'qg',sqg,istatus);ier=istatus+ier
     call gsi_bundlegetpointer(sval,'qh',sqh,istatus);ier=istatus+ier
     call gsi_bundlegetpointer(rval,'qi',rqi,istatus);ier=istatus+ier
     call gsi_bundlegetpointer(rval,'qs',rqs,istatus);ier=istatus+ier
     call gsi_bundlegetpointer(rval,'qg',rqg,istatus);ier=istatus+ier
     call gsi_bundlegetpointer(rval,'qh',rqh,istatus);ier=istatus+ier

  endif ! l_wcp_cwm

  if(ier/=0)return

  !swcpptr => swcphead
  swcpptr => swcpnode_typecast(swcphead)
  do while (associated(swcpptr))

     do k=1,nsig
        i1(k)=swcpptr%ij(1,k)
        i2(k)=swcpptr%ij(2,k)
        i3(k)=swcpptr%ij(3,k)
        i4(k)=swcpptr%ij(4,k)
     enddo
     w1=swcpptr%wij(1)
     w2=swcpptr%wij(2)
     w3=swcpptr%wij(3)
     w4=swcpptr%wij(4)
     
     val=zero

!    Forward model (linear operator)

     if (.not.l_wcp_cwm) then
        do k=1,nsig
           t_tl=w1* st(i1(k))+w2* st(i2(k))+w3* st(i3(k))+w4* st(i4(k))
           p_tl=w1* sp(i1(k))+w2* sp(i2(k))+w3* sp(i3(k))+w4* sp(i4(k))
           q_tl=w1* sq(i1(k))+w2* sq(i2(k))+w3* sq(i3(k))+w4* sq(i4(k))
           val = val + ( t_tl*swcpptr%jac_t(k) + &
                         p_tl*swcpptr%jac_p(k) + & 
                         q_tl*swcpptr%jac_q(k) ) ! tpwcon*r10*(piges(k)-piges(k+1)) already did in setupswcp.f90
        end do
     else
        do k=1,nsig
           qi_tl=w1* sqi(i1(k))+w2* sqi(i2(k))+w3* sqi(i3(k))+w4* sqi(i4(k))
           qs_tl=w1* sqs(i1(k))+w2* sqs(i2(k))+w3* sqs(i3(k))+w4* sqs(i4(k))
           qg_tl=w1* sqg(i1(k))+w2* sqg(i2(k))+w3* sqg(i3(k))+w4* sqg(i4(k))
           qh_tl=w1* sqh(i1(k))+w2* sqh(i2(k))+w3* sqh(i3(k))+w4* sqh(i4(k))
           val = val + ( qi_tl*swcpptr%jac_qi(k) + & 
                         qs_tl*swcpptr%jac_qs(k) + &
                         qg_tl*swcpptr%jac_qg(k) + & 
                         qh_tl*swcpptr%jac_qh(k) ) ! tpwcon*r10*(piges(k)-piges(k+1)) already did in setupswcp.f90
        end do
     endif ! l_wcp_cwm

     if(luse_obsdiag)then
        if (lsaveobsens) then
           grad = val*swcpptr%raterr2*swcpptr%err2
           !-- swcpptr%diags%obssen(jiter) = grad
           call obsdiagnode_set(swcpptr%diags,jiter=jiter,obssen=grad)
        else
           !-- if (swcpptr%luse) swcpptr%diags%tldepart(jiter)=val
           if (swcpptr%luse) call obsdiagnode_set(swcpptr%diags,jiter=jiter,tldepart=val)
        endif
     end if

     if (l_do_adjoint) then
        if (.not. lsaveobsens) then
!          Difference from observation
           if( .not. ladtest_obs)   val=val-swcpptr%res
!          needed for gradient of nonlinear qc operator
           if (nlnqc_iter .and. swcpptr%pg > tiny_r_kind .and.  &
                                swcpptr%b  > tiny_r_kind) then
              pg_swcp=swcpptr%pg*varqc_iter
              cg_swcp=cg_term/swcpptr%b
              wnotgross= one-pg_swcp
              wgross = pg_swcp*cg_swcp/wnotgross
              p0   = wgross/(wgross+exp(-half*swcpptr%err2*val**2))
              val = val*(one-p0)
           endif
           if( ladtest_obs) then
              grad = val
           else
              grad = val*swcpptr%raterr2*swcpptr%err2
           end if
        endif

!       Adjoint
        if (.not.l_wcp_cwm) then
           do k=1,nsig
              t_ad = grad*swcpptr%jac_t(k)
              rt(i1(k))=rt(i1(k))+w1*t_ad
              rt(i2(k))=rt(i2(k))+w2*t_ad
              rt(i3(k))=rt(i3(k))+w3*t_ad
              rt(i4(k))=rt(i4(k))+w4*t_ad
              p_ad = grad*swcpptr%jac_p(k)
              rp(i1(k))=rp(i1(k))+w1*p_ad
              rp(i2(k))=rp(i2(k))+w2*p_ad
              rp(i3(k))=rp(i3(k))+w3*p_ad
              rp(i4(k))=rp(i4(k))+w4*p_ad
              q_ad = grad*swcpptr%jac_q(k)
              rq(i1(k))=rq(i1(k))+w1*q_ad
              rq(i2(k))=rq(i2(k))+w2*q_ad
              rq(i3(k))=rq(i3(k))+w3*q_ad
              rq(i4(k))=rq(i4(k))+w4*q_ad
           enddo
        else
           do k=1,nsig
              qi_ad = grad*swcpptr%jac_qi(k)
              rqi(i1(k))=rqi(i1(k))+w1*qi_ad
              rqi(i2(k))=rqi(i2(k))+w2*qi_ad
              rqi(i3(k))=rqi(i3(k))+w3*qi_ad
              rqi(i4(k))=rqi(i4(k))+w4*qi_ad
              qs_ad = grad*swcpptr%jac_qs(k)
              rqs(i1(k))=rqs(i1(k))+w1*qs_ad
              rqs(i2(k))=rqs(i2(k))+w2*qs_ad
              rqs(i3(k))=rqs(i3(k))+w3*qs_ad
              rqs(i4(k))=rqs(i4(k))+w4*qs_ad
              qg_ad = grad*swcpptr%jac_qg(k)
              rqg(i1(k))=rqg(i1(k))+w1*qg_ad
              rqg(i2(k))=rqg(i2(k))+w2*qg_ad
              rqg(i3(k))=rqg(i3(k))+w3*qg_ad
              rqg(i4(k))=rqg(i4(k))+w4*qg_ad
              qh_ad = grad*swcpptr%jac_qh(k)
              rqh(i1(k))=rqh(i1(k))+w1*qh_ad
              rqh(i2(k))=rqh(i2(k))+w2*qh_ad
              rqh(i3(k))=rqh(i3(k))+w3*qh_ad
              rqh(i4(k))=rqh(i4(k))+w4*qh_ad
           enddo
        endif ! l_wcp_cwm
     endif ! l_do_adjoint


     !swcpptr => swcpptr%llpoint
     swcpptr => swcpnode_nextcast(swcpptr)

  end do

  return

end subroutine intswcp_

end module intswcpmod
