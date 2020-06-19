module intlightmod

!$$$ module documentation block
!           .      .    .                                       .
! module:   intlightmod    int module for the observation operator for lightning flash rate (lfr)
!                 
!   prgmmr: k apodaca <karina.apodaca@colostate.edu>
!      org: CSU/CIRA, Data Assimilation Group
!     date: 2016-05-04
!
! abstract: module for the tangent linear (flashrate_tl) and adjoint models (flashrate_ad)
!           of lfr
!
! program history log:
!   2016-05-04  apodaca  - implement tl and ad of the lfr observation operator  
!   2018-02-08  apodaca  - replaced ob_type with polymorphic obsnode through type casting
!   2019-03-01  j guo    - encapsulated access to obsdiagnode through obsdiagnode_set()
!
! subroutines included:
!   sub intlight_
!
! variable definitions:
!
! attributes:
!   language: Fortran 90 and/or above
!   machine:
!
!$$$ end documentation block
use m_obsnode, only: obsnode
use m_lightnode, only: lightnode
use m_lightnode, only: lightnode_typecast
use m_lightnode, only: lightnode_nextcast
use m_obsdiagnode, only: obsdiagnode_set
implicit none

private
public intlight

interface intlight; module procedure &
   intlight_
end interface

contains

subroutine intlight_(lighthead,rval,sval)

!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    intlight      tl and subsequent ad of the forward observation operator for lfr
!     prgmmr:    k apodaca          
!        org:    CSU/CIRA, Data Assimilation Group             
!       date:    2016-05-04
!
! abstract: In this program, the tangent linear and adjoint models of a 
!           lightning flash rate observation operator are calculated 
!           using a 12-point horizontal grid to calculate finite-difference 
!           derivatives and for interpolation in specific quadrants.
! 
!           The tangent linear equations represent a way to map 
!           the perturbation vectors for the control variables 
!           q, qi, qs, qg, t, u, and v.
!
!           Moreover, the adjoint equations map the sensitivity gradient 
!           vectors for the control variables (q, qi, qs, qg, t, u, v), 
!           thus providing a first order aproximation or linear projection
!           of the sesitivity (impact) of observations.
!
! program history log:
!     2018-01-18 k apodaca revision of ad code
!     2018-08-18 k apodaca add a the tl and ad of second oservation operator for lightning 
!                          observations suitable for non-hydrostatic, cloud-resolving models
!                          with additional ice-phase hydrometeor control variables
!
!   input argument list:
!     lighthead   - obs type pointer to obs structure
!     sq          - q  increment in grid space
!     sqi         - qi increment in grid space
!     sqs         - qs increment in grid space
!     sqg         - qg increment in grid space
!     st          - t  increment in grid space
!     su          - u  increment in grid space
!     sv          - v  increment in grid space
!
!   output argument list:
!     rq, rqi, rqs, rqg       - control variabble updates resulting from 
!     rt, ru, rv                the assimilation of lightning flash rate 
!                               observations 
!
!   comments:
!
! attributes:
!   language: Fortran 90 and/or above
!   machine: 
!
!$$$ end subprogram documentation block

  use kinds,         only: r_kind,i_kind
  use obsmod,        only: lsaveobsens,l_do_adjoint,luse_obsdiag
  use gridmod,       only: nsig
  use gridmod,       only: wrf_mass_regional,regional
  use qcmod,         only: nlnqc_iter,varqc_iter
  use constants,     only: zero,fv,one,half,two,tiny_r_kind,cg_term
  use jfunc,         only: jiter
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundlegetpointer
  use gsi_4dvar,     only: ladtest_obs
  implicit none

! Declare passed variables
  class(obsnode), pointer, intent(in   ) :: lighthead
  type(gsi_bundle),        intent(in   ) :: sval
  type(gsi_bundle),        intent(inout) :: rval

! Declare local variables
  integer(i_kind) k,ier,istatus
  integer(i_kind),dimension(nsig)           :: i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12
  real(r_kind) val,w1,w2,w3,w4
  real(r_kind) cg_light,grad,p0,wnotgross,wgross,pg_light
  real(r_kind),pointer,dimension(:)         :: sq,sqi,sqs,sqg,su,sv,st
  real(r_kind),pointer,dimension(:)         :: rq,rqi,rqs,rqg,ru,rv,rt
  type(lightnode),  pointer             :: lightptr

! Variables for tl and ad of lightning flash rate
  real(r_kind),dimension(1:nsig)            :: z_tl 
  real(r_kind),dimension(1:nsig)            :: horiz_adv_tl
  real(r_kind),dimension(1:nsig)            :: vert_adv_tl
  real(r_kind),dimension(1:nsig)            :: w_tl
  real(r_kind)                              :: wmaxi1_tl,wmaxi2_tl,wmaxi3_tl,wmaxi4_tl
  real(r_kind)                              :: flashrate_tl,flashratei1_tl,flashratei2_tl
  real(r_kind)                              :: flashratei3_tl, flashratei4_tl
  real(r_kind)                              :: h1i1_tl,h1i2_tl,h1i3_tl,h1i4_tl
  real(r_kind)                              :: h2i1_tl,h2i2_tl,h2i3_tl,h2i4_tl
  real(r_kind)                              :: totice_colinti1_tl,totice_colinti2_tl
  real(r_kind)                              :: totice_colinti3_tl,totice_colinti4_tl
  real(r_kind)                              :: htot_tl,htoti1_tl,htoti2_tl,htoti3_tl,htoti4_tl
  real(r_kind)                              :: flashrate_ad,flashratei1_ad,flashratei2_ad
  real(r_kind)                              :: flashratei3_ad,flashratei4_ad
  real(r_kind)                              :: wmaxi1_ad,wmaxi2_ad,wmaxi3_ad,wmaxi4_ad
  real(r_kind)                              :: h1i1_ad,h1i2_ad,h1i3_ad,h1i4_ad
  real(r_kind)                              :: h2i1_ad,h2i2_ad,h2i3_ad,h2i4_ad
  real(r_kind)                              :: totice_colinti1_ad,totice_colinti2_ad
  real(r_kind)                              :: totice_colinti3_ad,totice_colinti4_ad
  real(r_kind)                              :: htot_ad,htoti1_ad,htoti2_ad,htoti3_ad,htoti4_ad     
  real(r_kind),dimension(1:nsig)            :: z_ad
  real(r_kind),dimension(1:nsig)            :: w_ad
  real(r_kind),dimension(1:nsig)            :: vert_adv_ad,horiz_adv_ad
  real(r_kind),dimension(1:nsig)            :: diffq
  real(r_kind),dimension(1:nsig)            :: difft
  real(r_kind),dimension(1:nsig)            :: diffz
!  wmax variables for lightning flash rate
  real(r_kind)                             :: wmax   
  real(r_kind),parameter                   :: k3=0.95_r_kind
  
!  Output files
!  character :: tlh_file*40
     

!  If no light data return
  if(.not. associated(lighthead))return
! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'q',sq,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'q',rq,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(sval,'tsen',st,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'tsen',rt,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(sval,'u',su,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'u',ru,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(sval,'v',sv,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'v',rv,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(sval,'qi',sqi,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'qi',rqi,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(sval,'qg',sqg,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'qg',rqg,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(sval,'qs',sqs,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'qs',rqs,istatus);ier=istatus+ier
  if(ier/=0)return



  lightptr => lightnode_typecast(lighthead)
  do while (associated(lightptr))

! Load location information into local variables

     w1=lightptr%wij(1)
     w2=lightptr%wij(2)
     w3=lightptr%wij(3)
     w4=lightptr%wij(4)

     do k=1,nsig
        i1(k)=lightptr%ij(1,k)
        i2(k)=lightptr%ij(2,k)
        i3(k)=lightptr%ij(3,k)
        i4(k)=lightptr%ij(4,k)
        i5(k)=lightptr%ij(5,k)
        i6(k)=lightptr%ij(6,k)
        i7(k)=lightptr%ij(7,k)
        i8(k)=lightptr%ij(8,k)
        i9(k)=lightptr%ij(9,k)
        i10(k)=lightptr%ij(10,k)
        i11(k)=lightptr%ij(11,k)
        i12(k)=lightptr%ij(12,k)  
     end do 
     
!                .      .    .                                       .

! In the case of lightning observations (e.g. goes/glm), the schematic shown below is
! used for bi-linear interpolation of background fields to the location of an observation 
! (+) and for the finite-difference derivation method used in the calculation of the tl of
! the observation operator for lightning flash rate. Calculations are done
! at each quadrant, i.e., central, north, south, east, and west.
!
!         i6-------i8
!          |       |
!          |       |
! i10-----i2-------i4-----i12
!  |       |       |       |
!  |       |     + |       |
! i9------i1-------i3-----i11
!          |       |
!          |       |
!         i5-------i7
!

!                .      .    .                                       .
     
! In the following section, the tangent linear of the lightning flash rate observation  
! operator is calculated by being broken into parts.                                               

! Tangent linear of height (z)

     z_tl(:)=zero
     horiz_adv_tl(:)=zero

     do k=2,nsig-1

        z_tl(i1(1))=lightptr%jac_z0i1
        z_tl(i2(1))=lightptr%jac_z0i2
        z_tl(i3(1))=lightptr%jac_z0i3
        z_tl(i4(1))=lightptr%jac_z0i4
        z_tl(i5(1))=lightptr%jac_z0i5
        z_tl(i6(1))=lightptr%jac_z0i6
        z_tl(i7(1))=lightptr%jac_z0i7
        z_tl(i8(1))=lightptr%jac_z0i8
        z_tl(i9(1))=lightptr%jac_z0i9
        z_tl(i10(1))=lightptr%jac_z0i10
        z_tl(i11(1))=lightptr%jac_z0i11
        z_tl(i12(1))=lightptr%jac_z0i12


        z_tl(i1(k))=z_tl(i1(k-1))+lightptr%jac_vertti1(k)*st(i1(k))          &
                   +lightptr%jac_vertqi1(k)*sq(i1(k))

        z_tl(i2(k))=z_tl(i2(k-1))+lightptr%jac_vertti2(k)*st(i2(k))          &
                   +lightptr%jac_vertqi2(k)*sq(i2(k)) 

        z_tl(i3(k))=z_tl(i3(k-1))+lightptr%jac_vertti3(k)*st(i3(k))          &
                   +lightptr%jac_vertqi3(k)*sq(i3(k))
   
        z_tl(i4(k))=z_tl(i4(k-1))+lightptr%jac_vertti4(k)*st(i4(k))          &
                   +lightptr%jac_vertqi4(k)*sq(i4(k))

        z_tl(i5(k))=z_tl(i5(k-1))+lightptr%jac_vertti5(k)*st(i5(k))          &
                   +lightptr%jac_vertqi5(k)*sq(i5(k))

        z_tl(i6(k))=z_tl(i6(k-1))+lightptr%jac_vertti6(k)*st(i6(k))          &
                   +lightptr%jac_vertqi6(k)*sq(i6(k))

        z_tl(i7(k))=z_tl(i7(k-1))+lightptr%jac_vertti7(k)*st(i7(k))          &
                   +lightptr%jac_vertqi7(k)*sq(i7(k))

        z_tl(i8(k))=z_tl(i8(k-1))+lightptr%jac_vertti8(k)*st(i8(k))          &
                   +lightptr%jac_vertqi8(k)*sq(i8(k))

        z_tl(i9(k))=z_tl(i9(k-1))+lightptr%jac_vertti9(k)*st(i9(k))          &
                   +lightptr%jac_vertqi9(k)*sq(i9(k))

        z_tl(i10(k))=z_tl(i10(k-1))+lightptr%jac_vertti10(k)*st(i10(k))      &
                    +lightptr%jac_vertqi10(k)*sq(i10(k))

        z_tl(i11(k))=z_tl(i11(k-1))+lightptr%jac_vertti11(k)*st(i11(k))      & 
                  +lightptr%jac_vertqi11(k)*sq(i11(k))

        z_tl(i12(k))=z_tl(i12(k-1))+lightptr%jac_vertti12(k)*st(i12(k))      &
                  +lightptr%jac_vertqi12(k)*sq(i12(k))


! Tangent linear of the horizontal advection section

 
        horiz_adv_tl(i1(k))=lightptr%jac_zdxi1(k)*su(i1(k))                  &
                           +lightptr%jac_zdyi1(k)*sv(i1(k))                  &
                           +lightptr%jac_udxi1(k)*(z_tl(i3(k))-z_tl(i9(k)))  &
                           +lightptr%jac_vdyi1(k)*(z_tl(i2(k))-z_tl(i5(k)))
        horiz_adv_tl(i2(k))=lightptr%jac_zdxi2(k)*su(i2(k))                  &
                           +lightptr%jac_zdyi2(k)*sv(i2(k))                  &
                           +lightptr%jac_udxi2(k)*(z_tl(i4(k))-z_tl(i10(k))) &
                           +lightptr%jac_vdyi2(k)*(z_tl(i6(k))-z_tl(i1 (k)))

        horiz_adv_tl(i3(k))=lightptr%jac_zdxi3(k)*su(i3(k))                  &
                           +lightptr%jac_zdyi3(k)*sv(i3(k))                  &
                           +lightptr%jac_udxi3(k)*(z_tl(i11(k))-z_tl(i1(k))) &
                           +lightptr%jac_vdyi3(k)*(z_tl(i4 (k))-z_tl(i7(k)))

        horiz_adv_tl(i4(k))=lightptr%jac_zdxi4(k)*su(i4(k))                  &
                           +lightptr%jac_zdyi4(k)*sv(i4(k))                  &
                           +lightptr%jac_udxi4(k)*(z_tl(i12(k))-z_tl(i2(k))) &
                           +lightptr%jac_vdyi4(k)*(z_tl(i8 (k))-z_tl(i3(k)))

     enddo ! do k=2,nsig-1

! Tangent linear of the vertical advection section

     vert_adv_tl(:)=zero
     w_tl(:)=zero

     do k=1,nsig-1

        vert_adv_tl(i1(k))=-lightptr%jac_vert(k)*lightptr%jac_sigdoti1(k)*    &
                          (((one+fv*lightptr%jac_qi1(k))*st(i1(k)))           & 
                          +(lightptr%jac_ti1(k)*fv*sq(i1(k))))

        vert_adv_tl(i2(k))=-lightptr%jac_vert(k)*lightptr%jac_sigdoti2(k)*    &
                          (((one+fv*lightptr%jac_qi2(k))*st(i2(k)))           &
                          +(lightptr%jac_ti2(k)*fv*sq(i2(k))))

        vert_adv_tl(i3(k))=-lightptr%jac_vert(k)*lightptr%jac_sigdoti3(k)*    &
                          (((one+fv*lightptr%jac_qi3(k))*st(i3(k)))           &
                          +(lightptr%jac_ti3(k)*fv*sq(i3(k))))

        vert_adv_tl(i4(k))=-lightptr%jac_vert(k)*lightptr%jac_sigdoti4(k)*    &
                          (((one+fv*lightptr%jac_qi4(k))*st(i4(k)))           &
                          +(lightptr%jac_ti4(k)*fv*sq(i4(k))))
      


! Tangent linear of vertical velocity


        w_tl(i1(k))=horiz_adv_tl(i1(k))+vert_adv_tl(i1(k))
        w_tl(i2(k))=horiz_adv_tl(i2(k))+vert_adv_tl(i2(k))
        w_tl(i3(k))=horiz_adv_tl(i3(k))+vert_adv_tl(i3(k))
        w_tl(i4(k))=horiz_adv_tl(i4(k))+vert_adv_tl(i4(k))

     enddo !do k=1,nsig-1
 
!                .      .    .                                       .
! Tangent linear of lightning flash rate

!                .      .    .                                       .
! Regional

     if (regional) then

!-- wrf-arw

        if (wrf_mass_regional) then

! Tangent linear - lightning flash rate as a function of
! vertical graupel flux within the mixed-phase region
! (-15 deg C)


           if (lightptr%kboti1 > zero) then
              h1i1_tl=lightptr%jac_qgmai1(lightptr%kboti1)*sqg(i1(lightptr%kboti1))+&
                      lightptr%jac_qgmbi1(lightptr%kboti1)*&
                      (half*(w_tl(i1(lightptr%kboti1))+w_tl(i1(lightptr%kboti1+1))))
              h1i1_tl=h1i1_tl/(abs(h1i1_tl))
           else
              h1i1_tl=zero
           endif

           if (lightptr%kboti2 > zero) then
              h1i2_tl=lightptr%jac_qgmai2(lightptr%kboti2)*sqg(i2(lightptr%kboti2))+&
                      lightptr%jac_qgmbi2(lightptr%kboti2)*&
                      (half*(w_tl(i2(lightptr%kboti2))+w_tl(i2(lightptr%kboti2+1))))
              h1i2_tl=h1i2_tl/(abs(h1i2_tl))
           else
              h1i2_tl=zero
           endif

           if (lightptr%kboti3 > zero) then
              h1i3_tl=lightptr%jac_qgmai3(lightptr%kboti3)*sqg(i3(lightptr%kboti3))+&
                      lightptr%jac_qgmbi3(lightptr%kboti3)*&
                      (half*(w_tl(i3(lightptr%kboti3))+w_tl(i3(lightptr%kboti3+1))))
              h1i3_tl=h1i3_tl/(abs(h1i3_tl))
           else
              h1i3_tl=zero
           endif
        
           if (lightptr%kboti4 > zero) then
              h1i4_tl=lightptr%jac_qgmai4(lightptr%kboti4)*sqg(i4(lightptr%kboti4))+&
                      lightptr%jac_qgmbi4(lightptr%kboti4)*&
                      (half*(w_tl(i4(lightptr%kboti4))+w_tl(i4(lightptr%kboti4+1))))
              h1i4_tl=h1i4_tl/(abs(h1i4_tl))
           else
              h1i4_tl=zero
           endif


! Tangent linear - lightning flash rate as a function of total column-integrated
! ice-phase hydrometeors

           totice_colinti1_tl=zero
           totice_colinti2_tl=zero
           totice_colinti3_tl=zero
           totice_colinti4_tl=zero

           do k=1,nsig-1

              totice_colinti1_tl = totice_colinti1_tl+lightptr%jac_icei1(k) * &
                                  (sqi(i1(k))+sqs(i1(k))+sqg(i1(k)))+&
                                   lightptr%jac_zicei1(k)*z_tl(i1(k))

              totice_colinti2_tl = totice_colinti2_tl+lightptr%jac_icei2(k) * &
                                  (sqi(i2(k))+sqs(i2(k))+sqg(i2(k)))+&
                                   lightptr%jac_zicei2(k)*z_tl(i2(k))

              totice_colinti3_tl = totice_colinti3_tl+lightptr%jac_icei3(k) * &
                                  (sqi(i3(k))+sqs(i3(k))+sqg(i3(k)))+&
                                   lightptr%jac_zicei3(k)*z_tl(i3(k))

              totice_colinti4_tl = totice_colinti4_tl+lightptr%jac_icei4(k) * &
                                  (sqi(i4(k))+sqs(i4(k))+sqg(i4(k)))+&
                                   lightptr%jac_zicei4(k)*z_tl(i4(k))

           enddo !do k=1,nsig-1

           h2i1_tl=(1-k3)*totice_colinti1_tl
           h2i2_tl=(1-k3)*totice_colinti2_tl
           h2i3_tl=(1-k3)*totice_colinti3_tl
           h2i4_tl=(1-k3)*totice_colinti4_tl


           htoti1_tl= h1i1_tl+h2i1_tl
           htoti2_tl= h1i2_tl+h2i2_tl
           htoti3_tl= h1i3_tl+h2i3_tl
           htoti4_tl= h1i4_tl+h2i4_tl

!  Interpolation of lightning flash rate to observation location (2d field)
!  forward model

           htot_tl = (w1*htoti1_tl + w2*htoti2_tl + &
                      w3*htoti3_tl + w4*htoti4_tl)
           val = htot_tl

        endif ! wrf_mass_regional      

     endif !if (regional) then
   
!                .      .    .                                       .
! Global 
   
     if (.not. regional) then ! Global

! Cloud mask

! If clouds are present, find the maximum value of vertical velocity
! (wmax_tl) at four points sorounding an observation (+)
! and amongst all vertical levels, otherwise set wmax_tl to zero.

        wmaxi1_tl=zero
        wmaxi2_tl=zero
        wmaxi3_tl=zero
        wmaxi4_tl=zero

        if (lightptr%jac_wmaxflagi1) then
           wmax=-1.e+10_r_kind
           do k=1,nsig-1
              if (w_tl(i1(k)) > wmax) then
                 lightptr%jac_kverti1=k
                 wmaxi1_tl=w_tl(i1(lightptr%jac_kverti1))
              endif
              if (wmaxi1_tl < zero) then
                 wmaxi1_tl=zero
              endif
           enddo ! k loop
        endif

        if (lightptr%jac_wmaxflagi2) then
           wmax=-1.e+10_r_kind
           do k=1,nsig-1
              if (w_tl(i2(k)) > wmax) then
                 lightptr%jac_kverti2=k
                 wmaxi2_tl=w_tl(i2(lightptr%jac_kverti2))
              endif
              if (wmaxi2_tl <  zero) then
                 wmaxi2_tl=zero
              endif
           enddo ! k loop
        endif

        if (lightptr%jac_wmaxflagi3) then
           wmax=-1.e+10_r_kind
           do k=1,nsig-1
              if (w_tl(i3(k)) > wmax) then
                 lightptr%jac_kverti3=k
                 wmaxi3_tl=w_tl(i3(lightptr%jac_kverti3))
              endif
              if (wmaxi3_tl <  zero) then
                 wmaxi3_tl=zero
              endif
           enddo ! k loop
        endif

        if (lightptr%jac_wmaxflagi4) then
           wmax=-1.e+10_r_kind
           do k=1,nsig-1
              if (w_tl(i4(k)) > wmax) then
                 lightptr%jac_kverti4=k
                 wmaxi4_tl=w_tl(i4(lightptr%jac_kverti4))
              endif
              if (wmaxi4_tl < zero) then
                 wmaxi4_tl=zero
              endif
           enddo ! k loop
        endif

! Tangent linear of lightning flash rate
    
        flashratei1_tl=lightptr%jac_fratei1*wmaxi1_tl
        flashratei2_tl=lightptr%jac_fratei1*wmaxi2_tl
        flashratei3_tl=lightptr%jac_fratei1*wmaxi3_tl
        flashratei4_tl=lightptr%jac_fratei1*wmaxi4_tl

!  Interpolation of lightning flash rate to observation location (2d field)
!  forward model

        flashrate_tl = (w1*flashratei1_tl + w2*flashratei2_tl + & 
                     w3*flashratei3_tl + w4*flashratei4_tl)
        val =  flashrate_tl

     end if ! global block

     if (luse_obsdiag)then
        if (lsaveobsens) then
           grad = val*lightptr%raterr2*lightptr%err2
           !-- lightptr%diags%obssen(jiter) = grad
           call obsdiagnode_set(lightptr%diags,jiter=jiter,obssen=grad)
        else
           !-- if (lightptr%luse) lightptr%diags%tldepart(jiter)=val
           if (lightptr%luse) call obsdiagnode_set(lightptr%diags,jiter=jiter,tldepart=val)
        endif
     end if 
   

!                .      .    .                                       .

! Adjoint test
 
     if (l_do_adjoint) then
! Difference from observation
        if (.not. lsaveobsens) then
           if (.not. ladtest_obs)  val=val-lightptr%res

!          needed for gradient of nonlinear qc operator
           if (nlnqc_iter .and. lightptr%pg > tiny_r_kind .and.  &
                                lightptr%b  > tiny_r_kind) then
              pg_light=lightptr%pg*varqc_iter
              cg_light=cg_term/lightptr%b
              wnotgross= one-pg_light
              wgross = pg_light*cg_light/wnotgross
              p0   = wgross/(wgross+exp(-half*lightptr%err2*val**2))
              val = val*(one-p0)
           endif

           if( ladtest_obs) then
              grad = val
           else
              grad = val*lightptr%raterr2*lightptr%err2
           end if
        endif


!                .      .    .                                       .

! Adjoint of the lightning flash rate observation operator   

!                .      .    .                                       .
! Variable initialization

        z_ad(:)=zero
        w_ad(:)=zero

! Regional

        if (regional) then

!-- wrf-arw

           if (wrf_mass_regional) then

              htot_ad=grad

! Adjoint - total lightning flash rate

              htoti1_ad=htoti1_ad+w1*htot_ad
              htoti2_ad=htoti2_ad+w1*htot_ad
              htoti3_ad=htoti3_ad+w1*htot_ad
              htoti4_ad=htoti4_ad+w1*htot_ad

              h1i1_ad=h1i1_ad+htoti1_ad
              h2i1_ad=h2i1_ad+htoti1_ad

              h1i2_ad=h1i2_ad+htoti2_ad
              h2i2_ad=h2i2_ad+htoti2_ad

              h1i3_ad=h1i3_ad+htoti3_ad
              h2i3_ad=h2i3_ad+htoti3_ad

              h1i4_ad=h1i4_ad+htoti4_ad
              h2i4_ad=h2i4_ad+htoti4_ad

              totice_colinti1_ad=totice_colinti1_ad+(1-k3)*h2i1_ad
              totice_colinti2_ad=totice_colinti2_ad+(1-k3)*h2i2_ad
              totice_colinti3_ad=totice_colinti3_ad+(1-k3)*h2i3_ad
              totice_colinti4_ad=totice_colinti4_ad+(1-k3)*h2i4_ad

! Adjoint - lightning flash rate as a function of total column-integrated
! ice-phase hydrometeors


              do k=nsig-1,1,-1

                 z_ad(i1(k))=z_ad(i1(k))+lightptr%jac_zicei1(k)*totice_colinti1_ad
                 rqi(i1(k))=rqi(i1(k))+lightptr%jac_icei1(k)*totice_colinti1_ad
                 rqs(i1(k))=rqs(i1(k))+lightptr%jac_icei1(k)*totice_colinti1_ad
                 rqg(i1(k))=rqg(i1(k))+lightptr%jac_icei1(k)*totice_colinti1_ad
                 totice_colinti1_ad=two*totice_colinti1_ad

                 z_ad(i2(k))=z_ad(i2(k))+lightptr%jac_zicei2(k)*totice_colinti2_ad
                 rqi(i2(k))=rqi(i2(k))+lightptr%jac_icei2(k)*totice_colinti2_ad
                 rqs(i2(k))=rqs(i2(k))+lightptr%jac_icei2(k)*totice_colinti2_ad
                 rqg(i2(k))=rqg(i2(k))+lightptr%jac_icei2(k)*totice_colinti2_ad
                 totice_colinti2_ad=two*totice_colinti2_ad

                 z_ad(i3(k))=z_ad(i3(k))+lightptr%jac_zicei3(k)*totice_colinti3_ad
                 rqi(i3(k))=rqi(i3(k))+lightptr%jac_icei3(k)*totice_colinti3_ad
                 rqs(i3(k))=rqs(i3(k))+lightptr%jac_icei3(k)*totice_colinti3_ad
                 rqg(i3(k))=rqg(i3(k))+lightptr%jac_icei3(k)*totice_colinti3_ad
                 totice_colinti3_ad=two*totice_colinti3_ad

                 z_ad(i4(k))=z_ad(i4(k))+lightptr%jac_zicei4(k)*totice_colinti4_ad
                 rqi(i4(k))=rqi(i4(k))+lightptr%jac_icei4(k)*totice_colinti4_ad
                 rqs(i4(k))=rqs(i4(k))+lightptr%jac_icei4(k)*totice_colinti4_ad
                 rqg(i4(k))=rqg(i4(k))+lightptr%jac_icei4(k)*totice_colinti4_ad
                 totice_colinti4_ad=two*totice_colinti4_ad

! Adjoint - lightning flash rate as a function of
! vertical graupel flux within the mixed-phase region
! (-15 deg c)

                 if (lightptr%kboti1 > zero) then
                    h1i1_ad=h1i1_ad+(h1i1_ad/abs(h1i1_ad))
                    rqg(i1(lightptr%kboti1))=rqg(i1(lightptr%kboti1))+&
                                             lightptr%jac_qgmai1(lightptr%kboti1)*h1i1_ad
                    w_ad(i1(lightptr%kboti1))=w_ad(i1(lightptr%kboti1))+&
                                              half*lightptr%jac_qgmbi1(lightptr%kboti1)*h1i1_ad
                    w_ad(i1(lightptr%kboti1+1))=w_ad(i1(lightptr%kboti1+1))+&
                                                half*lightptr%jac_qgmbi1(lightptr%kboti1)*h1i1_ad
                 else
                    h1i1_ad=zero
                    rqg(i1(lightptr%kboti1))=zero
                    w_ad(i1(lightptr%kboti1))=zero
                    w_ad(i1(lightptr%kboti1+1))=zero
                 endif

                 if (lightptr%kboti2 > zero) then
                    h1i2_ad=h1i2_ad+(h1i2_ad/abs(h1i2_ad))
                    rqg(i2(lightptr%kboti2))=rqg(i2(lightptr%kboti2))+&
                                             lightptr%jac_qgmai2(lightptr%kboti2)*h1i2_ad
                    w_ad(i2(lightptr%kboti2))=w_ad(i2(lightptr%kboti2))+&
                                              half*lightptr%jac_qgmbi2(lightptr%kboti2)*h1i2_ad
                    w_ad(i2(lightptr%kboti2+1))=w_ad(i2(lightptr%kboti2+1))+&
                                                half*lightptr%jac_qgmbi2(lightptr%kboti2)*h1i2_ad
                 else
                    h1i2_ad=zero
                    rqg(i2(lightptr%kboti2))=zero
                    w_ad(i2(lightptr%kboti2+1))=zero
                 endif

                 if (lightptr%kboti3 > zero) then
                    h1i3_ad=h1i3_ad+(h1i3_ad/abs(h1i3_ad))
                    rqg(i3(lightptr%kboti3))=rqg(i3(lightptr%kboti3))+&
                                             lightptr%jac_qgmai3(lightptr%kboti3)*h1i3_ad
                    w_ad(i3(lightptr%kboti3))=w_ad(i3(lightptr%kboti3))+&
                                              half*lightptr%jac_qgmbi3(lightptr%kboti3)*h1i3_ad
                    w_ad(i3(lightptr%kboti3+1))=w_ad(i3(lightptr%kboti3+1))+&
                                                half*lightptr%jac_qgmbi3(lightptr%kboti3)*h1i3_ad
                 else
                    h1i3_ad=zero
                    rqg(i3(lightptr%kboti3))=zero
                    w_ad(i3(lightptr%kboti3+1))=zero
                 endif

                 if (lightptr%kboti4 > zero) then
                    h1i4_ad=h1i4_ad+(h1i4_ad/abs(h1i4_ad))
                    rqg(i4(lightptr%kboti4))=rqg(i4(lightptr%kboti4))+&
                                             lightptr%jac_qgmai4(lightptr%kboti4)*h1i4_ad
                    w_ad(i4(lightptr%kboti4))=w_ad(i4(lightptr%kboti4))+&
                                              half*lightptr%jac_qgmbi4(lightptr%kboti4)*h1i4_ad
                    w_ad(i4(lightptr%kboti4+1))=w_ad(i4(lightptr%kboti4+1))+&
                                                half*lightptr%jac_qgmbi4(lightptr%kboti4)*h1i4_ad
                 else
                    h1i4_ad=zero
                    rqg(i4(lightptr%kboti4))=zero
                    w_ad(i4(lightptr%kboti4+1))=zero
                 endif


              enddo !do k=nsig-1,1,-1

           endif ! wrf_mass_regional

        endif !if (regional) then

!                .      .    .                                       .
! Global

        if (.not. regional) then

           flashrate_ad=grad
 

           flashratei1_ad=flashratei1_ad+w1*flashrate_ad
           flashratei2_ad=flashratei2_ad+w2*flashrate_ad 
           flashratei3_ad=flashratei3_ad+w3*flashrate_ad 
           flashratei4_ad=flashratei4_ad+w4*flashrate_ad

! Adjoint of maximum vertical velocity 

           wmaxi1_ad=wmaxi1_ad+lightptr%jac_fratei1*flashratei1_ad
           wmaxi2_ad=wmaxi2_ad+lightptr%jac_fratei2*flashratei2_ad
           wmaxi3_ad=wmaxi3_ad+lightptr%jac_fratei3*flashratei3_ad
           wmaxi4_ad=wmaxi3_ad+lightptr%jac_fratei4*flashratei4_ad

           if (lightptr%jac_wmaxflagi1) then
              wmax=-1.e+10_r_kind
              do k=nsig-1,1,-1
                 if (wmaxi1_ad <  zero) then
                    wmaxi1_ad=zero
                 endif
                 if (wmaxi1_ad > wmax) then
                    lightptr%jac_kverti1=k
                    w_ad(i1(lightptr%jac_kverti1))=w_ad(i1(lightptr%jac_kverti1))+wmaxi1_ad
                 endif
              enddo 
           endif
                
           if (lightptr%jac_wmaxflagi2) then
              wmax=-1.e+10_r_kind
              do k=nsig-1,1,-1
                 if (wmaxi2_ad <  zero) then
                    wmaxi2_ad=zero
                 endif
                 if (wmaxi2_ad > wmax) then
                    lightptr%jac_kverti2=k
                    w_ad(i2(lightptr%jac_kverti2))=w_ad(i2(lightptr%jac_kverti2))+wmaxi2_ad
                 endif
              enddo
           endif

           if (lightptr%jac_wmaxflagi3) then
              wmax=-1.e+10_r_kind
              do k=nsig-1,1,-1
                 if (wmaxi3_ad <  zero) then
                    wmaxi3_ad=zero
                 endif
                 if (wmaxi3_ad > wmax) then
                    lightptr%jac_kverti3=k
                    w_ad(i3(lightptr%jac_kverti3))=w_ad(i3(lightptr%jac_kverti3))+wmaxi3_ad
                 endif
              enddo
           endif

           if (lightptr%jac_wmaxflagi4) then
              wmax=-1.e+10_r_kind
              do k=nsig-1,1,-1
                 if (wmaxi4_ad <  zero) then
                    wmaxi4_ad=zero
                 endif
                 if (wmaxi4_ad > wmax) then
                    lightptr%jac_kverti4=k
                    w_ad(i4(lightptr%jac_kverti4))=w_ad(i4(lightptr%jac_kverti4))+wmaxi4_ad
                 endif
              enddo
           endif


        endif !  end global block
!                .      .    .                                       .

! Adjoint of vertical velocity (from vertical and horizontal advection)

        vert_adv_ad(:)=zero
 
        do k=nsig-1,1,-1

           vert_adv_ad(i4(k))=vert_adv_ad(i4(k))+w_ad(i4(k))
           vert_adv_ad(i3(k))=vert_adv_ad(i3(k))+w_ad(i3(k))
           vert_adv_ad(i2(k))=vert_adv_ad(i2(k))+w_ad(i2(k))
           vert_adv_ad(i1(k))=vert_adv_ad(i1(k))+w_ad(i1(k))

        enddo

        horiz_adv_ad(:)=zero

        do k=nsig-1,2,-1

           horiz_adv_ad(i4(k))=horiz_adv_ad(i4(k))+w_ad(i4(k))
           horiz_adv_ad(i4(k))=horiz_adv_ad(i3(k))+w_ad(i3(k))
           horiz_adv_ad(i2(k))=horiz_adv_ad(i2(k))+w_ad(i2(k))
           horiz_adv_ad(i1(k))=horiz_adv_ad(i1(k))+w_ad(i1(k))

        enddo

! Adjoint of q and t from the vertical advection section

        diffq(:)=zero
        difft(:)=zero

        do k=nsig-1,1,-1

           diffq(i1(k))=-(lightptr%jac_ti1(k)*fv*lightptr%jac_vert(K)        &
                        *lightptr%jac_sigdoti1(k))*vert_adv_ad(i1(k))
           difft(i1(k))=-((one+fv*lightptr%jac_qi1(k))*lightptr%jac_vert(k)  &
                        *lightptr%jac_sigdoti1(k))*vert_adv_ad(i1(k))
           diffq(i2(k))=-(lightptr%jac_ti2(k)*fv*lightptr%jac_vert(k)        &
                        *lightptr%jac_sigdoti2(k))*vert_adv_ad(i2(k))
           difft(i2(k))=-((one+fv*lightptr%jac_qi2(k))*lightptr%jac_vert(k)  &
                        *lightptr%jac_sigdoti2(k))*vert_adv_ad(i2(k))
           diffq(i3(k))=-(lightptr%jac_ti3(k)*fv*lightptr%jac_vert(k)        &
                        *lightptr%jac_sigdoti3(k))*vert_adv_ad(i3(k))
           difft(i3(k))=-((one+fv*lightptr%jac_qi3(k))*lightptr%jac_vert(k)  &
                        *lightptr%jac_sigdoti3(k))*vert_adv_ad(i3(k))
           diffq(i4(k))=-(lightptr%jac_ti4(k)*fv*lightptr%jac_vert(k)        &
                        *lightptr%jac_sigdoti4(k))*vert_adv_ad(i4(k))
           difft(i4(k))=-((one+fv*lightptr%jac_qi4(k))*lightptr%jac_vert(k)  &
                        *lightptr%jac_sigdoti4(k))*vert_adv_ad(i4(k))

           rq(i1(k))=rq(i1(k))+diffq(i1(k))

           rt(i1(k))=rt(i1(k))+difft(i1(k))

           rq(i2(k))=rq(i2(k))+diffq(i2(k))

           rq(i3(k))=rq(i3(k))+diffq(i3(k))

           rt(i3(k))=rt(i3(k))+difft(i3(k))

           rq(i4(k))=rq(i4(k))+diffq(i4(k))

           rt(i4(k))=rt(i4(k))+difft(i4(k))

        enddo


 
! Adjoint of z, u, and v from the horizontal advection section
   
        diffz(:)=zero
        z_ad(:)=zero

        do k=nsig-1,2,-1

           diffz(i5(k))=-lightptr%jac_vdyi1(k)*horiz_adv_ad(i1(k))
           diffz(i9(k))=-lightptr%jac_udxi1(k)*horiz_adv_ad(i1(k))

           z_ad(i5(k))=z_ad(i5(k))+diffz(i5(k))
           z_ad(i2(k))=z_ad(i2(k))+(lightptr%jac_vdyi1(k)*horiz_adv_ad(i1(k)))
           z_ad(i9(k))=z_ad(i9(k))+(diffz(i9(k)))
           z_ad(i3(k))=z_ad(i3(k))+(lightptr%jac_udxi1(k)*horiz_adv_ad(i1(k)))
  
           rv(i1(k))=rv(i1(k))+(lightptr%jac_zdyi1(k)*horiz_adv_ad(i1(k)))
           ru(i1(k))=ru(i1(k))+(lightptr%jac_zdxi1(k)*horiz_adv_ad(i1(k)))

           diffz(i1(k)) =-lightptr%jac_vdyi2(k)*horiz_adv_ad(i2(k))
           diffz(i10(k))=-lightptr%jac_udxi2(k)*horiz_adv_ad(i2(k))

           z_ad(i1(k))=z_ad(i1(k))+(diffz(i1(k)))
           z_ad(i6(k))=z_ad(i6(k))+(lightptr%jac_vdyi2(k)*horiz_adv_ad(i2(k)))
           z_ad(i10(k))=z_ad(i10(k))+(diffz(i10(k)))
           z_ad(i4(k))=z_ad(i4(k))+(lightptr%jac_udxi2(k)*horiz_adv_ad(i2(k)))
           rv(i2(k))=rv(i2(k))+(lightptr%jac_zdyi2(k)*horiz_adv_ad(i2(k)))
           ru(i2(k))=ru(i2(k))+(lightptr%jac_zdxi2(k)*horiz_adv_ad(i2(k)))

           diffz(i7(k))= -lightptr%jac_vdyi3(k)*horiz_adv_ad(i3(k))
           diffz(i1(k))= -lightptr%jac_udxi3(k)*horiz_adv_ad(i3(k))

           z_ad(i7(k)) =  z_ad(i7(k))+diffz(i7(k))
           z_ad(i4(k)) =  z_ad(i4(k))+(lightptr%jac_vdyi3(k)*horiz_adv_ad(i3(k)))
           z_ad(i1(k)) =  z_ad(i1(k))+diffz(i1(k))
           z_ad(i11(k))=  z_ad(i11(k))+(lightptr%jac_udxi3(k)*horiz_adv_ad(i3(k)))
           rv(i3(k)) =  rv(i3(k))+(lightptr%jac_zdyi3(k)*horiz_adv_ad(i3(k)))
           ru(i3(k)) =  ru(i3(k))+(lightptr%jac_zdxi3(k)*horiz_adv_ad(i3(k)))

           diffz(i3(k))=-lightptr%jac_vdyi4(k)*horiz_adv_ad(i4(k))
           diffz(i2(k))=-z_tl(i2(k))-lightptr%jac_udxi4(k)*horiz_adv_ad(i4(k))

           z_ad(i3(k)) =  z_ad(i3(k))+diffz(i3(k))
           z_ad(i8(k)) =  z_ad(i8(k))+(lightptr%jac_vdyi4(k)*horiz_adv_ad(i4(k)))
           z_ad(i2(k)) = z_tl(i2(k))+diffz(i2(k))
           z_ad(i12(k))= z_ad(i12(k))+(lightptr%jac_udxi4(k)*horiz_adv_ad(i4(k)))
           rv(i4(k)) = rv(i4(k))+(lightptr%jac_zdyi4(k)*horiz_adv_ad(i4(k)))
           ru(i4(k)) = ru(i4(k))+(lightptr%jac_zdxi4(k)*horiz_adv_ad(i4(k)))

        enddo

! Adjoint of q and t from the calculation of height (z)
  
        do k=nsig-1,2,-1

           rq(i1(k))=rq(i1(k))+lightptr%jac_vertqi1(k)*z_ad(i1(k))
           rt(i1(k))=rt(i1(k))+lightptr%jac_vertti1(k)*z_ad(i1(k)) 
           z_ad(i1(k-1))=z_ad(i1(k-1))+z_ad(i1(k))
           z_ad(i1(k))=zero

           rq(i2(k))=rq(i2(k))+lightptr%jac_vertqi2(k)*z_ad(i2(k)) 
           rt(i2(k))=rt(i2(k))+lightptr%jac_vertti12(k)*z_ad(i2(k))
           z_ad(i2(k-1))=z_ad(i2(k-1))+z_ad(i2(k))
           z_ad(i2(k))=zero

           rq(i3(k))=rq(i3(k))+lightptr%jac_vertqi3(k)*z_ad(i3(k))
           rt(i3(k))=rt(i3(k))+lightptr%jac_vertti3(k)*z_ad(i3(k))
           z_ad(i3(k-1))=z_ad(i3(k-1))+z_ad(i3(k))
           z_ad(i3(k))=zero

           rq(i4(k))=rq(i4(k))+lightptr%jac_vertqi4(k)*z_ad(i4(k))
           rt(i4(k))=rt(i4(k))+lightptr%jac_vertti4(k)*z_ad(i4(k))
           z_ad(i4(k-1))=z_ad(i4(k-1))+z_ad(i4(k))
           z_ad(i4(k))=zero

           rq(i5(k))=rq(i5(k))+lightptr%jac_vertqi5(k)*z_ad(i5(k))
           rt(i5(k))=rt(i5(k))+lightptr%jac_vertti5(k)*z_ad(i5(k))
           z_ad(i5(k-1))=z_ad(i5(k-1))+z_ad(i5(k))
           z_ad(i5(k))=zero

           rq(i6(k))=rq(i6(k))+lightptr%jac_vertqi6(k)*z_ad(i6(k))
           rt(i6(k))=rt(i6(k))+lightptr%jac_vertti6(k)*z_ad(i6(k))
           z_ad(i6(k-1))=z_ad(i6(k-1))+z_ad(i6(k))
           z_ad(i6(k))=zero

           rq(i7(k))=rq(i7(k))+lightptr%jac_vertqi7(k)*z_ad(i7(k))
           rt(i7(k))=rt(i7(k))+lightptr%jac_vertti7(k)*z_ad(i7(k))
           z_ad(i7(k-1))=z_ad(i7(k-1))+z_ad(i7(k))
           z_ad(i7(k))=zero

           rq(i8(k))=rq(i8(k))+lightptr%jac_vertqi8(k)*z_ad(i8(k))
           rt(i8(k))=rt(i8(k))+lightptr%jac_vertti8(k)*z_ad(i8(k))
           z_ad(i8(k-1))=z_ad(i8(k-1))+z_ad(i8(k))
           z_ad(i8(k))=zero

           rq(i9(k))=rq(i9(k))+lightptr%jac_vertqi9(k)*z_ad(i9(k))
           rt(i9(k))=rt(i9(k))+lightptr%jac_vertti9(k)*z_ad(i9(k))
           z_ad(i9(k-1))=z_ad(i9(k-1))+z_ad(i9(k))
           z_ad(i9(k))=zero

           rq(i10(k))=rq(i10(k))+lightptr%jac_vertqi10(k)*z_ad(i10(k))
           rt(i10(k))=rt(i10(k))+lightptr%jac_vertti10(k)*z_ad(i10(k))
           z_ad(i10(k-1))=z_ad(i10(k-1))+z_ad(i10(k))
           z_ad(i10(k))=zero

           rq(i11(k))=rq(i11(k))+lightptr%jac_vertqi11(k)*z_ad(i11(k))   
           rt(i11(k))=rt(i11(k))+lightptr%jac_vertti11(k)*z_ad(i11(k))
           z_ad(i11(k-1))=z_ad(i11(k-1))+z_ad(i11(k))
           z_ad(i11(k))=zero

           rq(i12(k))=rq(i12(k))+lightptr%jac_vertqi12(k)*z_ad(i12(k))
           rt(i12(k))=rt(i12(k))+lightptr%jac_vertti12(k)*z_ad(i12(k))
           z_ad(i12(k-1))=z_ad(i12(k-1))+z_ad(i12(k))
           z_ad(i12(k))=zero

        enddo
 

     endif !Adjoint

     lightptr => lightnode_nextcast(lightptr)

  enddo  ! do while (associated(lightptr))
    
  
  return
end subroutine intlight_

end module intlightmod
