module intjomod
!$$$ module documentation block
!           .      .    .                                       .
! module:   intjo    module for intjo
!   prgmmr:
!
! abstract: module for h'r^{-1}h
!
! program history log:
!   2008-12-01  Todling - wrap in module
!   2009-08-13  lueken - update documentation
!   2015-09-03  guo     - obsmod::obs_handle has been replaced with obsheadbundle,
!                         defined by m_obsheadbundle.
!   2016-08-29  J Guo   - Separated calls to intozlay() and intozlev()
!   2018-08-10  guo     - a new generic intjo() implementation replaced type
!                         specific intxyz() calls with polymorphic %intjo().
!
! subroutines included:
!   sub intjo_
!
! variable definitions:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

use gsi_obopertypemanager, only: oboper_count
use gsi_obopertypemanager, only: oboper_typeinfo
use gsi_obopertypemanager, only: &
  ioboper_t,          ioboper_pw,         ioboper_q,                                                    &
                                          ioboper_cldtot,     ioboper_w,          ioboper_dw,           &
  ioboper_rw,         ioboper_dbz,                                                                      &
                      ioboper_spd,        ioboper_oz,         ioboper_o3l,        ioboper_colvk,        &
  ioboper_pm2_5,      ioboper_pm10,       ioboper_ps,         ioboper_tcp,        ioboper_sst,          &
  ioboper_gpsbend,    ioboper_gpsref,                                                                   &
                      ioboper_rad,        ioboper_pcp,        ioboper_aero,       ioboper_gust,         &
  ioboper_vis,        ioboper_pblh,       ioboper_wspd10m,    ioboper_td2m,       ioboper_mxtm,         &
  ioboper_mitm,       ioboper_pmsl,       ioboper_howv,       ioboper_tcamt,      ioboper_lcbas,        &
  ioboper_cldch,      ioboper_uwnd10m,    ioboper_vwnd10m,    ioboper_swcp,       ioboper_lwcp,         &
  ioboper_light
use kinds, only: i_kind
use mpeu_util, only: perr,die

implicit none

private
public intjo

interface intjo; module procedure &
   intjo_, intjo_reduced_
end interface

! This is a mapping to the exact %intjo() calling sequence in the earlier
! non-polymorphic implementation, to reproduce the exactly same summation
! ordering for rval and qpred.  It is not necessary, and can be removed once
! some non-zero-diff modifications are introduced.

integer(i_kind),parameter,dimension(oboper_count):: ix_obtype = (/ &
  ioboper_t,          ioboper_pw,         ioboper_q,                                                    &
                                          ioboper_cldtot,     ioboper_w,          ioboper_dw,           &
  ioboper_rw,         ioboper_dbz,                                                                      &
                      ioboper_spd,        ioboper_oz,         ioboper_o3l,        ioboper_colvk,        &
  ioboper_pm2_5,      ioboper_pm10,       ioboper_ps,         ioboper_tcp,        ioboper_sst,          &
  ioboper_gpsbend,    ioboper_gpsref,                                                                   &
                      ioboper_rad,        ioboper_pcp,        ioboper_aero,       ioboper_gust,         &
  ioboper_vis,        ioboper_pblh,       ioboper_wspd10m,    ioboper_td2m,       ioboper_mxtm,         &
  ioboper_mitm,       ioboper_pmsl,       ioboper_howv,       ioboper_tcamt,      ioboper_lcbas,        &
  ioboper_cldch,      ioboper_uwnd10m,    ioboper_vwnd10m,    ioboper_swcp,       ioboper_lwcp,         &
  ioboper_light                                                                                         /)
!...|....1....|....2....|....3....|....4....|....5....|....6....|....7....|....8....|....9....|....0

character(len=*),parameter:: myname="intjomod"

contains

subroutine intjo_(rval,qpred,sval,sbias)

!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    intjo      calculate rhs for analysis equation
!   prgmmr: derber           org: np23                date: 2003-12-18
!
! abstract: calculate rhs for all variables (nonlinear qc version)
!
!    A description of nonlinear qc follows:
!
!    The observation penalty jo is defined as
!
!          jo =  - (sum over obs) 2*log(po)
!
!      where,
!
!          po = wnotgross*exp(-.5*(hn(x+xb) - yo)**2 ) + wgross
!            with
!                hn = the forward model (possibly non-linear) normalized by 
!                     observation error
!                x  = the current estimate of the analysis increment
!                xb = the background state
!                yo = the observation normalized by observation error
!
!            Note:  The factor 2 in definition of jo is present because the 
!                   penalty jo as used in this code is 2*(usual definition 
!                   of penalty)
!
!          wgross = pgross*cg
!
!          wnotgross = 1 - wgross
!
!          pgross = probability of gross error for observation (assumed
!                   here to have uniform distribution over the possible
!                   range of values)
!
!          cg = sqrt(2*pi)/2b
!
!          b = possible range of variable for gross errors, normalized by 
!              observation error
!
!    The values for the above parameters that bill collins used in the
!    eta 3dvar are:
!
!          cg = cg_term/b, where cg_term = sqrt(2*pi)/2 
!
!          b = 10.        ! range for gross errors, normalized by obs error
!
!          pg_q=.002      ! probability of gross error for specific humidity
!          pg_pw=.002     ! probability of gross error for precipitable water
!          pg_p=.002      ! probability of gross error for pressure
!          pg_w=.005      ! probability of gross error for wind
!          pg_t=.007      ! probability of gross error for temperature
!          pg_rad=.002    ! probability of gross error for radiances
!
!
!    Given the above jo, the gradient of jo is as follows:
!
!                                             t
!        gradx(jo) = - (sum over observations) 2*h (hn(x+xb)-yo)*(po - wgross)/po
!
!      where, 
!
!          h = tangent linear model of hn about x+xb
!
! 
!    Note that if pgross = 0.0, then wnotgross=1.0 and wgross=0.0.  That is,
!    the code runs as though nonlinear quality control were not present
!    (which is indeed the case since the gross error probability is 0).  
!
!    As a result the same int* routines may be used for use with or without
!    nonlinear quality control.
!    
!
! program history log:
!   2003-12-18  derber
!   2004-07-23  derber  - modify to include conventional sst
!   2004-07-28  treadon - add only to module use, add intent in/out
!   2004-10-06  parrish - add nonlinear qc option
!   2004-10-06  kleist  - separate control vector for u,v, & convert int
!                         for wind components into int for st,vp
!   2004-11-30  treadon - add brightness temperatures to nonlinear 
!                         quality control
!   2004-12-03  treadon - replace mpe_iallreduce (ibm extension) with
!                         standard mpi_allreduce
!   2005-01-20  okamoto - add u,v to intrad
!   2005-02-23  wu      - changes related to normalized rh option
!   2005-04-11  treadon - rename intall_qc as intall
!   2005-05-18  yanqiu zhu - add 'use int*mod',and modify call interfaces for using these modules
!   2005-05-24  pondeca - take into consideration that npred=npredp=0
!                         for 2dvar only surface analysis option
!   2005-06-03  parrish - add horizontal derivatives
!   2005-07-10  kleist  - add dynamic constraint term
!   2005-09-29  kleist  - expand jc term, include time derivatives vector
!   2005-11-21  kleist  - separate tendencies from jc term, add call to calctends adjoint
!   2005-12-01  cucurull - add code for gps local bending angle, add use obsmod for ref_obs
!   2005-12-20  parrish - add arguments to call to intt to allow for option of using boundary
!                         layer forward tlm.
!   2006-02-03  derber  - modify to increase reproducibility
!   2006-03-17  park    - correct error in call to intt--rval,sval --> rvaluv,svaluv
!                          in order to correctly pass wind variables.
!   2006-04-06  kleist  - include both jc formulations
!   2006-07-26  parrish - correct inconsistency in computation of space and time derivatives of q
!                          currently, if derivatives computed, for q it is normalized q, but
!                          should be mixing ratio.
!   2006-07-26  parrish - add strong constraint initialization option
!   2007-03-19  tremolet - binning of observations
!   2007-04-13  tremolet - split jo from other components of intall
!   2007-06-04  derber  - use quad precision to get reproducibility over number of processors
!   2008-11-27  todling  - add tendencies for foto support and new interface to int's
!   2009-01-08  todling  - remove reference to ozohead
!   2009-03-23  meunier  - Add call to intlag (lagrangian observations)
!   2010-01-11  zhang,b  - Bug fix: bias predictors need to be accumulated over nbins
!   2010-03-24  zhu      - change the interfaces of intt,intrad,intpcp for generalizing control variable
!   2010-05-13  todling  - harmonized interfaces to int* routines when it comes to state_vector (add only's)
!   2010-06-13  todling  - add intco call
!   2010-10-15  pagowski - add intpm2_5 call
!   2010-10-20  hclin    - added aod
!   2011-02-20  zhu      - add intgust,intvis,intpblh calls
!   2013-05-20  zhu      - add codes related to aircraft temperature bias correction
!   2014-06-18  carley/zhu - add lcbas and tcamt 
!   2014-03-19  pondeca  - add intwspd10m
!   2014-04-10  pondeca  - add inttd2m,intmxtm,intmitm,intpmsl
!   2014-05-07  pondeca  - add inthowv
!   2015-07-10  pondeca  - add intcldch
!   2016-03-07  pondeca  - add intuwnd10m,intvwnd10m
!
!   input argument list:
!     ibin
!     yobs
!     sval     - solution on grid
!     sbias
!     rval
!     qpred
!
!   output argument list:      
!     rval     - rhs on grid
!     qpred
!
! remarks:
!     1) if strong initialization, then svalt, svalp, svaluv
!         are all grid fields after strong initialization.
!
!     2) The two interfaces to the int-routines should be temporary.
!        In the framework of the 4dvar-code, foto can be re-implemented as 
!        an approximate m and m' to the model matrices in 4dvar. Once that
!        is done, the int-routines should no longer need the time derivatives.
!        (todling)
!     3) Notice that now (2010-05-13) int routines handle non-essential
!        variables internally; also, when pointers non-existent, int routines 
!        simply return (todling).
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: i_kind,r_quad
  use gsi_bundlemod, only: gsi_bundle
  use bias_predictors, only: predictors
  use m_obsdiags, only: oboper_create
  use m_obsdiags, only: oboper_destroy
  use gsi_oboper, only: oboper

  use intradmod, only: setrad

  implicit none

! Declare passed variables
  type(gsi_bundle), dimension(  :), intent(inout) :: rval         ! (nobs_bins)
  type(gsi_bundle), dimension(  :), intent(in   ) :: sval         ! (nobs_bins)
  real(r_quad    ), dimension(:,:), intent(inout) :: qpred        ! (:,nobs_bins)
  type(predictors),                 intent(in   ) :: sbias

  character(len=*),parameter:: myname_=myname//"::intjo_"

! Declare local variables
  integer(i_kind):: ibin,it,ix
  class(oboper),pointer:: it_oboper

!******************************************************************************
  call setrad(sval(1))

! "rhs for jo", as it was labeled in intall().
!$omp parallel do  schedule(dynamic,1) private(ibin,it,ix,it_oboper)
  do ibin=1,size(sval)
     do it=1,oboper_count
        !ix=ix_obtype(it)  ! Use this line to ensure the same jo summartion
                           ! sequence as intjo was in its early implementation,
                           ! for reproducibility.

        ix=it     ! Using this line, jo summation sequence is not the same as
                  ! it used to be, nor the same if someone chooses to change
                  ! enumration sequence of obopers in gsi_obopertypemanager.F90.
                  ! But it would make this code more portable to new oboper
                  ! extensions.

        it_oboper => oboper_create(ix)

        if(.not.associated(it_oboper)) then
           call perr(myname_,'unexpected obOper, associated(it_obOper) =',associated(it_oboper))
           call perr(myname_,'                  obOper_typeInfo(ioper) =',oboper_typeinfo(ix))
           call perr(myname_,'                                   ioper =',ix)
           call perr(myname_,'                                      it =',it)
           call perr(myname_,'                            obOper_count =',oboper_count)
           call perr(myname_,'                                    ibin =',ibin)
           call  die(myname_)
        endif

        if(.not.associated(it_oboper%obsll)) then
           call perr(myname_,'unexpected component, associated(%obsLL) =',associated(it_oboper%obsll))
           call perr(myname_,'                  obOper_typeInfo(ioper) =',oboper_typeinfo(ix))
           call perr(myname_,'                                   ioper =',ix)
           call perr(myname_,'                                      it =',it)
           call perr(myname_,'                            obOper_count =',oboper_count)
           call perr(myname_,'                                    ibin =',ibin)
           call  die(myname_)
        endif

        call it_oboper%intjo(ibin,rval(ibin),sval(ibin),qpred(:,ibin),sbias)
        call oboper_destroy(it_oboper)
     enddo
  end do

  return
end subroutine intjo_

subroutine intjo_reduced_(rval,qpred,sval,sbias)
  use kinds, only: i_kind,r_quad
  use gsi_bundlemod, only: gsi_bundle
  use bias_predictors, only: predictors
  use constants, only: zero_quad
  implicit none
! Declare passed variables
  type(gsi_bundle), dimension(:), intent(inout) :: rval
  real(r_quad    ), dimension(:), intent(inout) :: qpred

  type(gsi_bundle), dimension(:), intent(in   ) :: sval
  type(predictors),               intent(in   ) :: sbias

!----------------------------------------
  real(kind(qpred)),allocatable,dimension(:,:):: qpred_bin
  integer(i_kind):: ibin

  allocate(qpred_bin(size(qpred),size(rval)))
  qpred_bin=zero_quad

  call intjo_(rval,qpred_bin,sval,sbias)

  do ibin=1,size(rval)
     qpred(:)=qpred(:)+qpred_bin(:,ibin)
  enddo

  deallocate(qpred_bin)

end subroutine intjo_reduced_

end module intjomod
