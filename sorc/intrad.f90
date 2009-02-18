subroutine intrad(rt,rq,roz,ru,rv,rst,st,sq,soz,su,sv,sst,rpred,spred,   &
                  drt,drq,droz,dru,drv,dst,dsq,dsoz,dsu,dsv)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    intrad      sat radiance nonlin qc obs operator
!   prgmmr: parrish          org: np22                date: 1990-10-11
!
! abstract: apply satellite radiance operator and adjoint with
!            addition of nonlinear qc operator.
!
! program history log:
!   1990-10-11  parrish
!   1992-07-21
!   1995-07-17  derber
!   1997-03-10  wu
!   1997-12-22  weiyu yang
!   1999-08-24  derber, j., treadon, r., yang, w., first frozen mpp version
!   2004-08-02  treadon - add only to module use, add intent in/out
!   2004-10-07  parrish - add nonlinear qc option
!   2005-01-20  okamoto - add wind components
!   2005-04-11  treadon - merge intrad and intrad_qc into single routine
!   2005-09-28  derber  - modify var qc and change location and weight arrays
!   2006-04-03  derber  - clean up code
!   2006-07-28  derber  - modify to use new inner loop obs data structure
!                       - unify NL qc
!   2007-06-04  derber  - use quad precision to get reproducability over number of processors
!
!   input argument list:
!     st       - input temperature correction field
!     sq       - input q correction field
!     soz      - input ozone correction field
!     su       - input u correction field
!     sv       - input v correction field
!     spred    - input predictor values
!     sst      - input skin temp. vector
!     dst      - input time derivative of temperature correction field
!     dsq      - input time derivative of q correction field
!     dsoz     - input time derivative of ozone correction field
!     dsu      - input time derivative of u correction field
!     dsv      - input time derivative of v correction field
!
!   output argument list:
!     rt       - output t vector after inclusion of radiance info.
!     rq       - output q vector after inclusion of radiance info.
!     ro       - output ozone vector after inclusion of radiance info.
!     ru       - output u vector after inclusion of radiance info.
!     rv       - output v vector after inclusion of radiance info.
!     rpred    - output predictor vector after inclusion of radiance info.
!     rst      - output skin temp. vector after inclusion of radiance info.
!     drt      - output time derivative of t vector after inclusion of radiance info.
!     drq      - output time derivative of q vector after inclusion of radiance info.
!     dro      - output time derivative of ozone vector after inclusion of radiance info.
!     dru      - output time derivative of u vector after inclusion of radiance info.
!     drv      - output time derivative of v vector after inclusion of radiance info.
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind,r_quad
  use radinfo, only: npred,npred1,jpch_rad,pg_rad,b_rad
  use obsmod, only: radptr,radhead
  use gridmod, only: latlon11,latlon1n,nsig,nsig2,nsig3,nsig4,&
       nsig3p1,nsig3p2,nsig3p3
  use qcmod, only: nlnqc_iter
  use constants, only: zero,half,one,tiny_r_kind,cg_term
  implicit none

! Declare passed variables
  real(r_kind),dimension(latlon1n),intent(in):: st,sq,soz,su,sv
  real(r_kind),dimension(latlon1n),intent(in):: dst,dsq,dsoz,dsu,dsv
  real(r_kind),dimension(latlon11),intent(in):: sst
  real(r_kind),dimension(jpch_rad,npred),intent(in):: spred
  real(r_kind),dimension(latlon1n),intent(inout):: rt,rq,roz,ru,rv
  real(r_kind),dimension(latlon1n),intent(inout):: drt,drq,droz,dru,drv
  real(r_kind),dimension(latlon11),intent(inout):: rst
  real(r_quad),dimension(jpch_rad,npred),intent(inout):: rpred

! Declare local variables
  integer(i_kind) i,j1,j2,j3,j4,i1,i2,i3,i4,n,n_1,n_2,k,ic,nn
  real(r_kind) valx,tlap,val,tlap2,w1,w2,w3,w4
! real(r_kind) penalty,p1
  real(r_kind),dimension(nsig3p3):: tval,tdir
  real(r_kind) cg_rad,p0,wnotgross,wgross,time_rad,y,temp
  integer(i_kind),dimension(nsig) :: i1n,i2n,i3n,i4n

  radptr => radhead
  do while (associated(radptr))
     j1=radptr%ij(1)
     j2=radptr%ij(2)
     j3=radptr%ij(3)
     j4=radptr%ij(4)
     w1=radptr%wij(1)
     w2=radptr%wij(2)
     w3=radptr%wij(3)
     w4=radptr%wij(4)

!    p1=zero
     do k=1,nsig3p3
        tval(k)=zero
     end do

     time_rad=radptr%time
!  Begin Forward model
!  calculate temperature, q, ozone, sst vector at observation location
     i1n(1) = j1
     i2n(1) = j2
     i3n(1) = j3
     i4n(1) = j4
     do k=2,nsig
      i1n(k) = i1n(k-1)+latlon11
      i2n(k) = i2n(k-1)+latlon11
      i3n(k) = i3n(k-1)+latlon11
      i4n(k) = i4n(k-1)+latlon11
     enddo
!$omp parallel do private(k,i1,i2,i3,i4)
    do k=1,nsig
        i1 = i1n(k)
        i2 = i2n(k)
        i3 = i3n(k)
        i4 = i4n(k)
        tdir(k)=      w1*  st(i1)+w2*  st(i2)+ &
                      w3*  st(i3)+w4*  st(i4)+ &
                     (w1* dst(i1)+w2* dst(i2)+ &
                      w3* dst(i3)+w4* dst(i4))*time_rad
        tdir(nsig+k)= w1*  sq(i1)+w2*  sq(i2)+ &
                      w3*  sq(i3)+w4*  sq(i4)+ &
                     (w1* dsq(i1)+w2* dsq(i2)+ &
                      w3* dsq(i3)+w4* dsq(i4))*time_rad
        tdir(nsig2+k)=w1* soz(i1)+w2* soz(i2)+ &
                      w3* soz(i3)+w4* soz(i4)+ &
                     (w1*dsoz(i1)+w2*dsoz(i2)+ &
                      w3*dsoz(i3)+w4*dsoz(i4))*time_rad
     end do
!$omp end parallel do
     tdir(nsig3p1)=   w1* su(j1) +w2* su(j2)+ &
                      w3* su(j3) +w4* su(j4)+ &
                     (w1*dsu(j1) +w2*dsu(j2)+ &
                      w3*dsu(j3) +w4*dsu(j4))*time_rad
     tdir(nsig3p2)=   w1* sv(j1) +w2* sv(j2)+ &
                      w3* sv(j3) +w4* sv(j4)+ &
                     (w1*dsv(j1) +w2*dsv(j2)+ &
                      w3*dsv(j3) +w4*dsv(j4))*time_rad
     tdir(nsig3p3)=   w1*sst(j1) +w2*sst(j2)+ &
                      w3*sst(j3) +w4*sst(j4)
!  begin channel specific calculations
     do nn=1,radptr%nchan
        ic=radptr%icx(nn)

!       include observation increment and lapse rate contributions to bias correction
        tlap=radptr%pred2(nn)
        tlap2=tlap*tlap
        valx=-radptr%res(nn)+spred(ic,npred)*tlap+spred(ic,npred1)*tlap2

!       Include contributions from atmospheric jacobian
        do k=1,nsig3p3
           valx=valx+tdir(k)*radptr%dtb_dvar(k,nn)
        end do

!       Include contributions from remaining bias correction terms
        do n=1,npred-2
           valx=valx+spred(ic,n)*radptr%pred1(n)
        end do

!       Multiply by variance.
        if (nlnqc_iter .and. pg_rad(ic) > tiny_r_kind .and. &
                             b_rad(ic)  > tiny_r_kind) then
           cg_rad=cg_term/b_rad(ic)
           wnotgross= one-pg_rad(ic)
           wgross = pg_rad(ic)*cg_rad/wnotgross
           p0   = wgross/(wgross+exp(-half*radptr%err2(nn)*valx*valx))
           valx = valx*(one-p0)
        endif

        val      = valx*radptr%err2(nn)
        val      = val*radptr%raterr2(nn)

!       Begin adjoint

!       Extract contributions from atmospheric jacobian
        do k=1,nsig3p3
           tval(k)=tval(k)+radptr%dtb_dvar(k,nn)*val
        end do

!       Extract contributions from bias correction terms
!       use compensated summation
        if(radptr%luse)then
          do n=1,npred-2
             rpred(ic,n)=rpred(ic,n)+radptr%pred1(n)*val
          end do
          rpred(ic,npred) =rpred(ic,npred) +val*tlap
          rpred(ic,npred1)=rpred(ic,npred1)+val*tlap2
        end if
     end do


!    Distribute adjoint contributions over surrounding grid points
     do k=1,nsig
        n_1=k+nsig
        n_2=k+nsig2
        i1 = i1n(k)
        i2 = i2n(k)
        i3 = i3n(k)
        i4 = i4n(k)

        rt(i1)=rt(i1)+w1*tval(k)
        rt(i2)=rt(i2)+w2*tval(k)
        rt(i3)=rt(i3)+w3*tval(k)
        rt(i4)=rt(i4)+w4*tval(k)
        drt(i1)=drt(i1)+w1*tval(k)*time_rad
        drt(i2)=drt(i2)+w2*tval(k)*time_rad
        drt(i3)=drt(i3)+w3*tval(k)*time_rad
        drt(i4)=drt(i4)+w4*tval(k)*time_rad

        rq(i1)=rq(i1)+w1*tval(n_1)
        rq(i2)=rq(i2)+w2*tval(n_1)
        rq(i3)=rq(i3)+w3*tval(n_1)
        rq(i4)=rq(i4)+w4*tval(n_1)
        drq(i1)=drq(i1)+w1*tval(n_1)*time_rad
        drq(i2)=drq(i2)+w2*tval(n_1)*time_rad
        drq(i3)=drq(i3)+w3*tval(n_1)*time_rad
        drq(i4)=drq(i4)+w4*tval(n_1)*time_rad

        roz(i1)=roz(i1)+w1*tval(n_2)
        roz(i2)=roz(i2)+w2*tval(n_2)
        roz(i3)=roz(i3)+w3*tval(n_2)
        roz(i4)=roz(i4)+w4*tval(n_2)
        droz(i1)=droz(i1)+w1*tval(n_2)*time_rad
        droz(i2)=droz(i2)+w2*tval(n_2)*time_rad
        droz(i3)=droz(i3)+w3*tval(n_2)*time_rad
        droz(i4)=droz(i4)+w4*tval(n_2)*time_rad

     end do

     ru(j1)=ru(j1)+w1*tval(nsig3p1)
     ru(j2)=ru(j2)+w2*tval(nsig3p1)
     ru(j3)=ru(j3)+w3*tval(nsig3p1)
     ru(j4)=ru(j4)+w4*tval(nsig3p1)
     dru(j1)=dru(j1)+w1*tval(nsig3p1)*time_rad
     dru(j2)=dru(j2)+w2*tval(nsig3p1)*time_rad
     dru(j3)=dru(j3)+w3*tval(nsig3p1)*time_rad
     dru(j4)=dru(j4)+w4*tval(nsig3p1)*time_rad
 
     rv(j1)=rv(j1)+w1*tval(nsig3p2)
     rv(j2)=rv(j2)+w2*tval(nsig3p2)
     rv(j3)=rv(j3)+w3*tval(nsig3p2)
     rv(j4)=rv(j4)+w4*tval(nsig3p2)
     drv(j1)=drv(j1)+w1*tval(nsig3p2)*time_rad
     drv(j2)=drv(j2)+w2*tval(nsig3p2)*time_rad
     drv(j3)=drv(j3)+w3*tval(nsig3p2)*time_rad
     drv(j4)=drv(j4)+w4*tval(nsig3p2)*time_rad

     rst(j1)=rst(j1)+w1*tval(nsig3p3)
     rst(j2)=rst(j2)+w2*tval(nsig3p3)
     rst(j3)=rst(j3)+w3*tval(nsig3p3)
     rst(j4)=rst(j4)+w4*tval(nsig3p3)

     radptr => radptr%llpoint
  end do
  return
end subroutine intrad
