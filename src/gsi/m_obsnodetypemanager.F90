module m_obsnodetypemanager
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module m_obsnodetypemanager
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2015-08-13
!
! abstract: obsnode type manager, as an enumerated type molder.
!
! program history log:
!   2015-08-13  j guo   - added this document block.
!   2016-05-18  j guo   - finished its initial polymorphic implementation,
!                         with total 33 obs-types.
!   2018-01-23  k apodaca - add a new observation type i.e. lightning (light)
!                           suitable for the goes/glm instrument
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

  use m_psnode   , only:    psnode
  use m_tnode    , only:     tnode
  use m_wnode    , only:     wnode
  use m_qnode    , only:     qnode
  use m_spdnode  , only:   spdnode
  use m_rwnode   , only:    rwnode
  use m_dwnode   , only:    dwnode
  use m_sstnode  , only:   sstnode
  use m_pwnode   , only:    pwnode
  use m_pcpnode  , only:   pcpnode
  use m_oznode   , only:    oznode
  use m_o3lnode  , only:   o3lnode
  use m_gpsnode  , only:   gpsnode
  use m_radnode  , only:   radnode
  use m_tcpnode  , only:   tcpnode
  use m_lagnode  , only:   lagnode
  use m_colvknode, only: colvknode
  use m_aeronode , only:  aeronode
  use m_aerolnode, only: aerolnode
  use m_pm2_5node, only: pm2_5node
  use m_gustnode , only:  gustnode
  use m_visnode  , only:   visnode
  use m_pblhnode , only:  pblhnode

  use m_wspd10mnode, only: wspd10mnode
  use m_uwnd10mnode, only: uwnd10mnode
  use m_vwnd10mnode, only: vwnd10mnode

  use m_td2mnode , only:  td2mnode
  use m_mxtmnode , only:  mxtmnode
  use m_mitmnode , only:  mitmnode
  use m_pmslnode , only:  pmslnode
  use m_howvnode , only:  howvnode
  use m_tcamtnode, only: tcamtnode
  use m_lcbasnode, only: lcbasnode
  use m_pm10node , only:  pm10node
  use m_cldchnode, only: cldchnode

  use m_swcpnode , only:  swcpnode
  use m_lwcpnode , only:  lwcpnode

  use m_lightnode, only: lightnode
  use m_dbznode  , only:   dbznode

  use kinds, only: i_kind
  use m_obsnode, only: obsnode
  use mpeu_util, only: perr,die

  implicit none
  private       ! except

  public:: obsnodetype_undef
  public:: obsnodetype_lbound
  public:: obsnodetype_ubound
  public:: obsnodetype_count

  public:: iobsnode_kind
  public:: iobsnode_ps
  public:: iobsnode_t
  public:: iobsnode_w
  public:: iobsnode_q
  public:: iobsnode_spd
  public:: iobsnode_rw
  public:: iobsnode_dw
  public:: iobsnode_sst
  public:: iobsnode_pw
  public:: iobsnode_pcp
  public:: iobsnode_oz
  public:: iobsnode_o3l
  public:: iobsnode_gps
  public:: iobsnode_rad
  public:: iobsnode_tcp
  public:: iobsnode_lag
  public:: iobsnode_colvk
  public:: iobsnode_aero
  public:: iobsnode_aerol
  public:: iobsnode_pm2_5
  public:: iobsnode_gust
  public:: iobsnode_vis
  public:: iobsnode_pblh
  public:: iobsnode_wspd10m
  public:: iobsnode_uwnd10m
  public:: iobsnode_vwnd10m
  public:: iobsnode_td2m
  public:: iobsnode_mxtm
  public:: iobsnode_mitm
  public:: iobsnode_pmsl
  public:: iobsnode_howv
  public:: iobsnode_tcamt
  public:: iobsnode_lcbas
  public:: iobsnode_pm10
  public:: iobsnode_cldch
  public:: iobsnode_swcp
  public:: iobsnode_lwcp

  public:: iobsnode_light
  public:: iobsnode_dbz

  public :: obsnode_typemold
  public :: obsnode_typeindex

  interface obsnode_typemold; module procedure &
     index2vmold_, &
     vname2vmold_
  end interface
  interface obsnode_typeindex; module procedure &
     vmold2index_, &
     vname2index_
  end interface

  type(psnode   ), target, save::    ps_mold
  type(tnode    ), target, save::     t_mold
  type(wnode    ), target, save::     w_mold
  type(qnode    ), target, save::     q_mold
  type(spdnode  ), target, save::   spd_mold
  type(rwnode   ), target, save::    rw_mold
  type(dwnode   ), target, save::    dw_mold
  type(sstnode  ), target, save::   sst_mold
  type(pwnode   ), target, save::    pw_mold
  type(pcpnode  ), target, save::   pcp_mold
  type(oznode   ), target, save::    oz_mold
  type(o3lnode  ), target, save::   o3l_mold
  type(gpsnode  ), target, save::   gps_mold
  type(radnode  ), target, save::   rad_mold
  type(tcpnode  ), target, save::   tcp_mold
  type(lagnode  ), target, save::   lag_mold
  type(colvknode), target, save:: colvk_mold
  type(aeronode ), target, save::  aero_mold
  type(aerolnode), target, save:: aerol_mold
  type(pm2_5node), target, save:: pm2_5_mold
  type(gustnode ), target, save::  gust_mold
  type(visnode  ), target, save::   vis_mold
  type(pblhnode ), target, save::  pblh_mold

  type(wspd10mnode), target, save:: wspd10m_mold
  type(uwnd10mnode), target, save:: uwnd10m_mold
  type(vwnd10mnode), target, save:: vwnd10m_mold

  type(   td2mnode), target, save::    td2m_mold
  type(   mxtmnode), target, save::    mxtm_mold
  type(   mitmnode), target, save::    mitm_mold
  type(   pmslnode), target, save::    pmsl_mold
  type(   howvnode), target, save::    howv_mold
  type(  tcamtnode), target, save::   tcamt_mold
  type(  lcbasnode), target, save::   lcbas_mold
  type(   pm10node), target, save::    pm10_mold
  type(  cldchnode), target, save::   cldch_mold

  type(   swcpnode), target, save::    swcp_mold
  type(   lwcpnode), target, save::    lwcp_mold
  type(  lightnode), target, save::   light_mold
  type(  dbznode),   target, save::     dbz_mold
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname='m_obsnodetypemanager'

! UseCase 1: configuration of a single mold
!
!   use m_obsnodetypemanager, only: obsnode_typemold
!   use m_psnode, only: i_psnode
!   ...
!   allocate(psllist%mold, source=obsnode_typemold(i_psnode))
! or, for Fortran 2008 allocate() with mold= specifier
!   allocate(psllist%mold,   mold=obsnode_typemold(i_psnode))
!
! UseCase 2: configuration of molds in an array
!
!   use m_obsllist, only: obsllist_moldconfig
!   use m_obsnodetypemanager, only: obsnode_typemold
!   ...
!   do jtype=lbound(obsdiags,2),ubound(obsdiags,2)
!      do ibin=lbound(obsdiags,1),ubound(obsdiags,1)
!         call obsllist_moldconfig(obsdiags(ibin,jtype),mold=obsnode_typemold(jtype))
!      enddo
!   enddo
!

  enum, bind(C)
     enumerator:: iobsnode_zero_   = 0

     enumerator:: iobsnode_ps
     enumerator:: iobsnode_t
     enumerator:: iobsnode_w
     enumerator:: iobsnode_q
     enumerator:: iobsnode_spd
     enumerator:: iobsnode_rw
     enumerator:: iobsnode_dw
     enumerator:: iobsnode_sst
     enumerator:: iobsnode_pw
     enumerator:: iobsnode_pcp
     enumerator:: iobsnode_oz
     enumerator:: iobsnode_o3l
     enumerator:: iobsnode_gps
     enumerator:: iobsnode_rad
     enumerator:: iobsnode_tcp
     enumerator:: iobsnode_lag
     enumerator:: iobsnode_colvk
     enumerator:: iobsnode_aero
     enumerator:: iobsnode_aerol
     enumerator:: iobsnode_pm2_5
     enumerator:: iobsnode_gust
     enumerator:: iobsnode_vis
     enumerator:: iobsnode_pblh
     enumerator:: iobsnode_wspd10m
     enumerator:: iobsnode_uwnd10m
     enumerator:: iobsnode_vwnd10m
     enumerator:: iobsnode_td2m
     enumerator:: iobsnode_mxtm
     enumerator:: iobsnode_mitm
     enumerator:: iobsnode_pmsl
     enumerator:: iobsnode_howv
     enumerator:: iobsnode_tcamt
     enumerator:: iobsnode_lcbas
     enumerator:: iobsnode_pm10
     enumerator:: iobsnode_cldch
     enumerator:: iobsnode_swcp
     enumerator:: iobsnode_lwcp
     enumerator:: iobsnode_light
     enumerator:: iobsnode_dbz

     enumerator:: iobsnode_extra_
  end enum
  
  integer(i_kind),parameter:: iobsnode_kind = kind(iobsnode_zero_)

  integer(iobsnode_kind),parameter:: obsnodetype_undef  = -1_iobsnode_kind
  integer(iobsnode_kind),parameter:: obsnodetype_lbound = iobsnode_zero_ +1
  integer(iobsnode_kind),parameter:: obsnodetype_ubound = iobsnode_extra_-1
  integer(iobsnode_kind),parameter:: obsnodetype_count  = obsnodetype_ubound-obsnodetype_lbound+1

contains
function vname2index_(vname) result(index_)
  use mpeu_util, only: lowercase
  implicit none
  integer(i_kind):: index_
  character(len=*),intent(in):: vname
  character(len=len(vname)):: vname_
  vname_=lowercase(vname)

  index_=0      ! a default return value, if the given name is unknown.
  select case(vname_)
     case("ps"   ,   "[psnode]"); index_ = iobsnode_ps
     case("t"    ,    "[tnode]"); index_ = iobsnode_t
     case("w"    ,    "[wnode]"); index_ = iobsnode_w
     case("q"    ,    "[qnode]"); index_ = iobsnode_q
     case("spd"  ,  "[spdnode]"); index_ = iobsnode_spd
     case("rw"   ,   "[rwnode]"); index_ = iobsnode_rw
     case("dw"   ,   "[dwnode]"); index_ = iobsnode_dw
     case("sst"  ,  "[sstnode]"); index_ = iobsnode_sst
     case("pw"   ,   "[pwnode]"); index_ = iobsnode_pw
     case("pcp"  ,  "[pcpnode]"); index_ = iobsnode_pcp
     case("oz"   ,   "[oznode]"); index_ = iobsnode_oz
     case("o3l"  ,  "[o3lnode]"); index_ = iobsnode_o3l
     case("gps"  ,  "[gpsnode]"); index_ = iobsnode_gps
     case("rad"  ,  "[radnode]"); index_ = iobsnode_rad
     case("tcp"  ,  "[tcpnode]"); index_ = iobsnode_tcp
     case("lag"  ,  "[lagnode]"); index_ = iobsnode_lag
     case("colvk","[colvknode]"); index_ = iobsnode_colvk
     case("aero" , "[aeronode]"); index_ = iobsnode_aero
     case("aerol","[aerolnode]"); index_ = iobsnode_aerol
     case("pm2_5","[pm2_5node]"); index_ = iobsnode_pm2_5
     case("gust" , "[gustnode]"); index_ = iobsnode_gust
     case("vis"  ,  "[visnode]"); index_ = iobsnode_vis
     case("pblh" , "[pblhnode]"); index_ = iobsnode_pblh

     case("wspd10m", &
                "[wspd10mnode]"); index_ = iobsnode_wspd10m
     case("uwnd10m", &
                "[uwnd10mnode]"); index_ = iobsnode_uwnd10m
     case("vwnd10m", &
                "[vwnd10mnode]"); index_ = iobsnode_vwnd10m

     case("td2m" , "[td2mnode]"); index_ = iobsnode_td2m
     case("mxtm" , "[mxtmnode]"); index_ = iobsnode_mxtm
     case("mitm" , "[mitmnode]"); index_ = iobsnode_mitm
     case("pmsl" , "[pmslnode]"); index_ = iobsnode_pmsl
     case("howv" , "[howvnode]"); index_ = iobsnode_howv
     case("tcamt","[tcamtnode]"); index_ = iobsnode_tcamt
     case("lcbas","[lcbasnode]"); index_ = iobsnode_lcbas

     case("pm10" , "[pm10node]"); index_ = iobsnode_pm10
     case("cldch","[cldchnode]"); index_ = iobsnode_cldch

     case("swcp" , "[swcpnode]"); index_ = iobsnode_swcp
     case("lwcp" , "[lwcpnode]"); index_ = iobsnode_lwcp

     case("light","[lightnode]"); index_ = iobsnode_light
     case("dbz"  ,  "[dbznode]"); index_ = iobsnode_dbz

  end select
end function vname2index_

function vmold2index_(mold) result(index_)
  implicit none
  integer(i_kind):: index_
  class(obsnode),target,intent(in):: mold

  index_=vname2index_(mold%mytype())
end function vmold2index_

function vmold2index_select_(mold) result(index_)
  implicit none
  integer(i_kind):: index_
  class(obsnode),target,intent(in):: mold

  index_=0
  select type(mold)
     type is(   psnode); index_ = iobsnode_ps
     type is(    tnode); index_ = iobsnode_t
     type is(    wnode); index_ = iobsnode_w
     type is(    qnode); index_ = iobsnode_q
     type is(  spdnode); index_ = iobsnode_spd
     type is(   rwnode); index_ = iobsnode_rw
     type is(   dwnode); index_ = iobsnode_dw
     type is(  sstnode); index_ = iobsnode_sst
     type is(   pwnode); index_ = iobsnode_pw
     type is(  pcpnode); index_ = iobsnode_pcp
     type is(   oznode); index_ = iobsnode_oz
     type is(  o3lnode); index_ = iobsnode_o3l
     type is(  gpsnode); index_ = iobsnode_gps
     type is(  radnode); index_ = iobsnode_rad
     type is(  tcpnode); index_ = iobsnode_tcp
     type is(  lagnode); index_ = iobsnode_lag
     type is(colvknode); index_ = iobsnode_colvk
     type is( aeronode); index_ = iobsnode_aero
     type is(aerolnode); index_ = iobsnode_aerol
     type is(pm2_5node); index_ = iobsnode_pm2_5
     type is( gustnode); index_ = iobsnode_gust
     type is(  visnode); index_ = iobsnode_vis
     type is( pblhnode); index_ = iobsnode_pblh

     type is(wspd10mnode); index_ = iobsnode_wspd10m
     type is(uwnd10mnode); index_ = iobsnode_uwnd10m
     type is(vwnd10mnode); index_ = iobsnode_vwnd10m
 
     type is( td2mnode); index_ = iobsnode_td2m
     type is( mxtmnode); index_ = iobsnode_mxtm
     type is( mitmnode); index_ = iobsnode_mitm
     type is( pmslnode); index_ = iobsnode_pmsl
     type is( howvnode); index_ = iobsnode_howv
     type is(tcamtnode); index_ = iobsnode_tcamt
     type is(lcbasnode); index_ = iobsnode_lcbas

     type is( pm10node); index_ = iobsnode_pm10
     type is(cldchnode); index_ = iobsnode_cldch
 
     type is( swcpnode); index_ = iobsnode_swcp
     type is( lwcpnode); index_ = iobsnode_lwcp

     type is(lightnode); index_ = iobsnode_light
     type is(  dbznode); index_ = iobsnode_dbz

  end select
end function vmold2index_select_

function index2vmold_(i_obtype) result(obsmold_)
  implicit none
  class(obsnode),pointer:: obsmold_
  integer(kind=i_kind),intent(in):: i_obtype

  character(len=*),parameter:: myname_=myname//"::index2vmold_"

  obsmold_ => null()
  select case(i_obtype)
     case(iobsnode_ps   ); obsmold_ =>    ps_mold
     case(iobsnode_t    ); obsmold_ =>     t_mold
     case(iobsnode_w    ); obsmold_ =>     w_mold
     case(iobsnode_q    ); obsmold_ =>     q_mold
     case(iobsnode_spd  ); obsmold_ =>   spd_mold
     case(iobsnode_rw   ); obsmold_ =>    rw_mold
     case(iobsnode_dw   ); obsmold_ =>    dw_mold
     case(iobsnode_sst  ); obsmold_ =>   sst_mold
     case(iobsnode_pw   ); obsmold_ =>    pw_mold
     case(iobsnode_pcp  ); obsmold_ =>   pcp_mold
     case(iobsnode_oz   ); obsmold_ =>    oz_mold
     case(iobsnode_o3l  ); obsmold_ =>   o3l_mold
     case(iobsnode_gps  ); obsmold_ =>   gps_mold
     case(iobsnode_rad  ); obsmold_ =>   rad_mold
     case(iobsnode_tcp  ); obsmold_ =>   tcp_mold
     case(iobsnode_lag  ); obsmold_ =>   lag_mold
     case(iobsnode_colvk); obsmold_ => colvk_mold
     case(iobsnode_aero ); obsmold_ =>  aero_mold
     case(iobsnode_aerol); obsmold_ => aerol_mold
     case(iobsnode_pm2_5); obsmold_ => pm2_5_mold
     case(iobsnode_gust ); obsmold_ =>  gust_mold
     case(iobsnode_vis  ); obsmold_ =>   vis_mold
     case(iobsnode_pblh ); obsmold_ =>  pblh_mold
 
     case(iobsnode_wspd10m); obsmold_ => wspd10m_mold
     case(iobsnode_uwnd10m); obsmold_ => uwnd10m_mold
     case(iobsnode_vwnd10m); obsmold_ => vwnd10m_mold
 
     case(iobsnode_td2m ); obsmold_ =>    td2m_mold
     case(iobsnode_mxtm ); obsmold_ =>    mxtm_mold
     case(iobsnode_mitm ); obsmold_ =>    mitm_mold
     case(iobsnode_pmsl ); obsmold_ =>    pmsl_mold
     case(iobsnode_howv ); obsmold_ =>    howv_mold
     case(iobsnode_tcamt); obsmold_ =>   tcamt_mold
     case(iobsnode_lcbas); obsmold_ =>   lcbas_mold
 
     case(iobsnode_pm10 ); obsmold_ =>    pm10_mold
     case(iobsnode_cldch); obsmold_ =>   cldch_mold

     case(iobsnode_swcp ); obsmold_ =>    swcp_mold
     case(iobsnode_lwcp ); obsmold_ =>    lwcp_mold

     case(iobsnode_light); obsmold_ =>   light_mold
     case(iobsnode_dbz);   obsmold_ =>     dbz_mold

  end select
end function index2vmold_

function vname2vmold_(vname) result(obsmold_)
  implicit none
  class(obsnode),pointer:: obsmold_
  character(len=*),intent(in):: vname

  character(len=*),parameter:: myname_=myname//"::vname2vmold_"
  integer(kind=i_kind):: i_obtype

  i_obtype=vname2index_(vname)
  obsmold_ => index2vmold_(i_obtype)
end function vname2vmold_

end module m_obsnodetypemanager
