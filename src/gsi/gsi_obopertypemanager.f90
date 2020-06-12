module gsi_obopertypemanager
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    module gsi_obopertypemanager
!   prgmmr:      j guo <jguo@nasa.gov>
!      org:      NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:      2018-07-12
!
! abstract: GSI observation operator (oboper) type manager
!
! program history log:
!   2018-07-12  j guo   - a type-manager for all oboper extensions.
!                       - an enum mapping of obsinput::dtype(:) to oboper type
!                         extensions.
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

  use gsi_oboper, only: oboper

  use gsi_aerooper    , only: aerooper
  use gsi_cldchoper   , only: cldchoper
  use gsi_colvkoper   , only: colvkoper
  use gsi_dwoper      , only: dwoper
  use gsi_gpsbendoper , only: gpsbendoper
  use gsi_gpsrefoper  , only: gpsrefoper
  use gsi_gustoper    , only: gustoper
  use gsi_howvoper    , only: howvoper
  use gsi_lcbasoper   , only: lcbasoper
  use gsi_lwcpoper    , only: lwcpoper
  use gsi_mitmoper    , only: mitmoper
  use gsi_mxtmoper    , only: mxtmoper
  use gsi_o3loper     , only: o3loper
  use gsi_ozoper      , only: ozoper
  use gsi_pblhoper    , only: pblhoper
  use gsi_pcpoper     , only: pcpoper
  use gsi_pm10oper    , only: pm10oper
  use gsi_pm2_5oper   , only: pm2_5oper
  use gsi_pmsloper    , only: pmsloper
  use gsi_psoper      , only: psoper
  use gsi_pwoper      , only: pwoper
  use gsi_qoper       , only: qoper
  use gsi_radoper     , only: radoper
  use gsi_rwoper      , only: rwoper
  use gsi_spdoper     , only: spdoper
  use gsi_sstoper     , only: sstoper
  use gsi_swcpoper    , only: swcpoper
  use gsi_tcamtoper   , only: tcamtoper
  use gsi_tcpoper     , only: tcpoper
  use gsi_td2moper    , only: td2moper
  use gsi_toper       , only: toper
  use gsi_uwnd10moper , only: uwnd10moper
  use gsi_visoper     , only: visoper
  use gsi_vwnd10moper , only: vwnd10moper
  use gsi_woper       , only: woper
  use gsi_wspd10moper , only: wspd10moper

  use gsi_lightoper   , only: lightoper
  use gsi_dbzoper     , only: dbzoper
  use gsi_cldtotoper  , only: cldtotoper

  use kinds     , only: i_kind
  use mpeu_util , only: perr,die
  implicit none
  private       ! except

  public:: oboper_typeMold
  public:: oboper_typeIndex
  public:: oboper_typeInfo
  interface oboper_typeMold; module procedure &
     dtype2vmold_,   &
     index2vmold_    ; end interface
  interface oboper_typeIndex; module procedure &
     vmold2index_,   &
     dtype2index_    ; end interface
  interface oboper_typeInfo; module procedure &
     vmold2tinfo_,   &
     index2tinfo_    ; end interface

  !public:: oboper_config
  !      interface oboper_config; module procedure config_; end interface

  public:: oboper_undef
  public:: oboper_lbound
  public:: oboper_ubound
  !public:: oboper_size
  public:: oboper_count

  public:: ioboper_kind
  public:: ioboper_ps
  public:: ioboper_t
  public:: ioboper_w
  public:: ioboper_q
  public:: ioboper_spd
  public:: ioboper_rw
  public:: ioboper_dw
  public:: ioboper_sst
  public:: ioboper_pw
  public:: ioboper_pcp
  public:: ioboper_oz
  public:: ioboper_o3l
  public:: ioboper_gpsbend
  public:: ioboper_gpsref
  public:: ioboper_rad
  public:: ioboper_tcp
 !public:: ioboper_lag
  public:: ioboper_colvk
  public:: ioboper_aero
 !public:: ioboper_aerol
  public:: ioboper_pm2_5
  public:: ioboper_gust
  public:: ioboper_vis
  public:: ioboper_pblh
  public:: ioboper_wspd10m
  public:: ioboper_td2m
  public:: ioboper_mxtm
  public:: ioboper_mitm
  public:: ioboper_pmsl
  public:: ioboper_howv
  public:: ioboper_tcamt
  public:: ioboper_lcbas
  public:: ioboper_pm10
  public:: ioboper_cldch
  public:: ioboper_uwnd10m
  public:: ioboper_vwnd10m
  public:: ioboper_swcp
  public:: ioboper_lwcp
  public:: ioboper_light
  public:: ioboper_dbz
  public:: ioboper_cldtot

  enum, bind(c)
     enumerator:: ioboper_zero_   = 0

     enumerator:: ioboper_ps
     enumerator:: ioboper_t
     enumerator:: ioboper_w
     enumerator:: ioboper_q
     enumerator:: ioboper_spd
     enumerator:: ioboper_rw
     enumerator:: ioboper_dw
     enumerator:: ioboper_sst
     enumerator:: ioboper_pw
     enumerator:: ioboper_pcp
     enumerator:: ioboper_oz
     enumerator:: ioboper_o3l
     enumerator:: ioboper_gpsbend
     enumerator:: ioboper_gpsref
     enumerator:: ioboper_rad
     enumerator:: ioboper_tcp
    !enumerator:: ioboper_lag
     enumerator:: ioboper_colvk
     enumerator:: ioboper_aero
    !enumerator:: ioboper_aerol
     enumerator:: ioboper_pm2_5
     enumerator:: ioboper_gust
     enumerator:: ioboper_vis
     enumerator:: ioboper_pblh
     enumerator:: ioboper_wspd10m
     enumerator:: ioboper_td2m
     enumerator:: ioboper_mxtm
     enumerator:: ioboper_mitm
     enumerator:: ioboper_pmsl
     enumerator:: ioboper_howv
     enumerator:: ioboper_tcamt
     enumerator:: ioboper_lcbas
     enumerator:: ioboper_pm10
     enumerator:: ioboper_cldch
     enumerator:: ioboper_uwnd10m
     enumerator:: ioboper_vwnd10m
     enumerator:: ioboper_swcp
     enumerator:: ioboper_lwcp
     enumerator:: ioboper_light
     enumerator:: ioboper_dbz
     enumerator:: ioboper_cldtot
 
     enumerator:: ioboper_extra_
  end enum

  integer(i_kind),parameter:: enum_kind = kind(ioboper_zero_)
  integer(i_kind),parameter:: ioboper_kind = enum_kind

  integer(enum_kind),parameter:: oboper_undef  = -1_enum_kind
  integer(enum_kind),parameter:: oboper_lbound = ioboper_zero_ +1
  integer(enum_kind),parameter:: oboper_ubound = ioboper_extra_-1
  integer(enum_kind),parameter:: oboper_size   = oboper_ubound-oboper_lbound+1
  integer(enum_kind),parameter:: oboper_count  = oboper_size

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname='gsi_obopertypemanager'
  logical,save:: oboper_configured_ = .false.

  character(len=20),dimension(oboper_lbound:oboper_ubound):: cobstype
  logical,save:: cobstype_configured_=.false.

  type(     psoper), target, save::       psoper_mold
  type(      toper), target, save::        toper_mold
  type(      woper), target, save::        woper_mold
  type(      qoper), target, save::        qoper_mold
  type(    spdoper), target, save::      spdoper_mold
  type(     rwoper), target, save::       rwoper_mold
  type(     dwoper), target, save::       dwoper_mold
  type(    sstoper), target, save::      sstoper_mold
  type(     pwoper), target, save::       pwoper_mold
  type(    pcpoper), target, save::      pcpoper_mold
  type(     ozoper), target, save::       ozoper_mold
  type(    o3loper), target, save::      o3loper_mold
  type(gpsbendoper), target, save::  gpsbendoper_mold
  type( gpsrefoper), target, save::   gpsrefoper_mold
  type(    radoper), target, save::      radoper_mold
  type(    tcpoper), target, save::      tcpoper_mold
 !type(    lagoper), target, save::      lagoper_mold
  type(  colvkoper), target, save::    colvkoper_mold
  type(   aerooper), target, save::     aerooper_mold
 !type(  aeroloper), target, save::    aeroloper_mold
  type(  pm2_5oper), target, save::    pm2_5oper_mold
  type(   gustoper), target, save::     gustoper_mold
  type(    visoper), target, save::      visoper_mold
  type(   pblhoper), target, save::     pblhoper_mold
  type(wspd10moper), target, save::  wspd10moper_mold
  type(   td2moper), target, save::     td2moper_mold
  type(   mxtmoper), target, save::     mxtmoper_mold
  type(   mitmoper), target, save::     mitmoper_mold
  type(   pmsloper), target, save::     pmsloper_mold
  type(   howvoper), target, save::     howvoper_mold
  type(  tcamtoper), target, save::    tcamtoper_mold
  type(  lcbasoper), target, save::    lcbasoper_mold
  type(   pm10oper), target, save::     pm10oper_mold
  type(  cldchoper), target, save::    cldchoper_mold
  type(uwnd10moper), target, save::  uwnd10moper_mold
  type(vwnd10moper), target, save::  vwnd10moper_mold
  type(   swcpoper), target, save::     swcpoper_mold
  type(   lwcpoper), target, save::     lwcpoper_mold
  type(  lightoper), target, save::    lightoper_mold
  type(    dbzoper), target, save::      dbzoper_mold
  type( cldtotoper), target, save::   cldtotoper_mold

contains
function dtype2index_(dtype) result(index_)
  use mpeu_util, only: lowercase
  implicit none
  integer(i_kind):: index_
  character(len=*),intent(in):: dtype

  select case(lowercase(dtype))
     case("ps"     ,"[psoper]"     ); index_= ioboper_ps
     case("t"      ,"[toper]"      ); index_= ioboper_t

     case("w"      ,"[woper]"      ); index_= ioboper_w
     case("uv"     ); index_= ioboper_w

     case("q"      ,"[qoper]"      ); index_= ioboper_q
     case("spd"    ,"[spdoper]"    ); index_= ioboper_spd
     case("rw"     ,"[rwoper]"     ); index_= ioboper_rw
     case("dw"     ,"[dwoper]"     ); index_= ioboper_dw
     case("sst"    ,"[sstoper]"    ); index_= ioboper_sst
     case("pw"     ,"[pwoper]"     ); index_= ioboper_pw

     case("pcp"    ,"[pcpoper]"    ); index_= ioboper_pcp
     case("pcp_ssmi"); index_= ioboper_pcp
     case("pcp_tmi" ); index_= ioboper_pcp

     case("oz"     ,"[ozoper]"     ); index_= ioboper_oz
     case("sbuv2"  ); index_= ioboper_oz
     case("omi"    ); index_= ioboper_oz
     case("gome"   ); index_= ioboper_oz
     case("ompstc8"); index_= ioboper_oz
     case("ompsnp" ); index_= ioboper_oz
     case("ompsnm" ); index_= ioboper_oz

     case("o3l"    ,"[o3loper]"    ); index_= ioboper_o3l
     case("o3lev"    ); index_= ioboper_o3l
     case("mls20"    ); index_= ioboper_o3l
     case("mls22"    ); index_= ioboper_o3l
     case("mls30"    ); index_= ioboper_o3l
     case("mls55"    ); index_= ioboper_o3l
     case("omieff"   ); index_= ioboper_o3l
     case("tomseff"  ); index_= ioboper_o3l
     case("ompslpuv" ); index_= ioboper_o3l
     case("ompslpvis"); index_= ioboper_o3l
     case("ompslp"   ); index_= ioboper_o3l

     case("gpsbend","[gpsbendoper]"); index_= ioboper_gpsbend
     case("gps_bnd"); index_= ioboper_gpsbend

     case("gpsref" ,"[gpsrefoper]" ); index_= ioboper_gpsref
     case("gps_ref"); index_= ioboper_gpsref

     case("rad"    ,"[radoper]"    ); index_= ioboper_rad
        !
     case("abi"    ); index_= ioboper_rad
        !
     case("amsua"  ); index_= ioboper_rad
     case("amsub"  ); index_= ioboper_rad
     case("msu"    ); index_= ioboper_rad
     case("mhs"    ); index_= ioboper_rad
     case("hirs2"  ); index_= ioboper_rad
     case("hirs3"  ); index_= ioboper_rad
     case("hirs4"  ); index_= ioboper_rad
     case("ssu"    ); index_= ioboper_rad
        !
     case("atms"   ); index_= ioboper_rad
     case("saphir" ); index_= ioboper_rad
         !
     case("airs"   ); index_= ioboper_rad
     case("hsb"    ); index_= ioboper_rad
        !
     case("iasi"   ); index_= ioboper_rad
     case("cris"   ); index_= ioboper_rad
     case("cris-fsr"  ); index_= ioboper_rad
        !
     case("sndr"   ); index_= ioboper_rad
     case("sndrd1" ); index_= ioboper_rad
     case("sndrd2" ); index_= ioboper_rad
     case("sndrd3" ); index_= ioboper_rad
     case("sndrd4" ); index_= ioboper_rad
        !
     case("ssmi"   ); index_= ioboper_rad
        !
     case("amsre"  ); index_= ioboper_rad
     case("amsre_low"); index_= ioboper_rad
     case("amsre_mid"); index_= ioboper_rad
     case("amsre_hig"); index_= ioboper_rad
        !
     case("ssmis"  ); index_= ioboper_rad
     case("ssmis_las"); index_= ioboper_rad
     case("ssmis_uas"); index_= ioboper_rad
     case("ssmis_img"); index_= ioboper_rad
     case("ssmis_env"); index_= ioboper_rad
        !
     case("amsr2"  ); index_= ioboper_rad
     case("goes_img"); index_= ioboper_rad
     case("gmi"    ); index_= ioboper_rad
     case("seviri" ); index_= ioboper_rad
     case("ahi"    ); index_= ioboper_rad
        !
     case("avhrr_navy"); index_= ioboper_rad
     case("avhrr"  ); index_= ioboper_rad
 
     case("tcp"    ,"[tcpoper]"    ); index_= ioboper_tcp

!    case("lag"    ,"[lagoper]"    ); index_= ioboper_lag

     case("colvk"  ,"[colvkoper]"  ); index_= ioboper_colvk
     case("mopitt" ); index_= ioboper_colvk

     case("aero"   ,"[aerooper]"   ); index_= ioboper_aero
     case("aod"      ); index_= ioboper_aero
     case("modis_aod"); index_= ioboper_aero

!    case("aerol"  ,"[aeroloper]"  ); index_= ioboper_aerol

     case("pm2_5"  ,"[pm2_5oper]"  ); index_= ioboper_pm2_5
     case("gust"   ,"[gustoper]"   ); index_= ioboper_gust
     case("vis"    ,"[visoper]"    ); index_= ioboper_vis
     case("pblh"   ,"[pblhoper]"   ); index_= ioboper_pblh

     case("wspd10m","[wspd10moper]"); index_= ioboper_wspd10m
     case("uwnd10m","[uwnd10moper]"); index_= ioboper_uwnd10m
     case("vwnd10m","[vwnd10moper]"); index_= ioboper_vwnd10m

     case("td2m"   ,"[td2moper]"   ); index_= ioboper_td2m
     case("mxtm"   ,"[mxtmoper]"   ); index_= ioboper_mxtm
     case("mitm"   ,"[mitmoper]"   ); index_= ioboper_mitm
     case("pmsl"   ,"[pmsloper]"   ); index_= ioboper_pmsl
     case("howv"   ,"[howvoper]"   ); index_= ioboper_howv
     case("tcamt"  ,"[tcamtoper]"  ); index_= ioboper_tcamt
     case("lcbas"  ,"[lcbasoper]"  ); index_= ioboper_lcbas
 
     case("pm10"   ,"[pm10oper]"   ); index_= ioboper_pm10
     case("cldch"  ,"[cldchoper]"  ); index_= ioboper_cldch

     case("swcp"   ,"[swcpoper]"   ); index_= ioboper_swcp
     case("lwcp"   ,"[lwcpoper]"   ); index_= ioboper_lwcp
 
     case("light"  ,"[lightoper]"  ); index_= ioboper_light
     case("goes_glm" ); index_= ioboper_light

     case("dbz"    ,"[dbzoper]"    ); index_= ioboper_dbz

     case("cldtot" ,"[cldtotoper]" ); index_= ioboper_cldtot
     case("mta_cld"  ); index_= ioboper_cldtot

        ! Known dtype values, but no known oboper type defined
     case("gos_ctp"); index_= oboper_undef
     case("rad_ref"); index_= oboper_undef
     case("lghtn"  ); index_= oboper_undef
     case("larccld"); index_= oboper_undef
     case("larcglb"); index_= oboper_undef

        ! A catch all case
     case default   ; index_= oboper_undef
  end select
end function dtype2index_

function vmold2index_(mold) result(index_)
  implicit none
  integer(i_kind):: index_
  class(oboper),target,intent(in):: mold

  character(len=*),parameter:: myname_=myname//"::vmold2index_"
  class(oboper),pointer:: ptr_
  ptr_ => mold
  if(.not.associated(ptr_)) call die(myname_,'not assoicated, argument mold')
  nullify(ptr_)

  index_=dtype2index_(mold%mytype())

  ! An alternative implementation to cache a managed ioboper value inside each
  ! oboper class.  This implementation requires two new tbps, %myinfo_get() and
  ! %myinfo_set().
  !
  ! call mold%myinfo_get(ioboper=index_)
  ! if(index_<oboper_lbound .or. index_>oboper_ubound) then
  !    index_=dtype2index_(mold%mytype())
  !    call mold%myinfo_set(ioboper_=index_)
  ! endif

end function vmold2index_

function dtype2vmold_(dtype) result(vmold_)
  implicit none
  class(oboper),pointer:: vmold_
  character(len=*),intent(in):: dtype

  integer(i_kind):: ioboper_
  ioboper_ = dtype2index_(dtype)
  vmold_ => index2vmold_(ioboper_)
end function dtype2vmold_

function index2vmold_(ioboper) result(vmold_)
  implicit none
  class(oboper),pointer:: vmold_
  integer(i_kind),intent(in):: ioboper
  select case(ioboper)

     case(ioboper_ps       ); vmold_ =>      psoper_mold
     case(ioboper_t        ); vmold_ =>       toper_mold
     case(ioboper_w        ); vmold_ =>       woper_mold
     case(ioboper_q        ); vmold_ =>       qoper_mold
     case(ioboper_spd      ); vmold_ =>     spdoper_mold
     case(ioboper_rw       ); vmold_ =>      rwoper_mold
     case(ioboper_dw       ); vmold_ =>      dwoper_mold
     case(ioboper_sst      ); vmold_ =>     sstoper_mold
     case(ioboper_pw       ); vmold_ =>      pwoper_mold
     case(ioboper_pcp      ); vmold_ =>     pcpoper_mold
     case(ioboper_oz       ); vmold_ =>      ozoper_mold
     case(ioboper_o3l      ); vmold_ =>     o3loper_mold
     case(ioboper_gpsbend  ); vmold_ => gpsbendoper_mold
     case(ioboper_gpsref   ); vmold_ =>  gpsrefoper_mold
     case(ioboper_rad      ); vmold_ =>     radoper_mold
     case(ioboper_tcp      ); vmold_ =>     tcpoper_mold
    !case(ioboper_lag      ); vmold_ =>     lagoper_mold
     case(ioboper_colvk    ); vmold_ =>   colvkoper_mold
     case(ioboper_aero     ); vmold_ =>    aerooper_mold
    !case(ioboper_aerol    ); vmold_ =>   aeroloper_mold
     case(ioboper_pm2_5    ); vmold_ =>   pm2_5oper_mold
     case(ioboper_gust     ); vmold_ =>    gustoper_mold
     case(ioboper_vis      ); vmold_ =>     visoper_mold
     case(ioboper_pblh     ); vmold_ =>    pblhoper_mold
     case(ioboper_wspd10m  ); vmold_ => wspd10moper_mold
     case(ioboper_td2m     ); vmold_ =>    td2moper_mold
     case(ioboper_mxtm     ); vmold_ =>    mxtmoper_mold
     case(ioboper_mitm     ); vmold_ =>    mitmoper_mold
     case(ioboper_pmsl     ); vmold_ =>    pmsloper_mold
     case(ioboper_howv     ); vmold_ =>    howvoper_mold
     case(ioboper_tcamt    ); vmold_ =>   tcamtoper_mold
     case(ioboper_lcbas    ); vmold_ =>   lcbasoper_mold
     case(ioboper_pm10     ); vmold_ =>    pm10oper_mold
     case(ioboper_cldch    ); vmold_ =>   cldchoper_mold
     case(ioboper_uwnd10m  ); vmold_ => uwnd10moper_mold
     case(ioboper_vwnd10m  ); vmold_ => vwnd10moper_mold
     case(ioboper_swcp     ); vmold_ =>    swcpoper_mold
     case(ioboper_lwcp     ); vmold_ =>    lwcpoper_mold
     case(ioboper_light    ); vmold_ =>   lightoper_mold
     case(ioboper_dbz      ); vmold_ =>     dbzoper_mold
     case(ioboper_cldtot   ); vmold_ =>  cldtotoper_mold

     case( oboper_undef    ); vmold_ => null()
     case default           ; vmold_ => null()
  end select
end function index2vmold_

function vmold2tinfo_(mold) result(info_)
!>> Simply mold%info(), but just in case one needs some indirection, with
!>> multiple oboper classes.
  implicit none
  character(len=:),allocatable:: info_
  class(oboper),target,intent(in):: mold

  character(len=*),parameter:: myname_=myname//"::vmold2tinfo_"
  class(oboper),pointer:: vmold__
  vmold__ => mold

  if(.not.associated(vmold__)) call die(myname_,'not assoicated, argument mold')
  nullify(vmold__)

  info_=index2tinfo_(vmold2index_(mold))
end function vmold2tinfo_

function index2tinfo_(ioboper) result(info_)
!>>
  implicit none
  character(len=:),allocatable:: info_
  integer(i_kind),intent(in):: ioboper

  if(.not.cobstype_configured_) call cobstype_config_()
  info_=""
  if(ioboper>=oboper_lbound .and. &
     ioboper<=oboper_ubound) info_=cobstype(ioboper)
end function index2tinfo_

subroutine config_()
  implicit none
  character(len=*),parameter:: myname_=myname//"::config_"
  class(oboper),pointer:: vmold_
  integer(i_kind):: iset_,iget_
  logical:: good_

  good_=.true.
  do iset_ = oboper_lbound, oboper_ubound
     vmold_ => index2vmold_(iset_)
     if(.not.associated(vmold_)) then
        call perr(myname_,'unexpected index, iset_ =',iset_)
        call perr(myname_,'          obOper_lbound =',oboper_lbound)
        call perr(myname_,'          obOper_ubound =',oboper_ubound)
        call  die(myname_)
     endif

     iget_=iset_         ! for additional test.
    !call vmold_%myinfo_set(ioboper=iset_)
    !call vmold_%myinfo_get(ioboper=iget_)
     if(iget_/=iset_) then
        call perr(myname_,'unexpected return, %myinfo_get(iobOper) =',iget_)
        call perr(myname_,'                   %myinfo_set(iobOper) =',iset_)
        call perr(myname_,'                              %mytype() =',vmold_%mytype())
        good_=.false.
     endif

     vmold_ => null()
  enddo
  if(.not.good_) call die(myname_)

  oboper_configured_ = .true.
end subroutine config_

subroutine cobstype_config_()
!>> Should this information be provided by individual oboper extensions, or
!>> be provided by this manager?  There are pros and cons in either approach.

  implicit none
    cobstype(ioboper_ps         )  ="surface pressure    " ! ps_ob_type
    cobstype(ioboper_t          )  ="temperature         " ! t_ob_type
    cobstype(ioboper_w          )  ="wind                " ! w_ob_type
    cobstype(ioboper_q          )  ="moisture            " ! q_ob_type
    cobstype(ioboper_spd        )  ="wind speed          " ! spd_ob_type
    cobstype(ioboper_rw         )  ="radial wind         " ! rw_ob_type
    cobstype(ioboper_dw         )  ="doppler wind        " ! dw_ob_type
    cobstype(ioboper_sst        )  ="sst                 " ! sst_ob_type
    cobstype(ioboper_pw         )  ="precipitable water  " ! pw_ob_type
    cobstype(ioboper_pcp        )  ="precipitation       " ! pcp_ob_type
    cobstype(ioboper_oz         )  ="ozone               " ! oz_ob_type
    cobstype(ioboper_o3l        )  ="level ozone         " ! o3l_ob_type
    cobstype(ioboper_gpsbend    )  ="gps bending angle   " ! using gps_ob_type
    cobstype(ioboper_gpsref     )  ="gps refractivity    " ! using gps_ob_type
    cobstype(ioboper_rad        )  ="radiance            " ! rad_ob_type
    cobstype(ioboper_tcp        )  ="tcp (tropic cyclone)" ! tcp_ob_type
    !cobstype(ioboper_lag        )  ="lagrangian tracer   " ! lag_ob_type
    cobstype(ioboper_colvk      )  ="carbon monoxide     " ! colvk_ob_type
    cobstype(ioboper_aero       )  ="aerosol aod         " ! aero_ob_type
    !cobstype(ioboper_aerol      )  ="level aero aod      " ! aerol_ob_type
    cobstype(ioboper_pm2_5      )  ="in-situ pm2_5 obs   " ! pm2_5_ob_type
    cobstype(ioboper_pm10       )  ="in-situ pm10 obs    " ! pm10_ob_type
    cobstype(ioboper_gust       )  ="gust                " ! gust_ob_type
    cobstype(ioboper_vis        )  ="vis                 " ! vis_ob_type
    cobstype(ioboper_pblh       )  ="pblh                " ! pblh_ob_type
    cobstype(ioboper_wspd10m    )  ="wspd10m             " ! wspd10m_ob_type
    cobstype(ioboper_td2m       )  ="td2m                " ! td2m_ob_type
    cobstype(ioboper_mxtm       )  ="mxtm                " ! mxtm_ob_type
    cobstype(ioboper_mitm       )  ="mitm                " ! mitm_ob_type
    cobstype(ioboper_pmsl       )  ="pmsl                " ! pmsl_ob_type
    cobstype(ioboper_howv       )  ="howv                " ! howv_ob_type
    cobstype(ioboper_tcamt      )  ="tcamt               " ! tcamt_ob_type
    cobstype(ioboper_lcbas      )  ="lcbas               " ! lcbas_ob_type
    cobstype(ioboper_cldch      )  ="cldch               " ! cldch_ob_type
    cobstype(ioboper_uwnd10m    )  ="uwnd10m             " ! uwnd10m_ob_type
    cobstype(ioboper_vwnd10m    )  ="vwnd10m             " ! vwnd10m_ob_type
    cobstype(ioboper_swcp       )  ="swcp                " ! swcp_ob_type
    cobstype(ioboper_lwcp       )  ="lwcp                " ! lwcp_ob_type
    cobstype(ioboper_light      )  ="light               " ! light_ob_type
    cobstype(ioboper_dbz        )  ="dbz                 " ! dbz_ob_type
    cobstype(ioboper_cldtot     )  ="cldtot              " ! using q_ob_type

  cobstype_configured_=.true.
end subroutine cobstype_config_

end module gsi_obopertypemanager
