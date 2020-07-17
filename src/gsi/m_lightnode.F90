module m_lightnode
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    module m_lightnode
!   prgmmr:      k apodaca <karina.apodaca@colostate.edu>
!      org:      CSU/CIRA, Data Assimilation Group
!     date:      2018-01-18
!
! abstract: class-module of obs-type lightnode (lightning)
!
! program history log:
!   2018-01-18  k apodaca   - add this document block for the implementation of
!                             variational lightning data assimilation.
!   2018-08-26  k apodaca   - add coefficients relaed to a second observaion operator
!                             for lighning flash rate, suitable for non-hydrostatic, 
!                             cloud-resolving models.
!                             

!                .      .    .                                       .

! In the case of lightning observations (e.g. goes/glm), the schematic shown below is
! used for the interpolation of background fields to the location of an observation (+)
! and for the finite-difference derivation method used in the calculation of the tl of
! the observation operator for lightning flash rate. Calculations are done at each
! quadrant (i.e. central, north, south, east, and west).
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

!   2019-03-01  j guo   - Merged in some cleaning up changes as in other
!                         obsnode types:
!                       . Added a type specific subroutine appendto_(), to avoid
!                         unnecessary type generalization between a generic
!                         append() and user code.
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
  use m_obsdiagnode, only: obs_diag
  use m_obsdiagnode, only: obs_diags

  use kinds , only: i_kind,r_kind
  use mpeu_util, only: assert_,die,perr,warn,tell
  use m_obsnode, only: obsnode
  implicit none
  private

  public:: lightnode

  type,extends(obsnode):: lightnode
     !type(light_ob_type),pointer :: llpoint => null()
     type(obs_diag), pointer :: diags => null()
     real(r_kind)    :: res    =0._r_kind    !  light residual
     real(r_kind)    :: err2   =0._r_kind    !  light error squared
     real(r_kind)    :: raterr2=0._r_kind    !  square of ratio of final obs error
                                      !  to original obs error
     !real(r_kind)    :: time          !  observation time in sec
     real(r_kind)    :: b      =0._r_kind    !  variational quality control parameter
     real(r_kind)    :: pg     =0._r_kind    !  variational quality control parameter
     real(r_kind)    :: wij(4) =0._r_kind    !  horizontal interpolation weights

! Central quadrant
     real(r_kind),pointer               :: jac_z0i1  => null() ! surface z at i1
     real(r_kind),pointer               :: jac_z0i2  => null() ! surface z at i2
     real(r_kind),pointer               :: jac_z0i3  => null() ! surface z at i3
     real(r_kind),pointer               :: jac_z0i4  => null() ! surface z at i4
     real(r_kind),dimension(:),pointer  :: jac_vertqi1 => null()
     real(r_kind),dimension(:),pointer  :: jac_vertqi2 => null()
     real(r_kind),dimension(:),pointer  :: jac_vertqi3 => null()
     real(r_kind),dimension(:),pointer  :: jac_vertqi4 => null()
     real(r_kind),dimension(:),pointer  :: jac_vertti1 => null()
     real(r_kind),dimension(:),pointer  :: jac_vertti2 => null()
     real(r_kind),dimension(:),pointer  :: jac_vertti3 => null()
     real(r_kind),dimension(:),pointer  :: jac_vertti4 => null()
     real(r_kind),dimension(:),pointer  :: jac_zdxi1 => null()
     real(r_kind),dimension(:),pointer  :: jac_zdxi2 => null()
     real(r_kind),dimension(:),pointer  :: jac_zdxi3 => null()
     real(r_kind),dimension(:),pointer  :: jac_zdxi4 => null()
     real(r_kind),dimension(:),pointer  :: jac_zdyi1 => null()
     real(r_kind),dimension(:),pointer  :: jac_zdyi2 => null()
     real(r_kind),dimension(:),pointer  :: jac_zdyi3 => null()
     real(r_kind),dimension(:),pointer  :: jac_zdyi4 => null()
     real(r_kind),dimension(:),pointer  :: jac_udxi1 => null()
     real(r_kind),dimension(:),pointer  :: jac_udxi2 => null()
     real(r_kind),dimension(:),pointer  :: jac_udxi3 => null()
     real(r_kind),dimension(:),pointer  :: jac_udxi4 => null()
     real(r_kind),dimension(:),pointer  :: jac_vdyi1 => null()
     real(r_kind),dimension(:),pointer  :: jac_vdyi2 => null()
     real(r_kind),dimension(:),pointer  :: jac_vdyi3 => null()
     real(r_kind),dimension(:),pointer  :: jac_vdyi4 => null()
     real(r_kind),dimension(:),pointer  :: jac_vert => null()
     real(r_kind),dimension(:),pointer  :: jac_sigdoti1 => null()
     real(r_kind),dimension(:),pointer  :: jac_sigdoti2 => null()
     real(r_kind),dimension(:),pointer  :: jac_sigdoti3 => null()
     real(r_kind),dimension(:),pointer  :: jac_sigdoti4 => null()
     real(r_kind),dimension(:),pointer  :: jac_qi1 => null()
     real(r_kind),dimension(:),pointer  :: jac_qi2 => null()
     real(r_kind),dimension(:),pointer  :: jac_qi3 => null()
     real(r_kind),dimension(:),pointer  :: jac_qi4 => null()
     real(r_kind),dimension(:),pointer  :: jac_ti1 => null()
     real(r_kind),dimension(:),pointer  :: jac_ti2 => null()
     real(r_kind),dimension(:),pointer  :: jac_ti3 => null()
     real(r_kind),dimension(:),pointer  :: jac_ti4 => null()
     real(r_kind),dimension(:),pointer  :: jac_qgmai1 => null()
     real(r_kind),dimension(:),pointer  :: jac_qgmai2 => null()
     real(r_kind),dimension(:),pointer  :: jac_qgmai3 => null()
     real(r_kind),dimension(:),pointer  :: jac_qgmai4 => null()
     real(r_kind),dimension(:),pointer  :: jac_qgmbi1 => null()
     real(r_kind),dimension(:),pointer  :: jac_qgmbi2 => null()
     real(r_kind),dimension(:),pointer  :: jac_qgmbi3 => null()
     real(r_kind),dimension(:),pointer  :: jac_qgmbi4 => null()
     real(r_kind),dimension(:),pointer  :: jac_icei1 => null()
     real(r_kind),dimension(:),pointer  :: jac_icei2 => null()
     real(r_kind),dimension(:),pointer  :: jac_icei3 => null()
     real(r_kind),dimension(:),pointer  :: jac_icei4 => null()
     real(r_kind),dimension(:),pointer  :: jac_zicei1 => null()
     real(r_kind),dimension(:),pointer  :: jac_zicei2 => null()
     real(r_kind),dimension(:),pointer  :: jac_zicei3 => null()
     real(r_kind),dimension(:),pointer  :: jac_zicei4 => null()
     real(r_kind),pointer               :: kboti1 => null()
     real(r_kind),pointer               :: kboti2 => null()
     real(r_kind),pointer               :: kboti3 => null()
     real(r_kind),pointer               :: kboti4 => null()
     real(r_kind),pointer               :: jac_kverti1 => null()
     real(r_kind),pointer               :: jac_kverti2 => null()
     real(r_kind),pointer               :: jac_kverti3 => null()
     real(r_kind),pointer               :: jac_kverti4 => null()
     real(r_kind),pointer               :: jac_fratei1 => null()
     real(r_kind),pointer               :: jac_fratei2 => null()
     real(r_kind),pointer               :: jac_fratei3 => null()
     real(r_kind),pointer               :: jac_fratei4 => null()
     logical,pointer                    :: jac_wmaxflagi1 => null() ! wmax flag at i1
     logical,pointer                    :: jac_wmaxflagi2 => null() ! wmax flag at i2
     logical,pointer                    :: jac_wmaxflagi3 => null() ! wmax flag at i3
     logical,pointer                    :: jac_wmaxflagi4 => null() ! wmax flag at i4

! South quadrant
     real(r_kind),pointer               :: jac_z0i5 => null()
     real(r_kind),pointer               :: jac_z0i7 => null()
     real(r_kind),dimension(:),pointer  :: jac_vertqi5 => null()
     real(r_kind),dimension(:),pointer  :: jac_vertqi7 => null()
     real(r_kind),dimension(:),pointer  :: jac_vertti5 => null()
     real(r_kind),dimension(:),pointer  :: jac_vertti7 => null()

! North quadrant
     real(r_kind),pointer               :: jac_z0i6 => null()
     real(r_kind),pointer               :: jac_z0i8 => null()
     real(r_kind),dimension(:),pointer  :: jac_vertqi6 => null()
     real(r_kind),dimension(:),pointer  :: jac_vertqi8 => null()
     real(r_kind),dimension(:),pointer  :: jac_vertti6 => null()
     real(r_kind),dimension(:),pointer  :: jac_vertti8 => null()

! West quadrant
     real(r_kind),pointer               :: jac_z0i9  => null() ! surface z at i9
     real(r_kind),pointer               :: jac_z0i10 => null() ! surface z at i10
     real(r_kind),dimension(:),pointer  :: jac_vertqi9  => null()
     real(r_kind),dimension(:),pointer  :: jac_vertqi10 => null()
     real(r_kind),dimension(:),pointer  :: jac_vertti9  => null()
     real(r_kind),dimension(:),pointer  :: jac_vertti10 => null()

! East quadrant
     real(r_kind),pointer               :: jac_z0i11 => null() ! surface z at i11
     real(r_kind),pointer               :: jac_z0i12 => null() ! surface z at i12
     real(r_kind),dimension(:),pointer  :: jac_vertqi11 => null()
     real(r_kind),dimension(:),pointer  :: jac_vertqi12 => null()
     real(r_kind),dimension(:),pointer  :: jac_vertti11 => null()
     real(r_kind),dimension(:),pointer  :: jac_vertti12 => null()

     integer(i_kind),dimension(:,:),pointer :: ij  => null()
     !logical         :: luse          !  flag indicating if ob is used in pen.

     !integer(i_kind) :: idv,iob              ! device id and obs index for sorting
     !real   (r_kind) :: elat, elon      ! earth lat-lon for redistribution
     !real   (r_kind) :: dlat, dlon      ! earth lat-lon for redistribution
  contains
     procedure,nopass::  mytype
     procedure::  sethop => obsnode_sethop_
     procedure::   xread => obsnode_xread_
     procedure::  xwrite => obsnode_xwrite_
     procedure:: isvalid => obsnode_isvalid_
     procedure::  gettlddp => gettlddp_

     procedure, nopass:: headerread  => obsheader_read_
     procedure, nopass:: headerwrite => obsheader_write_
     procedure:: init  => obsnode_init_
     procedure:: clean => obsnode_clean_
  end type lightnode

  public:: lightnode_typecast
  public:: lightnode_nextcast
     interface lightnode_typecast; module procedure typecast_ ; end interface
     interface lightnode_nextcast; module procedure nextcast_ ; end interface

  public:: lightnode_appendto
     interface lightnode_appendto; module procedure appendto_ ; end interface

  character(len=*),parameter:: myname="m_lightnode"

#include "myassert.H"
#include "mytrace.H"
contains
function typecast_(anode) result(ptr_)
!-- cast a class(obsnode) to a type(lightnode)
  use m_obsnode, only: obsnode
  implicit none
  type(lightnode),pointer:: ptr_
  class(obsnode ),pointer,intent(in):: anode

  ptr_ => null()
  if(.not.associated(anode)) return
        ! logically, typecast of a null-reference is a null pointer.
  select type(anode)
     type is(lightnode)
        ptr_ => anode
  end select
  return
end function typecast_

function nextcast_(anode) result(ptr_)
!-- cast an obsnode_next(obsnode) to a type(lightnode)
  use m_obsnode, only: obsnode,obsnode_next
  implicit none
  type(lightnode),pointer:: ptr_
  class(obsnode ),target ,intent(in):: anode

  class(obsnode),pointer:: inode_
  inode_ => obsnode_next(anode)
  ptr_ => typecast_(inode_)
  return
end function nextcast_

subroutine appendto_(anode,oll)
!-- append anode to linked-list oll
  use m_obsnode , only: obsnode
  use m_obsllist, only: obsllist,obsllist_appendnode
  implicit none
  type(lightnode),pointer,intent(in):: anode
  type(obsllist),intent(inout):: oll

  class(obsnode),pointer:: inode_
  inode_ => anode
  call obsllist_appendnode(oll,inode_)
  inode_ => null()
end subroutine appendto_

! obsnode implementations

function mytype()
  implicit none
  character(len=:),allocatable:: mytype
  mytype="[lightnode]"
end function mytype

subroutine obsHeader_read_(iunit,mobs,jread,istat)
  use gridmod, only: nsig
  implicit none
  integer(i_kind),intent(in ):: iunit
  integer(i_kind),intent(out):: mobs
  integer(i_kind),intent(out):: jread
  integer(i_kind),intent(out):: istat

  character(len=*),parameter:: myname_=myname//"::obsHeader_read"
  integer(i_kind):: msig
_ENTRY_(myname_)

  msig=-1
  read(iunit,iostat=istat) mobs,jread, msig
  if(istat==0 .and. msig/=nsig) then
     call perr(myname_,'unmatched dimension information, expecting nsig =',nsig)
     call perr(myname_,'                                  but read msig =',msig)
     call  die(myname_)
  endif
_EXIT_(myname_)
  return
end subroutine obsHeader_read_

subroutine obsHeader_write_(junit,mobs,jwrite,jstat)
  use gridmod, only: nsig
  implicit none
  integer(i_kind),intent(in ):: junit
  integer(i_kind),intent(in ):: mobs
  integer(i_kind),intent(in ):: jwrite
  integer(i_kind),intent(out):: jstat

  character(len=*),parameter:: myname_=myname//"::obsHeader_write"
_ENTRY_(myname_)
  write(junit,iostat=jstat) mobs,jwrite, nsig
_EXIT_(myname_)
  return
end subroutine obsHeader_write_

subroutine obsnode_init_(anode)
  use gridmod, only: nsig
  implicit none
  class(lightnode),intent(out):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_init_'
_ENTRY_(myname_)
  allocate(anode%jac_z0i1,          anode%jac_z0i2,           & 
           anode%jac_z0i3,          anode%jac_z0i4,           &
           anode%jac_z0i5,          anode%jac_z0i6,           &
           anode%jac_z0i7,          anode%jac_z0i8,           &
           anode%jac_z0i9,          anode%jac_z0i10,          &
           anode%jac_z0i11,         anode%jac_z0i12           )

  allocate(anode%jac_vertqi1(nsig), anode%jac_vertqi2(nsig),  &
           anode%jac_vertqi3(nsig), anode%jac_vertqi4(nsig),  &
           anode%jac_vertqi5(nsig), anode%jac_vertqi6(nsig),  &
           anode%jac_vertqi7(nsig), anode%jac_vertqi8(nsig),  &
           anode%jac_vertqi9(nsig), anode%jac_vertqi10(nsig), &
           anode%jac_vertqi11(nsig),anode%jac_vertqi12(nsig)  )


  allocate(anode%jac_vertti1(nsig), anode%jac_vertti2(nsig),  &
           anode%jac_vertti3(nsig), anode%jac_vertti4(nsig),  &
           anode%jac_vertti5(nsig), anode%jac_vertti6(nsig),  &
           anode%jac_vertti7(nsig), anode%jac_vertti8(nsig),  &
           anode%jac_vertti9(nsig), anode%jac_vertti10(nsig), &
           anode%jac_vertti11(nsig),anode%jac_vertti12(nsig)  )

  allocate(anode%jac_zdxi1(nsig),   anode%jac_zdxi2(nsig),    &
           anode%jac_zdxi3(nsig),   anode%jac_zdxi4(nsig)     )

  allocate(anode%jac_zdyi1(nsig),   anode%jac_zdyi2(nsig),    &
           anode%jac_zdyi3(nsig),   anode%jac_zdyi4(nsig)     )

  allocate(anode%jac_udxi1(nsig),   anode%jac_udxi2(nsig),    &
           anode%jac_udxi3(nsig),   anode%jac_udxi4(nsig)     )

  allocate(anode%jac_vdyi1(nsig),   anode%jac_vdyi2(nsig),    &
           anode%jac_vdyi3(nsig),   anode%jac_vdyi4(nsig)     )

  allocate(anode%jac_vert(nsig)                               ) 

  allocate(anode%jac_sigdoti1(nsig),anode%jac_sigdoti2(nsig), &
           anode%jac_sigdoti3(nsig),anode%jac_sigdoti4(nsig)  )

  allocate(anode%jac_qi1(nsig),     anode%jac_qi2(nsig),      &
           anode%jac_qi3(nsig),     anode%jac_qi4(nsig)       )

  allocate(anode%jac_ti1(nsig),     anode%jac_ti2(nsig),      &
           anode%jac_ti3(nsig),     anode%jac_ti4(nsig)       ) 

  allocate(anode%jac_qgmai1(nsig),  anode%jac_qgmai2(nsig),   &
           anode%jac_qgmai3(nsig),  anode%jac_qgmai4(nsig)    )

  allocate(anode%jac_qgmbi1(nsig),  anode%jac_qgmbi2(nsig),   &
           anode%jac_qgmbi3(nsig),  anode%jac_qgmbi4(nsig)    )

  allocate(anode%jac_icei1(nsig),   anode%jac_icei2(nsig),    &
           anode%jac_icei3(nsig),   anode%jac_icei4(nsig)     )

  allocate(anode%jac_zicei1(nsig),  anode%jac_zicei2(nsig),   &
           anode%jac_zicei3(nsig),  anode%jac_zicei4(nsig)    )

  allocate(anode%kboti1,            anode%kboti2,             &
           anode%kboti3,            anode%kboti4              )

  allocate(anode%jac_kverti1,       anode%jac_kverti2,        &
           anode%jac_kverti3,       anode%jac_kverti4         )

  allocate(anode%jac_fratei1,       anode%jac_fratei2,        &
           anode%jac_fratei3,       anode%jac_fratei4         )

  allocate(anode%jac_wmaxflagi1,    anode%jac_wmaxflagi2,     &
           anode%jac_wmaxflagi3,    anode%jac_wmaxflagi4,     &
           anode%ij(12,nsig)                                  )
_EXIT_(myname_)
  return
end subroutine obsnode_init_

subroutine obsnode_clean_(anode)
  implicit none
  class(lightnode),intent(inout):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_clean_'
_ENTRY_(myname_)
!_TRACEV_(myname_,'%mytype() =',anode%mytype())
  if(associated(anode%jac_z0i1 ))  deallocate(anode%jac_z0i1 )
  if(associated(anode%jac_z0i2 ))  deallocate(anode%jac_z0i2 )
  if(associated(anode%jac_z0i3 ))  deallocate(anode%jac_z0i3 )
  if(associated(anode%jac_z0i4 ))  deallocate(anode%jac_z0i4 )
  if(associated(anode%jac_z0i5 ))  deallocate(anode%jac_z0i5 )
  if(associated(anode%jac_z0i6 ))  deallocate(anode%jac_z0i6 )
  if(associated(anode%jac_z0i7 ))  deallocate(anode%jac_z0i7 )
  if(associated(anode%jac_z0i8 ))  deallocate(anode%jac_z0i8 )
  if(associated(anode%jac_z0i9 ))  deallocate(anode%jac_z0i9 )
  if(associated(anode%jac_z0i10))  deallocate(anode%jac_z0i10)
  if(associated(anode%jac_z0i11))  deallocate(anode%jac_z0i11)
  if(associated(anode%jac_z0i12))  deallocate(anode%jac_z0i12)

  if(associated(anode%jac_vertqi1 )) deallocate(anode%jac_vertqi1 )
  if(associated(anode%jac_vertqi2 )) deallocate(anode%jac_vertqi2 )
  if(associated(anode%jac_vertqi3 )) deallocate(anode%jac_vertqi3 )
  if(associated(anode%jac_vertqi4 )) deallocate(anode%jac_vertqi4 )
  if(associated(anode%jac_vertqi5 )) deallocate(anode%jac_vertqi5 )
  if(associated(anode%jac_vertqi6 )) deallocate(anode%jac_vertqi6 )
  if(associated(anode%jac_vertqi7 )) deallocate(anode%jac_vertqi7 )
  if(associated(anode%jac_vertqi8 )) deallocate(anode%jac_vertqi8 )
  if(associated(anode%jac_vertqi9 )) deallocate(anode%jac_vertqi9 )
  if(associated(anode%jac_vertqi10)) deallocate(anode%jac_vertqi10)
  if(associated(anode%jac_vertqi11)) deallocate(anode%jac_vertqi11)
  if(associated(anode%jac_vertqi12)) deallocate(anode%jac_vertqi12)

  if(associated(anode%jac_vertti1 )) deallocate(anode%jac_vertti1 )
  if(associated(anode%jac_vertti2 )) deallocate(anode%jac_vertti2 )
  if(associated(anode%jac_vertti3 )) deallocate(anode%jac_vertti3 )
  if(associated(anode%jac_vertti4 )) deallocate(anode%jac_vertti4 )
  if(associated(anode%jac_vertti5 )) deallocate(anode%jac_vertti5 )
  if(associated(anode%jac_vertti6 )) deallocate(anode%jac_vertti6 )
  if(associated(anode%jac_vertti7 )) deallocate(anode%jac_vertti7 )
  if(associated(anode%jac_vertti8 )) deallocate(anode%jac_vertti8 )
  if(associated(anode%jac_vertti9 )) deallocate(anode%jac_vertti9 )
  if(associated(anode%jac_vertti10)) deallocate(anode%jac_vertti10)
  if(associated(anode%jac_vertti11)) deallocate(anode%jac_vertti11)
  if(associated(anode%jac_vertti12)) deallocate(anode%jac_vertti12)

  if(associated(anode%jac_zdxi1)) deallocate(anode%jac_zdxi1)
  if(associated(anode%jac_zdxi2)) deallocate(anode%jac_zdxi2)
  if(associated(anode%jac_zdxi3)) deallocate(anode%jac_zdxi3)
  if(associated(anode%jac_zdxi4)) deallocate(anode%jac_zdxi4)

  if(associated(anode%jac_zdyi1)) deallocate(anode%jac_zdyi1)
  if(associated(anode%jac_zdyi2)) deallocate(anode%jac_zdyi2)
  if(associated(anode%jac_zdyi3)) deallocate(anode%jac_zdyi3)
  if(associated(anode%jac_zdyi4)) deallocate(anode%jac_zdyi4)

  if(associated(anode%jac_udxi1)) deallocate(anode%jac_udxi1)
  if(associated(anode%jac_udxi2)) deallocate(anode%jac_udxi2)
  if(associated(anode%jac_udxi3)) deallocate(anode%jac_udxi3)
  if(associated(anode%jac_udxi4)) deallocate(anode%jac_udxi4)

  if(associated(anode%jac_vdyi1)) deallocate(anode%jac_vdyi1)
  if(associated(anode%jac_vdyi2)) deallocate(anode%jac_vdyi2)
  if(associated(anode%jac_vdyi3)) deallocate(anode%jac_vdyi3)
  if(associated(anode%jac_vdyi4)) deallocate(anode%jac_vdyi4)

  if(associated(anode%jac_vert)) deallocate(anode%jac_vert)

  if(associated(anode%jac_sigdoti1)) deallocate(anode%jac_sigdoti1)
  if(associated(anode%jac_sigdoti2)) deallocate(anode%jac_sigdoti2)
  if(associated(anode%jac_sigdoti3)) deallocate(anode%jac_sigdoti3)
  if(associated(anode%jac_sigdoti1)) deallocate(anode%jac_sigdoti4)

  if(associated(anode%jac_qi1 )) deallocate(anode%jac_qi1 )
  if(associated(anode%jac_qi2 )) deallocate(anode%jac_qi2 )
  if(associated(anode%jac_qi3 )) deallocate(anode%jac_qi3 )
  if(associated(anode%jac_qi4 )) deallocate(anode%jac_qi4 )

  if(associated(anode%jac_ti1 )) deallocate(anode%jac_ti1 )
  if(associated(anode%jac_ti2 )) deallocate(anode%jac_ti2 )
  if(associated(anode%jac_ti3 )) deallocate(anode%jac_ti3 )
  if(associated(anode%jac_ti4 )) deallocate(anode%jac_ti4 )

  if(associated(anode%jac_qgmai1)) deallocate(anode%jac_qgmai1)
  if(associated(anode%jac_qgmai2)) deallocate(anode%jac_qgmai2)
  if(associated(anode%jac_qgmai3)) deallocate(anode%jac_qgmai3)
  if(associated(anode%jac_qgmai4)) deallocate(anode%jac_qgmai4)

  if(associated(anode%jac_qgmbi1)) deallocate(anode%jac_qgmbi1)
  if(associated(anode%jac_qgmbi2)) deallocate(anode%jac_qgmbi2)
  if(associated(anode%jac_qgmbi3)) deallocate(anode%jac_qgmbi3)
  if(associated(anode%jac_qgmbi4)) deallocate(anode%jac_qgmbi4)

  if(associated(anode%jac_icei1)) deallocate(anode%jac_icei1)
  if(associated(anode%jac_icei2)) deallocate(anode%jac_icei2)
  if(associated(anode%jac_icei3)) deallocate(anode%jac_icei3)
  if(associated(anode%jac_icei4)) deallocate(anode%jac_icei4)
 
  if(associated(anode%jac_zicei1)) deallocate(anode%jac_zicei1)
  if(associated(anode%jac_zicei2)) deallocate(anode%jac_zicei2)
  if(associated(anode%jac_zicei3)) deallocate(anode%jac_zicei3)
  if(associated(anode%jac_zicei4)) deallocate(anode%jac_zicei4)

  if(associated(anode%kboti1)) deallocate(anode%kboti1)
  if(associated(anode%kboti2)) deallocate(anode%kboti2)
  if(associated(anode%kboti3)) deallocate(anode%kboti3)
  if(associated(anode%kboti4)) deallocate(anode%kboti4)

  if(associated(anode%jac_kverti1)) deallocate(anode%jac_kverti1)
  if(associated(anode%jac_kverti2)) deallocate(anode%jac_kverti2)
  if(associated(anode%jac_kverti3)) deallocate(anode%jac_kverti3)
  if(associated(anode%jac_kverti4)) deallocate(anode%jac_kverti4)

  if(associated(anode%jac_fratei1)) deallocate(anode%jac_fratei1)
  if(associated(anode%jac_fratei2)) deallocate(anode%jac_fratei2)
  if(associated(anode%jac_fratei3)) deallocate(anode%jac_fratei3)
  if(associated(anode%jac_fratei4)) deallocate(anode%jac_fratei4)

  if(associated(anode%jac_wmaxflagi1)) deallocate(anode%jac_wmaxflagi1)
  if(associated(anode%jac_wmaxflagi2)) deallocate(anode%jac_wmaxflagi2)
  if(associated(anode%jac_wmaxflagi3)) deallocate(anode%jac_wmaxflagi3)
  if(associated(anode%jac_wmaxflagi4)) deallocate(anode%jac_wmaxflagi4)
  if(associated(anode%ij            )) deallocate(anode%ij            )    
_EXIT_(myname_)
  return
end subroutine obsnode_clean_

subroutine obsnode_xread_(anode,iunit,istat,diaglookup,skip)
  use m_obsdiagnode, only: obsdiaglookup_locate
  implicit none
  class(lightnode) , intent(inout):: anode
  integer(i_kind)  , intent(in   ):: iunit
  integer(i_kind)  , intent(  out):: istat
  type(obs_diags)  , intent(in   ):: diaglookup
  logical,optional , intent(in   ):: skip

  character(len=*),parameter:: myname_=myname//'::obsnode_xread_'
  logical:: skip_
_ENTRY_(myname_)
  skip_=.false.
  if(present(skip)) skip_=skip

  istat=0
  if(skip_) then
     read(iunit,iostat=istat)
     if(istat/=0) then
        call perr(myname_,'skipping read(%(res,err2,...)), iostat =',istat)
        _EXIT_(myname_)
        return
     endif

  else
     read(iunit,iostat=istat)    anode%res           , &
                                 anode%err2          , &
                                 anode%raterr2       , &
                                 anode%b             , &
                                 anode%pg            , & 
                                 anode%jac_z0i1      , &
                                 anode%jac_z0i2      , &
                                 anode%jac_z0i3      , &
                                 anode%jac_z0i4      , &
                                 anode%jac_z0i5      , &
                                 anode%jac_z0i6      , &
                                 anode%jac_z0i7      , &
                                 anode%jac_z0i8      , &
                                 anode%jac_z0i9      , &
                                 anode%jac_z0i10     , &
                                 anode%jac_z0i11     , &
                                 anode%jac_z0i12     , & 
                                 anode%jac_vertqi1   , & !(   nsig)
                                 anode%jac_vertqi2   , & !(   nsig)
                                 anode%jac_vertqi3   , & !(   nsig)
                                 anode%jac_vertqi4   , & !(   nsig)
                                 anode%jac_vertqi5   , & !(   nsig)
                                 anode%jac_vertqi6   , & !(   nsig)
                                 anode%jac_vertqi7   , & !(   nsig)
                                 anode%jac_vertqi8   , & !(   nsig)
                                 anode%jac_vertqi9   , & !(   nsig)
                                 anode%jac_vertqi10  , & !(   nsig)
                                 anode%jac_vertqi11  , & !(   nsig)
                                 anode%jac_vertqi12  , & !(   nsig)
                                 anode%jac_vertti1   , & !(   nsig)
                                 anode%jac_vertti2   , & !(   nsig)
                                 anode%jac_vertti3   , & !(   nsig)
                                 anode%jac_vertti4   , & !(   nsig)
                                 anode%jac_vertti5   , & !(   nsig)
                                 anode%jac_vertti6   , & !(   nsig)
                                 anode%jac_vertti7   , & !(   nsig)
                                 anode%jac_vertti8   , & !(   nsig)
                                 anode%jac_vertti9   , & !(   nsig)
                                 anode%jac_vertti10  , & !(   nsig)
                                 anode%jac_vertti11  , & !(   nsig)
                                 anode%jac_vertti12  , & !(   nsig)
                                 anode%jac_zdxi1     , & !(   nsig)
                                 anode%jac_zdxi2     , & !(   nsig)
                                 anode%jac_zdxi3     , & !(   nsig)
                                 anode%jac_zdxi4     , & !(   nsig)
                                 anode%jac_zdyi1     , & !(   nsig)
                                 anode%jac_zdyi2     , & !(   nsig)
                                 anode%jac_zdyi3     , & !(   nsig)
                                 anode%jac_zdyi4     , & !(   nsig)
                                 anode%jac_udxi1     , & !(   nsig)
                                 anode%jac_udxi2     , & !(   nsig)
                                 anode%jac_udxi3     , & !(   nsig)
                                 anode%jac_udxi4     , & !(   nsig)
                                 anode%jac_vdyi1     , & !(   nsig)
                                 anode%jac_vdyi2     , & !(   nsig)
                                 anode%jac_vdyi3     , & !(   nsig)
                                 anode%jac_vdyi4     , & !(   nsig)
                                 anode%jac_vert      , & !(   nsig)
                                 anode%jac_sigdoti1  , & !(   nsig)
                                 anode%jac_sigdoti2  , & !(   nsig)
                                 anode%jac_sigdoti3  , & !(   nsig)
                                 anode%jac_sigdoti4  , & !(   nsig)
                                 anode%jac_qi1       , & !(   nsig)
                                 anode%jac_qi2       , & !(   nsig)
                                 anode%jac_qi3       , & !(   nsig)
                                 anode%jac_qi4       , & !(   nsig)
                                 anode%jac_ti1       , & !(   nsig)
                                 anode%jac_ti2       , & !(   nsig)
                                 anode%jac_ti3       , & !(   nsig)
                                 anode%jac_ti4       , & !(   nsig)
                                 anode%jac_qgmai1    , & !(   nsig)
                                 anode%jac_qgmai2    , & !(   nsig)
                                 anode%jac_qgmai3    , & !(   nsig)
                                 anode%jac_qgmai4    , & !(   nsig)
                                 anode%jac_qgmbi1    , & !(   nsig)
                                 anode%jac_qgmbi2    , & !(   nsig)
                                 anode%jac_qgmbi3    , & !(   nsig)
                                 anode%jac_qgmbi4    , & !(   nsig)
                                 anode%jac_icei1     , & !(   nsig)
                                 anode%jac_icei2     , & !(   nsig)
                                 anode%jac_icei3     , & !(   nsig)
                                 anode%jac_icei4     , & !(   nsig)
                                 anode%jac_zicei1    , & !(   nsig)
                                 anode%jac_zicei2    , & !(   nsig)
                                 anode%jac_zicei3    , & !(   nsig)
                                 anode%jac_zicei4    , & !(   nsig)
                                 anode%kboti1        , & 
                                 anode%kboti2        , & 
                                 anode%kboti3        , &
                                 anode%kboti4        , & 
                                 anode%jac_kverti1   , &
                                 anode%jac_kverti2   , &
                                 anode%jac_kverti3   , &
                                 anode%jac_kverti4   , &
                                 anode%jac_fratei1   , &
                                 anode%jac_fratei2   , &
                                 anode%jac_fratei3   , &
                                 anode%jac_fratei4   , &
                                 anode%jac_wmaxflagi1, &
                                 anode%jac_wmaxflagi2, &
                                 anode%jac_wmaxflagi3, &
                                 anode%jac_wmaxflagi4, &
                                 anode%wij           , & !(4)
                                 anode%ij                !(12,nsig)
     if (istat/=0) then
        call perr(myname_,'read(%(res,err2,...)), iostat =',istat)
        _EXIT_(myname_)
        return
     end if

     anode%diags => obsdiaglookup_locate(diaglookup,anode%idv,anode%iob,1)
     if(.not.associated(anode%diags)) then
        call perr(myname_,'obsdiaglookup_locate(), %idv =',anode%idv)
        call perr(myname_,'                        %iob =',anode%iob)
        call  die(myname_)
     endif
  endif
_EXIT_(myname_)
  return
end subroutine obsnode_xread_

subroutine obsnode_xwrite_(anode,junit,jstat)
  implicit none
  class(lightnode),intent(in):: anode
  integer(i_kind),intent(in   ):: junit
  integer(i_kind),intent(  out):: jstat

  character(len=*),parameter:: myname_=myname//'::obsnode_xwrite_'
_ENTRY_(myname_)

  jstat=0
  write(junit,iostat=jstat)     anode%res           , &
                                anode%err2          , &
                                anode%raterr2       , &
                                anode%b             , &
                                anode%pg            , &
                                anode%jac_z0i1      , &
                                anode%jac_z0i2      , &
                                anode%jac_z0i3      , &
                                anode%jac_z0i4      , &
                                anode%jac_z0i5      , &
                                anode%jac_z0i6      , &
                                anode%jac_z0i7      , &
                                anode%jac_z0i8      , &
                                anode%jac_z0i9      , &
                                anode%jac_z0i10     , &
                                anode%jac_z0i11     , &
                                anode%jac_z0i12     , &
                                anode%jac_vertqi1   , & !(   nsig)
                                anode%jac_vertqi2   , & !(   nsig)
                                anode%jac_vertqi3   , & !(   nsig)
                                anode%jac_vertqi4   , & !(   nsig)
                                anode%jac_vertqi5   , & !(   nsig)
                                anode%jac_vertqi6   , & !(   nsig)
                                anode%jac_vertqi7   , & !(   nsig)
                                anode%jac_vertqi8   , & !(   nsig)
                                anode%jac_vertqi9   , & !(   nsig)
                                anode%jac_vertqi10  , & !(   nsig)
                                anode%jac_vertqi11  , & !(   nsig)
                                anode%jac_vertqi12  , & !(   nsig)
                                anode%jac_vertti1   , & !(   nsig)
                                anode%jac_vertti2   , & !(   nsig)
                                anode%jac_vertti3   , & !(   nsig)
                                anode%jac_vertti4   , & !(   nsig)
                                anode%jac_vertti5   , & !(   nsig)
                                anode%jac_vertti6   , & !(   nsig)
                                anode%jac_vertti7   , & !(   nsig)
                                anode%jac_vertti8   , & !(   nsig)
                                anode%jac_vertti9   , & !(   nsig)
                                anode%jac_vertti10  , & !(   nsig)
                                anode%jac_vertti11  , & !(   nsig)
                                anode%jac_vertti12  , & !(   nsig)
                                anode%jac_zdxi1     , & !(   nsig)
                                anode%jac_zdxi2     , & !(   nsig)
                                anode%jac_zdxi3     , & !(   nsig)
                                anode%jac_zdxi4     , & !(   nsig)
                                anode%jac_zdyi1     , & !(   nsig)
                                anode%jac_zdyi2     , & !(   nsig)
                                anode%jac_zdyi3     , & !(   nsig)
                                anode%jac_zdyi4     , & !(   nsig)
                                anode%jac_udxi1     , & !(   nsig)
                                anode%jac_udxi2     , & !(   nsig)
                                anode%jac_udxi3     , & !(   nsig)
                                anode%jac_udxi4     , & !(   nsig)
                                anode%jac_vdyi1     , & !(   nsig)
                                anode%jac_vdyi2     , & !(   nsig)
                                anode%jac_vdyi3     , & !(   nsig)
                                anode%jac_vdyi4     , & !(   nsig)
                                anode%jac_vert      , & !(   nsig)
                                anode%jac_sigdoti1  , & !(   nsig)
                                anode%jac_sigdoti2  , & !(   nsig)
                                anode%jac_sigdoti3  , & !(   nsig)
                                anode%jac_sigdoti4  , & !(   nsig)
                                anode%jac_qi1       , & !(   nsig)
                                anode%jac_qi2       , & !(   nsig)
                                anode%jac_qi3       , & !(   nsig)
                                anode%jac_qi4       , & !(   nsig)
                                anode%jac_ti1       , & !(   nsig)
                                anode%jac_ti2       , & !(   nsig)
                                anode%jac_ti3       , & !(   nsig)
                                anode%jac_ti4       , & !(   nsig)
                                anode%jac_qgmai1    , & !(   nsig)
                                anode%jac_qgmai2    , & !(   nsig)
                                anode%jac_qgmai3    , & !(   nsig)
                                anode%jac_qgmai4    , & !(   nsig)
                                anode%jac_qgmbi1    , & !(   nsig)
                                anode%jac_qgmbi2    , & !(   nsig)
                                anode%jac_qgmbi3    , & !(   nsig)
                                anode%jac_qgmbi4    , & !(   nsig)
                                anode%jac_icei1     , & !(   nsig)
                                anode%jac_icei2     , & !(   nsig)
                                anode%jac_icei3     , & !(   nsig)
                                anode%jac_icei4     , & !(   nsig)
                                anode%jac_zicei1    , & !(   nsig)
                                anode%jac_zicei2    , & !(   nsig)
                                anode%jac_zicei3    , & !(   nsig)
                                anode%jac_zicei4    , & !(   nsig)
                                anode%kboti1        , & 
                                anode%kboti2        , & 
                                anode%kboti3        , & 
                                anode%kboti4        , & 
                                anode%jac_kverti1   , &
                                anode%jac_kverti2   , &
                                anode%jac_kverti3   , &
                                anode%jac_kverti4   , &
                                anode%jac_fratei1   , &
                                anode%jac_fratei2   , &
                                anode%jac_fratei3   , &
                                anode%jac_fratei4   , &
                                anode%jac_wmaxflagi1, &
                                anode%jac_wmaxflagi2, &
                                anode%jac_wmaxflagi3, &
                                anode%jac_wmaxflagi4, &
                                anode%wij           , & !(4)
                                anode%ij                !(12,nsig)
  if (jstat/=0) then
     call perr(myname_,'write(%(res,err2,...)), iostat =',jstat)
     _EXIT_(myname_)
     return
  end if
_EXIT_(myname_)
  return
end subroutine obsnode_xwrite_

subroutine obsnode_sethop_(anode)
  use m_cvgridlookup, only: cvgridlookup_getiw
  use gridmod, only: nsig,latlon11
  implicit none
  class(lightnode),intent(inout):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_sethop_'
  integer(i_kind):: k
_ENTRY_(myname_)

  ASSERT(size(anode%ij,2)==nsig)
  ASSERT(nsig>0)

  call cvgridlookup_getiw(anode%elat,anode%elon,anode%ij(:,1),anode%wij)
  do k=2,nsig
     anode%ij(:,k) = anode%ij(:,1)+(k-1)*latlon11
  enddo
_EXIT_(myname_)
  return
end subroutine obsnode_sethop_

function obsnode_isvalid_(anode) result(isvalid_)
  implicit none
  logical:: isvalid_
  class(lightnode),intent(in):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_isvalid_'
_ENTRY_(myname_)
  isvalid_=associated(anode%diags)
_EXIT_(myname_)
end function obsnode_isvalid_

pure subroutine gettlddp_(anode,jiter,tlddp,nob)
  use kinds, only: r_kind
  implicit none
  class(lightnode), intent(in):: anode
  integer(kind=i_kind),intent(in):: jiter
  real(kind=r_kind),intent(inout):: tlddp
  integer(kind=i_kind),optional,intent(inout):: nob

  tlddp = tlddp + anode%diags%tldepart(jiter)*anode%diags%tldepart(jiter)
  if(present(nob)) nob=nob+1
  return
end subroutine gettlddp_

end module m_lightnode
