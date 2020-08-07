module m_radnode
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module m_radnode
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2016-05-18
!
! abstract: class-module of obs-type radnode (radiances)
!
! program history log:
!   2016-05-18  j guo   - added this document block for the initial polymorphic
!                         implementation.
!   2016-07-19  kbathmann - add rsqrtinv and use_corr_obs to rad_ob_type
!   2019-04-22  kbathmann - change rsqrtinv to rpred
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
  use m_obsdiagnode, only: obs_diag,aofp_obs_diag => fptr_obsdiagnode
  use m_obsdiagnode, only: obs_diags
  use kinds , only: i_kind,r_kind
  use mpeu_util, only: assert_,die,perr,warn,tell
  use m_obsnode, only: obsnode
  implicit none
  private

  public:: radnode

  type,extends(obsnode):: radnode
     type(aofp_obs_diag), dimension(:), pointer :: diags => null()
     real(r_kind),dimension(:),pointer :: res => null()
                                      !  obs-guess residual (nchan)
     real(r_kind),dimension(:),pointer :: err2 => null()
                                      !  error variances squared (nchan)
     real(r_kind),dimension(:),pointer :: raterr2 => null()
                                      !  ratio of error variances squared (nchan)
     real(r_kind)    :: wij(4)        !  horizontal interpolation weights
     real(r_kind),dimension(:,:),pointer :: pred => null()
                                      !  predictors (npred,nchan)
     real(r_kind),dimension(:,:),pointer :: dtb_dvar => null()
                                      !  radiance jacobian (nsigradjac,nchan)

     real(r_kind),dimension(:,:),pointer :: rpred => null()
                                      !  square root of inverse of Rr multiplied
                                      !  by bias predictor jacobian
                                      !  only used if using correlated obs
     real(r_kind),dimension(:  ),pointer :: rsqrtinv => null()
                                      !  square root of inverse of r, only used
                                      !  if using correlated obs

     integer(i_kind),dimension(:),pointer :: icx => null()
     integer(i_kind),dimension(:),pointer :: ich => null()
     integer(i_kind) :: nchan         !  number of channels for this profile
     integer(i_kind) :: ij(4)         !  horizontal locations
     logical         :: use_corr_obs  = .false. !  to indicate if correlated obs is implemented
     integer(i_kind) :: iuse_predoper_type = 0 !  indicate which type of correlated predictor operator is implemented
                                        ! 0 uses none (diagonal)
                                        ! 1 uses rsqrtinv
                                        ! 2 uses rpred

!!! Is %isis or %isfctype ever being assigned somewhere in the code?
!!! They are used in intrad().
!!!
!!! Now, they are not written to an obsdiags file, nor read from one.

     character(20) :: isis            ! sensor/instrument/satellite id, e.g. amsua_n15
     !integer(i_kind) :: isfctype      ! surf mask: ocean=0,land=1,ice=2,snow=3,mixed=4
     character(80) :: covtype      ! surf mask: ocean=0,land=1,ice=2,snow=3,mixed=4
  contains
     procedure,nopass::  mytype
     procedure::  sethop => obsnode_sethop_
     procedure::   xread => obsnode_xread_
     procedure::  xwrite => obsnode_xwrite_
     procedure:: isvalid => obsnode_isvalid_
     procedure::  gettlddp => gettlddp_

     procedure, nopass:: headerread  => obsheader_read_
     procedure, nopass:: headerwrite => obsheader_write_
     ! procedure:: init  => obsnode_init_
     procedure:: clean => obsnode_clean_
  end type radnode

  public:: radnode_typecast
  public:: radnode_nextcast
     interface radnode_typecast; module procedure typecast_ ; end interface
     interface radnode_nextcast; module procedure nextcast_ ; end interface

  public:: radnode_appendto
     interface radnode_appendto; module procedure appendto_ ; end interface

  character(len=*),parameter:: myname="m_radnode"

#include "myassert.H"
#include "mytrace.H"
contains
function typecast_(anode) result(ptr_)
!-- cast a class(obsnode) to a type(radnode)
  use m_obsnode, only: obsnode
  implicit none
  type(radnode ),pointer:: ptr_
  class(obsnode),pointer,intent(in):: anode
  ptr_ => null()
  if(.not.associated(anode)) return
        ! logically, typecast of a null-reference is a null pointer.
  select type(anode)
     type is(radnode)
        ptr_ => anode
  end select
  return
end function typecast_

function nextcast_(anode) result(ptr_)
!-- cast an obsnode_next(obsnode) to a type(radnode)
  use m_obsnode, only: obsnode,obsnode_next
  implicit none
  type(radnode ),pointer:: ptr_
  class(obsnode),target ,intent(in):: anode

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
  type(radnode),pointer,intent(in):: anode
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
  mytype="[radnode]"
end function mytype

subroutine obsheader_read_(iunit,mobs,jread,istat)
  use radinfo, only: npred,nsigradjac
  implicit none
  integer(i_kind),intent(in ):: iunit
  integer(i_kind),intent(out):: mobs
  integer(i_kind),intent(out):: jread
  integer(i_kind),intent(out):: istat
  
  character(len=*),parameter:: myname_=myname//'.obsheader_read_'
  integer(i_kind):: mpred,msigradjac
_ENTRY_(myname_)
  
  read(iunit,iostat=istat) mobs,jread, mpred,msigradjac
  if(istat==0 .and. (npred/=mpred .or. nsigradjac/=msigradjac)) then
     call perr(myname_,'unmatched dimension information, npred or nsigradjac')
     if(npred/=mpred) then
        call perr(myname_,'     expecting npred =',npred)
        call perr(myname_,'      but read mpred =',mpred)
     endif
     if(nsigradjac/=msigradjac) then
        call perr(myname_,'expecting nsigradjac =',nsigradjac)
        call perr(myname_,' but read msigradjac =',msigradjac)
     endif
     call die(myname_)
  endif
_EXIT_(myname_)
  return
end subroutine obsheader_read_

subroutine obsheader_write_(junit,mobs,jwrite,jstat)
  use radinfo, only: npred,nsigradjac
  implicit none
  integer(i_kind),intent(in ):: junit
  integer(i_kind),intent(in ):: mobs
  integer(i_kind),intent(in ):: jwrite
  integer(i_kind),intent(out):: jstat
  
  character(len=*),parameter:: myname_=myname//'.obsheader_write_'
_ENTRY_(myname_)
  write(junit,iostat=jstat) mobs,jwrite, npred,nsigradjac
_EXIT_(myname_)
  return
end subroutine obsheader_write_

subroutine obsnode_clean_(anode)
  implicit none
  class(radnode),intent(inout):: anode

  character(len=*),parameter:: myname_=myname//'.obsnode_clean_'
_ENTRY_(myname_)
!_TRACEV_(myname_,'%mytype() =',anode%mytype())
  if(associated(anode%diags   )) deallocate(anode%diags   )
  if(associated(anode%ich     )) deallocate(anode%ich     )
  if(associated(anode%res     )) deallocate(anode%res     )
  if(associated(anode%err2    )) deallocate(anode%err2    )
  if(associated(anode%raterr2 )) deallocate(anode%raterr2 )
  if(associated(anode%pred    )) deallocate(anode%pred    )
  if(associated(anode%dtb_dvar)) deallocate(anode%dtb_dvar)
  if(associated(anode%rpred   )) deallocate(anode%rpred   )
  if(associated(anode%rsqrtinv)) deallocate(anode%rsqrtinv)
  if(associated(anode%icx     )) deallocate(anode%icx     )
_EXIT_(myname_)
  return
end subroutine obsnode_clean_

subroutine obsnode_xread_(anode,iunit,istat,diaglookup,skip)
  use m_obsdiagnode, only: obsdiaglookup_locate
  use radinfo, only: npred,nsigradjac
  implicit none
  class(radnode),intent(inout):: anode
  integer(i_kind),intent(in   ):: iunit
  integer(i_kind),intent(  out):: istat
  type(obs_diags),intent(in   ):: diaglookup
  logical,optional,intent(in   ):: skip

  character(len=*),parameter:: myname_=myname//'.obsnode_xread_'
  integer(i_kind):: k,nchan
  logical:: skip_
_ENTRY_(myname_)
  skip_=.false.
  if(present(skip)) skip_=skip

  istat=0
  if(skip_) then
     read(iunit,iostat=istat)
     if (istat/=0) then
        call perr(myname_,'skipping read(%(nchan,iuse_predoper_type)), iostat =',istat)
        _EXIT_(myname_)
        return
     end if

     read(iunit,iostat=istat)
     if(istat/=0) then
        call perr(myname_,'skipping read(%(res,err2,...)), iostat =',istat)
        _EXIT_(myname_)
        return
     endif

     read(iunit,iostat=istat)
     if(istat/=0) then
        call perr(myname_,'skipping read(%(rpred||rsqrtinv)), iostat =',istat)
        _EXIT_(myname_)
        return
     endif

  else
     read(iunit,iostat=istat) anode%nchan,anode%use_corr_obs,anode%iuse_predoper_type
     if (istat/=0) then
        call perr(myname_,'read(%(nchan,use_corr_obs,iuse_predoper_type)), iostat =',istat)
        _EXIT_(myname_)
        return
     end if

     if(associated(anode%diags   )) deallocate(anode%diags   )
     if(associated(anode%ich     )) deallocate(anode%ich     )
     if(associated(anode%res     )) deallocate(anode%res     )
     if(associated(anode%err2    )) deallocate(anode%err2    )
     if(associated(anode%raterr2 )) deallocate(anode%raterr2 )
     if(associated(anode%pred    )) deallocate(anode%pred    )
     if(associated(anode%dtb_dvar)) deallocate(anode%dtb_dvar)
     if(associated(anode%rpred   )) deallocate(anode%rpred)
     if(associated(anode%rsqrtinv)) deallocate(anode%rsqrtinv)
     if(associated(anode%icx     )) deallocate(anode%icx     )

     nchan=anode%nchan
     allocate( anode%diags(nchan), &
               anode%res  (nchan), &
               anode%err2 (nchan), &
               anode%raterr2            (nchan), &
               anode%pred         (npred,nchan), &
               anode%dtb_dvar(nsigradjac,nchan), &
               anode%ich  (nchan), &
               anode%icx  (nchan)  )

     read(iunit,iostat=istat)    anode%ich     , &
                                 anode%res     , &
                                 anode%err2    , &
                                 anode%raterr2 , &
                                 anode%pred    , &
                                 anode%icx     , &
                                 anode%dtb_dvar, &
                                 anode%wij     , &
                                 anode%ij
     if (istat/=0) then
        call perr(myname_,'read(%(res,err2,...)), iostat =',istat)
        _EXIT_(myname_)
        return
     end if

     if(.not.anode%use_corr_obs) anode%iuse_predoper_type=0
     select case(anode%iuse_predoper_type)
        case(1)
           allocate(anode%rsqrtinv(((nchan+1)*nchan)/2))
           read(iunit,iostat=istat) anode%rsqrtinv
           if (istat/=0) then
              call perr(myname_,'read(%rsqrtinv), iostat =',istat)
              _EXIT_(myname_)
              return
           end if
        case(2)
           allocate(anode%rpred(((nchan+1)*nchan)/2,npred))
           read(iunit,iostat=istat) anode%rpred
           if (istat/=0) then
              call perr(myname_,'read(%rpred), iostat =',istat)
              _EXIT_(myname_)
              return
           end if

        case default
           read(iunit,iostat=istat)
     end select

     do k=1,nchan
        anode%diags(k)%ptr => obsdiaglookup_locate(diaglookup,anode%idv,anode%iob,anode%ich(k))
        if(.not.associated(anode%diags(k)%ptr)) then
           call perr(myname_,'obsdiaglookup_locate(k), k =',k)
           call perr(myname_,'                      %idv =',anode%idv)
           call perr(myname_,'                      %iob =',anode%iob)
           call perr(myname_,'                   %ich(k) =',anode%ich(k))
           call  die(myname_)
        endif
     enddo
  endif
_EXIT_(myname_)
  return
end subroutine obsnode_xread_

subroutine obsnode_xwrite_(anode,junit,jstat)
  use radinfo, only: npred
  implicit none
  class(radnode),intent(in):: anode
  integer(i_kind),intent(in   ):: junit
  integer(i_kind),intent(  out):: jstat

  character(len=*),parameter:: myname_=myname//'.obsnode_xwrite_'
  integer(i_kind):: k
  integer(i_kind):: iuse_predoper_type
_ENTRY_(myname_)

  jstat=0
  iuse_predoper_type=0
  if(anode%use_corr_obs) iuse_predoper_type=anode%iuse_predoper_type
  write(junit,iostat=jstat) anode%nchan,anode%use_corr_obs,iuse_predoper_type
  if (jstat/=0) then
     call perr(myname_,'write(%(nchan,use_corr_obs, etc.)), iostat =',jstat)
     _EXIT_(myname_)
     return
  end if

  write(junit,iostat=jstat) (/ (anode%ich(k),k=1,anode%nchan) /), &
                                anode%res     , &
                                anode%err2    , &
                                anode%raterr2 , &
                                anode%pred    , &
                                anode%icx     , &
                                anode%dtb_dvar, &
                                anode%wij     , &
                                anode%ij
  if (jstat/=0) then
     call perr(myname_,'write(%(ich,res,err2,...)), iostat =',jstat)
     _EXIT_(myname_)
     return
  end if

  select case(iuse_predoper_type)
     case(1)
        ASSERT(size(anode%rsqrtinv)==((anode%nchan+1)*anode%nchan)/2)
        write(junit,iostat=jstat) anode%rsqrtinv
        if (jstat/=0) then
           call perr(myname_,'write(%rsqrtinv), iostat =',jstat)
           _EXIT_(myname_)
           return
        end if
     case(2)
        ASSERT(size(anode%rpred,1)==((anode%nchan+1)*anode%nchan)/2)
        ASSERT(size(anode%rpred,2)==npred)
        write(junit,iostat=jstat) anode%rpred
        if (jstat/=0) then
           call perr(myname_,'write(%rpred), iostat =',jstat)
           _EXIT_(myname_)
           return
        end if

     case default
        write(junit,iostat=jstat)
        if (jstat/=0) then
           call perr(myname_,'write as skip record, iostat =',jstat)
           _EXIT_(myname_)
           return
        end if
  end select

_EXIT_(myname_)
  return
end subroutine obsnode_xwrite_

subroutine obsnode_sethop_(anode)
  use m_cvgridlookup, only: cvgridlookup_getiw
  implicit none
  class(radnode),intent(inout):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_sethop_'
_ENTRY_(myname_)
  call cvgridlookup_getiw(anode%elat,anode%elon,anode%ij,anode%wij)
_EXIT_(myname_)
  return
end subroutine obsnode_sethop_

function obsnode_isvalid_(anode) result(isvalid_)
  implicit none
  logical:: isvalid_
  class(radnode),intent(in):: anode

  character(len=*),parameter:: myname_=myname//'::obsnode_isvalid_'
  integer(i_kind):: k
_ENTRY_(myname_)
  isvalid_=all( (/ (associated(anode%diags(k)%ptr),k=1,anode%nchan) /) )
_EXIT_(myname_)
  return
end function obsnode_isvalid_

pure subroutine gettlddp_(anode,jiter,tlddp,nob)
  use kinds, only: r_kind
  implicit none
  class(radnode), intent(in):: anode
  integer(kind=i_kind),intent(in):: jiter
  real(kind=r_kind),intent(inout):: tlddp
  integer(kind=i_kind),optional,intent(inout):: nob

  integer(kind=i_kind):: k
  do k=1,anode%nchan
    tlddp = tlddp + anode%diags(k)%ptr%tldepart(jiter)*anode%diags(k)%ptr%tldepart(jiter)
  enddo
  if(present(nob)) nob=nob+anode%nchan
  return
end subroutine gettlddp_

end module m_radnode
