module m_abstracttimer
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module m_abstracttimer
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 610.3
!     date:	 2017-06-30
!
! abstract: an abstract multi-timer replacing stub_timermod.f90 with m_stubtimer
!
! program history log:
!   2017-06-30  j guo   - Replaced stub_timermod with this module and module
!                         m_stubtimer, in the same file m_stubtimer.f90.
!                       . With abstracttimer type and stubtimer type, this
!                         implementation is extensible either from abstracttimer
!                         or from default stubtimer.
!
!   input argument list: see Fortran inline document below
!
!   output argument list: see Fortran inline document below

! attributes:
!   language: Fortran 2003/2008 and/or above
!   machine:
!
!$$$  end subprogram documentation block

! module interface:
  use kinds    , only: i_kind
  use mpeu_util, only: tell
  implicit none
  private
  public:: abstracttimer
  public:: abstracttimer_typename

     interface abstracttimer_typename; module procedure typename_; end interface

  type, abstract:: abstracttimer
     private
     contains
     procedure(mytype  ),nopass,deferred:: mytype        ! typename inquiry
     procedure(on      ),       deferred:: on            ! turn on a single named timer
     procedure(off     ),       deferred:: off           ! turn off a single named timer
     procedure(reset   ),       deferred:: reset         ! reset all timers
     procedure(flush   ),       deferred:: flush         ! summerize all local timers
     procedure(allflush),       deferred:: allflush      ! reduce-summarize distributed timers
  end type abstracttimer

  abstract interface
     function mytype() result(type_)
        implicit none
        character(:),allocatable:: type_
     end function mytype
  end interface

  abstract interface
     subroutine on(tm,name)
        import abstracttimer
        implicit none
        class(abstracttimer), intent(inout):: tm
        character(len=*)    , intent(in   ):: name
     end subroutine on
  end interface

  abstract interface
     subroutine off(tm,name)
        import abstracttimer
        implicit none
        class(abstracttimer), intent(inout):: tm
        character(len=*)    , intent(in   ):: name
     end subroutine off
  end interface

  abstract interface
     subroutine reset(tm)
        import abstracttimer
        implicit none
        class(abstracttimer), intent(inout):: tm
     end subroutine reset
  end interface

  abstract interface
     subroutine flush(tm,lu)
        import abstracttimer
        import i_kind
        implicit none
        class(abstracttimer), intent(in):: tm
        integer(kind=i_kind), intent(in):: lu
     end subroutine flush
  end interface

  abstract interface
     subroutine allflush(tm,lu,comm,root)
        import abstracttimer
        import i_kind
        implicit none
        class(abstracttimer), intent(in):: tm
        integer(kind=i_kind), intent(in):: lu
        integer(kind=i_kind), intent(in):: comm
        integer(kind=i_kind), intent(in):: root
     end subroutine allflush
  end interface

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname='m_abstracttimer'

contains
function typename_() result(typename)
!-- Return the type name.
  implicit none
  character(len=:),allocatable:: typename
  typename="[abstracttimer]"
end function typename_

end module m_abstracttimer

module m_stubtimer
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 module m_abstracttimer
!   prgmmr:      todling          org: gmao                date: 2007-10-01
!
! abstract: a do-nothing multi-timer
!
! program history log:
!   2007-10-01  todling - Original stub_timermod
!   2009-02-26  todling - if-def from gmao_fvgsi to geos_pert
!   2009-08-13  lueken - update documentation
!   2010-06-16  guo - separated stub implementation with implicit interfaces
!		      from module implementation with explicit interfaces.
!   2011-08-01  lueken  - replaced F90 with f90 (no machine logic)
!   2017-06-30  j guo   - replaced stub_timermod.f90 with this module and module
!                         m_stubtimer, in the same file m_stubtimer.f90.
!                       . With abstracttimer type and stubtimer type, this
!                         implementation is extensible either from abstracttimer
!                         or from default stubtimer.
!
!   input argument list: see Fortran inline document below
!
!   output argument list: see Fortran inline document below

! attributes:
!   language: Fortran 2003/2008 and/or above
!   machine:
!
!$$$  end subprogram documentation block

  use m_abstracttimer, only: abstracttimer
  use kinds    , only: i_kind
  use mpeu_util, only: tell,die
  implicit none
  private
  public:: timer
  public:: timer_typemold

  type, extends(abstracttimer):: timer
     private
     contains
        ! see m_abstracttimer for more information
     procedure,nopass:: mytype
     procedure:: on
     procedure:: off
     procedure:: reset
     procedure:: flush
     procedure:: allflush
  end type timer

  character(len=*),parameter:: myname  ="m_stubtimer"
  type(timer),target:: typemold_

  logical,parameter:: verbose=.false.
  !logical,parameter:: verbose=.true.

contains

function timer_typemold() result(typemold)
!-- return a mold of timer
  implicit none
  type(timer),pointer:: typemold
  typemold => typemold_
end function timer_typemold

!--------------------------------------------------
! type-bound-procedures.  See type(abstrcttimer) in module
! m_abstracttimer for specifications.
function mytype()
  implicit none
  character(len=:), allocatable:: mytype
  mytype="["//myname//"::timer]"
end function mytype

subroutine on(tm,name)
  implicit none
  class(timer), intent(inout):: tm
  character(len=*), intent(in):: name
  if(verbose) call tell(tm%mytype()//'%on','timer = ',trim(name))
end subroutine on

subroutine off(tm,name)
  implicit none
  class(timer), intent(inout):: tm
  character(len=*), intent(in):: name
  if(verbose) call tell(tm%mytype()//'%off','timer = ',trim(name))
end subroutine off

subroutine reset(tm)
  implicit none
  class(timer), intent(inout):: tm
  if(verbose) call tell(tm%mytype()//'%reset','no action taken')
end subroutine reset

subroutine flush(tm,lu)
  implicit none
  class(timer)        , intent(in):: tm
  integer(kind=i_kind), intent(in):: lu
  if(verbose) call tell(tm%mytype()//'%flush','no action taken, lu =',lu)
end subroutine flush

subroutine allflush(tm,lu,comm,root)
  use mpeu_mpif,only: mpi_ikind
  implicit none
  class(timer)        , intent(in):: tm         ! a handle to this timer
  integer(kind=i_kind), intent(in):: lu         ! output logic unit
  integer(kind=i_kind), intent(in):: comm       ! communicator
  integer(kind=i_kind), intent(in):: root       ! root pe

  character(len=*),parameter:: myname_=myname//'::allflush'
  integer(kind=mpi_ikind):: mype,ier

  call mpi_comm_rank(comm,mype,ier)
  if(ier/=0) call die(myname_,'mpi_comm_rank(), ierror =',ier)
  if(verbose.and.mype==root) call tell(tm%mytype()//'%allflush','no action taken, lu =',lu)
end subroutine allflush

end module m_stubtimer
