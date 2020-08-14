module m_sortind

!$$$ module documentation block
!           .      .    .                                       .
! module:   m_sortind   finds indices to sort an array in ascending order 
!   prgmmr: eliu
!
! abstract: module to find indices to sort an array in ascending order 
! assimilation
!
! program history log:
!   1996-10-01  Joiner/Karki - initial coding from nasa/gmao
!   2012-02-15  eliu         - reformat to use in gsi
!
! subroutines included:
!
! variable definitions:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

  use kinds,only : i_kind, r_kind
  interface sortind
     module procedure r_sortind
     module procedure i_sortind
  end interface

  contains

  function r_sortind(arr) result(arr2)
  implicit none

  !input parameters:
  real(r_kind), dimension(:) :: arr ! input vector to sort
  !output parameters:
  integer(i_kind), dimension(size(arr)) :: arr2

  call indexx(size(arr),arr, arr2)

  end function r_sortind

  function i_sortind(arr) result(arr2)

  implicit none

  integer(i_kind), dimension(:) :: arr
  integer(i_kind), dimension(size(arr)) :: arr2

  call iindexx(size(arr),arr, arr2)

  end function i_sortind

  subroutine indexx(n,arr,indx)

  implicit none

  integer(i_kind):: n,indx(n),m,nstack
  real(r_kind)   ::  arr(n)
  parameter (m=7,nstack=50)
  integer(i_kind)::  i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
  real(r_kind)   ::  a
   do 11 j=1,n
      indx(j)=j
11 continue
   jstack=0
   l=1
   ir=n
   loop2: do 
      if(ir-l<m)then
         loop: do j=l+1,ir
            indxt=indx(j)
            a=arr(indxt)
            do 12 i=j-1,1,-1
               if(arr(indx(i))<=a)then
                  indx(i+1)=indxt
                  cycle loop 
               end if
               indx(i+1)=indx(i)
12          continue
            indx(1)=indxt
         end do loop
         if(jstack==0)return
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2
         itemp=indx(k)
         indx(k)=indx(l+1)
         indx(l+1)=itemp
         if(arr(indx(l+1)) > arr(indx(ir)))then
            itemp=indx(l+1)
            indx(l+1)=indx(ir)
            indx(ir)=itemp
         endif
         if(arr(indx(l)) > arr(indx(ir)))then
            itemp=indx(l)
            indx(l)=indx(ir)
            indx(ir)=itemp
         endif
         if(arr(indx(l+1)) >  arr(indx(l)))then
            itemp=indx(l+1)
            indx(l+1)=indx(l)
            indx(l)=itemp
         endif
         i=l+1
         j=ir
         indxt=indx(l)
         a=arr(indxt)
         loop1: do 
            i=i+1
            if(arr(indx(i))<a)cycle loop1
            do 
4              continue
               j=j-1
               if(arr(indx(j))<=a)exit
            end do
            if(j<i)exit loop1
            itemp=indx(i)
            indx(i)=indx(j)
            indx(j)=itemp
         end do loop1
5        indx(l)=indx(j)
         indx(j)=indxt
         jstack=jstack+2
         if(jstack>nstack)pause 'NSTACK too small in indexx'
         if(ir-i+1>=j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         endif
      endif
   end do loop2
  end subroutine indexx

  subroutine iindexx(n,arr,indx)

  implicit none

  integer(i_kind):: n,indx(n),m,nstack
  integer(i_kind):: arr(n)
  parameter (m=7,nstack=50)
  integer(i_kind):: i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
  integer(i_kind):: a
   do 11 j=1,n
      indx(j)=j
11 continue
   jstack=0
   l=1
   ir=n
   loop2: do 
      if(ir-l<m)then
         loop: do j=l+1,ir
            indxt=indx(j)
            a=arr(indxt)
            do 12 i=j-1,1,-1
               if(arr(indx(i))<=a)then
                  indx(i+1)=indxt
                  cycle loop
               end if
               indx(i+1)=indx(i)
12          continue
            indx(1)=indxt
         end do loop
         if(jstack==0)return
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2
         itemp=indx(k)
         indx(k)=indx(l+1)
         indx(l+1)=itemp
         if(arr(indx(l+1))>arr(indx(ir)))then
            itemp=indx(l+1)
            indx(l+1)=indx(ir)
            indx(ir)=itemp
         endif
         if(arr(indx(l))>arr(indx(ir)))then
            itemp=indx(l)
            indx(l)=indx(ir)
            indx(ir)=itemp
         endif
         if(arr(indx(l+1))>arr(indx(l)))then
            itemp=indx(l+1)
            indx(l+1)=indx(l)
            indx(l)=itemp
         endif
         i=l+1
         j=ir
         indxt=indx(l)
         a=arr(indxt)
         loop1: do 
            do 
               i=i+1
               if(arr(indx(i))>=a)exit
            end do
            do 
               j=j-1
               if(arr(indx(j))<=a)exit
            end do
            if(j<i)exit loop1
            itemp=indx(i)
            indx(i)=indx(j)
            indx(j)=itemp
         end do loop1
         indx(l)=indx(j)
         indx(j)=indxt
         jstack=jstack+2
         if(jstack>nstack)pause 'NSTACK too small in indexx'
         if(ir-i+1>=j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         endif
      endif
   end do loop2
  end subroutine iindexx

end module m_sortind
