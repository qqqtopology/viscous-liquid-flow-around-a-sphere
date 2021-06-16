!******************************************************************************
!  Program: (B)oundary (E)lements (M)ethod
!  Module: MatVec_lib
!  Copyright (C): Denis V. Esipov (2009 -- 2010)
!******************************************************************************
module MatVec_lib
use RP_lib
implicit none

contains


subroutine Mult_Av(w,A,v,n)
!******************************************************************************
! Subroutine to multiply Matrix A and Vector v
!******************************************************************************
implicit none

! Input-Output variables
integer :: n
real(RP) :: w(:), A(:,:), v(:)
integer :: i, j

w = 0._RP

do j=1,n
  do i=1,n
    w(i) = w(i) + A(i,j)*v(j)
  end do
end do

return
end subroutine Mult_Av


function Dot_vv(v1,v2,n)
!******************************************************************************
! Function to compute the dot product of two vectors
!******************************************************************************
implicit none

! Input-Output variables
real(RP) :: Dot_vv                                                             ! Function result
real(RP) :: v1(:), v2(:)                                                       ! Vectors
integer :: n                                                                   ! Dimension

! Local variables
real(RP) :: w                                                                  ! Partial sum
integer :: i                                                                   ! Loop index

w = 0._RP

do i=1,n
  w = w + v1(i)*v2(i)
end do

Dot_vv = w

return
end function Dot_vv


subroutine Mult_AA(C,A,B,n)
!******************************************************************************
! Subroutine to multiply Matrix A and Vector v
!******************************************************************************
implicit none

! Input-Output variables
integer :: n
real(RP) :: A(:,:), B(:,:), C(:,:)
integer :: i, j, k

do j=1,n
  do i=1,n
    C(i,j) = 0._RP
    do k=1,n
      C(i,j) = C(i,j) + A(i,k)*B(k,j)
    end do
  end do
end do

return
end subroutine Mult_AA


end module MatVec_lib
