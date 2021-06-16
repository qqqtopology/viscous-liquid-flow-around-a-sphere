!******************************************************************************
!  Program: (B)oundary (E)lements (M)ethod
!  Module: Lin_solve_lib
!  Copyright (C): Denis V. Esipov (2009 -- 2010)
!******************************************************************************
module Lin_solve_lib
use RP_lib
use MatVec_lib
implicit none

contains

!> Soubroutine to solve Ax=b problem by simple (Gauss) ellimination method
subroutine Gauss(A,B,X,n)
implicit none

! Input-output variables
real(RP) :: A(:,:), B(:), X(:)
integer :: n

! Local variables
integer :: i, j, k
real(RP) :: c

do i=1,(N-1)

  ! Interchange rows to get non zero diagonal coefficient
  if (A(i,i) == 0._RP) then
    k = i
    do j=i+1,N
      if (A(j,i) /= 0._RP) then
        k = j
        exit
      end if
    end do
    A(i,i:N) = A(i,i:N) + A(k,i:N)
    B(i)     = B(i)     + B(k)
  end if

  ! Divide row by diagonal coefficient
  c = 1._RP/A(i,i)
  A(i,(i+1):N) = A(i,(i+1):N)*c
  B(i)         = B(i)*c

  ! Eliminate unknown X(k) from row i
  do j=(i+1),N
    A(j,(i+1):N) = A(j,(i+1):N) - A(i,(i+1):N)*A(j,i)
    B(j)         = B(j) - B(i)*A(j,i)
    A(j,i) = 0._RP
  end do
end do
B(N)   = B(N)/A(N,N)
A(N,N) = 1._RP

! Backsubstitute to compute remaining unknowns
do i=N,2,-1
  B(1:(i-1)) = B(1:(i-1)) - A(1:(i-1),i)*B(i)
end do

X = B

return
end subroutine Gauss

!> Soubroutine to solve Ax=b problem by (G)auss-(Seidel) (L)eft method
subroutine GSeidelL(A,B,X,n,tol)
implicit none

! Input-output variables
real(RP) :: A(:,:), B(:), X(:)
real(RP) :: tol
integer :: n

! Local variables
integer :: i, j
real(RP) :: c, xp(1:n), xpnorm

X = B
do i=1,(2*n)
  xp = X
  xpnorm = DSqrt(Dot_product(xp(1:n),xp(1:n)))
  do j=1,n
    c = B(j) + A(j,j)*X(j) - Dot_product(A(j,1:n),X(1:n))
    X(j) = c/A(j,j)
  end do
  if (DSqrt(Dot_product(X(1:n) - xp(1:n),X(1:n) - xp(1:n))) < tol*xpnorm) exit
end do

return
end subroutine GSeidelL

!> Soubroutine to solve Ax=b problem by (G)auss-(Seidel) (R)ight method
subroutine GSeidelR(A,B,X,n,tol)
implicit none

! Input-output variables
real(RP) :: A(:,:), B(:), X(:)
real(RP) :: tol
integer :: n

! Local variables
integer :: i, j
real(RP) :: c, xp(1:n), xpnorm

X = B
do i=1,(2*n)
  xp = X
  xpnorm = DSqrt(Dot_product(xp(1:n),xp(1:n)))
  do j=n,1,-1
    c = B(j) + A(j,j)*X(j) - Dot_product(A(j,1:n),X(1:n))
    X(j) = c/A(j,j)
  end do
  if (DSqrt(Dot_product(X(1:n) - xp(1:n),X(1:n) - xp(1:n))) < tol*xpnorm) exit
end do

return
end subroutine GSeidelR

!> Soubroutine to solve Ax=b problem by
!>                                    (G)eneralized (M)inimal (RES)idual method
subroutine GMRESp(A,B,X,n,tol,rho,ni,gemm)
implicit none

! Input-output variables
real(RP) :: A(:,:), B(:), X(:)                                                 ! AX = B
real(RP) :: tol, rho                                                           ! Tolerance
integer :: n, ni                                                               !
integer :: gemm

! Local variables
integer :: iter, iter1, i, j
real(RP) :: prod
real(RP) :: wnorm, bnorm
real(RP) :: q(1:(2*ni)),y(1:(n+1))
real(RP), allocatable :: v(:,:),hes(:,:)
integer :: ios

allocate(v(1:n,1:(ni+1)),hes(1:(ni+1),1:ni),stat=ios)
if (ios /= 0) then
!  call error_msg(20017,err_type_error,'Cannot allocate array!')
  return
end if

prod = 1._RP

! Restarting loop
do j=1,(2*n/ni + 1)
  hes = 0._RP

  ! Compute first residual and first vector v
  if (gemm == 0) then
    call Mult_Av(v(1:n,1),A,X(1:n),n)
  else
    !call dmatvecmul(A,X(1:n),v(1:n,1),n)
  end if
  v(1:n,1) = B(1:n) - v(1:n,1)
  bnorm = Dsqrt(Dot_product(v(1:n,1),v(1:n,1)))
  if (bnorm == 0._RP) return
  v(1:n,1) = v(1:n,1)/bnorm

  ! Iteration loop
  do iter=1,ni
    iter1 = iter + 1

    ! Compute Hessenberg matrix and next vector v
    if (gemm == 0) then
      call Mult_Av(v(1:n,iter1),A,v(1:n,iter),n)
    else
      !call dmatvecmul(A,v(1:n,iter),v(1:n,iter1),n)
    end if
    do i=1,iter
      hes(i,iter) = Dot_product(v(1:n,i),v(1:n,iter1))
      v(1:n,iter1) = v(1:n,iter1) - hes(i,iter)*v(1:n,i)
    end do
    hes(iter1,iter) = DSqrt(Dot_product(v(1:n,iter1),v(1:n,iter1)))
    v(1:n,iter1) = v(1:n,iter1)/hes(iter1,iter)

    ! Using rotations transform Hessenberg matrix to R
    call dhesrot(hes,iter,q,ni)

    prod = prod*q(2*iter)
    rho = Dabs(prod)

    !write(*,*) j, iter, rho
    if (rho <= tol) exit
  end do
  if (iter > ni) iter = iter - 1

  ! Compute vector y
  y = 0._RP
  y(1) = bnorm
  call dhesy(hes,iter,q,y,ni)

  ! Compute the solution
  do i=1,iter
    X(1:n) = X(1:n) + y(i)*v(1:n,i)
  end do

  if (rho <= tol) exit
end do

if (allocated(v)) deallocate(v,stat=ios)
! if (ios /= 0) call error_msg(10006,err_type_warning,'Cannot deallocate array!')
if (allocated(hes)) deallocate(hes,stat=ios)
! if (ios /= 0) call error_msg(10007,err_type_warning,'Cannot deallocate array!')

return
end subroutine GMRESp

!> Soubroutine to make Arnoldi orthogonalization
subroutine dorth(hes,v,iter,iter1,n)
implicit none

! Input-Output variables
integer :: iter, iter1                                                       ! Iteration number
real(RP) :: hes(:,:), v(:,:)                                                 ! Hessenberg matrix and vectors v
integer :: n                                                                 ! Dimension

! Local variables
integer :: i                                                                 ! Loop index

do i=1,iter
  hes(i,iter) = Dot_product(v(1:n,i),v(1:n,iter1))
  v(1:n,iter1) = v(1:n,iter1) - hes(i,iter)*v(1:n,i)
end do
hes(iter1,iter) = DSqrt(Dot_product(v(1:n,iter1),v(1:n,iter1)))
v(1:n,iter1) = v(1:n,iter1)/hes(iter1,iter)

return
end subroutine dorth

!> Soubroutine to
subroutine dhesrot(hes,N,q,ni)
implicit none

integer :: i
integer :: N, ni
real(RP) :: q(1:(2*ni))
real(RP) :: T,T1,T2
real(RP) :: cos, sin
real(RP) :: hes(:,:)

! Main loop of rotations
do i=1,(N-1)
  T1 = hes(i,N)
  T2 = hes(i+1,N)
  cos = q(2*i-1)
  sin = q(2*i)
  hes(i,N) = cos*T1 - sin*T2
  hes(i+1,N) = sin*T1 + cos*T2
end do

! Rotation to the last element
T1 = hes(N,N)
T2 = hes(N+1,N)
if (T2 == 0._RP) then
  cos = 1._RP
  sin = 0._RP
elseif (Dabs(T2) >= Dabs(T1)) then
  T = T1/T2
  sin = - 1._RP/DSqrt(1.0_RP + T**2)
  cos = - sin*T
else
  T = T2/T1
  cos = 1._RP/DSqrt(1.0_RP + T**2)
  sin = - cos*T
end if

! Save rotations
q(2*N-1) = cos
q(2*N)   = sin
hes(N,N) = cos*T1 - sin*T2

return
end subroutine dhesrot

!> Soubroutine to
subroutine dhesy(hes,N,Q,y,ni)
implicit none

integer :: i, N, ni
real(RP) :: hes(:,:)
real(RP) :: y(1:(ni+1)), q(1:(2*ni))
real(RP) :: cos, sin, T, T1, T2

do i=1,N
  cos = q(2*i - 1)
  sin = q(2*i)
  T1 = y(i)
  T2 = y(i+1)
  y(i) = cos*T1 - sin*T2
  y(i+1) = sin*T1 + cos*T2
end do
do i=N,1,-1
  y(i) = y(i)/hes(i,i)
  T = - y(i)
  y(1:(i-1)) = y(1:(i-1)) + T*hes(1:(i-1),i)
end do

return
end subroutine dhesy

!> Subroutine to make LU decomposition
subroutine LU_Decomposition(A,n,IPiv)
implicit none

! Input-Output variables
real(RP), intent(INOUT) :: A(:,:)                                              ! Matrix
integer, intent(IN) :: n                                                       ! Dimension
integer, intent(OUT) :: IPiv(:)                                                ! Pivoting vector

! Local variables
integer :: i, j, k                                                             ! Loop indices
real(RP) :: summ, Amax, A0, v(1:n)                                             ! Some useful variables
integer :: Imax                                                                ! Index of max element

! Find max element in columns
do i=1,n
  Amax = 0._RP
  do j=1,n
    A0 = DAbs(A(i,j))
    if (A0 > Amax) Amax = A0
  end do
  v(i) = 1._RP / Amax
end do

! Foward substitution
do j=1,n

  ! Compute matrix U
  do i=1,j-1
    summ = A(i,j) - sum(A(i,1:i-1)*A(1:i-1,j))
    A(i,j) = summ
  end do

  ! Compute matrix L
  Amax = 0._RP
  do i=j,n
    summ = A(i,j) - sum(A(i,1:j-1)*A(1:j-1,j))
    A(i,j) = summ
    A0 = v(i)*DAbs(summ)
    if (A0 >= Amax) then
      Imax = i
      Amax = A0
    end if
  end do

  ! Change columns
  if (j /= Imax) then
    do k=1,n
      A0 = A(Imax,k)
      A(Imax,k) = A(j,k)
      A(j,k) = A0
    end do
    v(Imax) = v(j)
  end if

  ! Normalize
  IPiv(j) = Imax
  if (j /= n) then
    A0 = 1._RP / A(j,j)
    A(j+1:n,j) = A(j+1:n,j)*A0
  end if
end do

return
end subroutine LU_Decomposition

!> Subroutine to make backward substitution process to LU decomposed matrix
subroutine LU_BackSubstitution(A,n,IPiv,B)
implicit none

! Input-Output variables
real(RP), intent(INOUT) :: A(:,:), B(:)                                        ! Matrix and vector X
integer, intent(IN) :: n                                                       ! Dimension
integer, intent(IN) :: IPiv(:)                                                 ! Pivoting vector

! Local variables
real(RP) :: summ                                                               ! Sum variable
integer :: i, j, k                                                             ! Loop indices

! Back substitution to L
k = 0
do i=1,n
  summ = B(IPiv(i))
  B(IPiv(i)) = B(i)
  if (k /= 0) then
    summ = summ - sum(A(i,k:i-1)*B(k:i-1))
  else if (summ /= 0._RP) then
    k = i
  end if
  B(i) = summ
end do

! Back substitution to U
do i=n,1,-1
  summ = B(i)
  if (i < n) then
    summ = summ - sum(A(i,(i+1):n)*B((i+1):n))
  end if
  B(i) = summ / A(i,i)
end do

return
end subroutine LU_BackSubstitution

!> Subroutine to inverse Matrix A by Gaussian elimination
subroutine Inverse_A(A1,A,n)
implicit none

! Input-output variables
integer :: n                                                                   ! Dimension
real(RP) :: A(:,:), A1(:,:)                                                    ! Main and inverse matrices

! Local variables
integer :: i, j, k                                                             ! Loop indices
real(RP) :: c                                                                  ! Useful constant
real(RP) :: tol                                                                ! Zero tolerance
integer :: i1, j1

! Define tolerance
tol = 0.0000001_RP

! Define identity matrix
A1 = 0._RP
do i=1,n
  A1(i,i) = 1._RP
end do

do i=1,(n-1)

  ! Find row with non-zero pivot
  if (Dabs(A(i,i)) < tol) then
    do j=i+1,n
      if (Dabs(A(j,i)) > tol) then
        k = j
        exit
      end if
    end do

    ! Add to main row
    if (k <= n) then
      A(i,i:n)  = A(i,i:n)  + A(k,i:n)
      A1(i,1:n) = A1(i,1:n) + A1(k,1:n)
    end if
  end if

  ! Devise row to pivot
  c = 1._RP/A(i,i)
  A(i,i:n)  = A(i,i:n)*c
  A1(i,1:n) = A1(i,1:n)*c

  ! Eliminate column
  do j=(i+1),n
    c = A(j,i)
    A(j,i:n)  = A(j,i:n)  - A(i,i:n)*c
    A1(j,1:n) = A1(j,1:n) - A1(i,1:n)*c
  end do
end do

! Last row
c = 1._RP/A(n,n)
A(n,n) = 1._RP
A1(n,1:n) = A1(n,1:n)*c

! Back elimination
do i=n,2,-1
  do j=(i-1),1,-1
    c = A(j,i)
    A1(j,1:n) = A1(j,1:n) - A1(i,1:n)*c
  end do
end do

return
end subroutine Inverse_A


end module Lin_solve_lib
