!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Module to implement subroutines to solve boundary value problem for Helmholtz
! equation in the rectangular domain:
! - div(grad p) + q*p = rhs
! The equation is approximated by the second order difference schema:
! (q*I - Lambda_xx - Lambda_yy)p = rhs
! and obtained SLAE is solved by the Successive over-relaxation (SOR) method or
! subroutines from MKL (Pardiso solver). Poisson subroutine from MKL has order
! and type restrictions in approximation of boundary conditions.
! To start main subroutine 6 arrayes should be set:
!   p(0:nx+1,0:ny+1) - initial guess
!   ax(0:ny+1,1:3), bx(0:ny+1,1:3) - boundary conditions (left and right)
!   ay(0:nx+1,1:3), by(0:nx+1,1:3) - boundary conditions (bottom and top)
!   rhs(0:nx+1,0:ny+1) - right hand side function
! The boundary conditions realize in the following manner:
!   p(0,i)    + ax(i,1)*p(1,i)  + ax(i,2)*p(2,i)    = ax(i,3)
!   p(nx+1,i) + bx(i,1)*p(nx,i) + bx(i,2)*p(nx-1,i) = bx(i,3)
!   p(i,0)    + ay(i,1)*p(i,1)  + ay(i,2)*p(i,2)    = ay(i,3)
!   p(i,ny+1) + by(i,1)*p(i,ny) + by(i,2)*p(i,ny-1) = by(i,3)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program: ibm2d_cyl -- SIMPLE + IBM method for incomressible viscous Navier -
! Stockes equations in 2D statement for flow past a cylinder
! Copyright (C): Denis V. Esipov, Grigorii R. Loskutov (2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module helmholtz
use types
use strings
use mkl_poisson
use mkl_pardiso
implicit none

  public helmholtz_solve, helmholtz_solve_sor, helmholtz_solve_pardiso
  private helmholtz_solve_mkl
  private ipar, dpar, xhandle

  ! Variables to start MKL routines
  integer :: ipar(1:128)
  real(8), allocatable :: dpar(:)
  type(DFTI_DESCRIPTOR), pointer :: xhandle

contains

!> subroutine to solve Helmholtz
subroutine helmholtz_solve(p,ax,bx,ay,by,rhs,bc,q,nx,ny,hx,hy,nit,eps,deb)
implicit none

real(DP), allocatable :: p(:,:), ax(:,:), bx(:,:), ay(:,:), by(:,:), rhs(:,:)   ! Arrayes of variables, boundary conditions and right hand sides
character(4) :: bc                                                              ! Types of booundary conditions
real(DP) :: q                                                                   ! Coefficient of u term
integer :: nx, ny                                                               ! Numbers of intervals in x and y directions
real(DP):: hx, hy                                                               ! Space steps
integer :: nit                                                                  ! Number of iterations
real(DP) :: eps                                                                 ! Value for stop condition
integer :: deb                                                                  ! Debug flag

real(DP), allocatable :: ax_mkl(:), bx_mkl(:), ay_mkl(:), by_mkl(:)             ! Arrayes of values of boundary conditions for MKL solver
real(DP), allocatable :: rhs_mkl(:,:)                                           ! Right hand side for MKL solver
real(DP) :: pc                                                                  ! Correction constant

  if (mkl == 0) then ! SOR method
    call helmholtz_solve_sor(p,ax,bx,ay,by,rhs,q,nx,ny,hx,hy,nit,eps,deb)
  else if (mkl == 1) then ! MKL solver with first order approximation of boundary conditions

    ! Allocate arrayes
    allocate(ax_mkl(0:ny+1))
    allocate(bx_mkl(0:ny+1))
    allocate(ay_mkl(0:nx+1))
    allocate(by_mkl(0:nx+1))

    ! Treat boundary conditions
    ax_mkl(0:ny+1) = ax(0:ny+1,3)
    bx_mkl(0:ny+1) = bx(0:ny+1,3)
    ay_mkl(0:nx+1) = ay(0:nx+1,3)
    by_mkl(0:nx+1) = by(0:nx+1,3)

    ! Call main subroutine
    call helmholtz_solve_mkl(p,ax_mkl,bx_mkl,ay_mkl,by_mkl,rhs,bc,q,nx,ny,hx,hy)

    ! Deallocate arrayes
    deallocate(ax_mkl)
    deallocate(bx_mkl)
    deallocate(ay_mkl)
    deallocate(by_mkl)

  else if (mkl == 2) then ! Hybrid variant MKL+SOR with second order approximation of boundary conditions

    ! Allocate arrayes
    allocate(ax_mkl(0:ny+1))
    allocate(bx_mkl(0:ny+1))
    allocate(ay_mkl(0:nx+1))
    allocate(by_mkl(0:nx+1))
    allocate(rhs_mkl(0:nx+1,0:ny+1))

    ! Treat boundary conditions
    ax_mkl(0:ny+1) = ax(0:ny+1,3)
    bx_mkl(0:ny+1) = bx(0:ny+1,3)
    ay_mkl(0:nx+1) = ay(0:nx+1,3)
    by_mkl(0:nx+1) = by(0:nx+1,3)

    ! Copy rhs array
    rhs_mkl = rhs

    ! Call main subroutine
    call helmholtz_solve_mkl(p,ax_mkl,bx_mkl,ay_mkl,by_mkl,rhs_mkl,bc,q,nx,ny,hx,hy)

    ! Deallocate arrayes
    deallocate(ax_mkl)
    deallocate(bx_mkl)
    deallocate(ay_mkl)
    deallocate(by_mkl)
    deallocate(rhs_mkl)

    ! Apply correction and improve solution by SOR method
    pc = -.5_DP*p(nx,ny/2)
    p = p + pc
    call helmholtz_solve_sor(p,ax,bx,ay,by,rhs,q,nx,ny,hx,hy,nit,eps,deb)

  else if (mkl == 3) then !Pardiso
    call helmholtz_solve_pardiso(p,ax,bx,ay,by,rhs,q,nx,ny,hx,hy)
  end if

end subroutine helmholtz_solve

!> subroutine to solve Helmholtz problem by the Successive
!                                                   over-relaxation (SOR) method
subroutine helmholtz_solve_sor(p,ax,bx,ay,by,rhs,q,nx,ny,hx,hy,nit,eps,deb,omg)
implicit none

real(DP), allocatable :: p(:,:), ax(:,:), bx(:,:), ay(:,:), by(:,:), rhs(:,:)   ! Arrayes of variables, boundary conditions and right hand sides
real(DP) :: q                                                                   ! Coefficient of u term
integer :: nx, ny                                                               ! Numbers of intervals in x and y directions
real(DP):: hx, hy                                                               ! Space steps
integer :: nit                                                                  ! Number of iterations
real(DP) :: eps                                                                 ! Value for stop condition
integer :: deb                                                                  ! Debug flag
real(DP), optional :: omg                                                       ! Relaxation parameter

integer :: i, j, k                                                              ! Loop indices
real(DP) :: st, pmax, c                                                         ! Stop condition, maximal value and usefull variable
real(DP) :: omega, omega1                                                       ! Relaxation parameters
real(DP) :: rxx1, ryy1, rp                                                      ! Step ratios

  ! Check the arrayes boundary
  if (deb > 0) then
    if (lbound(p,1) > 0 .or. ubound(p,1) < nx+1 .or. lbound(p,2) > 0 .or. ubound(p,2) < ny+1 .or. & 
        lbound(rhs,1) > 0 .or. ubound(rhs,1) < nx+1 .or. lbound(rhs,2) > 0 .or. ubound(rhs,2) < ny+1 .or. & 
        lbound(ax,1) > 0 .or. ubound(ax,1) < ny+1 .or. lbound(ay,1) > 0 .or. ubound(ay,1) < nx+1 .or. & 
        lbound(bx,1) > 0 .or. ubound(bx,1) < ny+1 .or. lbound(by,1) > 0 .or. ubound(by,1) < nx+1) then
      write(*,'(A)') 'Boundaries of arrayes are inconsistent in subroutine helmholtz_solve_sor!'
      stop
    end if
  end if

  ! Set relaxation parameters
  if (present(omg)) then
    omega = omg
  else
    omega = 2._DP/(1._DP + sin(pi/min(nx,ny))) ! Close to optimal value
    ! See page 135 in "Khakimzyanov G.S. Cherny S.G. Method of computations.
    ! Part 3. Numerical methods to solve problems for parabolic and elliptic
    ! problems. Novosibirsk: NSU, 2007. 160 p.".
  end if
  omega1 = 1._DP - omega

  ! Compute ratios
  rxx1 = 1._DP/(hx*hx)
  ryy1 = 1._DP/(hy*hy)
  rp   = 1._DP/(2._DP*(rxx1 + ryy1) + q)

  ! Main loop of SOR iterations
  do k = 1,nit

    ! Apply boundary conditions
    do i=0,ny+1
      p(0,i)    = ax(i,3) - ax(i,1)*p(1,i)  - ax(i,2)*p(2,i)
      p(nx+1,i) = bx(i,3) - bx(i,1)*p(nx,i) - bx(i,2)*p(nx-1,i)
    end do
    do i=0,nx+1
      p(i,0)    = ay(i,3) - ay(i,1)*p(i,1)  - ay(i,2)*p(i,2)
      p(i,ny+1) = by(i,3) - by(i,1)*p(i,ny) - by(i,2)*p(i,ny-1)
    end do

    ! One SOR iteration
    st = 0._DP
    pmax = 0._DP

    do i=1,nx
      do j=1,ny
        c = (rxx1*(p(i-1,j) + p(i+1,j)) + ryy1*(p(i,j-1) + p(i,j+1)) + rhs(i,j))*rp
        st = max(st,abs(p(i,j) - c))
        p(i,j) = omega1*p(i,j) + omega*c
        pmax = max(pmax,abs(p(i,j)))
      end do
    end do
    if (pmax < eps) pmax = 1._DP ! To avoid divide by zero

    !  Check stop condition
    if (st*omega/pmax < eps) then
      exit
    end if

    ! Write history of convergence
    if (deb > 1) then
      write(*,'(A,E12.5)') '  SOR - dp^'// trim(int2str(k)) //' = ', st*omega/pmax
    end if

  end do

  ! Write message on the screen if at the same time iteration process is stopped and stop condition is not fullfilled
  if (k == nit + 1) then
    write(*,'(A)') 'Stop condition in the Successive over-relaxation (SOR) method is not fulfilled!'
    write(*,'(A,E12.5)') ' ||p^'// trim(int2str(nit)) //' - p^'// trim(int2str(nit-1)) //'|| = ', st
  end if

  ! Check the solution
  if (deb > 0) then
    do i=1,nx
      do j=1,ny
        c = rxx1*(p(i-1,j) - 2._DP*p(i,j) + p(i+1,j)) + ryy1*(p(i,j-1) - 2._DP*p(i,j) + p(i,j+1)) - q*p(i,j) + rhs(i,j)
        st = max(st,abs(c))
      end do
    end do
    do i=1,ny
      c = p(0,i) + ax(i,1)*p(1,i) + ax(i,2)*p(2,i) - ax(i,3)
      st = max(st,abs(c))
      c = p(nx+1,i) + bx(i,1)*p(nx,i) + bx(i,2)*p(nx-1,i) - bx(i,3)
      st = max(st,abs(c))
    end do
    do i=1,nx
      c = p(i,0) + ay(i,1)*p(i,1) + ay(i,2)*p(i,2) - ay(i,3)
      st = max(st,abs(c))
      c = p(i,ny+1) + by(i,1)*p(i,ny) + by(i,2)*p(i,ny-1) - by(i,3)
      st = max(st,abs(c))
    end do
    write(*,'(A,E12.5)') ' Residual of the SOR method: res^'// trim(int2str(k)) // ' = ', st
  end if

end subroutine helmholtz_solve_sor

!> subroutine to solve Helmholtz problem by Pardiso subroutine
subroutine helmholtz_solve_pardiso(p,ax,bx,ay,by,rhs,q,nx,ny,hx,hy)
implicit none

real(DP), allocatable :: p(:,:), ax(:,:), bx(:,:), ay(:,:), by(:,:), rhs(:,:)   ! Arrayes of variables, boundary conditions and right hand sides
real(DP) :: q                                                                   ! Coefficient of u term
integer :: nx, ny                                                               ! Numbers of intervals in x and y directions
real(DP):: hx, hy                                                               ! Space steps

integer :: i, j, k, l, s                                                        ! Loop indices
real(DP) :: rxx1, ryy1, rxy                                                     ! Step ratios
real(DP), allocatable :: a(:), f(:), x(:)                                       ! a is the matrix in CSR format; f is right hand side vector; x is the computed solution
integer, allocatable :: ja(:), ia(:)                                            ! ja carries the corresponding column indices; ia indicate where every row starts
integer :: n                                                                    ! Number of unknowns
integer, allocatable :: iparm(:), perm(:)                                       ! Job arrayes for pardiso subroutine
integer :: error                                                                ! Job number for pardiso subroutine
type(MKL_PARDISO_HANDLE), allocatable :: pt(:)                                  ! MKL Pointer

  ! Compute ratios
  rxx1 = -1._DP/(hx*hx)
  ryy1 = -1._DP/(hy*hy)
  rxy  =  2._DP/(hx*hx) + 2._DP/(hy*hy)

  !Allocate arrayes for SLAE
  n = nx*ny
  allocate(a(5*n))
  allocate(ja(5*n))
  allocate(ia(n+1))
  allocate(f(n))
  allocate(x(n))
  allocate(perm(n))
  allocate(iparm(64))
  allocate(pt(64))

  ! Fill a matrix in CSR format (a,ia,ja) corresponding to a constant coefficient
  ! five-point stencil
  k = 0
  l = 0
  ia(1) = 1
  do i = 1,nx
    do j = 1,ny
      k = k + 1
      l = l + 1
      s = l
      a(l)  = rxy + q
      ja(l) = k
      f(k)  = rhs(i,j)
      if (j < ny) then
        l = l + 1
        a(l) = ryy1
        ja(l) = k + 1
      else
        a(s)   = a(s)   - ryy1*by(i,1)
        f(k)   = f(k)   - ryy1*by(i,3)
        a(s-1) = a(s-1) - ryy1*by(i,2)
      end if
      if (i < nx) then
        l = l + 1
        a(l) = rxx1
        ja(l) = k + ny
      else
        a(s)   = a(s)   - rxx1*bx(j,1)
        f(k)   = f(k)   - rxx1*bx(j,3)
        a(s-1) = a(s-1) - rxx1*bx(j,2)
      end if
      if (j > 1) then
        l = l + 1
        a(l) = ryy1
        ja(l) = k - 1
      else
        a(s)   = a(s)   - ryy1*ay(i,1)
        f(k)   = f(k)   - ryy1*ay(i,3)
        a(s+1) = a(s+1) - ryy1*ay(i,2)
      end if
      if (i > 1) then
        l = l + 1
        a(l) = rxx1
        ja(l) = k - ny
      else
        a(s)   = a(s)   - rxx1*ax(j,1)
        f(k)   = f(k)   - rxx1*ax(j,3)
        a(s+1) = a(s+1) - rxx1*ax(j,2)
      end if
      ia(k+1) = l+1
    end do
  end do

  ! Set parameters for Pardiso subroutine
  pt(1:64)%DUMMY = 0
  error = 0
  iparm(1:64) = 0

  ! Analysis, numerical factorization, solving, and iterative refinement
  call pardiso_d(pt, 1, 1, 1, 13, n, a, ia, ja, perm, 1, iparm, 0, f, x, error)

  ! Store solution into 2D array
  do i = 1,nx
    do j = 1,ny
      p(i,j) = x(ny*(i-1)+j)
    end do
  end do

  ! Compute boundary values
  do i=0,ny+1
    p(0,i)    = ax(i,3) - ax(i,1)*p(1,i)  - ax(i,2)*p(2,i)
    p(nx+1,i) = bx(i,3) - bx(i,1)*p(nx,i) - bx(i,2)*p(nx-1,i)
  end do
  do i=0,nx+1
    p(i,0)    = ay(i,3) - ay(i,1)*p(i,1)  - ay(i,2)*p(i,2)
    p(i,ny+1) = by(i,3) - by(i,1)*p(i,ny) - by(i,2)*p(i,ny-1)
  end do

  !Release all internal memory for all matrices
  call pardiso_d(pt, 1, 1, 1, -1, n, a, ia, ja, perm, 1, iparm, 0, f, x, error)

end subroutine helmholtz_solve_pardiso

!> subroutine to solve Helmholtz problem by MKL method
subroutine helmholtz_solve_mkl(p,ax,bx,ay,by,rhs,bc,q,nx,ny,hx,hy)
implicit none

real(DP) :: p(:,:), ax(:), bx(:), ay(:), by(:), rhs(:,:)                        ! Arrayes of variables, boundary conditions and right hand sides
character(4) :: bc                                                              ! Types of boundary conditions
real(DP) :: q                                                                   ! Coefficient of u term
integer :: nx, ny                                                               ! Numbers of intervals in x and y directions
real(DP):: hx, hy                                                               ! Space steps

integer :: stat                                                                 ! Status of execution
real(DP) :: lx, ly                                                              ! Lengths

  ! Compute lenghts
  lx = hx*(nx+1)
  ly = hy*(ny+1)

  ! Allocate temporary array of parameters
  allocate(dpar(1:5*(nx+1)/2+8)) ! It needs to add one to recomended value in MKL documentation

  ! Execute MKL subroitines
  call d_init_Helmholtz_2D (0._DP,lx,0._DP,ly,nx+1,ny+1,bc,q,ipar,dpar,stat)
  if (stat /= 0) write(*,*) 'MKL error!'
  call d_commit_Helmholtz_2D(rhs,ax,bx,ay,by,xhandle,ipar,dpar,stat)
  if (stat /= 0) write(*,*) 'MKL error!'
  call d_Helmholtz_2D(rhs,ax,bx,ay,by,xhandle,ipar,dpar,stat)
  if (stat /= 0) write(*,*) 'MKL error!'
  call free_Helmholtz_2D(xhandle,ipar,stat)
  if (stat /= 0) write(*,*) 'MKL error!'

  ! Store solution
  p = rhs

  ! Deallocate array
  deallocate(dpar)

end subroutine helmholtz_solve_mkl

end module helmholtz