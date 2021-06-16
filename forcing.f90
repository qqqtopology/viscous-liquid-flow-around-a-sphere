!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Module to implement forcing procedures
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program: ibm2d_cyl -- SIMPLE + IBM method for incomressible viscous Navier -
! Stockes equations in 2D statement for flow past a cylinder
! Copyright (C): Denis V. Esipov (2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module forcing
use types
use Lin_solve_lib
implicit none

  ! Constants for computation of regularized delta function
  public  force_comp, correct_vel_ibm
  private compose_w_u, interpolate_u, spread_u, interpolate_v, spread_v
  private delta, dh
  real(DP), parameter, private :: e1_3 = 1._DP/3._DP                            ! 1/3
  real(DP), parameter, private :: e1_6 = 1._DP/6._DP                            ! 1/6

contains

!> subroutine to compute force induced by immersed boundary
subroutine force_comp(u,v,fx,fy,xl,ul,vl,nx,ny,nl,hx,hy,hl,tau,fe)
implicit none

real(DP), allocatable :: u(:,:), v(:,:)                                         ! Arrayes of variables
real(DP), allocatable :: fx(:,:), fy(:,:)                                       ! Arrayes for force terms
real(DP), allocatable :: xl(:,:), ul(:), vl(:)                                  ! Arrayes of variables
integer :: nx, ny, nl                                                           ! Numbers of cells in x and y directions and Lagrangian nodes
real(DP):: hx, hy, hl, tau                                                      ! Space and time steps
real(DP) :: fe(1:2)                                                             ! Computed total force

integer :: i, j, k, ii, jj                                                      ! Loop indices
real(DP) :: x(1:2), xx(1:2)                                                     ! Points
integer :: i0, i1, j0, j1                                                       ! Secondary indices
real(DP), allocatable :: uel(:), vel(:)                                         ! Velocities of fluid in Lagrangian nodes

  ! Allocate arrayes
  allocate(uel(1:nl))
  allocate(vel(1:nl))

  ! Compute and correct horizontal velocity u
  call interpolate_u(u,xl,uel,nl,hx,hy)

  ! Correct horizontal velocity u in the Eulerian nodes
  fx = 0._DP
  do i=1,nl
    i0 = nint(xl(i,1)/hx - .5_DP) ! Set up rectangle
    i1 = nint(xl(i,1)/hx + 2.5_DP)
    j0 = nint(xl(i,2)/hy - 1._DP)
    j1 = nint(xl(i,2)/hy + 2._DP)
    do ii=i0,i1
      x(1) = (ii - 1)*hx
      do jj=j0,j1
        x(2) = (jj - .5_DP)*hy
        xx(1:2) = xl(i,1:2) - x(1:2)
        fx(ii,jj) = fx(ii,jj) + (uel(i) - ul(i))*delta(xx(1:2),hx,hy)*hl*min(hx,hy)
      end do
    end do
  end do

  ! Compute vertical velocity v in the Lagrangian nodes
  call interpolate_v(v,xl,vel,nl,hx,hy)

  ! Correct vertical velocity v in the Eulerian nodes
  fy = 0._DP
  do i=1,nl
    i0 = nint(xl(i,1)/hx - 1._DP) ! Set up rectangle
    i1 = nint(xl(i,1)/hx + 2._DP)
    j0 = nint(xl(i,2)/hy - .5_DP)
    j1 = nint(xl(i,2)/hy + 2.5_DP)
    do ii=i0,i1
      x(1) = (ii - .5_DP)*hx
      do jj=j0,j1
        x(2) = (jj - 1)*hy
        xx(1:2) = xl(i,1:2) - x(1:2)
        fy(ii,jj) = fy(ii,jj) + (vel(i) - vl(i))*delta(xx(1:2),hx,hy)*hl*min(hx,hy)
      end do
    end do
  end do

  ! Devide on tau
  fx = fx/tau
  fy = fy/tau
  fe(1) = sum(fx)*hx*hy
  fe(2) = sum(fy)*hx*hy

  ! Write information to temporary file
  write(8000,'(4D)') maxval(abs(uel(1:nl) - ul(1:nl))), maxval(abs(vel(1:nl) - vl(1:nl))), fe(1:2)

  ! Deallocate arrayes
  deallocate(uel)
  deallocate(vel)

end subroutine force_comp

!> subroutine to correct velocity field using the IBM method
subroutine correct_vel_ibm(u,v,xl,ul,vl,nx,ny,nl,hx,hy,hl,tau,fe)
implicit none

real(DP), allocatable :: u(:,:), v(:,:), xl(:,:), ul(:), vl(:)                  ! Arrayes of variables
integer :: nx, ny, nl                                                           ! Numbers of cells in x and y directions and Lagrangian nodes
real(DP):: hx, hy, hl, tau                                                      ! Space and time steps
real(DP) :: fe(1:2)                                                             ! Computed total force

integer :: i, j, k, ii, jj                                                      ! Loop indices
real(DP) :: x(1:2), xx(1:2)                                                     ! Points
integer :: i0, i1, j0, j1                                                       ! Secondary indices
real(DP), allocatable :: uel(:), vel(:)                                         ! Velocities of fluid in Lagrangian nodes
real(DP) :: f                                                                   ! Force in the node

real(DP), allocatable :: w(:,:)                                                 ! Array for operator U_d*u_d
real(DP), allocatable :: rhs(:)                                                 ! Right hand side vector
real(DP), allocatable :: ful(:), fvl(:)                                         ! Force in Lagrangian nodes

  ! Allocate arrayes
  allocate(uel(1:nl))
  allocate(vel(1:nl))
  allocate(ful(1:nl))
  allocate(fvl(1:nl))
  ful = 0._DP
  fvl = 0._DP
  fe = 0._DP

  if (md > 0) then ! Direct or multidirect forcing

    ! Forcing for horizontal velocity
    do k=1,md

      ! Compute and correct horizontal velocity u
      call interpolate_u(u,xl,uel,nl,hx,hy)
      ful = -1.95_DP*(uel - ul)
      call spread_u(xl,ful,u,nl,hx,hy,hl)
      fe(1) = fe(1) - sum(ful)

      ! Compute and correct vertical velocity v
      call interpolate_v(v,xl,vel,nl,hx,hy)
      fvl = -1.95_DP*(vel - vl)
      call spread_v(xl,fvl,v,nl,hx,hy,hl)
      fe(2) = fe(2) - sum(fvl)

    end do

  else ! Exact forcing

    ! Allocate array
    allocate(w(1:nl,1:nl))
    allocate(rhs(1:nl))

    ! Correct horizontal velocity u
    call interpolate_u(u,xl,uel,nl,hx,hy)
    call compose_w_u(w,xl,nx,ny,nl,hx,hy,hl)
    rhs = uel - ul
    call Gauss(w,rhs,ful,nl)
    ful = - ful
    call spread_u(xl,ful,u,nl,hx,hy,hl)
    fe(1) = fe(1) - sum(ful)

    ! Correct vertical velocity v
    call interpolate_v(v,xl,vel,nl,hx,hy)
    call compose_w_v(w,xl,nx,ny,nl,hx,hy,hl)
    rhs = vel - vl
    call Gauss(w,rhs,fvl,nl)
    fvl = - fvl
    call spread_v(xl,fvl,v,nl,hx,hy,hl)
    fe(2) = fe(2) - sum(fvl)

    ! Deallocate arrayes
    deallocate(w)
    deallocate(rhs)

  end if

  ! Devide on tau the value of force
  fe = fe*hl*min(hx,hy)/tau

  ! Write information to temporary file
  call interpolate_u(u,xl,uel,nl,hx,hy)
  call interpolate_v(v,xl,vel,nl,hx,hy)
  write(8000,'(4D)') maxval(abs(uel(1:nl) - ul(1:nl))), maxval(abs(vel(1:nl) - vl(1:nl))), fe(1:2)

  ! Write distribution of lagrangian force into file
  !write(8080,'(D)') ful
  !write(8081,'(D)') fvl

  ! Deallocate arrayes
  deallocate(uel)
  deallocate(vel)
  deallocate(ful)
  deallocate(fvl)

end subroutine correct_vel_ibm

!> subroutine to compose matrix for operator w = U_d*u_d for horizontal
!>                                                                    velocity u
subroutine compose_w_u(w,xl,nx,ny,nl,hx,hy,hl)
implicit none

real(DP), allocatable :: w(:,:), xl(:,:)                                        ! Arrayes of variables
integer :: nx, ny, nl                                                           ! Numbers of cells in x and y directions and Lagrangian nodes
real(DP):: hx, hy, hl                                                           ! Space steps

real(DP), allocatable :: u(:,:),  ul(:)                                         ! Arrayes of variables
integer :: i                                                                    ! Loop index

  ! Allocate array
  allocate(u(0:nx+2,0:ny+1))
  allocate(ul(1:nl))

  ! Main loop
  do i=1,nl

    ! Set unity vector
    ul = 0._DP
    ul(i) = 1._DP

    ! Compute the value of operator U_d*u_d
    u = 0._DP
    call spread_u(xl,ul,u,nl,hx,hy,hl)
    call interpolate_u(u,xl,ul,nl,hx,hy)

    ! Add column into the matrix
    w(1:nl,i) = ul(1:nl)

  end do

  ! Deallocate arrayes
  deallocate(u)
  deallocate(ul)

end subroutine compose_w_u

!> subroutine to compose matrix for operator w = U_d*u_d for vertical velocity v
subroutine compose_w_v(w,xl,nx,ny,nl,hx,hy,hl)
implicit none

real(DP), allocatable :: w(:,:), xl(:,:)                                        ! Arrayes of variables
integer :: nx, ny, nl                                                           ! Numbers of cells in x and y directions and Lagrangian nodes
real(DP):: hx, hy, hl                                                           ! Space steps

real(DP), allocatable :: v(:,:),  vl(:)                                         ! Arrayes of variables
integer :: i                                                                    ! Loop index

  ! Allocate array
  allocate(v(0:nx+1,0:ny+2))
  allocate(vl(1:nl))

  ! Main loop
  do i=1,nl

    ! Set unity vector
    vl = 0._DP
    vl(i) = 1._DP

    ! Compute the value of operator U_d*u_d
    v = 0._DP
    call spread_v(xl,vl,v,nl,hx,hy,hl)
    call interpolate_v(v,xl,vl,nl,hx,hy)

    ! Add column into the matrix
    w(1:nl,i) = vl(1:nl)

  end do

  ! Deallocate arrayes
  deallocate(v)
  deallocate(vl)

end subroutine compose_w_v

!> subroutine to interpolate values from Eulerian grid to Lagrangian mesh
!>                                                     for horizontal velocity u
subroutine interpolate_u(u,xl,ul,nl,hx,hy)
implicit none

real(DP), allocatable :: u(:,:), xl(:,:), ul(:)                                 ! Arrayes of variables
integer :: nl                                                                   ! Number of Lagrangian nodes
real(DP):: hx, hy                                                               ! Space steps

integer :: i, j, k, ii, jj                                                      ! Loop indices
real(DP) :: x(1:2), xx(1:2)                                                     ! Points
integer :: i0, i1, j0, j1                                                       ! Secondary indices

  ! Compute horizontal velocity u in the Lagrangian nodes
  ul = 0._DP
  do i=1,nl
    i0 = nint(xl(i,1)/hx - .5_DP) ! Set up rectangle
    i1 = nint(xl(i,1)/hx + 2.5_DP)
    j0 = nint(xl(i,2)/hy - 1._DP)
    j1 = nint(xl(i,2)/hy + 2._DP)
    do ii=i0,i1
      x(1) = (ii - 1)*hx
      do jj=j0,j1
        x(2) = (jj - .5_DP)*hy
        xx(1:2) = xl(i,1:2) - x(1:2)
        ul(i) = ul(i) + u(ii,jj)*delta(xx(1:2),hx,hy)*hx*hy
      end do
    end do
  end do

end subroutine interpolate_u

!> subroutine to spread values from Lagrangian mesh to Eulerian grid
!>                                                     for horizontal velocity u
subroutine spread_u(xl,ul,u,nl,hx,hy,hl)
implicit none

real(DP), allocatable :: u(:,:), xl(:,:), ul(:)                                 ! Arrayes of variables
integer :: nl                                                                   ! Number of Lagrangian nodes
real(DP):: hx, hy, hl                                                           ! Space steps

integer :: i, j, k, ii, jj                                                      ! Loop indices
real(DP) :: x(1:2), xx(1:2)                                                     ! Points
integer :: i0, i1, j0, j1                                                       ! Secondary indices
real(DP) :: f                                                                   ! Force in the node

  ! Correct horizontal velocity u in the Eulerian nodes
  do i=1,nl
    i0 = nint(xl(i,1)/hx - .5_DP) ! Set up a rectangle
    i1 = nint(xl(i,1)/hx + 2.5_DP)
    j0 = nint(xl(i,2)/hy - 1._DP)
    j1 = nint(xl(i,2)/hy + 2._DP)
    do ii=i0,i1
      x(1) = (ii - 1)*hx
      do jj=j0,j1
        x(2) = (jj - .5_DP)*hy
        xx(1:2) = xl(i,1:2) - x(1:2)
        u(ii,jj) = u(ii,jj) + ul(i)*delta(xx(1:2),hx,hy)*hl*min(hx,hy)
      end do
    end do
  end do

end subroutine spread_u

!> subroutine to interpolate values from Eulerian grid to Lagrangian mesh
!>                                                     for horizontal velocity v
subroutine interpolate_v(v,xl,vl,nl,hx,hy)
implicit none

real(DP), allocatable :: v(:,:), xl(:,:), vl(:)                                 ! Arrayes of variables
integer :: nl                                                                   ! Number of Lagrangian nodes
real(DP):: hx, hy                                                               ! Space steps

integer :: i, j, k, ii, jj                                                      ! Loop indices
real(DP) :: x(1:2), xx(1:2)                                                     ! Points
integer :: i0, i1, j0, j1                                                       ! Secondary indices

  ! Compute vertical velocity v in the Lagrangian nodes
  vl = 0._DP
  do i=1,nl
    i0 = nint(xl(i,1)/hx - 1._DP) ! Set up a rectangle
    i1 = nint(xl(i,1)/hx + 2._DP)
    j0 = nint(xl(i,2)/hy - .5_DP)
    j1 = nint(xl(i,2)/hy + 2.5_DP)
    do ii=i0,i1
      x(1) = (ii - .5_DP)*hx
      do jj=j0,j1
        x(2) = (jj - 1)*hy
        xx(1:2) = xl(i,1:2) - x(1:2)
        vl(i) = vl(i) + v(ii,jj)*delta(xx(1:2),hx,hy)*hx*hy
      end do
    end do
  end do

end subroutine interpolate_v

!> subroutine to spread values from Lagrangian mesh to Eulerian grid
!>                                                     for horizontal velocity v
subroutine spread_v(xl,vl,v,nl,hx,hy,hl)
implicit none

real(DP), allocatable :: v(:,:), xl(:,:), vl(:)                                 ! Arrayes of variables
integer :: nl                                                                   ! Number of Lagrangian nodes
real(DP):: hx, hy, hl                                                           ! Space steps

integer :: i, j, k, ii, jj                                                      ! Loop indices
real(DP) :: x(1:2), xx(1:2)                                                     ! Points
integer :: i0, i1, j0, j1                                                       ! Secondary indices

  ! Correct vertical velocity v in the Eulerian nodes
  do i=1,nl
    i0 = nint(xl(i,1)/hx - 1._DP) ! Set up rectangle
    i1 = nint(xl(i,1)/hx + 2._DP)
    j0 = nint(xl(i,2)/hy - .5_DP)
    j1 = nint(xl(i,2)/hy + 2.5_DP)
    do ii=i0,i1
      x(1) = (ii - .5_DP)*hx
      do jj=j0,j1
        x(2) = (jj - 1)*hy
        xx(1:2) = xl(i,1:2) - x(1:2)
        v(ii,jj) = v(ii,jj) + vl(i)*delta(xx(1:2),hx,hy)*hl*min(hx,hy)
      end do
    end do
  end do

end subroutine spread_v

!> function to compute discrete delta function
function delta(x,hx,hy)
implicit none

real(DP) :: x(1:2)                                                              ! Coordinates of point
real(DP) :: hx, hy                                                              ! Space steps
real(DP) :: delta                                                               ! Function result

  delta = dh(x(1)/hx)*dh(x(2)/hy)/(hx*hy)

end function delta

!> function to compute 1D approximation of delta function
function dh(r)
implicit none

real(DP) :: r, dh                                                               ! Argument and fucntion result

real(DP) :: ra                                                                  ! Absolute value of r

  ! Peskin variant
  ! from Peskin C.S. The immersed boundary method // Acta numerica. 2002.
  !      Vol. 11. P. 479-517
!  ra = abs(r)
!  if (ra < 1._DP) then
!    dh = .125_DP*(3._DP - 2._DP*ra + sqrt(1._DP + 4._DP*ra - 4._DP*ra**2))
!  else if (ra < 2._DP) then
!    dh = .125_DP*(5._DP - 2._DP*ra - sqrt(-7._DP + 12._DP*ra - 4._DP*ra**2))
!  else
!    dh = 0._DP
!  end if

  ! Roma variant
  ! from Roma A.M., Peskin C.S., Berger M.J. An adaptive version of the immersed
  !      boundary method // J. Computational Physics. 1999. Vol. 153, iss. 2.
  !      P. 509-534.
  ra = abs(r)
  if (ra < .5_DP) then
    dh = e1_3*(1._DP + sqrt(1 - 3._DP*ra**2))
  else if (ra < 1.5_DP) then
    dh = e1_6*(5._DP - 3._DP*ra - sqrt(1 - 3._DP*(1 - ra)**2))
  else
    dh = 0._DP
  end if

end function dh

end module forcing