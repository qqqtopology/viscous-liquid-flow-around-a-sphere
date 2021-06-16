!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Program to solve 2D viscous incompressible Navier-Stokes equations:
!   u_t + u*u_x + v*u_y = (1/Re)*(u_xx + u_yy) - p_x + fx
!   v_t + u*v_x + v*v_y = (1/Re)*(v_xx + v_yy) - p_y + fy
!   u_x + v_y = 0
! by SIMPLE method (3 variants distinguished by approximation of convective
! terms). Computation is performed on the squared staggered mesh with
! fictitious cells near boundaries to obtain second order of approximation of
! the whole problem. Computational domain is rectangle of size lx*ly.
! The flow past a cylinder placed at point (lx/3,ly/2) is set up.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (C): Denis V. Esipov (2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program ibm2d_test
use types
use helmholtz
use forcing
use inout
implicit none

! Variables
real(DP), allocatable :: u(:,:), v(:,:), p(:,:), u0(:,:), v0(:,:)               ! Arrayes of Eulerian variables
real(DP), allocatable :: u1(:,:), v1(:,:), p0(:,:)                              ! Arrayes of Eulerian variables on the previous steps
real(DP), allocatable :: ut(:,:), vt(:,:), pt(:,:)                              ! Arrayes of Eulerian variables for undisturbed flow
real(DP), allocatable :: fu(:,:), fv(:,:)                                       ! Arrayes for source terms in motion equation
real(DP), allocatable :: fx(:,:), fy(:,:)                                       ! Arrayes for force terms
real(DP), allocatable :: xl(:,:), ul(:), vl(:)                                  ! Arrayes of Lagrangian variables (coordinates and velocities)
real(DP), allocatable :: rhsu(:,:), rhsv(:,:), rhs(:,:), pc(:,:)                ! Arrayes of usefull variables (RHS and pressure correction fields)
real(DP) :: lx, ly                                                              ! Size of computational domain
integer :: nx, ny, nl, nt                                                       ! Numbers of cells in x and y directions, Lagrangian nodes and time steps
real(DP):: hx, hy, hl, tau                                                      ! Space and time steps
real(DP) :: Re                                                                  ! Reynolds number
integer :: it, strand, i, j, k, nk                                              ! Loop indices
type(sratios) :: r, r2                                                          ! Step ratios

real(DP) :: pvar                                                                ! Pressure variation on internal iterations
real(DP) :: x, y, t                                                             ! Coordinate and time
real(DP) :: t_end                                                               ! Computation stop time
real(DP) :: eps_pvar                                                            ! Stop condition for internal iterations (pressure variation)
integer :: deb                                                                  ! Debug flag
integer :: start                                                                ! Start flag (1 -- load data from 'start.dat', 0 -- nothing)
real(DP) :: fe(1:2)                                                             ! Computed total force

  ! Set up flags
  deb = 1
  start = 0

  ! Set up model parameters
  lx = 15._DP
  ly = 15._DP
  Re = 20._DP

  ! Set up mesh and time step (set of minimal parameters: 150,150,0.01,100)
  nx = 150
  ny = 150
  tau = 0.01_DP
  nl = 20
  t_end = 100._DP
  eps_pvar = 1.e-8_DP
  hx = lx/nx
  hy = ly/ny
  nt = t_end/tau
  ! phi = 2._DP*pi/nl ! Circular cylinder of diameter 1.
  hl = sin(pi/nl) ! Precise hl
  call step_ratios(r,hx,hy,tau,Re)
  call step_ratios(r2,hx,hy,.5_DP*tau,Re)

  ! Allocate arrayes of variables
  ! Scheme of indexation for every computational cell
  ! 1. Node for pressure p(i,j)                +-x-+
  ! 2. Node for horizontal velocity u(i,j)     |   | Computation cell
  ! 3. Node for vertical velocity v(i,j)       2 1 x
  !                                            |   |
  !                                            +-3-+
  allocate(u(0:nx+2,0:ny+1)) ! First index corrected by adding 0.5, second index is the same as for pressure
  allocate(v(0:nx+1,0:ny+2)) ! First index is the same as for pressure, second index corrected by adding 0.5
  allocate(u0(0:nx+2,0:ny+1)); allocate(v0(0:nx+1,0:ny+2))
  allocate(u1(0:nx+2,0:ny+1)); allocate(v1(0:nx+1,0:ny+2))
  allocate(ut(0:nx+2,0:ny+1)); allocate(vt(0:nx+1,0:ny+2)); allocate(pt(0:nx+1,0:ny+1));
  allocate(p(0:nx+1,0:ny+1));  allocate(p0(0:nx+1,0:ny+1))
  allocate(fu(0:nx+2,0:ny+1)); allocate(fv(0:nx+1,0:ny+2))
  allocate(fx(0:nx+2,0:ny+1)); allocate(fy(0:nx+1,0:ny+2))
  allocate(xl(1:nl,1:2));      allocate(ul(1:nl));          allocate(vl(1:nl))

  ! Allocate arrayes of usefull variables
  allocate(rhsu(1:nx+1,1:ny)); allocate(rhsv(1:nx,1:ny+1)); allocate(rhs(0:nx+1,0:ny+1))
  allocate(pc(0:nx+1,0:ny+1))

  ! Set up initial conditions
  ! Starting from the uniform flow situation
  u = 1._DP; v = 0._DP; u0 = 1._DP; v0 = 0._DP; p = 0._DP
  pt = p ! For met == 6

  ! Set up Lagrangian mesh for circular cylinder of diameter 1 and placed at point (lx/3,ly/2)
  do i=1,nl
    xl(i,1) = lx/3._DP + .5_DP*cos((2._DP*pi*i)/nl)
    xl(i,2) = .5_DP*ly + .5_DP*sin((2._DP*pi*i)/nl)
    ul(i) = 0._DP ! The cylinder is fixed
    vl(i) = 0._DP
  end do

  ! Start from old data
  if (start == 1) then
    call load_dat(u,v,p,u0,v0,xl,ul,vl,nx,ny,nl,hx,hy,hl,Re)
  end if

  ! Save initial data
  strand = 0
  call save_plt(u,v,p,xl,ul,vl,nx,ny,nl,0,hx,hy,tau,strand,Re)

  ! Main loop
  t = 0._DP
  pc = 0._DP
  fu = 0._DP
  fv = 0._DP
  fx = 0._DP ! Set zero for met == 1-4
  fy = 0._DP
  do it=1,nt

    ! Store variables from previous time step
    u1 = u0; v1 = v0; u0 = u; v0 = v; p0 = p

    ! Increase time
    t = it*tau

    ! Choose the maximal number of iterations
    !nk = 50000 ! The maximum is 50000 iterations
    nk = 1 ! One iteration method

    ! Compute force for met 6
    if (met == 6) then

      fx = 0._DP ! Set zero to avoid forcing for undisturbed velocity computation
      fy = 0._DP

      ! Compute RHS
      call make_uv_RHS(rhsu,rhsv,u0,v0,u1,v1,p0,fu,fv,nx,ny,r2)

      do k=1,nk

        ! Make prediction of velocities
        call predict_uv(ut,vt,u0,v0,u1,v1,pt,rhsu,rhsv,fx,fy,nx,ny,hx,hy,tau,Re,t,r2)

        ! Solve Poisson equation for pressure correction
        call solve_pcp(pc,rhs,ut,vt,pt,nx,ny,hx,hy,tau,Re,t,k,r2)

        ! Correct pressure and velocities
        pt = pt + pc
        call correct_vel(ut,vt,pc,nx,ny,hx,hy,tau,t,r2)

        ! Compute variation of pressure field at the intermediate iterations
        pvar = 0._DP
        do i=1,nx
          do j=1,ny
            pvar = max(pvar,abs(pc(i,j)))
          end do
        end do
        write(*,'(2I,D)') it, k, pvar

        ! Stop condition
        if (pvar < eps_pvar) exit

      end do

      ! Compute force
      call force_comp(ut,vt,fx,fy,xl,ul,vl,nx,ny,nl,hx,hy,hl,tau,fe)

    end if

    ! Compute RHS
    call make_uv_RHS(rhsu,rhsv,u0,v0,u1,v1,p0,fu,fv,nx,ny,r)

    ! Iterative SIMPLE method
    do k=1,nk

      ! Make prediction of velocities
      if (met == 1 .or. met == 2 .or. met == 3 .or. met == 4) then
        call predict_uv(u,v,u0,v0,u1,v1,p,rhsu,rhsv,fx,fy,nx,ny,hx,hy,tau,Re,t,r)
        call correct_vel_ibm(u,v,xl,ul,vl,nx,ny,nl,hx,hy,hl,tau,fe) ! IBM correction
      else if (met == 5) then
        fx = 0._DP ! Set zero for separate force computation
        fy = 0._DP
        call predict_uv(u,v,u0,v0,u1,v1,p,rhsu,rhsv,fx,fy,nx,ny,hx,hy,tau,Re,t,r)
        call force_comp(u,v,fx,fy,xl,ul,vl,nx,ny,nl,hx,hy,hl,tau,fe)
        call predict_uv(u,v,u0,v0,u1,v1,p,rhsu,rhsv,fx,fy,nx,ny,hx,hy,tau,Re,t,r) ! Prediction with IBM correction (adding force)
      else if (met == 6) then
        call predict_uv(u,v,u0,v0,u1,v1,p,rhsu,rhsv,fx,fy,nx,ny,hx,hy,tau,Re,t,r)
      end if

      ! Solve Poisson equation for pressure correction
      call solve_pcp(pc,rhs,u,v,p,nx,ny,hx,hy,tau,Re,t,k,r)

      ! Correct pressure and velocities
      p = p + pc
      call correct_vel(u,v,pc,nx,ny,hx,hy,tau,t,r)

      ! Compute variation of pressure field at the intermediate iterations
      pvar = 0._DP
      do i=1,nx
        do j=1,ny
          pvar = max(pvar,abs(pc(i,j)))
        end do
      end do
      write(*,'(2I,D,2F12.4)') it, k, pvar, 2._DP*fe
       call move_particles(particles,1,xl,ul,vl,fl,lx,ly,hx,hy,hl,tau,Re,Bu,gamma)
      ! Stop condition
      if (pvar < eps_pvar) exit

    end do

    ! Write message on the screen if at the same time iteration process is stopped and stop condition is not fullfilled
    if (k > nk) then
      write(*,'(A)') 'Stop condition for internal iterations is not fulfilled!'
      write(*,'(A,E12.5)') 'pvar = ', pvar
    end if

    ! Save results
    if (mod(it,200) == 0) then
      call save_dat(u,v,p,u0,v0,xl,ul,vl,nx,ny,nl,hx,hy,hl,it,tau,Re)
      call save_plt(u,v,p,xl,ul,vl,nx,ny,nl,it,hx,hy,tau,strand,Re)
    end if
  end do

  ! Deallocate arrayes
  deallocate(u);     deallocate(v);    deallocate(p)
  deallocate(xl);    deallocate(ul);   deallocate(vl)
  deallocate(u0);    deallocate(v0);   deallocate(p0)
  deallocate(u1);    deallocate(v1)
  deallocate(ut);    deallocate(vt);   deallocate(pt)
  deallocate(fu);    deallocate(fv)
  deallocate(fx);    deallocate(fy)
  deallocate(rhsu);  deallocate(rhsv); deallocate(rhs)
  deallocate(pc)

contains

!> subroutine to compute right hand sides of impulse equations (see predict_uv
!                                                                    subroutine)
subroutine make_uv_RHS(rhsu,rhsv,u0,v0,u1,v1,p0,fu,fv,nx,ny,r)
implicit none

real(DP), allocatable :: rhsu(:,:), rhsv(:,:), u0(:,:), v0(:,:)                 ! Arrayes of variables
real(DP), allocatable :: p0(:,:), u1(:,:), v1(:,:), fu(:,:), fv(:,:)            ! Arrayes of variables
integer :: nx, ny                                                               ! Numbers of cells in x and y directions
type(sratios) :: r                                                              ! Step ratios

integer :: i, j                                                                 ! Loop indices
real(DP) :: a1, a2                                                              ! Approximate advection velocities

  ! Compute rhs arrayes
  if (met == 1) then
    do i=1,nx+1
      do j=1,ny
        rhsu(i,j) = r%re2_t*u0(i,j) & 
                  + r%e_xx *(u0(i+1,j) - 2._DP*u0(i,j) + u0(i-1,j)) & 
                  + r%e_yy *(u0(i,j+1) - 2._DP*u0(i,j) + u0(i,j-1)) & 
                  + r%re2  *fu(i,j)
      end do
    end do
    do i=1,nx
      do j=1,ny+1
        rhsv(i,j) = r%re2_t*v0(i,j) & 
                  + r%e_xx *(v0(i+1,j) - 2._DP*v0(i,j) + v0(i-1,j)) & 
                  + r%e_yy *(v0(i,j+1) - 2._DP*v0(i,j) + v0(i,j-1)) & 
                  + r%re2  *fv(i,j)
      end do
    end do
  else if (met == 2) then
    do i=1,nx+1
      do j=1,ny
        a1 = 1.5_DP*u0(i,j) - .5_DP*u1(i,j)
        a2 = .375_DP*(v0(i-1,j) + v0(i-1,j+1) + v0(i,j) + v0(i,j+1)) - .125_DP*(v1(i-1,j) + v1(i-1,j+1) + v1(i,j) + v1(i,j+1))
        rhsu(i,j) = r%re2_t*u0(i,j) & 
                  - r%re_2x*a1*(u0(i+1,j) - u0(i-1,j)) & 
                  - r%re_2y*a2*(u0(i,j+1) - u0(i,j-1)) & 
                  + r%e_xx *(u0(i+1,j) - 2._DP*u0(i,j) + u0(i-1,j)) & 
                  + r%e_yy *(u0(i,j+1) - 2._DP*u0(i,j) + u0(i,j-1)) & 
                  + r%re2  *fu(i,j)
      end do
    end do
    do i=1,nx
      do j=1,ny+1
        a1 = .375_DP*(u0(i,j-1) + u0(i+1,j-1) + u0(i,j) + u0(i+1,j)) - .125_DP*(u1(i,j-1) + u1(i+1,j-1) + u1(i,j) + u1(i+1,j))
        a2 = 1.5_DP*v0(i,j) - .5_DP*v1(i,j)
        rhsv(i,j) = r%re2_t*v0(i,j) & 
                  - r%re_2x*a1*(v0(i+1,j) - v0(i-1,j)) & 
                  - r%re_2y*a2*(v0(i,j+1) - v0(i,j-1)) & 
                  + r%e_xx *(v0(i+1,j) - 2._DP*v0(i,j) + v0(i-1,j)) & 
                  + r%e_yy *(v0(i,j+1) - 2._DP*v0(i,j) + v0(i,j-1)) & 
                  + r%re2  *fv(i,j)
      end do
    end do
  else if (met == 3) then
    do i=1,nx+1
      do j=1,ny
        rhsu(i,j) = r%re2_t*u0(i,j) & 
                  - r%re_x *(1.5_DP*(u0(i+1,j)**2 - u0(i-1,j)**2) - .5_DP*(u1(i+1,j)**2 - u1(i-1,j)**2)) & 
                  - r%re_2y*(1.5_DP*((u0(i,j) + u0(i,j+1))*(v0(i,j+1) + v0(i-1,j+1)) - (u0(i,j) + u0(i,j-1))*(v0(i,j) + v0(i-1,j))) & 
                             -.5_DP*((u1(i,j) + u1(i,j+1))*(v1(i,j+1) + v1(i-1,j+1)) - (u1(i,j) + u1(i,j-1))*(v1(i,j) + v1(i-1,j)))) & 
                  + r%e_xx *(u0(i+1,j) - 2._DP*u0(i,j) + u0(i-1,j)) & 
                  + r%e_yy *(u0(i,j+1) - 2._DP*u0(i,j) + u0(i,j-1)) & 
                  + r%re2  *fu(i,j)
      end do
    end do
    do i=1,nx
      do j=1,ny+1
        rhsv(i,j) = r%re2_t*v0(i,j) & 
                  - r%re_2x*(1.5_DP*((u0(i+1,j) + u0(i+1,j-1))*(v0(i,j) + v0(i+1,j)) - (u0(i,j) + u0(i,j-1))*(v0(i,j) + v0(i-1,j))) & 
                             -.5_DP*((u1(i+1,j) + u1(i+1,j-1))*(v1(i,j) + v1(i+1,j)) - (u1(i,j) + u1(i,j-1))*(v1(i,j) + v1(i-1,j)))) & 
                  - r%re_y *(1.5_DP*(v0(i,j+1)**2 - v0(i,j-1)**2) - .5_DP*(v1(i,j+1)**2 - v1(i,j-1)**2)) & 
                  + r%e_xx *(v0(i+1,j) - 2._DP*v0(i,j) + v0(i-1,j)) & 
                  + r%e_yy *(v0(i,j+1) - 2._DP*v0(i,j) + v0(i,j-1)) & 
                  + r%re2  *fv(i,j)
      end do
    end do
  else if (met == 4) then
    do i=1,nx+1
      do j=1,ny
        rhsu(i,j) = r%re2_t*u0(i,j) & 
                  - r%re_x *(1.5_DP*(u0(i+1,j)**2 - u0(i-1,j)**2) - .5_DP*(u1(i+1,j)**2 - u1(i-1,j)**2)) & 
                  - r%re_2y*(1.5_DP*((u0(i,j) + u0(i,j+1))*(v0(i,j+1) + v0(i-1,j+1)) - (u0(i,j) + u0(i,j-1))*(v0(i,j) + v0(i-1,j))) & 
                             -.5_DP*((u1(i,j) + u1(i,j+1))*(v1(i,j+1) + v1(i-1,j+1)) - (u1(i,j) + u1(i,j-1))*(v1(i,j) + v1(i-1,j)))) & 
                  + r%e_xx *(u0(i+1,j) - 2._DP*u0(i,j) + u0(i-1,j)) & 
                  + r%e_yy *(u0(i,j+1) - 2._DP*u0(i,j) + u0(i,j-1)) & 
                  + r%re2  *fu(i,j)
      end do
    end do
    do i=1,nx
      do j=1,ny+1
        rhsv(i,j) = r%re2_t*v0(i,j) & 
                  - r%re_2x*(1.5_DP*((u0(i+1,j) + u0(i+1,j-1))*(v0(i,j) + v0(i+1,j)) - (u0(i,j) + u0(i,j-1))*(v0(i,j) + v0(i-1,j))) & 
                             -.5_DP*((u1(i+1,j) + u1(i+1,j-1))*(v1(i,j) + v1(i+1,j)) - (u1(i,j) + u1(i,j-1))*(v1(i,j) + v1(i-1,j)))) & 
                  - r%re_y *(1.5_DP*(v0(i,j+1)**2 - v0(i,j-1)**2) - .5_DP*(v1(i,j+1)**2 - v1(i,j-1)**2)) & 
                  + r%e_xx *(v0(i+1,j) - 2._DP*v0(i,j) + v0(i-1,j)) & 
                  + r%e_yy *(v0(i,j+1) - 2._DP*v0(i,j) + v0(i,j-1)) & 
                  + r%re2  *fv(i,j)
      end do
    end do
  else if (met == 5) then
    do i=1,nx+1
      do j=1,ny
        rhsu(i,j) = r%re2_t*u0(i,j) & 
                  - r%re_x *(1.5_DP*(u0(i+1,j)**2 - u0(i-1,j)**2) - .5_DP*(u1(i+1,j)**2 - u1(i-1,j)**2)) & 
                  - r%re_2y*(1.5_DP*((u0(i,j) + u0(i,j+1))*(v0(i,j+1) + v0(i-1,j+1)) - (u0(i,j) + u0(i,j-1))*(v0(i,j) + v0(i-1,j))) & 
                             -.5_DP*((u1(i,j) + u1(i,j+1))*(v1(i,j+1) + v1(i-1,j+1)) - (u1(i,j) + u1(i,j-1))*(v1(i,j) + v1(i-1,j)))) & 
                  + r%e_xx *(u0(i+1,j) - 2._DP*u0(i,j) + u0(i-1,j)) & 
                  + r%e_yy *(u0(i,j+1) - 2._DP*u0(i,j) + u0(i,j-1)) & 
                  + r%re2  *fu(i,j)
      end do
    end do
    do i=1,nx
      do j=1,ny+1
        rhsv(i,j) = r%re2_t*v0(i,j) & 
                  - r%re_2x*(1.5_DP*((u0(i+1,j) + u0(i+1,j-1))*(v0(i,j) + v0(i+1,j)) - (u0(i,j) + u0(i,j-1))*(v0(i,j) + v0(i-1,j))) & 
                             -.5_DP*((u1(i+1,j) + u1(i+1,j-1))*(v1(i,j) + v1(i+1,j)) - (u1(i,j) + u1(i,j-1))*(v1(i,j) + v1(i-1,j)))) & 
                  - r%re_y *(1.5_DP*(v0(i,j+1)**2 - v0(i,j-1)**2) - .5_DP*(v1(i,j+1)**2 - v1(i,j-1)**2)) & 
                  + r%e_xx *(v0(i+1,j) - 2._DP*v0(i,j) + v0(i-1,j)) & 
                  + r%e_yy *(v0(i,j+1) - 2._DP*v0(i,j) + v0(i,j-1)) & 
                  + r%re2  *fv(i,j)
      end do
    end do
  else if (met == 6) then
    do i=1,nx+1
      do j=1,ny
        rhsu(i,j) = r%re2_t*u0(i,j) & 
                  + r%e_xx *(u0(i+1,j) - 2._DP*u0(i,j) + u0(i-1,j)) & 
                  + r%e_yy *(u0(i,j+1) - 2._DP*u0(i,j) + u0(i,j-1)) & 
                  + r%re2  *fu(i,j)
      end do
    end do
    do i=1,nx
      do j=1,ny+1
        rhsv(i,j) = r%re2_t*v0(i,j) & 
                  + r%e_xx *(v0(i+1,j) - 2._DP*v0(i,j) + v0(i-1,j)) & 
                  + r%e_yy *(v0(i,j+1) - 2._DP*v0(i,j) + v0(i,j-1)) & 
                  + r%re2  *fv(i,j)
      end do
    end do
  else
    write(*,'(A,I1,A)') 'Wrong method (met = ', met,')!'
    stop
  end if

end subroutine make_uv_RHS

!> subroutine to predict velocities using the formulae:
! for met == 1 (for v relation is analogous)
! (u*(i,j) - u0(i,j))/tau + u(i,j)%*(u(i+1,j)% - u(i-1,j)%)/(2*hx)
!   + (v(i-1,j)% + v(i-1,j+1)% + v(i,j)% + v(i,j+1))*(u(i,j+1)% - u(i,j-1)%)/(8*hy)
!   = (1/Re)*(u(i+1,j)% - 2*u0(i,j)% + u0(i-1,j)%)/hx^2
!   + (1/Re)*(u(i,j+1)% - 2*u0(i,j)% + u0(i,j-1)%)/hy^2 - (p*(i,j) - p*(i-1,j))/hx
! sign "%" means averaging on two time steps: a% = 0.5a^(n+1) + 0.5a^n
! p everywhere is considered at t+1/2 time steps.
! for met == 2
! (u*(i,j) - u0(i,j))/tau + u0(i,j)'*(u(i+1,j)% - u0(i-1,j)%)/(2*hx)
!   + (v0(i-1,j)' + v0(i-1,j+1)' + v0(i,j)' + v0(i,j+1)')*(u0(i,j+1)% - u0(i,j-1)%)/(8*hy)
!   = (1/Re)*(u0(i+1,j)% - 2*u0(i,j)% + u0(i-1,j)%)/hx^2
!   + (1/Re)*(u0(i,j+1)% - 2*u0(i,j)% + u0(i,j-1)%)/hy^2 - (p*(i,j) - p*(i-1,j))/hx
! sign "'" means the Adams - Bashforth approximation from two previous time steps a' = 1.5a^n - 0.5a^(n-1)
! for met == 3
! (u*(i,j) - u0(i,j))/tau + [u(i,j)*(u(i+1,j) - u(i-1,j))/(2*hx)]'
!   + [(v0(i-1,j) + v0(i-1,j+1) + v0(i,j) + v0(i,j+1))*(u0(i,j+1) - u0(i,j-1))/(8*hy)]'
!   = (1/Re)*(u0(i+1,j)% - 2*u0(i,j)% + u0(i-1,j)%)/hx^2
!   + (1/Re)*(u0(i,j+1)% - 2*u0(i,j)% + u0(i,j-1)%)/hy^2 - (p*(i,j) - p*(i-1,j))/hx
! for met == 4
! (u*(i,j) - u0(i,j))/tau + [u(i,j)*(u(i+1,j) - u(i-1,j))/(2*hx)]'
!   + [(v0(i-1,j) + v0(i-1,j+1) + v0(i,j) + v0(i,j+1))*(u0(i,j+1) - u0(i,j-1))/(8*hy)]'
!   = (1/Re)*(u0(i+1,j)% - 2*u0(i,j)% + u0(i-1,j)%)/hx^2
!   + (1/Re)*(u0(i,j+1)% - 2*u0(i,j)% + u0(i,j-1)%)/hy^2 - (p*(i,j) - p*(i-1,j))/hx
! for met == 5
! (u*(i,j) - u0(i,j))/tau + [u(i,j)*(u(i+1,j) - u(i-1,j))/(2*hx)]'
!   + [(v0(i-1,j) + v0(i-1,j+1) + v0(i,j) + v0(i,j+1))*(u0(i,j+1) - u0(i,j-1))/(8*hy)]'
!   = (1/Re)*(u0(i+1,j)% - 2*u0(i,j)% + u0(i-1,j)%)/hx^2
!   + (1/Re)*(u0(i,j+1)% - 2*u0(i,j)% + u0(i,j-1)%)/hy^2 - (p*(i,j) - p*(i-1,j))/hx
! for met == 6
! it needs to be written
subroutine predict_uv(u,v,u0,v0,u1,v1,p,rhsu,rhsv,fx,fy,nx,ny,hx,hy,tau,Re,t,r)
implicit none

real(DP), allocatable :: u(:,:), v(:,:), u0(:,:), v0(:,:),  u1(:,:), v1(:,:)    ! Arrayes of variables
real(DP), allocatable :: p(:,:), rhsu(:,:), rhsv(:,:)                           ! Arrayes of variables
real(DP), allocatable :: fx(:,:), fy(:,:)                                       ! Arrayes for force terms
integer :: nx, ny                                                               ! Numbers of cells in x and y directions
real(DP):: hx, hy, tau                                                          ! Space and time steps
real(DP) :: Re                                                                  ! Reynolds number
real(DP) :: t                                                                   ! Time
type(sratios) :: r                                                              ! Step ratios

integer :: i, j                                                                 ! Loop indices
real(DP), allocatable :: fu(:,:), fv(:,:)                                       ! Right hand sides
real(DP) :: a1, a2                                                              ! Approximate advection velocities
real(DP), allocatable :: ax(:,:), bx(:,:), ay(:,:), by(:,:)                     ! Values of boundary conditions

  ! Allocate arrayes
  allocate(fu(0:nx+2,0:ny+1))
  allocate(fv(0:nx+1,0:ny+2))

  ! Compute right hand sides
  fu(0,0:ny+1)    = 0._DP ! These 8 lines added for consistency
  fu(nx+2,0:ny+1) = 0._DP
  fu(1:nx+1,0)    = 0._DP
  fu(1:nx+1,ny+1) = 0._DP
  fv(0,0:ny+2)    = 0._DP
  fv(nx+1,0:ny+2) = 0._DP
  fv(1:nx,0)      = 0._DP
  fv(1:nx,ny+2)   = 0._DP
  if (met == 1) then
    do i=1,nx+1
      do j=1,ny
        a1 = .5_DP*(u(i,j) + u0(i,j))
        a2 = .125_DP*(v(i-1,j) + v(i-1,j+1) + v(i,j) + v(i,j+1) + v0(i-1,j) + v0(i-1,j+1) + v0(i,j) + v0(i,j+1))
        fu(i,j) = rhsu(i,j) & 
                - r%re_2x*a1*(u(i+1,j) - u(i-1,j) + u0(i+1,j) - u0(i-1,j)) & 
                - r%re_2y*a2*(u(i,j+1) - u(i,j-1) + u0(i,j+1) - u0(i,j-1)) & 
                - r%re2_x *(p(i,j) - p(i-1,j))
      end do
    end do
    do i=1,nx
      do j=1,ny+1
        a1 = .125_DP*(u(i,j-1) + u(i+1,j-1) + u(i,j) + u(i+1,j) + u0(i,j-1) + u0(i+1,j-1) + u0(i,j) + u0(i+1,j))
        a2 = .5_DP*(v(i,j) + v0(i,j))
        fv(i,j) = rhsv(i,j) & 
                - r%re_2x*a1*(v(i+1,j) - v(i-1,j) + v0(i+1,j) - v0(i-1,j)) & 
                - r%re_2y*a2*(v(i,j+1) - v(i,j-1) + v0(i,j+1) - v0(i,j-1)) & 
                - r%re2_y *(p(i,j) - p(i,j-1))
      end do
    end do
  else if (met == 2) then
    do i=1,nx+1
      do j=1,ny
        a1 = 1.5_DP*u0(i,j) - .5_DP*u1(i,j)
        a2 = .375_DP*(v0(i-1,j) + v0(i-1,j+1) + v0(i,j) + v0(i,j+1)) - .125_DP*(v1(i-1,j) + v1(i-1,j+1) + v1(i,j) + v1(i,j+1))
        fu(i,j) = rhsu(i,j) & 
                - r%re_2x*a1*(u(i+1,j) - u(i-1,j)) & 
                - r%re_2y*a2*(u(i,j+1) - u(i,j-1)) & 
                - r%re2_x *(p(i,j) - p(i-1,j))
      end do
    end do
    do i=1,nx
      do j=1,ny+1
        a1 = .375_DP*(u0(i,j-1) + u0(i+1,j-1) + u0(i,j) + u0(i+1,j)) - .125_DP*(u1(i,j-1) + u1(i+1,j-1) + u1(i,j) + u1(i+1,j))
        a2 = 1.5_DP*v0(i,j) - .5_DP*v1(i,j)
        fv(i,j) = rhsv(i,j) & 
                - r%re_2x*a1*(v(i+1,j) - v(i-1,j)) & 
                - r%re_2y*a2*(v(i,j+1) - v(i,j-1)) & 
                - r%re2_y *(p(i,j) - p(i,j-1))
      end do
    end do
  else if (met == 3) then
    do i=1,nx+1
      do j=1,ny
        fu(i,j) = rhsu(i,j) & 
                - r%re2_x *(p(i,j) - p(i-1,j))
      end do
    end do
    do i=1,nx
      do j=1,ny+1
        fv(i,j) = rhsv(i,j) & 
                - r%re2_y *(p(i,j) - p(i,j-1))
      end do
    end do
  else if (met == 4) then
    do i=1,nx+1
      do j=1,ny
        fu(i,j) = rhsu(i,j) & 
                - r%re2_x *(p(i,j) - p(i-1,j))
      end do
    end do
    do i=1,nx
      do j=1,ny+1
        fv(i,j) = rhsv(i,j) & 
                - r%re2_y *(p(i,j) - p(i,j-1))
      end do
    end do
  else if (met == 5) then
    do i=1,nx+1
      do j=1,ny
        fu(i,j) = rhsu(i,j) & 
                - r%re2_x *(p(i,j) - p(i-1,j)) & 
                - fx(i,j)
      end do
    end do
    do i=1,nx
      do j=1,ny+1
        fv(i,j) = rhsv(i,j) & 
                - r%re2_y *(p(i,j) - p(i,j-1)) & 
                - fy(i,j)
      end do
    end do
  else if (met == 6) then
    do i=1,nx+1
      do j=1,ny
        a1 = .5_DP*(u(i,j) + u0(i,j))
        a2 = .125_DP*(v(i-1,j) + v(i-1,j+1) + v(i,j) + v(i,j+1) + v0(i-1,j) + v0(i-1,j+1) + v0(i,j) + v0(i,j+1))
        fu(i,j) = rhsu(i,j) & 
                - r%re_2x*a1*(u(i+1,j) - u(i-1,j) + u0(i+1,j) - u0(i-1,j)) & 
                - r%re_2y*a2*(u(i,j+1) - u(i,j-1) + u0(i,j+1) - u0(i,j-1)) & 
                - r%re2_x *(p(i,j) - p(i-1,j)) & 
                - fx(i,j)
      end do
    end do
    do i=1,nx
      do j=1,ny+1
        a1 = .125_DP*(u(i,j-1) + u(i+1,j-1) + u(i,j) + u(i+1,j) + u0(i,j-1) + u0(i+1,j-1) + u0(i,j) + u0(i+1,j))
        a2 = .5_DP*(v(i,j) + v0(i,j))
        fv(i,j) = rhsv(i,j) & 
                - r%re_2x*a1*(v(i+1,j) - v(i-1,j) + v0(i+1,j) - v0(i-1,j)) & 
                - r%re_2y*a2*(v(i,j+1) - v(i,j-1) + v0(i,j+1) - v0(i,j-1)) & 
                - r%re2_y *(p(i,j) - p(i,j-1)) & 
                - fy(i,j)
      end do
    end do
  else
    write(*,'(A,I1,A)') 'Wrong method (met = ', met,')!'
    stop
  end if

  ! Allocate arrayes
  allocate(ax(0:ny+1,1:3))
  allocate(bx(0:ny+1,1:3))
  allocate(ay(0:nx+2,1:3))
  allocate(by(0:nx+2,1:3))

  ! Set up boundary conditions
  do i=0,ny+1
    ax(i,1:3) = (/ 0._DP,  1._DP, 2._DP /) ! Inflow boundary condition: u = 1
    bx(i,1:3) = (/ 0._DP, -1._DP, 0._DP /) ! Outflow boundary condition: du/dn = 0 (soft)
  end do
  do i=0,nx+2
    ay(i,1:3) = (/ -1._DP, 0._DP, 0._DP /) ! Bottom boundary condition: du/dn = 0 (soft)
    by(i,1:3) = (/ -1._DP, 0._DP, 0._DP /) ! Top boundary condition: du/dn = 0 (soft)
  end do

  ! Solve Boundary value problem for Helmholtz equation by Gauss–Seidel method (it is faster here than SOR method)
  call helmholtz_solve_sor(u,ax,bx,ay,by,fu,r%re2_t,nx+1,ny,hx,hy,200000,1.E-11_DP,deb-1,1._DP)

  ! Deallocate arrayes
  deallocate(ax)
  deallocate(bx)
  deallocate(ay)
  deallocate(by)

  ! Allocate arrayes
  allocate(ax(0:ny+2,1:3))
  allocate(bx(0:ny+2,1:3))
  allocate(ay(0:nx+1,1:3))
  allocate(by(0:nx+1,1:3))

  ! Set up boundary conditions
  do i=0,ny+2
    ax(i,1:3) = (/ 1._DP, 0._DP, 0._DP /)  ! Inflow boundary condition: v = 0
    bx(i,1:3) = (/ -1._DP, 0._DP, 0._DP /) ! Outflow boundary condition: dv/dn = 0 (soft)
  end do
  do i=0,nx+1
    ay(i,1:3) = (/ 0._DP, 1._DP, 0._DP /) ! Bottom boundary condition: v = 0 (no-slip)
    by(i,1:3) = (/ 0._DP, 1._DP, 0._DP /) ! Top boundary condition: v = 0 (no-slip)
  end do

  ! Solve Boundary value problem for Helmholtz equation by Gauss–Seidel method (it is faster here than SOR method)
  call helmholtz_solve_sor(v,ax,bx,ay,by,fv,r%re2_t,nx,ny+1,hx,hy,200000,1.E-11_DP,deb-1,1._DP)

  ! Deallocate arrayes
  deallocate(ax)
  deallocate(bx)
  deallocate(ay)
  deallocate(by)

  ! Deallocate arrayes
  deallocate(fu)
  deallocate(fv)

end subroutine predict_uv

!> subroutine to compute prediction velocity with added force terms
subroutine force_add(u,v,rhsu,rhsv,p,fx,fy,nx,ny,hx,hy,r)
implicit none

real(DP), allocatable :: u(:,:), v(:,:), p(:,:)                                 ! Arrayes of variables
real(DP), allocatable :: rhsu(:,:), rhsv(:,:)                                   ! Arrayes of variables
real(DP), allocatable :: fx(:,:), fy(:,:)                                       ! Arrayes for force terms
integer :: nx, ny                                                               ! Numbers of cells in x and y directions
real(DP):: hx, hy                                                               ! Space and time steps
type(sratios) :: r                                                              ! Step ratios

integer :: i, j                                                                 ! Loop indices
real(DP), allocatable :: fu(:,:), fv(:,:)                                       ! Right hand sides
real(DP), allocatable :: ax(:,:), bx(:,:), ay(:,:), by(:,:)                     ! Values of boundary conditions

  ! Allocate arrayes
  allocate(fu(0:nx+2,0:ny+1))
  allocate(fv(0:nx+1,0:ny+2))

  ! Compute right hand sides
  fu(0,0:ny+1)    = 0._DP ! These 8 lines added for consistency
  fu(nx+2,0:ny+1) = 0._DP
  fu(1:nx+1,0)    = 0._DP
  fu(1:nx+1,ny+1) = 0._DP
  fv(0,0:ny+2)    = 0._DP
  fv(nx+1,0:ny+2) = 0._DP
  fv(1:nx,0)      = 0._DP
  fv(1:nx,ny+2)   = 0._DP
  if (met == 4) then
    do i=1,nx+1
      do j=1,ny
        fu(i,j) = rhsu(i,j) & 
                - r%re2_x*(p(i,j) - p(i-1,j)) & 
                - r%re2  *fx(i,j)
      end do
    end do
    do i=1,nx
      do j=1,ny+1
        fv(i,j) = rhsv(i,j) & 
                - r%re2_y*(p(i,j) - p(i,j-1)) & 
                - r%re2  *fy(i,j)
      end do
    end do
  else
    write(*,'(A,I1,A)') 'Wrong method (met = ', met,')!'
    stop
  end if

  ! Allocate arrayes
  allocate(ax(0:ny+1,1:3))
  allocate(bx(0:ny+1,1:3))
  allocate(ay(0:nx+2,1:3))
  allocate(by(0:nx+2,1:3))

  ! Set up boundary conditions
  do i=0,ny+1
    ax(i,1:3) = (/ 0._DP,  1._DP, 2._DP /) ! Inflow boundary condition: u = 1
    bx(i,1:3) = (/ 0._DP, -1._DP, 0._DP /) ! Outflow boundary condition: du/dn = 0 (soft)
  end do
  do i=0,nx+2
    ay(i,1:3) = (/ -1._DP, 0._DP, 0._DP /) ! Bottom boundary condition: du/dn = 0 (soft)
    by(i,1:3) = (/ -1._DP, 0._DP, 0._DP /) ! Top boundary condition: du/dn = 0 (soft)
  end do

  ! Solve Boundary value problem for Helmholtz equation by Gauss–Seidel method (it is faster here than SOR method)
  call helmholtz_solve_sor(u,ax,bx,ay,by,fu,r%re2_t,nx+1,ny,hx,hy,200000,1.E-11_DP,deb-1,1._DP)

  ! Deallocate arrayes
  deallocate(ax)
  deallocate(bx)
  deallocate(ay)
  deallocate(by)

  ! Allocate arrayes
  allocate(ax(0:ny+2,1:3))
  allocate(bx(0:ny+2,1:3))
  allocate(ay(0:nx+1,1:3))
  allocate(by(0:nx+1,1:3))

  ! Set up boundary conditions
  do i=0,ny+2
    ax(i,1:3) = (/ 1._DP, 0._DP, 0._DP /)  ! Inflow boundary condition: v = 0
    bx(i,1:3) = (/ -1._DP, 0._DP, 0._DP /) ! Outflow boundary condition: dv/dn = 0 (soft)
  end do
  do i=0,nx+1
    ay(i,1:3) = (/ 0._DP, 1._DP, 0._DP /) ! Bottom boundary condition: v = 0 (no-slip)
    by(i,1:3) = (/ 0._DP, 1._DP, 0._DP /) ! Top boundary condition: v = 0 (no-slip)
  end do

  ! Solve Boundary value problem for Helmholtz equation by Gauss–Seidel method (it is faster here than SOR method)
  call helmholtz_solve_sor(v,ax,bx,ay,by,fv,r%re2_t,nx,ny+1,hx,hy,200000,1.E-11_DP,deb-1,1._DP)

  ! Deallocate arrayes
  deallocate(ax)
  deallocate(bx)
  deallocate(ay)
  deallocate(by)

  ! Deallocate arrayes
  deallocate(fu)
  deallocate(fv)

end subroutine force_add

!> subroutine to solve Poisson problem for pressure correction by Successive
!                                                   over-relaxation (SOR) method
subroutine solve_pcp(pc,rhs,u,v,p,nx,ny,hx,hy,tau,Re,t,s,r)
implicit none

real(DP), allocatable :: pc(:,:), rhs(:,:), u(:,:), v(:,:), p(:,:)              ! Arrayes of variables
integer :: nx, ny                                                               ! Numbers of cells in x and y directions
real(DP):: hx, hy, tau                                                          ! Space and time steps
real(DP) :: Re                                                                  ! Reynolds number
real(DP) :: t                                                                   ! Time
integer :: s                                                                    ! NUmber of iteration
type(sratios) :: r                                                              ! Step ratios

real(DP), allocatable :: ax(:,:), bx(:,:), ay(:,:), by(:,:)                     ! Values of boundary conditions
integer :: i, j                                                                 ! Loop indices
real(DP):: rx, ry                                                               ! Ratios
character(4) :: bc                                                              ! Types of boundary conditions

  ! Choose rx and ry depending on the method
  rx = r%e_xt
  ry = r%e_yt

  ! Allocate arrayes
  allocate(ax(0:ny+1,1:3))
  allocate(bx(0:ny+1,1:3))
  allocate(ay(0:nx+1,1:3))
  allocate(by(0:nx+1,1:3))

  ! Set up boundary conditions
  bc = 'NDNN'
!  bc = 'NNNN'
  do i=0,ny+1
    ax(i,1:3) = (/ -1._DP, 0._DP, 0._DP /) ! Inflow boundary condition: d(pc)/dn = 0 (dp/dn = 0)
    bx(i,1:3) = (/ 1._DP, 0._DP, 0._DP /)  ! Outflow boundary condition: pc = 0 (no correction)
!    bx(i,1:3) = (/ -1._DP, 0._DP, 0._DP /)  ! ***
  end do
  do i=0,nx+1
    ay(i,1:3) = (/ -1._DP, 0._DP, 0._DP /) ! Bottom boundary condition: d(pc)/dn = 0 (dp/dn = 0)
    by(i,1:3) = (/ -1._DP, 0._DP, 0._DP /) ! Top boundary condition: d(pc)/dn = 0 (dp/dn = 0)
  end do
!  bx(ny/2,1:3) = (/ 1._DP, 0._DP, 0._DP /)  ! ***

  ! Set up array of righ hand side
  rhs(0:nx+1,0)    = 0._DP ! These 4 lines added for consistency
  rhs(0:nx+1,ny+1) = 0._DP
  rhs(0,1:ny)      = 0._DP
  rhs(nx+1,1:ny)   = 0._DP
  do i=1,nx
    do j=1,ny
      rhs(i,j) = rx*(u(i,j) - u(i+1,j)) + ry*(v(i,j) - v(i,j+1)) ! Negative value due to the statement: - div(grad pc) = rhs(x,y)
    end do
  end do

  ! Solve Boundary value problem for Helmholtz equation
  call helmholtz_solve(pc,ax,bx,ay,by,rhs,bc,0._DP,nx,ny,hx,hy,200000,1.E-11_DP,deb-1)

  ! Deallocate arrayes
  deallocate(ax)
  deallocate(bx)
  deallocate(ay)
  deallocate(by)

end subroutine solve_pcp

!> subroutine to correct velocity field after pressure correction
subroutine correct_vel(u,v,pc,nx,ny,hx,hy,tau,t,r)
implicit none

real(DP), allocatable :: u(:,:), v(:,:), pc(:,:)                                ! Arrayes of variables
integer :: nx, ny                                                               ! Numbers of cells in x and y directions
real(DP):: hx, hy, tau                                                          ! Space and time steps
real(DP) :: t
type(sratios) :: r                                                              ! Step ratios

integer :: i, j                                                                 ! Loop indices
real(DP):: rx, ry                                                               ! Ratios

  ! Choose rx and ry depending on the method
  rx = r%t_x
  ry = r%t_y

  ! Correction of the main field of u
  do i=1,nx+1
    do j=1,ny
      u(i,j) = u(i,j) - rx*(pc(i,j) - pc(i-1,j))
    end do
  end do


  ! Apply boundary conditions
if (.false.) then
  do i=0,ny+1
    u(0,i)    = 2._DP - u(2,i)  ! Inflow boundary condition: u = 1
    u(nx+2,i) = u(nx,i)         ! Outflow boundary condition: du/dn = 0 (soft)
  end do
  do i=0,nx+2
    u(i,0)    = u(i,1)  ! Bottom boundary condition: du/dn = 0 (soft)
    u(i,ny+1) = u(i,ny) ! Top boundary condition: du/dn = 0 (soft)
  end do
end if

  ! Correction of the main field of v
  do i=1,nx
    do j=1,ny+1
      v(i,j) = v(i,j) - ry*(pc(i,j) - pc(i,j-1))
    end do
  end do

  ! Apply boundary conditions
if (.false.) then
  do i=0,ny+2
    v(0,i)    = - v(1,i)   ! Inflow boundary condition: v = 0
    v(nx+1,i) =   v(nx,i)  ! Outflow boundary condition: dv/dn = 0 (soft)
  end do
  do i=0,nx+1
    v(i,0)    = - v(i,2)  ! Bottom boundary condition: v = 0 (no-slip)
    v(i,ny+2) = - v(i,ny) ! Top boundary condition: v = 0 (no-slip)
  end do
end if

end subroutine correct_vel

!> function to check the rightness of the solution
function residuals(u,v,u0,v0,p,nx,ny,hx,hy,tau,Re)
implicit none

real(DP), allocatable :: u(:,:), v(:,:), u0(:,:), v0(:,:), p(:,:)               ! Arrayes of variables
integer :: nx, ny                                                               ! Numbers of cells in x and y directions
real(DP):: hx, hy, tau                                                          ! Space and time steps
real(DP) :: Re                                                                  ! Reynolds number
real(DP) :: residuals

real(DP) :: rt, rx, ry, rx_2, ry_2, rx_4, ry_4, rx_16, ry_16, rxx_2, ryy_2      ! Step ratios
real(DP) :: eq1, eq2, eq3, c                                                    ! Residuals
integer :: i                                                                    ! Loop index

  ! Compute step ratios
  rt    = 1._DP/tau
  rx    = 1._DP/hx
  ry    = 1._DP/hy
  rx_2  = .5_DP*rx
  ry_2  = .5_DP*ry
  rx_4  = .25_DP*rx
  ry_4  = .25_DP*ry
  rx_16 = .0625_DP*rx
  ry_16 = .0625_DP*ry
  rxx_2 = .5_DP/(hx*hx*Re)
  ryy_2 = .5_DP/(hy*hy*Re)

  ! Compute right hand sides
  eq1 = 0._DP
  do i=1,nx+1
    do j=1,ny
      c = rt*(u(i,j) - u0(i,j)) & 
        + rx_2 *(1.5_DP*(u0(i+1,j)**2 - u0(i-1,j)**2) - .5_DP*(u1(i+1,j)**2 - u1(i-1,j)**2)) & 
        + ry_4 *(1.5_DP*((u0(i,j) + u0(i,j+1))*(v0(i,j+1) + v0(i-1,j+1)) - (u0(i,j) + u0(i,j-1))*(v0(i,j) + v0(i-1,j))) & 
                 -.5_DP*((u1(i,j) + u1(i,j+1))*(v1(i,j+1) + v1(i-1,j+1)) - (u1(i,j) + u1(i,j-1))*(v1(i,j) + v1(i-1,j)))) & 
        - rxx_2*(u0(i+1,j) - 2._DP*u0(i,j) + u0(i-1,j)) & 
        - rxx_2*(u(i+1,j)  - 2._DP*u(i,j)  + u(i-1,j)) & 
        - ryy_2*(u0(i,j+1) - 2._DP*u0(i,j) + u0(i,j-1)) & 
        - ryy_2*(u(i,j+1)  - 2._DP*u(i,j)  + u(i,j-1)) & 
        + rx_2 *(p0(i,j) - p0(i-1,j)) & 
        + rx_2 *(p(i,j)  - p(i-1,j))
      eq1 = max(eq1,abs(c))
    end do
  end do

  eq2 = 0._DP
  do i=1,nx
    do j=1,ny+1
      c = rt*(v(i,j) - v0(i,j)) & 
        + rx_4 *(1.5_DP*((u0(i+1,j) + u0(i+1,j-1))*(v0(i,j) + v0(i+1,j)) - (u0(i,j) + u0(i,j-1))*(v0(i,j) + v0(i-1,j))) & 
                 -.5_DP*((u1(i+1,j) + u1(i+1,j-1))*(v1(i,j) + v1(i+1,j)) - (u1(i,j) + u1(i,j-1))*(v1(i,j) + v1(i-1,j)))) & 
        + ry_2 *(1.5_DP*(v0(i,j+1)**2 - v0(i,j-1)**2) - .5_DP*(v1(i,j+1)**2 - v1(i,j-1)**2)) & 
        - rxx_2*(v0(i+1,j) - 2._DP*v0(i,j) + v0(i-1,j)) & 
        - rxx_2*(v(i+1,j)  - 2._DP*v(i,j)  + v(i-1,j)) & 
        - ryy_2*(v0(i,j+1) - 2._DP*v0(i,j) + v0(i,j-1)) & 
        - ryy_2*(v(i,j+1)  - 2._DP*v(i,j)  + v(i,j-1)) & 
        + ry_2 *(p0(i,j) - p0(i,j-1)) & 
        + ry_2 *(p(i,j)  - p(i,j-1))
      eq2 = max(eq2,abs(c))
    end do
  end do

  eq3 = 0._DP
  do i=1,nx
    do j=1,ny
      c = rx*(u(i+1,j) - u(i,j)) & 
        + ry*(v(i,j+1) - v(i,j))
      eq3 = max(eq3,abs(c))
    end do
  end do

  residuals = max(abs(eq1),abs(eq2),abs(eq3))

end function residuals
subroutine move_particles(particles,np,xl,ul,vl,fl,lx,ly,hx,hy,hl,tau,Re,Bu,gamma)
implicit none

type(tparticle), allocatable :: particles(:)                                    ! Arrayes of particles
integer :: np                                                                   ! Number of particles
real(DP), allocatable :: xl(:,:), ul(:), vl(:)                                  ! Arrayes of Lagrangian variables (coordinates and velocities)
real(DP), allocatable :: fl(:,:)                                                ! Computed Lagrangian forces
real(DP) :: lx, ly                                                              ! Size of computational domain
real(DP):: hx, hy, hl, tau                                                      ! Space and time steps
real(DP) :: Re                                                                  ! Reynolds number
real(DP) :: Bu                                                                  ! Relative difference of particle density with the fluid one
real(DP) :: gamma                                                               ! Friction coefficient

integer :: i, j, k, it                                                          ! Loop indices
real(DP) :: r(1:2), rotv                                                        ! Vector and velocity induced by rotation
integer :: ndem                                                                 ! Number of intervals on which time step is splitted
real(DP) :: tau0, ttb, ttb_min                                                  ! Time step for DEM and for take back
real(DP) :: gap1, gap2                                                          ! Minimal distances between particles
real(DP) :: d, rc                                                               ! Distance and limiting distance for two particles
real(DP) :: norm(1:2)                                                           ! Normal direction
real(DP) :: veln_i, veln_j                                                      ! Normal velocities of the particles
real(DP) :: St, ed, mass                                                        ! Stokes number, dissipation of energy and mass
real(DP) :: xx(1:2)                                                             ! Vector between centers of particles
real(DP) :: ui, vi                                                              ! Velocity of i-th particle
real(DP), allocatable :: x0(:,:)                                                ! Coordinates of particles from previous time step
integer :: ctype                                                                ! Collision type
real(DP) :: omegai(1:2), omegaj(1:2)                                            ! Rotation velocities for friction
real(DP) :: beta, beta_r                                                        ! Updated friction coefficient

  ! Allocate arrayes and update friction coefficient
  if (.not. allocated(cont)) then
    allocate(cont(1:np,1:16))
    allocate(tcont(1:np,1:16))
    cont = 0
    tcont = 0
  end if
  allocate(x0(1:np,1:2))
  do i=1,np
    x0(i,1:2) = (/ particles(i)%x, particles(i)%y /)
  end do
  if (gamma < 0.99_DP) then
    beta = gamma/(1._DP - gamma)
  else
    beta = 200._DP
  end if

  ! Compute force and moment of force acting onto the particle
  call force_particles(particles,np,xl,fl,tau)

  ! Set main parameters of DEM
  gap1 = .1_DP*min(hx,hy) ! 0.1 from mesh steps
  gap2 = .2_DP*min(hx,hy) ! 0.2 from mesh steps
  ndem = int(5._DP*tau/(gap2 - gap1)) ! Five (Two) small time steps for distance between gaps
  tau0 = tau/ndem
  beta_r = .5_DP*tau0*beta

  ! Particle movement and collision loop
  do it=1,ndem
    do i=1,np

      ! Save velocities
      ui = particles(i)%u
      vi = particles(i)%v


      ! Update the position of particle
      particles(i)%x = particles(i)%x + tau0*ui
      particles(i)%y = particles(i)%y + tau0*vi

      ttb_min = - tau0

      ! Process particle-wall collisions
      if (particles(i)%y < particles(i)%r .or. particles(i)%y > 1._DP - particles(i)%r) then

        ! Calculate disipated energy and compute Stokes number [St = (rho_p + .5*rho_f)*v_app*d/(9*mu)]
        mass = particles(i)%Bu*(1._DP + 1._DP/Bu)
        ed = .5_DP*mass*particles(i)%v**2
        St = abs(particles(i)%v)*particles(i)%r*Re*(Bu + 1.5_DP)*2._DP/9._DP

        ! Change coordinate and velocity
        if (particles(i)%y < particles(i)%r) then ! Contact with bottom wall
          particles(i)%y = particles(i)%r
          particles(i)%v = 0._DP
          ctype = 0
        end if
        if (particles(i)%y > 1._DP - particles(i)%r) then ! Contact with top wall
          particles(i)%y = 1._DP - particles(i)%r
          particles(i)%v = 0._DP
          ctype = 1
        end if

        write(7051,'(A,I4,A,I1,A,D,X,F8.2,X,2(X,F9.5))') 'Wall[', i, ']<', ctype, '>:', ed, St, & 
                                                         particles(i)%x, particles(i)%y

      end if

      ! Process particle-particle collisions
      do k=1,16 ! Process existing particle collisions
        if (cont(i,k) > 0) then
          j = cont(i,k)
          rc = particles(i)%r + particles(j)%r + gap2
          xx(1:2) = (/ particles(j)%x - particles(i)%x, particles(j)%y - particles(i)%y /)
          d = dist(xx,lx,ctype)
          if (d < rc) then ! Check is particles in contact
            tcont(i,k) = 10
            rc = particles(i)%r + particles(j)%r + gap1
            if (d < rc) then

              ! Estimate approaching velocity
              norm = vect_xx(xx,lx,ctype)/d
              veln_i = particles(i)%u*norm(1) + particles(i)%v*norm(2)
              veln_j = particles(j)%u*norm(1) + particles(j)%v*norm(2)

              if (veln_i - veln_j > 0._DP) then ! Approaching velosity is positive

                ! Take back particle if it is possible
                if (veln_i > 0._DP) then
                  ttb = max((d - rc)/veln_i,ttb_min) ! Upper estimate for take back
                  particles(i)%x = particles(i)%x + ttb*ui
                  particles(i)%y = particles(i)%y + ttb*vi

                  ttb_min = min(ttb_min - ttb,0._DP)
                  xx(1:2) = (/ particles(j)%x - particles(i)%x, particles(j)%y - particles(i)%y /)
                  d = dist(xx,lx,ctype)

                  ! Correct approaching velocity
                  norm = vect_xx(xx,lx,ctype)/d
                  veln_i = particles(i)%u*norm(1) + particles(i)%v*norm(2)
                  veln_j = particles(j)%u*norm(1) + particles(j)%v*norm(2)

                end if

                ! Compute energy before collision
                ed = (particles(i)%u**2 + particles(i)%v**2) & 
                   + (particles(j)%u**2 + particles(j)%v**2)

                ! Correct velocities
                particles(i)%u = particles(i)%u - .5_DP*(veln_i - veln_j)*norm(1)
                particles(i)%v = particles(i)%v - .5_DP*(veln_i - veln_j)*norm(2)

                particles(j)%u = particles(j)%u - .5_DP*(veln_j - veln_i)*norm(1)
                particles(j)%v = particles(j)%v - .5_DP*(veln_j - veln_i)*norm(2)


                ! Account energy after collision and compute Stokes number [St = (rho_p + .5*rho_f)*v_app*d/(9*mu)]
                mass = particles(i)%Bu*(1._DP + 1._DP/Bu)
                ed = ed - (particles(i)%u**2 + particles(i)%v**2) & 
                        - (particles(j)%u**2 + particles(j)%v**2)
                ed = .5_DP*mass*ed
                St = (veln_i - veln_j)*particles(i)%r*Re*(Bu + 1.5_DP)*2._DP/9._DP
                write(7050,'(A,I4,A,I4,A,I1,A,D,X,F8.2,X,F8.5,3(X,3(X,F9.5)),I)') 'Overlap[', i, ',', j, ']<',ctype,'>:', ed, St, & 
                                                                                  d, norm(1:2), & 
                                                                                  particles(i)%x, particles(i)%y, & 
                                                                                  particles(j)%x, particles(j)%y
              end if
            end if
          else
            tcont = tcont(i,k) - 1
          end if
        end if
      end do

      ! Process tcont array
      k = 1
      do while (k <= 16)
        if (tcont(i,k) == 0 .and. cont(i,k) > 0) then
          cont(i,k:15)  = cont(i,(k+1):16)
          tcont(i,k:15) = tcont(i,(k+1):16)
          cont(i,16)  = 0
          tcont(i,16) = 0
        else
          k = k + 1
        end if
      end do

      ! Find new contacts with other particles
      do j=1,i-1
        if (minval(abs(cont(i,1:16) - j)) > 0) then
          rc = particles(i)%r + particles(j)%r + gap2
          xx(1:2) = (/ particles(j)%x - particles(i)%x, particles(j)%y - particles(i)%y /)
          d = dist(xx,lx,ctype)
          if (d < rc) then
            do k=1,16
              if (cont(i,k) == 0) then
                cont(i,k)  = j
                tcont(i,k) = 10
                exit
              end if
            end do
          end if
        end if
      end do
    end do
  end do

  ! Friction loop
  do i=1,np

    ! Particle-wall friction
    if (particles(i)%y < particles(i)%r + gap2) then
      if (beta < 100._DP) then
        particles(i)%omega(1) = particles(i)%omega(1)*(1 - beta_r)

      else
        particles(i)%omega(1) =  0! Absolute friction
      end if
    end if
    if (particles(i)%y > 1._DP - particles(i)%r - gap2) then
      if (beta < 100._DP) then
        particles(i)%omega(1) = particles(i)%omega(1)*(1 - beta_r)
      else
        particles(i)%omega(1) =  0 ! Absolute friction
      end if
    end if

    ! Particle-particle friction
    do k=1,16
      if (cont(i,k) > 0) then
        j = cont(i,k)

        ! Compute vector in the collision direction
        xx(1:2) = (/ particles(j)%x - particles(i)%x, particles(j)%y - particles(i)%y /)
        d = dist(xx,lx,ctype)
        norm = vect_xx(xx,lx,ctype)/d

        ! Compute rotation velocity for friction
        omegai = particles(i)%omega - dot_product(particles(i)%omega,norm)*norm
        omegaj = particles(j)%omega - dot_product(particles(j)%omega,norm)*norm

        ! Update rotation velocity
        if (beta < 100._DP) then
          particles(i)%omega = particles(i)%omega - beta_r*(omegai + omegaj)
          particles(j)%omega = particles(j)%omega - beta_r*(omegai + omegaj)
        else
          particles(i)%omega = particles(i)%omega - .5_DP*(omegai + omegaj) ! Absolute friction
          particles(j)%omega = particles(j)%omega - .5_DP*(omegai + omegaj)
        end if

      end if
    end do

  end do

  !$OMP PARALLEL DO PRIVATE(i,j,xx,r,rotv)
  do i=1,np

    ! Periodical transfer of the particle in x and z direction
    if (particles(i)%x >= lx) then ! Particle cross boundary x = lx
      particles(i)%x = particles(i)%x - lx
    end if
    if (particles(i)%x < 0._DP) then ! Particle cross boundary x = 0
      particles(i)%x = particles(i)%x + lx
    end if

    ! Update the position of nodes of the Lagrangian mesh
    xx(1:2) = (/ particles(i)%x, particles(i)%y/) - x0(i,1:2)
    do j=particles(i)%na,particles(i)%nb
      xl(j,1:2) = xl(j,1:2) + xx(1:2)
    end do

    ! Update velocities of nodes of the Lagrangian mesh
    xx(1:2) = (/ particles(i)%x, particles(i)%y/)
    do j=particles(i)%na,particles(i)%nb
      r = xl(j,1:2) - xx(1:2)
      rotv = X_product(r,particles(i)%omega)
      ul(j) = particles(i)%u + rotv
      vl(j) = particles(i)%v + rotv
    end do

  end do
  !$OMP END PARALLEL DO

  ! Deallocate array
  deallocate(x0)

end subroutine move_particles
subroutine force_particles(particles,np,xl,fl,tau)
implicit none

type(tparticle), allocatable :: particles(:)                                    ! Arrayes of particles
integer :: np                                                                   ! Number of particles
real(DP), allocatable :: xl(:,:)                                                ! Arrayes of Lagrangian coordinates
real(DP), allocatable :: fl(:,:)                                                ! Computed Lagrangian forces
real(DP):: tau                                                                  ! Time steps

integer :: i, j                                                                 ! Loop indices
real(DP) :: r(1:2)                                                              ! Vector and velocity induced by rotation
real(DP) :: f(1:2), mf(1:2)                                                     ! Force and moment of force
real(DP) :: xx(1:2)                                                             ! Vector between centers of particles

  ! Main loop
  !$OMP PARALLEL DO PRIVATE(i,j,r,f,mf,xx)
  do i=1,np

    f  = 0._DP
    mf = 0._DP
    xx(1:2) = (/ particles(i)%x, particles(i)%y /)
    do j=particles(i)%na,particles(i)%nb
      r = xx(1:2) - xl(j,1:2) !? TODO: Explain, why this [-r x f] is right.
      f = f + fl(j,1:2)
      mf = mf + X_product(r,fl(j,1:2))
    end do
    particles(i)%f(1:2) = f
    particles(i)%mf(1:2) = mf

    ! Update velocity and rotation velocity
    particles(i)%u = particles(i)%u + tau*particles(i)%f(1)/particles(i)%Bu
    particles(i)%v = particles(i)%v + tau*particles(i)%f(2)/particles(i)%Bu
    particles(i)%omega(1:2) = particles(i)%omega(1:2) + tau*particles(i)%mf(1:2)/(particles(i)%Bu*particles(i)%I)

  end do
  !$OMP END PARALLEL DO

end subroutine force_particles

end program ibm2d_test