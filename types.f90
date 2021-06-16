!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Module to implement basic types and constants
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program: ibm2d_cyl -- SIMPLE + IBM method for incomressible viscous Navier -
! Stockes equations in 2D statement for flow past a cylinder
! Copyright (C): Denis V. Esipov (2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module types
implicit none

  ! Global parameters
  integer, parameter :: DP = 8                                                  ! Precision of reals
  real(DP), parameter :: pi = 3.14159265358979323846264338327950288_DP          ! Pi number
  integer, parameter :: met = 1                                                 ! Method
                                                                                ! 1 - Original method like Crank–Nicolson method
                                                                                ! 2 - Adams-Bashforth approximation for advection term
                                                                                ! 3 - Adams-Bashforth approximation for convective term
                                                                                ! 4 - 3th method with the explicit addition of force's values
                                                                                ! 5 - 7th method with separate calculation of force
                                                                                ! 6 - experimental method (does not work)
  integer, parameter :: md = 1                                                  ! Forcing
                                                                                ! 0 - Exact forcing
                                                                                ! 1 - Direct forcing
                                                                                ! 2 - Multidirect forcing
  integer, parameter :: mkl = 1                                                 ! Use mkl solvers or not
                                                                                ! 0 - SOR method
                                                                                ! 1 - MKL solver
                                                                                ! 2 - Hybrid MKL+SOR
                                                                                ! 3 - Pardiso solver

  ! Structure to store all step ratios
  type sratios
    real(DP) :: re2_t, re2_x, re2_y, re_x, re_y, re_2x, re_2y
    real(DP) :: e_xx, e_yy, e_xt, e_yt, re2, t_x, t_y, t_2x, t_2y
  end type sratios

contains

!> subroutine to compute all needed step ratios
subroutine step_ratios(r,hx,hy,tau,Re)
implicit none

type(sratios) :: r                                                              ! Step ratios
real(DP):: hx, hy, tau                                                          ! Space and time steps
real(DP) :: Re                                                                  ! Reynolds number

  r%re2_t = 2._DP*Re/tau
  r%re2_x = 2._DP*Re/hx
  r%re2_y = 2._DP*Re/hy
  r%re_x  = Re/hx
  r%re_y  = Re/hy
  r%re_2x = Re/(2._DP*hx)
  r%re_2y = Re/(2._DP*hy)
  r%e_xx  = 1._DP/(hx*hx)
  r%e_yy  = 1._DP/(hy*hy)
  r%e_xt  = 1._DP/(tau*hx)
  r%e_yt  = 1._DP/(tau*hy)
  r%re2   = 2._DP*Re
  r%t_x   = tau/hx
  r%t_y   = tau/hy
  r%t_2x  = .5_DP*tau/hx
  r%t_2y  = .5_DP*tau/hy

end subroutine

end module types
