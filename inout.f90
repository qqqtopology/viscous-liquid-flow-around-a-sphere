!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Module to implement basic types and constants
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program: ibm2d_cyl -- SIMPLE + IBM method for incomressible viscous Navier -
! Stockes equations in 2D statement for flow past a cylinder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module inout
use types
implicit none

  public save_plt, save_dat, load_dat

contains

!> subroutine to write result into .plt file
subroutine save_plt(u,v,p,xl,ul,vl,nx,ny,nl,it,hx,hy,tau,strand,Re)
implicit none

real(DP), allocatable :: u(:,:), v(:,:), p(:,:), xl(:,:), ul(:), vl(:)          ! Arrayes of variables
integer :: nx, ny, nl                                                           ! Numbers of cells in x and y directions and Lagrangian nodes
integer :: it                                                                   ! Number of time step
real(DP):: hx, hy, tau                                                          ! Space and time steps
integer :: strand                                                               ! StrandID
real(DP) :: Re                                                                  ! Reynolds number

integer :: i, j                                                                 ! Loop indices

  ! Update StrandID
  strand = strand + 1

  ! Open file
  if (it == 0) then
    open(1001, file='time.plt', status='replace', action='write')
    write(1001,'(A,F9.2,A)') 'Title = " Re = ', Re, '"'
    write(1001,'(A)') 'Variables = "x", "y", "u", "v", "p"'
  else
    open(1001, file='time.plt', status='old', action='write', position='append')
  end if

  ! Write Eulerian mesh
  write(1001,'(A,I4,A,F9.5,A,I4,A,I4,A)') 'Zone T = "Flow", StrandID = ', strand, ', SolutionTime = ', it*tau, ', I = ', ny + 1, ', J = ', nx + 1, ', Datapacking = Point'
  do i=0,nx
    do j=0,ny
      write(1001,'(4(E12.5,X),E12.5)') i*hx, j*hy, & 
                                       .5_DP*(u(i+1,j) + u(i+1,j+1)), & ! Interpolate variables onto uniform mesh
                                       .5_DP*(v(i,j+1) + v(i+1,j+1)), & 
                                       .25_DP*(p(i,j) + p(i+1,j) + p(i,j+1) + p(i+1,j+1))
    end do
  end do

  ! Write Lagrangian mesh
  if (nl > 0) then
    write(1001,'(A,I4,A,F9.5,A,I4,A)') 'Zone T = "Body", StrandID = ', strand,', SolutionTime = ', it*tau, ', I = ', nl + 1, ', Datapacking = Point'
    do i=1,nl
      write(1001,'(4(E12.5,X),E12.5)') xl(i,1), xl(i,2), ul(i), vl(i), 0._DP
    end do
    write(1001,'(4(E12.5,X),E12.5)') xl(1,1), xl(1,2), ul(1), vl(1), 0._DP ! Close the line
  end if

  ! Close file
  close(1001)

end subroutine save_plt

!> subroutine to write last field
subroutine save_dat(u,v,p,u0,v0,xl,ul,vl,nx,ny,nl,hx,hy,hl,it,tau,Re)
implicit none

real(DP), allocatable :: u(:,:), v(:,:), p(:,:), u0(:,:), v0(:,:)               ! Arrayes of Eulerian variables
real(DP), allocatable :: xl(:,:), ul(:), vl(:)                                  ! Arrayes of Lagrangian variables
integer :: nx, ny, nl                                                           ! Numbers of cells in x and y directions and Lagrangian nodes
integer :: it                                                                   ! Number of time step
real(DP):: hx, hy, hl, tau                                                      ! Space and time steps
real(DP) :: Re                                                                  ! Reynolds number

integer :: i, j                                                                 ! Loop indices
real(DP) :: time

  ! Calculate time
  time = it*tau

  ! Open file
  open(1002, file='last.dat', status='replace', action='write')
  write(1002,'(F9.2,A)') Re,    '                  = Re'
  write(1002,'(I4,A)') nx, '                       = nx'
  write(1002,'(I4,A)') ny, '                       = ny'
  write(1002,'(I4,A)') nl, '                       = nl'
  write(1002,'(D,A)') hx,                       '  = hx'
  write(1002,'(D,A)') hy,                       '  = hy'
  write(1002,'(D,A)') hl,                       '  = hl'
  write(1002,'(D,A)') time,                   '  = time'

  ! Write Eulerian mesh
  do i=0,nx+2 ! Write velocities u, u0
    do j=0,ny+1
      write(1002,'(I4,X,I4,X,2D)') i, j, u(i,j), u0(i,j)
    end do
  end do
  do i=0,nx+1 ! Write velocities v, v0
    do j=0,ny+2
      write(1002,'(I4,X,I4,X,2D)') i, j, v(i,j), v0(i,j)
    end do
  end do
  do i=0,nx+1 ! Write pressure p
    do j=0,ny+1
      write(1002,'(I4,X,I4,X,D)') i, j, p(i,j)
    end do
  end do

  ! Write Lagrangian mesh
  do i=1,nl
    write(1002,'(I4,X,4D)') i, xl(i,1), xl(i,2), ul(i), vl(i)
  end do

  ! Close file
  close(1002)

end subroutine save_dat

!> subroutine to write last field
subroutine load_dat(u,v,p,u0,v0,xl,ul,vl,nx,ny,nl,hx,hy,hl,Re)
implicit none

real(DP), allocatable :: u(:,:), v(:,:), p(:,:), u0(:,:), v0(:,:)               ! Arrayes of Eulerian variables
real(DP), allocatable :: xl(:,:), ul(:), vl(:)                                  ! Arrayes of Lagrangian variables
integer :: nx, ny, nl                                                           ! Numbers of cells in x and y directions and Lagrangian nodes
real(DP):: hx, hy, hl                                                           ! Space steps
real(DP) :: Re                                                                  ! Reynolds number

integer :: i, j, i0, j0                                                         ! Loop indices

  ! Open file
  open(1003, file='start.dat')
  read(1003,*) Re
  read(1003,*) nx
  read(1003,*) ny
  read(1003,*) nl
  read(1003,*) hx
  read(1003,*) hy
  read(1003,*) hl
  read(1003,*)

  ! Read Eulerian mesh
  do i=0,nx+2 ! Read velocities u, u0
    do j=0,ny+1
      read(1003,*) i0, j0, u(i,j), u0(i,j)
    end do
  end do
  do i=0,nx+1 ! Read velocities v, v0
    do j=0,ny+2
      read(1003,*) i0, j0, v(i,j), v0(i,j)
    end do
  end do
  do i=0,nx+1 ! Read pressure p
    do j=0,ny+1
      read(1003,*) i0, j0, p(i,j)
    end do
  end do

  ! Read Lagrangian mesh
  do i=1,nl
    read(1003,*) i0, xl(i,1), xl(i,2), ul(i), vl(i)
  end do

  ! Close file
  close(1003)

end subroutine load_dat

end module inout
