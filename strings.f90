!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Module to implement some useful string functions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program: ibm2d_cyl -- SIMPLE + IBM method for incomressible viscous Navier -
! Stockes equations in 2D statement for flow past a cylinder
! Copyright (C): Denis V. Esipov (2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module strings
implicit none

  public int2str, double2str, str2int, str2double

contains

  !> function to write i integer to string
  !> example: i = 1234567890, int2str = '1234567890      '
  function int2str(i,fm)
  implicit none

  ! Input-Output variables
  character(16) :: int2str                                                      ! Function result
  integer :: i                                                                  ! Integer
  character(*), optional :: fm                                                  ! Format string
  ! Local variables
  integer :: ios                                                                ! Input-output status

    ! Convert number to string
    int2str = ''
    if (present(fm)) then
      write(int2str,fm,iostat=ios) i
    else
      write(int2str,'(I)',iostat=ios) i
      int2str = adjustl(int2str)
    end if

    ! Check errors
    if (ios /= 0) int2str = 'int2str error!'

  end function int2str

  !> function to write d double to string
  !> example d = -1234567890123456.0d84, double2str(1:24) = '-0.1234567890123456E+100'
  function double2str(d,fm)
  implicit none

  ! Input-Output variables
  character(64) :: double2str                                                   ! Function result
  double precision :: d                                                         ! Double
  character(*), optional :: fm                                                  ! Format string
  ! Local variables
  integer :: ios                                                                ! Input-output status

    ! Convert number to string
    double2str = ''
    if (present(fm)) then
      write(double2str,fm,iostat=ios) d
    else
      write(double2str,'(E)',iostat=ios) d
      double2str = adjustl(double2str)
    end if

    ! Check errors
    if (ios /= 0) double2str = 'double2str error!'

  end function double2str

  !> function to read integer from string str
  function str2int(str)
  implicit none

  ! Input-Output variables
  integer :: str2int                                                            ! Function result
  character(*) :: str                                                           ! String
  ! Local variables
  integer :: ios                                                                ! Input-output status

    read(str,*,iostat=ios) str2int
    if (ios /= 0) str2int = 0 ! If error, put zero into function result

  end function str2int

  !> function to read double from string str
  function str2double(str)
  implicit none

  ! Input-Output variables
  double precision :: str2double                                                ! Function result
  character(*) :: str                                                           ! String
  ! Local variables
  integer :: ios                                                                ! Input-output status

    read(str,*,iostat=ios) str2double
    if (ios /= 0) str2double = 0.0_8 ! If error, put zero into function result

  end function str2double

end module strings
