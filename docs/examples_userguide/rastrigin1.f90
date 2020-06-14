program rastrigin1

  implicit none

  integer, parameter :: dp = kind(1.0d0)

  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: a  = 10.0_dp
  real(dp), parameter :: b  = 2.0_dp * pi

  character(len=*), parameter :: pfile = 'params.txt'
  character(len=*), parameter :: ofile = 'out.txt'

  integer, parameter :: punit = 99
  integer, parameter :: ounit = 101

  real(dp), dimension(100) :: x ! parameters, up to 100 dimensions
  real(dp) :: out               ! output value
  integer  :: n                 ! number of dimensions

  integer  :: ios

  ! read parameters
  open(punit, file=pfile, status='old', action='read')
  ios = 0
  n = 1
  do while (ios==0)
     read(punit, fmt=*, iostat=ios) x(n)
     n = n + 1
  end do
  n = n - 2
  close(punit)

  ! calc function
  out = a * real(n,dp) + sum(x(1:n)**2 - a*cos(b*x(1:n)))

  ! write output file
  open(ounit, file=ofile)
  write(ounit,*) out
  close(ounit)

end program rastrigin1
