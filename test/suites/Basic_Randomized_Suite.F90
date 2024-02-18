module Basic_Randomized_Suite
  use, intrinsic :: iso_fortran_env, only: output_unit
  use Extrapolation_Integration_Kinds, only: sp, dp
  use Extrapolation_Integration, only: extrapolation
  use testdrive, only: new_unittest, unittest_type, error_type
  implicit none
  private

  real(sp), parameter :: tolsp = 1.0E1_sp
  real(dp), parameter :: toldp = 1.0E5_dp
  real(dp), parameter :: pi = acos(-1.0_dp)

  public :: Collect_Basic_Randomized_Tests

contains

  subroutine Collect_Basic_Randomized_Tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [new_unittest("Basic Test", test_basic), &
                 new_unittest("Randomized Test", test_randomized)]

  end subroutine Collect_Basic_Randomized_Tests

  subroutine test_basic(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: n = 129

    real(dp) :: x, res

    real(sp)    :: rsr, rsarrayvt(n + 1), rsarrayvn(1), rsarrayvtb(2), rsarray(n)
    complex(sp) :: csr, csarrayvt(n + 1), csarrayvn(1), csarrayvtb(2), csarray(n)
    real(dp)    :: rdr, rdarrayvt(n + 1), rdarrayvn(1), rdarrayvtb(2), rdarray(n)
    complex(dp) :: cdr, cdarrayvt(n + 1), cdarrayvn(1), cdarrayvtb(2), cdarray(n)

    integer :: i

    do i = 1, n + 1
      x = 1.0_dp*(real(i - 1, dp))/(real(n - 1, dp))
      rdarrayvt(i) = f(x, 1.0_dp)
      rsarrayvt(i) = real(rdarrayvt(i), sp)
      csarrayvt(i) = cmplx(rdarrayvt(i), 0.5_dp*rdarrayvt(i), sp)
      cdarrayvt(i) = cmplx(rdarrayvt(i), 0.5_dp*rdarrayvt(i), dp)
    enddo
    rsr = extrapolation(rsarrayvt)
    csr = extrapolation(csarrayvt)
    rdr = extrapolation(rdarrayvt)
    cdr = extrapolation(cdarrayvt)

    if (abs(rsr - 0.154382393_sp) > tolsp*epsilon(1.0_sp)) then
      allocate (error)
      return
    endif

    if (abs(abs(csr) - 0.172604769_sp) > tolsp*epsilon(1.0_sp)) then
      allocate (error)
      return
    endif

    if (abs(rdr - 0.15438247092706356_dp) > toldp*epsilon(1.0_dp)) then
      allocate (error)
      return
    endif

    if (abs(abs(cdr) - 0.17260484976364954_dp) > toldp*epsilon(1.0_dp)) then
      allocate (error)
      return
    endif

    rdarrayvn = f(0.5_dp, 1.0_dp)
    rsarrayvn = real(rdarrayvn, sp)
    csarrayvn = cmplx(rdarrayvn, 0.5_dp*rdarrayvn, sp)
    cdarrayvn = cmplx(rdarrayvn, 0.5_dp*rdarrayvn, dp)

    rsr = extrapolation(rsarrayvn)
    csr = extrapolation(csarrayvn)
    rdr = extrapolation(rdarrayvn)
    cdr = extrapolation(cdarrayvn)

    res = f(0.5_dp, 1.0_dp)

    if (abs(rsr - real(res)) > tolsp*epsilon(1.0_sp)) then
      allocate (error)
      return
    endif

    if (abs(abs(csr) - abs(real(res)*sqrt(5.0_sp)/2.0_sp)) > tolsp*epsilon(1.0_sp)) then
      allocate (error)
      return
    endif

    if (abs(rdr - res) > toldp*epsilon(1.0_dp)) then
      allocate (error)
      return
    endif

    if (abs(abs(cdr) - abs(res*sqrt(5.0_dp)/2.0_dp)) > toldp*epsilon(1.0_dp)) then
      allocate (error)
      return
    endif

    do i = 1, 2
      x = 1.0_dp*(real(i - 1, dp))/(real(2 - 1, dp))
      rdarrayvtb(i) = f(x, 1.0_dp)
      rsarrayvtb(i) = real(rdarrayvtb(i), sp)
      csarrayvtb(i) = cmplx(rdarrayvtb(i), 0.5_dp*rdarrayvtb(i), sp)
      cdarrayvtb(i) = cmplx(rdarrayvtb(i), 0.5_dp*rdarrayvtb(i), dp)
    enddo
    rsr = extrapolation(rsarrayvtb)
    csr = extrapolation(csarrayvtb)
    rdr = extrapolation(rdarrayvtb)
    cdr = extrapolation(cdarrayvtb)

    res = 0.5_dp*(f(0.0_dp, 1.0_dp) + f(1.0_dp, 1.0_dp))

    if (abs(rsr - real(res)) > tolsp*epsilon(1.0_sp)) then
      allocate (error)
      return
    endif

    if (abs(abs(csr) - abs(real(res)*sqrt(5.0_sp)/2.0_sp)) > tolsp*epsilon(1.0_sp)) then
      allocate (error)
      return
    endif

    if (abs(rdr - res) > toldp*epsilon(1.0_dp)) then
      allocate (error)
      return
    endif

    if (abs(abs(cdr) - abs(res*sqrt(5.0_dp)/2.0_dp)) > toldp*epsilon(1.0_dp)) then
      allocate (error)
      return
    endif

    do i = 1, n
      x = 1.0_dp*(real(i - 1, dp))/(real(n - 1, dp))
      rdarray(i) = f(x, 1.0_dp)
      rsarray(i) = real(rdarray(i), sp)
      csarray(i) = cmplx(rdarray(i), 0.5_dp*rdarray(i), sp)
      cdarray(i) = cmplx(rdarray(i), 0.5_dp*rdarray(i), dp)
    enddo
    rsr = extrapolation(rsarray)
    csr = extrapolation(csarray)
    rdr = extrapolation(rdarray)
    cdr = extrapolation(cdarray)

    res = (-Si(2*pi, 120) + Si(2*(exp(1.0_dp))*pi, 120))/1.0_dp

    if (abs(rsr - real(res)) > tolsp*epsilon(1.0_sp)) then
      allocate (error)
      return
    endif

    if (abs(abs(csr) - abs(real(res)*sqrt(5.0_sp)/2.0_sp)) > tolsp*epsilon(1.0_sp)) then
      allocate (error)
      return
    endif

    if (abs(rdr - res) > toldp*epsilon(1.0_dp)) then
      allocate (error)
      return
    endif

    if (abs(abs(cdr) - abs(res*sqrt(5.0_dp)/2.0_dp)) > toldp*epsilon(1.0_dp)) then
      allocate (error)
      return
    endif

  end subroutine test_basic

  subroutine test_randomized(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: omn = 1, omx = 11

    real(dp) :: x, res, or, ar, adpttol

    real(dp) :: rdr
    real(dp), allocatable :: rdarray(:)

    integer :: ord, n, i

    call random_seed()
    call random_number(or)

    ord = nint(real(omn, dp) + real(omx - omn, dp)*or)
    !ord = 4!DBG.
    n = 2**(ord) + 1

    write (output_unit, fmt="(A, i0, A)") &
      & "Selected number of samples: ", n, "."

    call random_number(ar)
    !ar = 1.0_dp !DBG.

    write (output_unit, fmt="(A, E14.8, A)") &
    & "Value of a: ", ar, "."

    adpttol = 5.0E10_dp*epsilon(1.0_dp) + 1.0E12_dp*(10.0_dp**(16.0_dp - real(4*ord, dp)))*epsilon(1.0_dp)
    if (adpttol > 1.0_dp) adpttol = 0.1_dp

    write (output_unit, fmt="(A, E14.8, A)") &
    & "Adaptive tolerance: ", adpttol, "."

    allocate (rdarray(n))

    do i = 1, n
      x = 1.0_dp*(real(i - 1, dp))/(real(n - 1, dp))
      rdarray(i) = f(x, ar)
    enddo
    rdr = extrapolation(rdarray)

    res = (-Si(2*pi, 120) + Si(2*(exp(ar))*pi, 120))/ar

    write (output_unit, fmt="(A, E14.8, A)") &
    & "True result: ", res, "."

    write (output_unit, fmt="(A, E14.8, A)") &
    & "Extrapolated result: ", rdr, "."

    write (output_unit, fmt="(A, E14.8, A)") &
    & "Relative error: ", 100*abs(res - rdr)/abs(res), "%."

    if (abs(res - rdr) > adpttol) then
      allocate (error)
      return
    endif

    deallocate (rdarray)

  end subroutine test_randomized

  function f(x, a)
    !Function to test.
    !Its integral is highly nonlinear depending on the
    !value of a.
    !The integral from 0 to 1 of f is given by
    !(-Si(2*pi) + Si(2*(e^a)*pi))/a.
    real(dp), intent(in) :: x, a
    real(dp) :: f
    f = sin(2*pi*exp(a*x))
  end function f

  function Si(x, n)
    !Sine integral.
    !Given in terms of its series expansion.
    real(dp), intent(in) :: x
    integer, intent(in) :: n

    real(dp) :: Si

    integer :: i, j
    real(dp) :: fct

    Si = 0.0_dp
    do i = 1, n
      fct = 1.0_dp
      do j = 1, 2*i - 1
        fct = fct*real(j, dp)
      enddo
      Si = Si + ((-1.0_dp)**(i - 1))*((x)**(2*i - 1))/((2*i - 1)*fct)
    enddo
  end function Si
end module Basic_Randomized_Suite
