program Integral_of_SINC
  use, intrinsic :: iso_fortran_env, only: output_unit
  use Extrapolation_Integration_Kinds, only: wp => dp
  use Extrapolation_Integration, only: extrapolation

  !In this example, we consider the integration of the
  !sinc function and demonstrate the extrapolation method.

  implicit none

  !Number of terms to be kept in the series expansion
  !of the sine integral.
  integer, parameter :: nterms = 100

  !Number of samples. For the extrapolation method to be
  !applied, this number must be expressible as 2^m + 1
  !for some natural number m.
  integer, parameter :: n = 65
  !Integration bounds.
  real(wp), parameter :: a = 0.0_wp, &
                         b = 3.0_wp

  real(wp) :: storage(n), result, true_result

  real(wp) :: x
  integer :: i

  do i = 1, n
    x = a + (b - a)*real(i - 1, wp)/real(n - 1, wp)
    storage(i) = sinc(x)
  enddo
  !"extrapolation" is an array reduction operation,
  !that is, given an array returns a scalar.
  result = extrapolation(storage)
  !The method assumes an array sampled in equally
  !spaced points in the [0, 1] interval. Thats why
  !we multiply by (b-a).
  result = result*(b - a)

  !We also compute the exact value of the integral,
  true_result = Si(b, nterms) - Si(a, nterms)

  write (output_unit, fmt="(A, F20.16, A)") &
  & "True result: ", true_result, "."
  write (output_unit, fmt="(A)") ""

  write (output_unit, fmt="(A, F20.16, A)") &
    & "Extrapolated result: ", result, "."

  write (output_unit, fmt="(A, F20.16, A)") &
    & "  -Relative error: ", 100*abs(true_result - result)/abs(true_result), "%."
  write (output_unit, fmt="(A)") ""
  !For comparison, consider also the rectangle approximation:

  write (output_unit, fmt="(A, F20.16, A)") &
  & "Rectangle approximation result: ", (b - a)*sum(storage)/n, "."

  write (output_unit, fmt="(A, F20.16, A)") &
    & "  -Relative error: ", 100*abs(true_result - (b - a)*sum(storage)/n)/abs(true_result), "%."

  !The user is now encouraged to change the number of samples
  !such that it is not expressible as 2^m + 1.

contains

  function sinc(x)
    !sinc function implementation.
    real(wp), intent(in) :: x

    real(wp) :: sinc

    if (abs(x) < 1.0E-12_wp) then
      sinc = 1.0_wp
    else
      sinc = sin(x)/x
    endif

  end function sinc

  function Si(x, n)
    !Integral of the sinc function from 0 to x,
    !a.k.a. the sine integral,
    !given in terms of its series expansion.
    real(wp), intent(in) :: x
    integer, intent(in) :: n

    real(wp) :: Si

    integer :: i, j
    real(wp) :: fct

    Si = 0.0_wp
    do i = 1, n
      fct = 1.0_wp
      do j = 1, 2*i - 1
        fct = fct*real(j, wp)
      enddo
      Si = Si + ((-1.0_wp)**(i - 1))*((x)**(2*i - 1))/((2*i - 1)*fct)
    enddo
  end function Si

end program Integral_of_SINC
