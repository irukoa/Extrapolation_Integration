#:set types = ['real(sp)', 'complex(sp)', 'real(dp)', 'complex(dp)']
#:set kinds = ['sp', 'sp', 'dp', 'dp']
#:set number_sets = ['R', 'C', 'R', 'C']
#:set suffixes = ['r', 'c', 'rdp', 'cdp']
#:set lst = list(zip(types, kinds, suffixes, number_sets))
module Extrapolation_Integration
  use Extrapolation_Integration_Kinds, only: sp, dp
  implicit none
  private

  real(dp), parameter :: tol = 1.0E3_dp

  interface extrapolation
  #: for type, kind, sfx, st in lst
    module procedure :: extrapol_${sfx}$
  #: endfor
  end interface extrapolation

  public :: extrapolation

contains

  #: for type, kind, sfx, st in lst
  function extrapol_${sfx}$(array)

    !Utility to integrate a discetization of a ${type}$
    !f(x) (${st}$->${st}$) function over the [0, 1] real interval using
    !the extrapolation method.
    !The data is stored in the input array
    !such that array(i) = f(x_i), x_i = (i-1)/(2^n), i \in [1, 2^n + 1].
    !size(array) must be expressible as 2^n + 1, else the trapecium
    !method approximation will be returned.

    ${type}$, intent(in) :: array(:)

    ${type}$ :: extrapol_${sfx}$

    integer :: n, sz
    ${type}$ :: ord_array(size(array))

    real(dp) :: nr

    sz = size(array)

    if (sz /= 1) then
      nr = log(real(sz - 1, dp))/log(2.0_dp)
      n = nint(nr)
      if ((abs(real(n, dp) - nr) > tol*epsilon(1.0_dp))) n = -2
    else
      n = -1
    endif

    select case (n)
    case (-2)
      !Extrapolation not applicable. Return trapecium approx.
      extrapol_${sfx}$ = trap_${sfx}$(array)
    case (-1)
      !Exception size(array) = 1, do not integrate.
      extrapol_${sfx}$ = array(1)
    case (0)
      !Special case size(array) = 2, use simple trapezium method.
      extrapol_${sfx}$ = sum(array)*0.5_${kind}$
    case default
      !size(array) = 3, 5, 9, ... proceed normally.
      ord_array = org_dt_${sfx}$(array, n)
      extrapol_${sfx}$ = extrapol_driver_${sfx}$(ord_array, n)
    end select

  contains

    function trap_${sfx}$(array)

      !Implementation of the composite trapezium method.
      !Assigns a weight of 0.5 to the endpoints and a weight of 1
      !to internal points. Later divides by the number of steps.
      !Data from array is assumed to be in the [0, 1] interval equally
      !spaced.

      ${type}$, intent(in) :: array(:)

      ${type}$ :: trap_${sfx}$

      integer :: i

      trap_${sfx}$ = 0.5_${kind}$*(array(1) + array(2))
      do i = 3, size(array)
        trap_${sfx}$ = trap_${sfx}$ + array(i)
      enddo
      trap_${sfx}$ = trap_${sfx}$/real(size(array) - 1, ${kind}$)

    end function trap_${sfx}$

    function org_dt_${sfx}$(array, n)

      !Implementation of the practical ordering map.
      !Passes data array form being canonically ordered to practically ordered.
      !The integer n is the "order" of the data. Relates to the quantity of data points
      !by: size(array) = 1 + 2^n.

      ${type}$, intent(in) :: array(:)
      integer, intent(in)  :: n

      ${type}$ :: org_dt_${sfx}$(size(array))

      integer :: i, j, k, l

      !Manual assignation to endpoints.
      org_dt_${sfx}$(1) = array(1)
      org_dt_${sfx}$(2) = array(size(array))

      !Intermediate points.
      l = 2
      do i = 2, 1 + n
        j = 2**(n - i + 1) + 1
        do k = 1, 2**(i - 2)
          l = l + 1
          org_dt_${sfx}$(l) = array(j + 2*(j - 1)*(k - 1))
        enddo
      enddo

    end function org_dt_${sfx}$

    function extrapol_driver_${sfx}$(array, n)

      !Implementation of the extrapolation scheme.
      !Calculates trapezium rule approximation for the data array and different
      !values of the step. Practical order is assumed in the data array.
      !The integer n is the "order" of the data. Relates to the quantity of data points
      !by: size(array) = 1 + 2^n.
      !The different I_i(1/2^j) are then calculated using extrapolation.
      !The output, extrapol_driver_${sfx}$, is the last extrapolated result.

      ${type}$, intent(in) :: array(:)
      integer, intent(in)  :: n

      ${type}$ :: extrapol_driver_${sfx}$

      ${type}$ :: workarray(n)
      integer :: i, j

      do i = 1, n !Trapezium rule. Implemented using practical order.
        workarray(i) = trap_${sfx}$(array(1:2**(i) + 1))
      enddo

      do i = 2, n !Extrapolation scheme. Calculates I_i(1/2^j).
        do j = 1, n - i + 1
          workarray(j) = ((4.0_${kind}$)**(i - 1))*workarray(j + 1) - workarray(j)
          workarray(j) = workarray(j)/(-1.0_${kind}$ + (4.0_${kind}$)**(i - 1))
        enddo
      enddo

      extrapol_driver_${sfx}$ = workarray(1)

    end function extrapol_driver_${sfx}$

  end function extrapol_${sfx}$

  #: endfor
end module Extrapolation_Integration
