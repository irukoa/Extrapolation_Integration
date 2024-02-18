module Extrapolation_Integration
  use Extrapolation_Integration_Kinds, only: sp, dp
  implicit none
  private

  real(dp), parameter :: tol = 1.0E3_dp

  interface extrapolation
    module procedure :: extrapol_r
    module procedure :: extrapol_c
    module procedure :: extrapol_rdp
    module procedure :: extrapol_cdp
  end interface extrapolation

  public :: extrapolation

contains

  function extrapol_r(array)

    !Utility to integrate a discetization of a real(sp)
    !f(x) (R->R) function over the [0, 1] real interval using
    !the extrapolation method.
    !The data is stored in the input array
    !such that array(i) = f(x_i), x_i = (i-1)/(2^n), i \in [1, 2^n + 1].
    !size(array) must be expressible as 2^n + 1, else the trapecium
    !method approximation will be returned.

    real(sp), intent(in) :: array(:)

    real(sp) :: extrapol_r

    integer :: n, sz
    real(sp) :: ord_array(size(array))

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
      extrapol_r = trap_r(array)
    case (-1)
      !Exception size(array) = 1, do not integrate.
      extrapol_r = array(1)
    case (0)
      !Special case size(array) = 2, use simple trapezium method.
      extrapol_r = sum(array)*0.5_sp
    case default
      !size(array) = 3, 5, 9, ... proceed normally.
      ord_array = org_dt_r(array, n)
      extrapol_r = extrapol_driver_r(ord_array, n)
    end select

  contains

    function trap_r(array)

      !Implementation of the composite trapezium method.
      !Assigns a weight of 0.5 to the endpoints and a weight of 1
      !to internal points. Later divides by the number of steps.
      !Data from array is assumed to be in the [0, 1] interval equally
      !spaced.

      real(sp), intent(in) :: array(:)

      real(sp) :: trap_r

      integer :: i

      trap_r = 0.5_sp*(array(1) + array(2))
      do i = 3, size(array)
        trap_r = trap_r + array(i)
      enddo
      trap_r = trap_r/real(size(array) - 1, sp)

    end function trap_r

    function org_dt_r(array, n)

      !Implementation of the practical ordering map.
      !Passes data array form being canonically ordered to practically ordered.
      !The integer n is the "order" of the data. Relates to the quantity of data points
      !by: size(array) = 1 + 2^n.

      real(sp), intent(in) :: array(:)
      integer, intent(in)  :: n

      real(sp) :: org_dt_r(size(array))

      integer :: i, j, k, l

      !Manual assignation to endpoints.
      org_dt_r(1) = array(1)
      org_dt_r(2) = array(size(array))

      !Intermediate points.
      l = 2
      do i = 2, 1 + n
        j = 2**(n - i + 1) + 1
        do k = 1, 2**(i - 2)
          l = l + 1
          org_dt_r(l) = array(j + 2*(j - 1)*(k - 1))
        enddo
      enddo

    end function org_dt_r

    function extrapol_driver_r(array, n)

      !Implementation of the extrapolation scheme.
      !Calculates trapezium rule approximation for the data array and different
      !values of the step. Practical order is assumed in the data array.
      !The integer n is the "order" of the data. Relates to the quantity of data points
      !by: size(array) = 1 + 2^n.
      !The different I_i(1/2^j) are then calculated using extrapolation.
      !The output, extrapol_driver_r, is the last extrapolated result.

      real(sp), intent(in) :: array(:)
      integer, intent(in)  :: n

      real(sp) :: extrapol_driver_r

      real(sp) :: workarray(n)
      integer :: i, j

      do i = 1, n !Trapezium rule. Implemented using practical order.
        workarray(i) = trap_r(array(1:2**(i) + 1))
      enddo

      do i = 2, n !Extrapolation scheme. Calculates I_i(1/2^j).
        do j = 1, n - i + 1
          workarray(j) = ((4.0_sp)**(i - 1))*workarray(j + 1) - workarray(j)
          workarray(j) = workarray(j)/(-1.0_sp + (4.0_sp)**(i - 1))
        enddo
      enddo

      extrapol_driver_r = workarray(1)

    end function extrapol_driver_r

  end function extrapol_r

  function extrapol_c(array)

    !Utility to integrate a discetization of a complex(sp)
    !f(x) (C->C) function over the [0, 1] real interval using
    !the extrapolation method.
    !The data is stored in the input array
    !such that array(i) = f(x_i), x_i = (i-1)/(2^n), i \in [1, 2^n + 1].
    !size(array) must be expressible as 2^n + 1, else the trapecium
    !method approximation will be returned.

    complex(sp), intent(in) :: array(:)

    complex(sp) :: extrapol_c

    integer :: n, sz
    complex(sp) :: ord_array(size(array))

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
      extrapol_c = trap_c(array)
    case (-1)
      !Exception size(array) = 1, do not integrate.
      extrapol_c = array(1)
    case (0)
      !Special case size(array) = 2, use simple trapezium method.
      extrapol_c = sum(array)*0.5_sp
    case default
      !size(array) = 3, 5, 9, ... proceed normally.
      ord_array = org_dt_c(array, n)
      extrapol_c = extrapol_driver_c(ord_array, n)
    end select

  contains

    function trap_c(array)

      !Implementation of the composite trapezium method.
      !Assigns a weight of 0.5 to the endpoints and a weight of 1
      !to internal points. Later divides by the number of steps.
      !Data from array is assumed to be in the [0, 1] interval equally
      !spaced.

      complex(sp), intent(in) :: array(:)

      complex(sp) :: trap_c

      integer :: i

      trap_c = 0.5_sp*(array(1) + array(2))
      do i = 3, size(array)
        trap_c = trap_c + array(i)
      enddo
      trap_c = trap_c/real(size(array) - 1, sp)

    end function trap_c

    function org_dt_c(array, n)

      !Implementation of the practical ordering map.
      !Passes data array form being canonically ordered to practically ordered.
      !The integer n is the "order" of the data. Relates to the quantity of data points
      !by: size(array) = 1 + 2^n.

      complex(sp), intent(in) :: array(:)
      integer, intent(in)  :: n

      complex(sp) :: org_dt_c(size(array))

      integer :: i, j, k, l

      !Manual assignation to endpoints.
      org_dt_c(1) = array(1)
      org_dt_c(2) = array(size(array))

      !Intermediate points.
      l = 2
      do i = 2, 1 + n
        j = 2**(n - i + 1) + 1
        do k = 1, 2**(i - 2)
          l = l + 1
          org_dt_c(l) = array(j + 2*(j - 1)*(k - 1))
        enddo
      enddo

    end function org_dt_c

    function extrapol_driver_c(array, n)

      !Implementation of the extrapolation scheme.
      !Calculates trapezium rule approximation for the data array and different
      !values of the step. Practical order is assumed in the data array.
      !The integer n is the "order" of the data. Relates to the quantity of data points
      !by: size(array) = 1 + 2^n.
      !The different I_i(1/2^j) are then calculated using extrapolation.
      !The output, extrapol_driver_c, is the last extrapolated result.

      complex(sp), intent(in) :: array(:)
      integer, intent(in)  :: n

      complex(sp) :: extrapol_driver_c

      complex(sp) :: workarray(n)
      integer :: i, j

      do i = 1, n !Trapezium rule. Implemented using practical order.
        workarray(i) = trap_c(array(1:2**(i) + 1))
      enddo

      do i = 2, n !Extrapolation scheme. Calculates I_i(1/2^j).
        do j = 1, n - i + 1
          workarray(j) = ((4.0_sp)**(i - 1))*workarray(j + 1) - workarray(j)
          workarray(j) = workarray(j)/(-1.0_sp + (4.0_sp)**(i - 1))
        enddo
      enddo

      extrapol_driver_c = workarray(1)

    end function extrapol_driver_c

  end function extrapol_c

  function extrapol_rdp(array)

    !Utility to integrate a discetization of a real(dp)
    !f(x) (R->R) function over the [0, 1] real interval using
    !the extrapolation method.
    !The data is stored in the input array
    !such that array(i) = f(x_i), x_i = (i-1)/(2^n), i \in [1, 2^n + 1].
    !size(array) must be expressible as 2^n + 1, else the trapecium
    !method approximation will be returned.

    real(dp), intent(in) :: array(:)

    real(dp) :: extrapol_rdp

    integer :: n, sz
    real(dp) :: ord_array(size(array))

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
      extrapol_rdp = trap_rdp(array)
    case (-1)
      !Exception size(array) = 1, do not integrate.
      extrapol_rdp = array(1)
    case (0)
      !Special case size(array) = 2, use simple trapezium method.
      extrapol_rdp = sum(array)*0.5_dp
    case default
      !size(array) = 3, 5, 9, ... proceed normally.
      ord_array = org_dt_rdp(array, n)
      extrapol_rdp = extrapol_driver_rdp(ord_array, n)
    end select

  contains

    function trap_rdp(array)

      !Implementation of the composite trapezium method.
      !Assigns a weight of 0.5 to the endpoints and a weight of 1
      !to internal points. Later divides by the number of steps.
      !Data from array is assumed to be in the [0, 1] interval equally
      !spaced.

      real(dp), intent(in) :: array(:)

      real(dp) :: trap_rdp

      integer :: i

      trap_rdp = 0.5_dp*(array(1) + array(2))
      do i = 3, size(array)
        trap_rdp = trap_rdp + array(i)
      enddo
      trap_rdp = trap_rdp/real(size(array) - 1, dp)

    end function trap_rdp

    function org_dt_rdp(array, n)

      !Implementation of the practical ordering map.
      !Passes data array form being canonically ordered to practically ordered.
      !The integer n is the "order" of the data. Relates to the quantity of data points
      !by: size(array) = 1 + 2^n.

      real(dp), intent(in) :: array(:)
      integer, intent(in)  :: n

      real(dp) :: org_dt_rdp(size(array))

      integer :: i, j, k, l

      !Manual assignation to endpoints.
      org_dt_rdp(1) = array(1)
      org_dt_rdp(2) = array(size(array))

      !Intermediate points.
      l = 2
      do i = 2, 1 + n
        j = 2**(n - i + 1) + 1
        do k = 1, 2**(i - 2)
          l = l + 1
          org_dt_rdp(l) = array(j + 2*(j - 1)*(k - 1))
        enddo
      enddo

    end function org_dt_rdp

    function extrapol_driver_rdp(array, n)

      !Implementation of the extrapolation scheme.
      !Calculates trapezium rule approximation for the data array and different
      !values of the step. Practical order is assumed in the data array.
      !The integer n is the "order" of the data. Relates to the quantity of data points
      !by: size(array) = 1 + 2^n.
      !The different I_i(1/2^j) are then calculated using extrapolation.
      !The output, extrapol_driver_rdp, is the last extrapolated result.

      real(dp), intent(in) :: array(:)
      integer, intent(in)  :: n

      real(dp) :: extrapol_driver_rdp

      real(dp) :: workarray(n)
      integer :: i, j

      do i = 1, n !Trapezium rule. Implemented using practical order.
        workarray(i) = trap_rdp(array(1:2**(i) + 1))
      enddo

      do i = 2, n !Extrapolation scheme. Calculates I_i(1/2^j).
        do j = 1, n - i + 1
          workarray(j) = ((4.0_dp)**(i - 1))*workarray(j + 1) - workarray(j)
          workarray(j) = workarray(j)/(-1.0_dp + (4.0_dp)**(i - 1))
        enddo
      enddo

      extrapol_driver_rdp = workarray(1)

    end function extrapol_driver_rdp

  end function extrapol_rdp

  function extrapol_cdp(array)

    !Utility to integrate a discetization of a complex(dp)
    !f(x) (C->C) function over the [0, 1] real interval using
    !the extrapolation method.
    !The data is stored in the input array
    !such that array(i) = f(x_i), x_i = (i-1)/(2^n), i \in [1, 2^n + 1].
    !size(array) must be expressible as 2^n + 1, else the trapecium
    !method approximation will be returned.

    complex(dp), intent(in) :: array(:)

    complex(dp) :: extrapol_cdp

    integer :: n, sz
    complex(dp) :: ord_array(size(array))

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
      extrapol_cdp = trap_cdp(array)
    case (-1)
      !Exception size(array) = 1, do not integrate.
      extrapol_cdp = array(1)
    case (0)
      !Special case size(array) = 2, use simple trapezium method.
      extrapol_cdp = sum(array)*0.5_dp
    case default
      !size(array) = 3, 5, 9, ... proceed normally.
      ord_array = org_dt_cdp(array, n)
      extrapol_cdp = extrapol_driver_cdp(ord_array, n)
    end select

  contains

    function trap_cdp(array)

      !Implementation of the composite trapezium method.
      !Assigns a weight of 0.5 to the endpoints and a weight of 1
      !to internal points. Later divides by the number of steps.
      !Data from array is assumed to be in the [0, 1] interval equally
      !spaced.

      complex(dp), intent(in) :: array(:)

      complex(dp) :: trap_cdp

      integer :: i

      trap_cdp = 0.5_dp*(array(1) + array(2))
      do i = 3, size(array)
        trap_cdp = trap_cdp + array(i)
      enddo
      trap_cdp = trap_cdp/real(size(array) - 1, dp)

    end function trap_cdp

    function org_dt_cdp(array, n)

      !Implementation of the practical ordering map.
      !Passes data array form being canonically ordered to practically ordered.
      !The integer n is the "order" of the data. Relates to the quantity of data points
      !by: size(array) = 1 + 2^n.

      complex(dp), intent(in) :: array(:)
      integer, intent(in)  :: n

      complex(dp) :: org_dt_cdp(size(array))

      integer :: i, j, k, l

      !Manual assignation to endpoints.
      org_dt_cdp(1) = array(1)
      org_dt_cdp(2) = array(size(array))

      !Intermediate points.
      l = 2
      do i = 2, 1 + n
        j = 2**(n - i + 1) + 1
        do k = 1, 2**(i - 2)
          l = l + 1
          org_dt_cdp(l) = array(j + 2*(j - 1)*(k - 1))
        enddo
      enddo

    end function org_dt_cdp

    function extrapol_driver_cdp(array, n)

      !Implementation of the extrapolation scheme.
      !Calculates trapezium rule approximation for the data array and different
      !values of the step. Practical order is assumed in the data array.
      !The integer n is the "order" of the data. Relates to the quantity of data points
      !by: size(array) = 1 + 2^n.
      !The different I_i(1/2^j) are then calculated using extrapolation.
      !The output, extrapol_driver_cdp, is the last extrapolated result.

      complex(dp), intent(in) :: array(:)
      integer, intent(in)  :: n

      complex(dp) :: extrapol_driver_cdp

      complex(dp) :: workarray(n)
      integer :: i, j

      do i = 1, n !Trapezium rule. Implemented using practical order.
        workarray(i) = trap_cdp(array(1:2**(i) + 1))
      enddo

      do i = 2, n !Extrapolation scheme. Calculates I_i(1/2^j).
        do j = 1, n - i + 1
          workarray(j) = ((4.0_dp)**(i - 1))*workarray(j + 1) - workarray(j)
          workarray(j) = workarray(j)/(-1.0_dp + (4.0_dp)**(i - 1))
        enddo
      enddo

      extrapol_driver_cdp = workarray(1)

    end function extrapol_driver_cdp

  end function extrapol_cdp

end module Extrapolation_Integration
