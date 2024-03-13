program main
  use iso_fortran_env, only: output_unit
  use Extrapolation_Integration_Kinds
  use Extrapolation_Integration

  implicit none

  character(len=10) :: ver = "1.0.0"

  write (unit=output_unit, fmt="(A)") ""
  write (unit=output_unit, fmt="(A)") "     ______     __                    __      __  _                "
  write (unit=output_unit, fmt="(A)") "    / ____/  __/ /_____ _____  ____  / /___ _/ /_(_)___  ____      "
  write (unit=output_unit, fmt="(A)") "   / __/ | |/_/ __/ __ `/ __ \/ __ \/ / __ `/ __/ / __ \/ __ \     "
  write (unit=output_unit, fmt="(A)") "  / /____>  </ /_/ /_/ / /_/ / /_/ / / /_/ / /_/ / /_/ / / / /     "
  write (unit=output_unit, fmt="(A)") " /_____/_/|_|\__/\__,_/ .___/\____/_/\__,_/\__/_/\____/_/ /_/      "
  write (unit=output_unit, fmt="(A)") "    /  _/___  / /____/_/___ __________ _/ /_(_)___  ____           "
  write (unit=output_unit, fmt="(A)") "    / // __ \/ __/ _ \/ __ `/ ___/ __ `/ __/ / __ \/ __ \          "
  write (unit=output_unit, fmt="(A)") "  _/ // / / / /_/  __/ /_/ / /  / /_/ / /_/ / /_/ / / / /          "
  write (unit=output_unit, fmt="(A)") " /___/_/ /_/\__/\___/\__, /_/   \__,_/\__/_/\____/_/ /_/           "
  write (unit=output_unit, fmt="(A)") "                    /____/                                         "
  write (unit=output_unit, fmt="(A)") ""
  write (unit=output_unit, fmt="(A)") " Extrapolation Integration v"//trim(adjustl(ver))//" built."

end program main
