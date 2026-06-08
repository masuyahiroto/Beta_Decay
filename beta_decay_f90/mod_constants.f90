module mod_constants
  implicit none

  real(8), parameter :: ALPHA    = 1.0d0 / 137.035999084d0
  real(8), parameter :: ME_MEV   = 0.51099895d0
  real(8), parameter :: G_A      = 1.2756d0
  real(8), parameter :: G_V      = 1.0d0
  real(8), parameter :: D_CONST  = 6289.0d0   ! [s]  t_{1/2}*f = D/(g_A^2 * B_GT)
  real(8), parameter :: PI       = 3.141592653589793d0

  integer, parameter :: MAX_SP   = 100         ! max single-particle states
  integer, parameter :: MAX_OBDM = 10000       ! max OBDM entries

end module mod_constants
