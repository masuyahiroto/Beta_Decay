module mod_phase_space
  use mod_constants, only: ALPHA, ME_MEV, PI
  use mod_shape_factor, only: ShapeFactor, eval_shape_factor
  implicit none

  integer, parameter :: N_QUAD = 200  ! quadrature points

contains

  ! Fermi function F(Z, W) for beta-minus decay (point nucleus approximation)
  ! eta = alpha * Z * W / p,  F = 2*pi*eta / (1 - exp(-2*pi*eta))
  function fermi_function(Z, W) result(F)
    integer, intent(in) :: Z
    real(8), intent(in) :: W
    real(8) :: F, p, eta, two_pi_eta

    if (W <= 1.0d0) then
      F = 0.0d0
      return
    end if

    p           = sqrt(W*W - 1.0d0)
    eta         = ALPHA * dble(Z) * W / p
    two_pi_eta  = 2.0d0 * PI * eta

    if (two_pi_eta > 500.0d0) then
      F = two_pi_eta
    else
      F = two_pi_eta / (1.0d0 - exp(-two_pi_eta))
    end if
  end function fermi_function

  ! Phase space integral using composite Simpson's rule:
  ! f = integral_1^W0  C(W) * F(Z,W) * p * W * (W0-W)^2  dW
  function phase_space_integral(Z, Q_MeV, sf) result(integral)
    integer,          intent(in) :: Z
    real(8),          intent(in) :: Q_MeV
    type(ShapeFactor),intent(in) :: sf
    real(8) :: integral

    real(8) :: W0, dW, We, ff_val, mom, shape_val
    integer :: i
    real(8) :: val(0:N_QUAD)

    integral = 0.0d0
    W0 = Q_MeV / ME_MEV + 1.0d0

    if (W0 <= 1.0d0) return

    dW = (W0 - 1.0d0) / dble(N_QUAD)

    do i = 0, N_QUAD
      We = 1.0d0 + dble(i) * dW
      ! avoid exact endpoints where mom=0 or (W0-We)=0
      if (i == 0)      We = 1.0d0 + 1.0d-8
      if (i == N_QUAD) We = W0    - 1.0d-8

      ff_val    = fermi_function(Z, We)
      mom       = sqrt(We*We - 1.0d0)
      shape_val = eval_shape_factor(sf, We)
      val(i) = shape_val * ff_val * mom * We * (W0 - We)**2
    end do

    ! Composite Simpson's rule (N_QUAD must be even)
    integral = val(0) + val(N_QUAD)
    do i = 1, N_QUAD - 1
      if (mod(i, 2) == 1) then
        integral = integral + 4.0d0 * val(i)
      else
        integral = integral + 2.0d0 * val(i)
      end if
    end do
    integral = integral * dW / 3.0d0
  end function phase_space_integral

  ! log10(f * t_{1/2})
  function log_ft(f, B, g_A_in) result(lft)
    real(8), intent(in) :: f, B, g_A_in
    real(8) :: lft
    real(8) :: t_half
    real(8), parameter :: D = 6289.0d0

    if (f <= 0.0d0 .or. B <= 0.0d0) then
      lft = huge(1.0d0)
      return
    end if
    t_half = D / (g_A_in**2 * f * B)
    lft    = log10(f * t_half)
  end function log_ft

end module mod_phase_space
