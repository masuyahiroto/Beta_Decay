module mod_quenching
  use mod_constants, only: G_A
  implicit none

  real(8), parameter :: Q_STANDARD = 0.74d0
  real(8), parameter :: Q_HEAVY    = 0.70d0

contains

  ! Effective axial coupling: g_A^eff = q * g_A
  function effective_ga(q) result(ga_eff)
    real(8), intent(in) :: q
    real(8) :: ga_eff
    ga_eff = q * G_A
  end function effective_ga

  ! Apply quenching to GT matrix element: M_eff = q * M_GT
  function apply_gt_quenching(M_GT, q) result(M_eff)
    real(8), intent(in) :: M_GT, q
    real(8) :: M_eff
    M_eff = q * M_GT
  end function apply_gt_quenching

  ! Apply quenching to FF matrix elements
  subroutine apply_ff_quenching(xi, zeta, q_gt, q_ff, xi_eff, zeta_eff)
    real(8), intent(in)  :: xi, zeta, q_gt, q_ff
    real(8), intent(out) :: xi_eff, zeta_eff
    xi_eff   = q_ff * xi
    zeta_eff = q_gt * zeta
  end subroutine apply_ff_quenching

  ! Mass-number-dependent quenching factor
  function quenching_from_mass(A) result(q)
    integer, intent(in) :: A
    real(8) :: q
    if (A <= 60) then
      q = Q_STANDARD
    else if (A <= 100) then
      q = 0.72d0
    else
      q = Q_HEAVY
    end if
  end function quenching_from_mass

end module mod_quenching
