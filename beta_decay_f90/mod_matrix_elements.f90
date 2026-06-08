module mod_matrix_elements
  use mod_wave_function
  implicit none

contains

  ! Single-particle reduced GT matrix element <jp || sigma || jn>
  ! Requires lp == ln (same orbital angular momentum)
  function sp_gt_reduced_me(p, n) result(me)
    type(SingleParticleState), intent(in) :: p, n
    real(8) :: me
    integer :: l, jp2, jn2

    me = 0.0d0
    if (p%l /= n%l) return

    l   = p%l
    jp2 = p%j2
    jn2 = n%j2

    if      (jp2 == 2*l+1 .and. jn2 == 2*l+1) then
      me =  sqrt(dble(2*(jp2+1)))
    else if (jp2 == 2*l+1 .and. jn2 == 2*l-1) then
      me = -sqrt(dble(jp2+1))
    else if (jp2 == 2*l-1 .and. jn2 == 2*l+1) then
      me = -sqrt(dble(jn2+1))
    else if (jp2 == 2*l-1 .and. jn2 == 2*l-1) then
      me = -sqrt(dble(2*(jp2+1)))
    end if
  end function sp_gt_reduced_me

  ! Raw GT matrix element M_GT = sum_{pn} <p||sigma||n> * rho_{pn}
  function gt_matrix_element(state) result(M_GT)
    type(NuclearState), intent(in) :: state
    real(8) :: M_GT
    integer :: k
    real(8) :: sp_me

    M_GT = 0.0d0
    do k = 1, state%n_obdm
      sp_me = sp_gt_reduced_me(state%sp(state%obdm(k)%p_idx), &
                                state%sp(state%obdm(k)%n_idx))
      M_GT = M_GT + sp_me * state%obdm(k)%val
    end do
  end function gt_matrix_element

  ! GT transition strength B(GT) = M_GT^2 / (2*Ji+1)
  function gamow_teller_strength(state, Ji2) result(B_GT)
    type(NuclearState), intent(in) :: state
    integer, intent(in) :: Ji2
    real(8) :: B_GT
    real(8) :: M_GT

    M_GT = gt_matrix_element(state)
    B_GT = M_GT**2 / dble(Ji2 + 1)
  end function gamow_teller_strength

  ! First-forbidden radial matrix element (approximation: <r> ~ r_mean if Δl=1)
  function sp_ff_r_me(p, n, r_mean) result(me)
    type(SingleParticleState), intent(in) :: p, n
    real(8), intent(in) :: r_mean
    real(8) :: me
    me = 0.0d0
    if (abs(p%l - n%l) == 1) me = r_mean
  end function sp_ff_r_me

  ! First-forbidden xi matrix element
  function ff_xi(state, r_mean) result(xi)
    type(NuclearState), intent(in) :: state
    real(8), intent(in) :: r_mean
    real(8) :: xi
    integer :: k
    xi = 0.0d0
    do k = 1, state%n_obdm
      xi = xi + sp_ff_r_me(state%sp(state%obdm(k)%p_idx), &
                            state%sp(state%obdm(k)%n_idx), r_mean) &
               * state%obdm(k)%val
    end do
  end function ff_xi

  ! First-forbidden zeta matrix element
  function ff_zeta(state, r_mean) result(zeta)
    type(NuclearState), intent(in) :: state
    real(8), intent(in) :: r_mean
    real(8) :: zeta
    integer :: k
    zeta = 0.0d0
    do k = 1, state%n_obdm
      zeta = zeta + sp_ff_r_me(state%sp(state%obdm(k)%p_idx), &
                                state%sp(state%obdm(k)%n_idx), r_mean) &
                   * state%obdm(k)%val
    end do
  end function ff_zeta

end module mod_matrix_elements
