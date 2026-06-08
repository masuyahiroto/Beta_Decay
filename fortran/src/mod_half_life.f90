module mod_half_life
  use mod_constants,    only: D_CONST, G_A
  use mod_shape_factor, only: ShapeFactor, make_allowed_sf, make_ff_sf
  use mod_phase_space,  only: phase_space_integral
  implicit none

  ! Single decay transition
  type :: Transition
    character(len=20) :: label    = ''
    real(8)           :: Q_MeV   = 0.0d0
    real(8)           :: B_GT    = 0.0d0
    real(8)           :: xi      = 0.0d0
    real(8)           :: zeta    = 0.0d0
    integer           :: Z_dau   = 0       ! daughter Z
    real(8)           :: ga_eff  = G_A
  end type Transition

  ! Result container
  type :: HalfLifeResult
    character(len=20) :: label
    real(8) :: f            ! phase space integral
    real(8) :: t_half_s     ! [seconds]
    real(8) :: t_half_ms    ! [milliseconds]
    real(8) :: log_ft_val
    real(8) :: exp_ms       ! experimental value (-1 = N/A)
  end type HalfLifeResult

contains

  ! Partial half-life for allowed GT transition [seconds]
  function partial_half_life_gt(Z, Q_MeV, B_GT, ga_eff) result(t)
    integer, intent(in) :: Z
    real(8), intent(in) :: Q_MeV, B_GT, ga_eff
    real(8) :: t
    type(ShapeFactor) :: sf
    real(8) :: f

    if (B_GT <= 0.0d0 .or. Q_MeV <= 0.0d0) then
      t = huge(1.0d0)
      return
    end if
    sf = make_allowed_sf(B_GT)
    f  = phase_space_integral(Z, Q_MeV, sf)
    if (f <= 0.0d0) then
      t = huge(1.0d0)
    else
      t = D_CONST / (ga_eff**2 * f)
    end if
  end function partial_half_life_gt

  ! Partial half-life for first-forbidden transition [seconds]
  function partial_half_life_ff(Z, Q_MeV, xi, zeta, ga_eff) result(t)
    integer, intent(in) :: Z
    real(8), intent(in) :: Q_MeV, xi, zeta, ga_eff
    real(8) :: t
    type(ShapeFactor) :: sf
    real(8) :: f

    if (Q_MeV <= 0.0d0) then
      t = huge(1.0d0)
      return
    end if
    sf = make_ff_sf(xi, zeta)
    f  = phase_space_integral(Z, Q_MeV, sf)
    if (f <= 0.0d0) then
      t = huge(1.0d0)
    else
      t = D_CONST / (ga_eff**2 * f)
    end if
  end function partial_half_life_ff

  ! Total half-life from multiple transitions (rates add)
  function total_half_life(transitions, n_tr) result(t_total)
    type(Transition), intent(in) :: transitions(:)
    integer,          intent(in) :: n_tr
    real(8) :: t_total
    integer :: i
    real(8) :: t_i, rate

    rate = 0.0d0
    do i = 1, n_tr
      if (transitions(i)%B_GT > 0.0d0) then
        t_i = partial_half_life_gt(transitions(i)%Z_dau, &
                                    transitions(i)%Q_MeV, &
                                    transitions(i)%B_GT,  &
                                    transitions(i)%ga_eff)
      else if (abs(transitions(i)%xi) > 0.0d0 .or. &
               abs(transitions(i)%zeta) > 0.0d0) then
        t_i = partial_half_life_ff(transitions(i)%Z_dau,  &
                                    transitions(i)%Q_MeV,  &
                                    transitions(i)%xi,     &
                                    transitions(i)%zeta,   &
                                    transitions(i)%ga_eff)
      else
        cycle
      end if
      if (t_i > 0.0d0 .and. t_i < huge(1.0d0)) rate = rate + 1.0d0 / t_i
    end do

    if (rate <= 0.0d0) then
      t_total = huge(1.0d0)
    else
      t_total = 1.0d0 / rate
    end if
  end function total_half_life

  ! Compute result for one transition and print comparison
  subroutine compute_and_print(tr, exp_ms)
    type(Transition), intent(in) :: tr
    real(8),          intent(in) :: exp_ms
    type(ShapeFactor) :: sf
    real(8) :: f, t_s, t_ms, lft, ratio

    if (tr%B_GT > 0.0d0) then
      sf = make_allowed_sf(tr%B_GT)
    else
      sf = make_ff_sf(tr%xi, tr%zeta)
    end if

    f    = phase_space_integral(tr%Z_dau, tr%Q_MeV, sf)
    t_s  = D_CONST / (tr%ga_eff**2 * f)
    t_ms = t_s * 1.0d3
    lft  = log10(f * t_s)

    write(*,'(A20,3F12.2,F10.4)') trim(tr%label), t_ms, exp_ms, &
                                   t_ms/exp_ms, lft
  end subroutine compute_and_print

  ! Print comparison table header
  subroutine print_table_header()
    write(*,'(A20,A12,A12,A10,A10)') 'Nucleus', 'Calc[ms]', 'Exp[ms]', 'Ratio', 'log(ft)'
    write(*,*) repeat('-', 65)
  end subroutine print_table_header

end module mod_half_life
