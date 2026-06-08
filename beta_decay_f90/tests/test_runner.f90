program test_runner
  use mod_constants
  use mod_wave_function
  use mod_matrix_elements
  use mod_quenching
  use mod_shape_factor
  use mod_phase_space
  use mod_half_life
  implicit none

  integer :: n_pass, n_fail
  n_pass = 0; n_fail = 0

  call test_wave_function(n_pass, n_fail)
  call test_matrix_elements(n_pass, n_fail)
  call test_quenching(n_pass, n_fail)
  call test_shape_factor(n_pass, n_fail)
  call test_phase_space(n_pass, n_fail)
  call test_half_life(n_pass, n_fail)

  write(*,*)
  write(*,'(A,I4,A,I4,A)') '=== Results: ', n_pass, ' passed, ', n_fail, ' failed ==='
  if (n_fail > 0) stop 1

contains

  subroutine assert_near(label, got, expected, tol, np, nf)
    character(*), intent(in) :: label
    real(8), intent(in) :: got, expected, tol
    integer, intent(inout) :: np, nf
    if (abs(got - expected) <= tol) then
      write(*,'(A,A40,A)') '  PASS  ', label, ''
      np = np + 1
    else
      write(*,'(A,A40,A,F15.8,A,F15.8)') '  FAIL  ', label, &
            '  got=', got, '  expected=', expected
      nf = nf + 1
    end if
  end subroutine assert_near

  subroutine assert_true(label, cond, np, nf)
    character(*), intent(in) :: label
    logical, intent(in) :: cond
    integer, intent(inout) :: np, nf
    if (cond) then
      write(*,'(A,A40)') '  PASS  ', label
      np = np + 1
    else
      write(*,'(A,A40)') '  FAIL  ', label
      nf = nf + 1
    end if
  end subroutine assert_true

  ! ----------------------------------------------------------------
  subroutine test_wave_function(np, nf)
    integer, intent(inout) :: np, nf
    type(NuclearState) :: st
    integer :: p_idx, n_idx
    real(8) :: val

    write(*,*) '--- test_wave_function ---'
    call init_state(st, J2=7, parity=-1, Z=20, N=21)
    call add_sp_state(st, 0, 3, 7, -1, p_idx)  ! f7/2 proton
    call add_sp_state(st, 0, 3, 7, +1, n_idx)  ! f7/2 neutron
    call set_obdm(st, p_idx, n_idx, 0.8d0)

    call assert_near('obdm set/get', get_obdm(st,p_idx,n_idx), 0.8d0, 1d-10, np, nf)
    call assert_near('missing obdm = 0', get_obdm(st,2,1), 0.0d0, 1d-10, np, nf)
    call assert_true('J2 correct', st%J2 == 7, np, nf)
    call assert_true('Z correct',  st%Z  == 20, np, nf)
    call assert_true('A correct',  st%Z + st%N == 41, np, nf)
    call assert_true('n_sp = 2',   st%n_sp == 2, np, nf)
  end subroutine test_wave_function

  ! ----------------------------------------------------------------
  subroutine test_matrix_elements(np, nf)
    integer, intent(inout) :: np, nf
    type(SingleParticleState) :: p, n_st
    type(NuclearState) :: state
    integer :: p_idx, n_idx
    real(8) :: me, B

    write(*,*) '--- test_matrix_elements ---'

    ! f7/2 -> f7/2: should be sqrt(2*(7+1)) = sqrt(16) = 4
    p%n=0; p%l=3; p%j2=7; p%tz2=-1
    n_st%n=0; n_st%l=3; n_st%j2=7; n_st%tz2=1
    me = sp_gt_reduced_me(p, n_st)
    call assert_near('f7/2->f7/2 sp_me', me, sqrt(16.0d0), 1d-10, np, nf)

    ! f7/2 -> f5/2: should be -sqrt(8)
    n_st%j2 = 5
    me = sp_gt_reduced_me(p, n_st)
    call assert_near('f7/2->f5/2 sp_me', me, -sqrt(8.0d0), 1d-10, np, nf)

    ! different l -> 0
    n_st%l=2; n_st%j2=5
    me = sp_gt_reduced_me(p, n_st)
    call assert_near('diff l -> 0', me, 0.0d0, 1d-10, np, nf)

    ! s1/2 -> s1/2: sqrt(2*2) = 2
    p%l=0; p%j2=1
    n_st%l=0; n_st%j2=1
    me = sp_gt_reduced_me(p, n_st)
    call assert_near('s1/2->s1/2 sp_me', me, 2.0d0, 1d-10, np, nf)

    ! B(GT) from OBDM
    call init_state(state, J2=7, parity=-1, Z=20, N=21)
    call add_sp_state(state, 0, 3, 7, -1, p_idx)
    call add_sp_state(state, 0, 3, 7, +1, n_idx)
    call set_obdm(state, p_idx, n_idx, 1.0d0)
    B = gamow_teller_strength(state, 7)
    call assert_near('B(GT) f7/2', B, 16.0d0/8.0d0, 1d-10, np, nf)
    call assert_true('B(GT) >= 0', B >= 0.0d0, np, nf)
  end subroutine test_matrix_elements

  ! ----------------------------------------------------------------
  subroutine test_quenching(np, nf)
    integer, intent(inout) :: np, nf
    real(8) :: ga_eff, M_eff, xi_eff, zeta_eff
    real(8) :: q_light, q_heavy

    write(*,*) '--- test_quenching ---'

    ga_eff = effective_ga(0.74d0)
    call assert_near('g_A^eff q=0.74', ga_eff, 0.74d0*G_A, 1d-10, np, nf)

    ga_eff = effective_ga(1.0d0)
    call assert_near('g_A^eff q=1', ga_eff, G_A, 1d-10, np, nf)

    M_eff = apply_gt_quenching(2.0d0, 0.74d0)
    call assert_near('GT quench 2.0', M_eff, 0.74d0*2.0d0, 1d-10, np, nf)

    call apply_ff_quenching(1.0d0, 1.0d0, 0.7d0, 0.9d0, xi_eff, zeta_eff)
    call assert_near('FF xi_eff', xi_eff, 0.9d0, 1d-10, np, nf)
    call assert_near('FF zeta_eff', zeta_eff, 0.7d0, 1d-10, np, nf)

    q_light = quenching_from_mass(40)
    q_heavy = quenching_from_mass(208)
    call assert_near('q light nucleus', q_light, Q_STANDARD, 1d-10, np, nf)
    call assert_near('q heavy nucleus', q_heavy, Q_HEAVY,    1d-10, np, nf)
    call assert_true('q heavy <= q light', q_heavy <= q_light, np, nf)
  end subroutine test_quenching

  ! ----------------------------------------------------------------
  subroutine test_shape_factor(np, nf)
    integer, intent(inout) :: np, nf
    type(ShapeFactor) :: sf
    real(8) :: C1, C2, C3

    write(*,*) '--- test_shape_factor ---'

    ! Allowed GT: constant
    sf = make_allowed_sf(2.5d0)
    C1 = eval_shape_factor(sf, 1.5d0)
    C2 = eval_shape_factor(sf, 5.0d0)
    call assert_near('allowed C(1.5)', C1, 2.5d0, 1d-10, np, nf)
    call assert_near('allowed C(5.0)', C2, 2.5d0, 1d-10, np, nf)

    ! Allowed GT = 0
    sf = make_allowed_sf(0.0d0)
    call assert_near('allowed C=0 for B=0', eval_shape_factor(sf, 2.0d0), 0.0d0, 1d-10, np, nf)

    ! FF: xi=1, zeta=0 -> C = xi^2
    sf = make_ff_sf(1.5d0, 0.0d0)
    call assert_near('FF xi-only', eval_shape_factor(sf, 3.0d0), 1.5d0**2, 1d-10, np, nf)

    ! FF formula at W=2: C = (xi - zeta*W)^2 + zeta^2*(W^2-1)/3
    sf = make_ff_sf(1.0d0, 0.5d0)
    C3 = (1.0d0 - 0.5d0*2.0d0)**2 + 0.5d0**2*(4.0d0-1.0d0)/3.0d0
    call assert_near('FF formula at W=2', eval_shape_factor(sf, 2.0d0), C3, 1d-10, np, nf)

    ! FF positive definite
    sf = make_ff_sf(1.0d0, 1.0d0)
    call assert_true('FF C >= 0 at W=1.5', eval_shape_factor(sf, 1.5d0) >= 0.0d0, np, nf)
    call assert_true('FF C >= 0 at W=5.0', eval_shape_factor(sf, 5.0d0) >= 0.0d0, np, nf)
  end subroutine test_shape_factor

  ! ----------------------------------------------------------------
  subroutine test_phase_space(np, nf)
    integer, intent(inout) :: np, nf
    type(ShapeFactor) :: sf
    real(8) :: F, f1, f2, f3, f4

    write(*,*) '--- test_phase_space ---'

    ! Fermi function: threshold
    F = fermi_function(50, 1.0d0)
    call assert_near('F(Z,W=1)=0', F, 0.0d0, 1d-10, np, nf)

    ! Fermi function: below threshold
    F = fermi_function(10, 0.9d0)
    call assert_near('F(Z,W<1)=0', F, 0.0d0, 1d-10, np, nf)

    ! Fermi function: positive for valid W
    F = fermi_function(20, 2.0d0)
    call assert_true('F > 0', F > 0.0d0, np, nf)

    ! Fermi function: high Z > low Z
    call assert_true('F high-Z > low-Z', &
         fermi_function(82, 2.0d0) > fermi_function(1, 2.0d0), np, nf)

    ! Phase space: Q=0 -> 0
    sf = make_allowed_sf(1.0d0)
    call assert_near('f(Q=0)=0', phase_space_integral(20, 0.0d0, sf), 0.0d0, 1d-10, np, nf)

    ! Phase space: Q<0 -> 0
    call assert_near('f(Q<0)=0', phase_space_integral(20, -1.0d0, sf), 0.0d0, 1d-10, np, nf)

    ! Phase space: positive for valid decay
    f1 = phase_space_integral(20, 5.0d0, sf)
    call assert_true('f > 0 valid', f1 > 0.0d0, np, nf)

    ! Phase space: increases with Q
    f1 = phase_space_integral(20, 2.0d0, sf)
    f2 = phase_space_integral(20, 5.0d0, sf)
    call assert_true('f increases with Q', f2 > f1, np, nf)

    ! Phase space: increases with Z
    f3 = phase_space_integral(1,  3.0d0, sf)
    f4 = phase_space_integral(82, 3.0d0, sf)
    call assert_true('f increases with Z', f4 > f3, np, nf)

    ! Neutron decay: f ~ 1.686 (standard value), check within 10%
    sf = make_fermi_sf(1.0d0)
    f1 = phase_space_integral(1, 0.7824d0, sf)
    call assert_true('neutron f in [1.4, 2.0]', f1 > 1.4d0 .and. f1 < 2.0d0, np, nf)

    ! B(GT) scales f linearly
    f1 = phase_space_integral(20, 5.0d0, make_allowed_sf(1.0d0))
    f2 = phase_space_integral(20, 5.0d0, make_allowed_sf(2.0d0))
    call assert_near('f scales with B_GT', f2/f1, 2.0d0, 1d-6, np, nf)
  end subroutine test_phase_space

  ! ----------------------------------------------------------------
  subroutine test_half_life(np, nf)
    integer, intent(inout) :: np, nf
    real(8) :: t1, t2, t3
    type(Transition) :: tr(3)

    write(*,*) '--- test_half_life ---'

    ! Valid decay -> finite positive
    t1 = partial_half_life_gt(20, 5.0d0, 1.0d0, G_A)
    call assert_true('t>0 valid GT', t1 > 0.0d0 .and. t1 < huge(1.0d0), np, nf)

    ! B=0 -> infinite
    t1 = partial_half_life_gt(20, 5.0d0, 0.0d0, G_A)
    call assert_true('t=inf for B=0', t1 >= huge(1.0d0)*0.9d0, np, nf)

    ! Q=0 -> infinite
    t1 = partial_half_life_gt(20, 0.0d0, 1.0d0, G_A)
    call assert_true('t=inf for Q=0', t1 >= huge(1.0d0)*0.9d0, np, nf)

    ! Larger B -> shorter t
    t1 = partial_half_life_gt(20, 5.0d0, 1.0d0, G_A)
    t2 = partial_half_life_gt(20, 5.0d0, 2.0d0, G_A)
    call assert_true('larger B -> shorter t', t2 < t1, np, nf)

    ! Larger Q -> shorter t
    t1 = partial_half_life_gt(20, 3.0d0, 1.0d0, G_A)
    t2 = partial_half_life_gt(20, 6.0d0, 1.0d0, G_A)
    call assert_true('larger Q -> shorter t', t2 < t1, np, nf)

    ! Quenching increases t
    t1 = partial_half_life_gt(20, 5.0d0, 1.0d0, G_A)
    t2 = partial_half_life_gt(20, 5.0d0, 1.0d0, 0.74d0*G_A)
    call assert_true('quenching -> longer t', t2 > t1, np, nf)

    ! Neutron-like decay: t in [1000, 8000] s (simplified GT only)
    t1 = partial_half_life_gt(1, 0.7824d0, 1.0d0, 1.0d0)
    call assert_true('neutron-like t [1000,8000]s', t1 > 1000.0d0 .and. t1 < 8000.0d0, np, nf)

    ! FF: valid decay
    t1 = partial_half_life_ff(82, 5.0d0, 1.0d0, 0.5d0, G_A)
    call assert_true('FF t > 0', t1 > 0.0d0 .and. t1 < huge(1.0d0), np, nf)

    ! Total half-life: 2 transitions -> shorter than either
    tr(1)%Q_MeV=5.0d0; tr(1)%B_GT=1.0d0; tr(1)%Z_dau=20; tr(1)%ga_eff=G_A
    tr(2)%Q_MeV=4.0d0; tr(2)%B_GT=0.5d0; tr(2)%Z_dau=20; tr(2)%ga_eff=G_A
    t1 = partial_half_life_gt(20, 5.0d0, 1.0d0, G_A)
    t2 = partial_half_life_gt(20, 4.0d0, 0.5d0, G_A)
    t3 = total_half_life(tr(1:2), 2)
    call assert_true('total < each', t3 < t1 .and. t3 < t2, np, nf)
  end subroutine test_half_life

end program test_runner
