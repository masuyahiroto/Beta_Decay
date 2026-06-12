! ============================================================
! debug_calc.f90  --  計算過程の中間値を全出力するデバッグ用プログラム
! ============================================================
program debug_calc
  use mod_constants,    only: D_CONST, ME_MEV, ALPHA, G_A, PI
  use mod_quenching,    only: effective_ga, Q_STANDARD
  use mod_shape_factor, only: ShapeFactor, make_allowed_sf, make_ff_sf, eval_shape_factor
  use mod_phase_space,  only: phase_space_integral, fermi_function
  implicit none

  integer, parameter :: NDATA = 8
  integer, parameter :: FUNIT = 20

  character(len=12) :: label(NDATA)
  integer  :: Z_par(NDATA), A_num(NDATA), Z_dau(NDATA)
  real(8)  :: Q_MeV(NDATA), B_GT(NDATA), xi(NDATA), zeta(NDATA), T_exp(NDATA)

  real(8)  :: q_factor, ga_eff
  real(8)  :: W0, f0, f_GT, T_calc, ratio, B_GT_needed
  real(8)  :: W_vals(5), F_vals(5)
  type(ShapeFactor) :: sf
  integer  :: i, j

  ! ------ 核種データ ------
  label = ['207Tl>207Pb ', '206Hg>206Tl ', '205Au>205Hg ', '204Pt>204Au ', &
           '206Tl>206Pb ', '205Hg>205Tl ', '208Tl>208Pb ', '207Hg>207Tl ']
  Z_par = [81, 80, 79, 78, 81, 80, 81, 80]
  A_num = [207, 206, 205, 204, 206, 205, 208, 207]
  Z_dau = [82, 81, 80, 79, 82, 81, 82, 81]
  Q_MeV = [1.427d0, 1.307d0, 2.660d0, 0.510d0, 1.534d0, 1.530d0, 4.990d0, 1.370d0]
  B_GT  = [0.0010d0, 0.0006d0, 0.0150d0, 0.0001d0, 0.0008d0, 0.0005d0, 0.1200d0, 0.0100d0]
  xi    = [0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.30d0, 0.10d0]
  zeta  = [0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.10d0, 0.00d0]
  T_exp = [286.2d0, 509.0d0, 186.0d0, -1.0d0, 248.0d0, 317.0d0, 183.2d0, 174.0d0]

  q_factor = Q_STANDARD
  ga_eff   = effective_ga(q_factor)

  open(unit=FUNIT, file='debug_output.txt', status='replace', action='write')

  ! ===== Section 1: 物理定数 =====
  write(FUNIT,'(A)') '=== SECTION 1: 物理定数 ==='
  write(FUNIT,'(A,F12.6)') 'kappa (D_CONST) [s]  = ', D_CONST
  write(FUNIT,'(A,F12.6)') 'g_A              = ', G_A
  write(FUNIT,'(A,F12.9)') 'm_e c^2  [MeV]   = ', ME_MEV
  write(FUNIT,'(A,F14.12)') 'alpha            = ', ALPHA

  ! ===== Section 2: クエンチング =====
  write(FUNIT,'(/,A)') '=== SECTION 2: クエンチング ==='
  write(FUNIT,'(A,F6.2)') 'q          = ', q_factor
  write(FUNIT,'(A,F8.4)') 'g_A^eff    = ', ga_eff
  write(FUNIT,'(A,F8.4)') 'g_A^eff^2  = ', ga_eff**2

  ! ===== Section 3: W0の計算 =====
  write(FUNIT,'(/,A)') '=== SECTION 3: W0 = Q / m_e + 1 ==='
  write(FUNIT,'(A12,A8,A10)') 'Nucleus', 'Q[MeV]', 'W0'
  do i = 1, NDATA
    W0 = Q_MeV(i) / ME_MEV + 1.0d0
    write(FUNIT,'(A12,F8.3,F10.4)') trim(label(i)), Q_MeV(i), W0
  end do

  ! ===== Section 4: Fermi関数 =====
  write(FUNIT,'(/,A)') '=== SECTION 4: Fermi関数 F(Z,W) ==='
  write(FUNIT,'(A)') '  eta = alpha*Z*W/p,  F = 2*pi*eta / (1 - exp(-2*pi*eta))'
  do i = 1, NDATA
    W0 = Q_MeV(i) / ME_MEV + 1.0d0
    W_vals(1) = 1.01d0
    W_vals(2) = 1.0d0 + (W0-1.0d0)*0.25d0
    W_vals(3) = 1.0d0 + (W0-1.0d0)*0.50d0
    W_vals(4) = 1.0d0 + (W0-1.0d0)*0.75d0
    W_vals(5) = W0 * 0.999d0
    write(FUNIT,'(/,A,A,A,I3,A,F6.3,A,F6.3,A)') &
      '  [', trim(label(i)), ']  Z_dau=', Z_dau(i), &
      '  W in [1, ', W0, ']  (W0=', W0, ')'
    write(FUNIT,'(A8,A12,A12,A12)') 'W', 'p=sqrt(W2-1)', 'eta', 'F(Z,W)'
    do j = 1, 5
      F_vals(j) = fermi_function(Z_dau(i), W_vals(j))
      write(FUNIT,'(F8.4,F12.6,F12.6,F12.4)') &
        W_vals(j), sqrt(max(W_vals(j)**2-1.0d0, 0.0d0)), &
        ALPHA * dble(Z_dau(i)) * W_vals(j) / max(sqrt(max(W_vals(j)**2-1.0d0,1d-10)),1d-10), &
        F_vals(j)
    end do
  end do

  ! ===== Section 5: 位相空間積分 f0 =====
  write(FUNIT,'(/,A)') '=== SECTION 5: 位相空間積分 f0 (C(W)=1) ==='
  write(FUNIT,'(A)') '  f0 = integral_1^W0 F(Z,W) * p * W * (W0-W)^2 dW  (Simpson 200点)'
  write(FUNIT,'(A12,A8,A10,A14)') 'Nucleus', 'Q[MeV]', 'W0', 'f0'
  do i = 1, NDATA
    sf = make_allowed_sf(1.0d0)
    f0 = phase_space_integral(Z_dau(i), Q_MeV(i), sf)
    write(FUNIT,'(A12,F8.3,F10.4,ES14.5)') trim(label(i)), Q_MeV(i), &
      Q_MeV(i)/ME_MEV+1.0d0, f0
  end do

  ! ===== Section 6: GT半減期の計算 =====
  write(FUNIT,'(/,A)') '=== SECTION 6: GT半減期  T = D / (g_A^eff^2 * B_GT * f0) ==='
  write(FUNIT,'(A12,A10,A12,A14,A14,A14,A10,A14)') &
    'Nucleus', 'B_GT', 'f0', 'f=B_GT*f0', 'T_calc[s]', 'T_exp[s]', 'ratio', 'B_GT_needed'
  do i = 1, NDATA
    sf = make_allowed_sf(1.0d0)
    f0 = phase_space_integral(Z_dau(i), Q_MeV(i), sf)
    f_GT  = B_GT(i) * f0
    T_calc = D_CONST / (ga_eff**2 * f_GT)
    if (T_exp(i) > 0.0d0) then
      ratio = T_calc / T_exp(i)
      B_GT_needed = D_CONST / (ga_eff**2 * f0 * T_exp(i))
    else
      ratio = -1.0d0
      B_GT_needed = -1.0d0
    end if
    write(FUNIT,'(A12,ES10.3,ES12.4,ES14.4,ES14.4,ES14.4,F10.3,ES14.4)') &
      trim(label(i)), B_GT(i), f0, f_GT, T_calc, T_exp(i), ratio, B_GT_needed
  end do

  ! ===== Section 7: B_GTの比較 =====
  write(FUNIT,'(/,A)') '=== SECTION 7: B_GT プレースホルダー vs 必要値 ==='
  write(FUNIT,'(A12,A12,A14,A10)') 'Nucleus', 'B_GT_input', 'B_GT_needed', 'ratio'
  do i = 1, NDATA
    if (T_exp(i) < 0.0d0) cycle
    sf = make_allowed_sf(1.0d0)
    f0 = phase_space_integral(Z_dau(i), Q_MeV(i), sf)
    B_GT_needed = D_CONST / (ga_eff**2 * f0 * T_exp(i))
    write(FUNIT,'(A12,ES12.4,ES14.4,F10.2)') &
      trim(label(i)), B_GT(i), B_GT_needed, B_GT_needed/B_GT(i)
  end do

  close(FUNIT)
  write(*,'(A)') 'debug_output.txt を生成しました'

end program debug_calc
