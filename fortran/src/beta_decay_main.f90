! ============================================================
! beta_decay_main.f90
! N=126, 125 核種の β⁻ 崩壊半減期計算
! Phys. Rev. C 109, 064319 (2024) の計算手順を実装
!
! 計算フロー:
!   1. 殻模型 (KSHELL) の B(GT), xi, zeta を入力
!   2. クエンチング (q = 0.74) を適用
!   3. 位相空間積分 f0 を数値計算 (シンプソン法, 200点)
!   4. 半減期 T1/2 = kappa / (g_A^eff^2 * B_GT * f0) を計算
!   5. 実験値と比較して log ft を出力
! ============================================================
program beta_decay_n126
  use mod_constants,    only: D_CONST
  use mod_quenching,    only: effective_ga, Q_STANDARD
  use mod_shape_factor, only: ShapeFactor, make_allowed_sf, make_ff_sf
  use mod_phase_space,  only: phase_space_integral
  use mod_half_life,    only: partial_half_life_gt, partial_half_life_ff, &
                               Transition, total_half_life
  implicit none

  ! ----------------------------------------------------------
  ! 核種データ構造
  ! ----------------------------------------------------------
  type :: NuclideDatum
    character(len=12) :: label     ! 核種名
    integer  :: Z_par              ! 親核の陽子数
    integer  :: A                  ! 質量数
    integer  :: Z_dau              ! 娘核の陽子数
    real(8)  :: Q_MeV             ! Q_β [MeV]  (AME2020)
    real(8)  :: B_GT              ! GT 遷移強度 B(GT) (クエンチング適用前)
    real(8)  :: xi                ! FF xi 行列要素 (クエンチング適用前)
    real(8)  :: zeta              ! FF zeta 行列要素 (クエンチング適用前)
    real(8)  :: T_exp_s           ! 実験半減期 [s] (不明=-1)
  end type NuclideDatum

  ! ----------------------------------------------------------
  ! 核種データ
  ! Q_MeV: AME2020 の値
  ! B_GT, xi, zeta: KSHELL 計算結果 (★ 現在はプレースホルダー ★)
  !   → 実際の KSHELL 出力に置き換えること
  ! T_exp_s: NUBASE2020 の実験値
  ! ----------------------------------------------------------
  integer, parameter :: NDATA = 8
  type(NuclideDatum) :: dat(NDATA)

  real(8)  :: ga_eff
  real(8)  :: q_factor
  real(8)  :: f0_gt             ! 純粋な位相空間因子 (GT)
  real(8)  :: f0_ff             ! 純粋な位相空間因子 (FF)
  real(8)  :: t_tot_s           ! 全半減期 [s]
  real(8)  :: lft_gt            ! log10(f0 * T1/2)  ← 正しい log ft
  real(8)  :: ratio             ! T_calc / T_exp
  type(ShapeFactor) :: sf_unit  ! C(W)=1 (位相空間のみ)
  type(ShapeFactor) :: sf_ff    ! FF シェイプファクター
  type(Transition)  :: tr(2)
  integer  :: i, n_tr
  integer, parameter :: FUNIT = 10
  character(len=*), parameter :: DATFILE = 'results.dat'

  ! ----------------------------------------------------------
  ! N=126 核種の設定
  ! ----------------------------------------------------------

  ! 207Tl (Z=81, N=126, A=207) → 207Pb  [N=126]
  dat(1) = NuclideDatum('207Tl>207Pb ', 81, 207, 82, &
                         1.427d0, 0.0010d0, 0.0d0, 0.0d0, 286.2d0)

  ! 206Hg (Z=80, N=126, A=206) → 206Tl  [N=126]
  dat(2) = NuclideDatum('206Hg>206Tl ', 80, 206, 81, &
                         1.307d0, 0.0006d0, 0.0d0, 0.0d0, 509.0d0)

  ! 205Au (Z=79, N=126, A=205) → 205Hg  [N=126]
  dat(3) = NuclideDatum('205Au>205Hg ', 79, 205, 80, &
                         2.660d0, 0.0150d0, 0.0d0, 0.0d0, 186.0d0)

  ! 204Pt (Z=78, N=126, A=204) → 204Au  [N=126]
  dat(4) = NuclideDatum('204Pt>204Au ', 78, 204, 79, &
                         0.510d0, 0.0001d0, 0.0d0, 0.0d0, -1.0d0)

  ! ----------------------------------------------------------
  ! N=125 核種の設定
  ! ----------------------------------------------------------

  ! 206Tl (Z=81, N=125, A=206) → 206Pb  [N=125]
  dat(5) = NuclideDatum('206Tl>206Pb ', 81, 206, 82, &
                         1.534d0, 0.0008d0, 0.0d0, 0.0d0, 248.0d0)

  ! 205Hg (Z=80, N=125, A=205) → 205Tl  [N=125]
  dat(6) = NuclideDatum('205Hg>205Tl ', 80, 205, 81, &
                         1.530d0, 0.0005d0, 0.0d0, 0.0d0, 317.0d0)

  ! ----------------------------------------------------------
  ! 論文の比較例 (N=127 核種も含む)
  ! ----------------------------------------------------------

  ! 208Tl (Z=81, N=127, A=208) → 208Pb  [GT + FF]
  dat(7) = NuclideDatum('208Tl>208Pb ', 81, 208, 82, &
                         4.990d0, 0.1200d0, 0.30d0, 0.10d0, 183.2d0)

  ! 207Hg (Z=80, N=127, A=207) → 207Tl  [GT + FF]
  dat(8) = NuclideDatum('207Hg>207Tl ', 80, 207, 81, &
                         1.370d0, 0.0100d0, 0.10d0, 0.00d0, 174.0d0)

  ! ----------------------------------------------------------
  ! クエンチング: q = 0.74 (論文 Phys.Rev.C 109, 064319)
  ! ----------------------------------------------------------
  q_factor = Q_STANDARD   ! 0.74
  ga_eff   = effective_ga(q_factor)

  ! ----------------------------------------------------------
  ! dat ファイルを開く
  ! ----------------------------------------------------------
  open(unit=FUNIT, file=DATFILE, status='replace', action='write')
  write(FUNIT,'(A)') '# label  A  Z  Q_MeV  T_calc_s  T_exp_s  ratio  logft'

  ! ----------------------------------------------------------
  ! ヘッダー
  ! ----------------------------------------------------------
  write(*,*)
  write(*,'(A)') repeat('=', 78)
  write(*,'(A)') '  β⁻ 崩壊半減期計算  (N=126, 125 核種)'
  write(*,'(A)') '  Phys. Rev. C 109, 064319 (2024) の計算手順'
  write(*,'(A,F5.2,A,F7.4)') &
        '  クエンチング q = ', q_factor, '   g_A^eff = ', ga_eff
  write(*,'(A)') repeat('=', 78)
  write(*,'(A12,A4,A5,A9,A13,A13,A9,A8)') &
        'Nucleus', 'Z', 'A', 'Q[MeV]', 'T_calc[s]', &
        'T_exp[s]', 'Calc/Exp', 'logft'
  write(*,'(A)') repeat('-', 78)

  ! ----------------------------------------------------------
  ! 核種ごとの計算
  ! ----------------------------------------------------------
  do i = 1, NDATA

    n_tr = 0

    ! ---- GT 遷移 ----
    if (dat(i)%B_GT > 0.0d0) then
      ! 位相空間因子 f0 を C(W)=1 で計算 (B_GT は別途掛ける)
      sf_unit = make_allowed_sf(1.0d0)
      f0_gt   = phase_space_integral(dat(i)%Z_dau, dat(i)%Q_MeV, sf_unit)

      n_tr = n_tr + 1
      tr(n_tr)%Q_MeV  = dat(i)%Q_MeV
      tr(n_tr)%B_GT   = dat(i)%B_GT
      tr(n_tr)%xi     = 0.0d0
      tr(n_tr)%zeta   = 0.0d0
      tr(n_tr)%Z_dau  = dat(i)%Z_dau
      tr(n_tr)%ga_eff = ga_eff
    else
      f0_gt = 0.0d0
    end if

    ! ---- FF 遷移 ----
    if (abs(dat(i)%xi) > 0.0d0 .or. abs(dat(i)%zeta) > 0.0d0) then
      sf_ff = make_ff_sf(dat(i)%xi, dat(i)%zeta)
      f0_ff = phase_space_integral(dat(i)%Z_dau, dat(i)%Q_MeV, sf_ff)

      n_tr = n_tr + 1
      tr(n_tr)%Q_MeV  = dat(i)%Q_MeV
      tr(n_tr)%B_GT   = 0.0d0
      tr(n_tr)%xi     = dat(i)%xi
      tr(n_tr)%zeta   = dat(i)%zeta
      tr(n_tr)%Z_dau  = dat(i)%Z_dau
      tr(n_tr)%ga_eff = ga_eff
    else
      f0_ff = 0.0d0
    end if

    ! ---- 全半減期 ----
    if (n_tr > 0) then
      t_tot_s = total_half_life(tr(1:n_tr), n_tr)
    else
      t_tot_s = huge(1.0d0)
    end if

    ! ---- log ft (GT 成分から計算) ----
    ! log ft = log10(kappa / (g_A^eff^2 * B_GT))
    if (dat(i)%B_GT > 0.0d0) then
      lft_gt = log10(D_CONST / (ga_eff**2 * dat(i)%B_GT))
    else
      lft_gt = huge(1.0d0)
    end if

    ! ---- 実験値との比 ----
    if (dat(i)%T_exp_s > 0.0d0 .and. t_tot_s < huge(1.0d0)*0.9d0) then
      ratio = t_tot_s / dat(i)%T_exp_s
    else
      ratio = -1.0d0
    end if

    ! ---- 標準出力 ----
    if (t_tot_s < huge(1.0d0)*0.9d0) then
      if (dat(i)%T_exp_s > 0.0d0) then
        write(*,'(A12,I4,I5,F9.3,F13.1,F13.1,F9.3,F8.2)') &
          trim(dat(i)%label), dat(i)%Z_par, dat(i)%A, &
          dat(i)%Q_MeV, t_tot_s, dat(i)%T_exp_s, ratio, lft_gt
      else
        write(*,'(A12,I4,I5,F9.3,F13.1,A13,A9,F8.2)') &
          trim(dat(i)%label), dat(i)%Z_par, dat(i)%A, &
          dat(i)%Q_MeV, t_tot_s, '    unknown  ', '   ---  ', lft_gt
      end if
    else
      write(*,'(A12,I4,I5,F9.3,A13,F13.1,A9,A8)') &
        trim(dat(i)%label), dat(i)%Z_par, dat(i)%A, &
        dat(i)%Q_MeV, '   (invalid) ', dat(i)%T_exp_s, '   ---  ', '  ---   '
    end if

    ! ---- dat ファイル出力 ----
    if (t_tot_s < huge(1.0d0)*0.9d0) then
      write(FUNIT,'(A12,I5,I4,F8.3,ES14.4,ES14.4,F9.3,F8.2)') &
        trim(dat(i)%label), dat(i)%A, dat(i)%Z_par, dat(i)%Q_MeV, &
        t_tot_s, dat(i)%T_exp_s, ratio, lft_gt
    else
      write(FUNIT,'(A12,I5,I4,F8.3,A14,ES14.4,A9,A8)') &
        trim(dat(i)%label), dat(i)%A, dat(i)%Z_par, dat(i)%Q_MeV, &
        '           NaN', dat(i)%T_exp_s, '      ---', '     ---'
    end if

  end do

  write(*,'(A)') repeat('=', 78)
  write(*,'(A)') '  ★ B_GT, xi, zeta は現在プレースホルダー値'
  write(*,'(A)') '    KSHELL 計算結果の OBDM から以下の手順で更新:'
  write(*,'(A)') '      B_GT = |M_GT|^2 / (2Ji+1)'
  write(*,'(A)') '      xi   = sum_pn <p||r C1||n> * rho_pn  (rank 0,2 成分)'
  write(*,'(A)') '      zeta = sum_pn <p||r C1 x sigma||n> * rho_pn  (rank 1 成分)'
  write(*,'(A)') repeat('=', 78)
  write(*,*)
  write(*,'(A,A)') '  dat ファイルを出力しました: ', DATFILE

  close(FUNIT)

end program beta_decay_n126
