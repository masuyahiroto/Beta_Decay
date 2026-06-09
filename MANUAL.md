# 取扱説明書 — β崩壊半減期計算プログラム

## 1. このプログラムでできること

N=126・125 近傍の核種について、β⁻ 崩壊の半減期を計算する。  
計算手順は Phys. Rev. C **109**, 064319 (2024) に基づく。

**入力**: 殻模型コード KSHELL が出力する行列要素（B_GT, xi, zeta）  
**出力**: 各核種の計算半減期・実験値との比・log ft（画面 + `results.dat`）

---

## 2. 必要な環境

| ソフトウェア | 確認コマンド | 用途 |
|---|---|---|
| gfortran（GCC 9以降） | `gfortran --version` | Fortran コンパイル |
| bash | `bash --version` | ビルドスクリプト実行 |
| Python 3（任意） | `python3 --version` | 結果の描画 |
| matplotlib, numpy（任意） | `pip install matplotlib numpy` | 描画ライブラリ |

---

## 3. ファイル構成

```
HF_beta/
├── MANUAL.md                ← このファイル
├── README.md                ← プロジェクト概要
├── plot_results.py          ← 描画スクリプト
├── results.pdf / results.png← 描画結果（実行後に生成）
├── fortran/
│   ├── build.sh             ← ビルドスクリプト（Makefileなしで使用可）
│   ├── Makefile             ← make が使える環境向け
│   ├── beta_decay_calc      ← 計算バイナリ（ビルド後に生成）
│   ├── results.dat          ← 計算結果（実行後に生成）
│   ├── src/                 ← Fortran ソースコード
│   │   ├── mod_constants.f90       物理定数
│   │   ├── mod_wave_function.f90   波動関数・OBDM
│   │   ├── mod_matrix_elements.f90 行列要素
│   │   ├── mod_quenching.f90       クエンチング
│   │   ├── mod_shape_factor.f90    シェイプファクター C(W)
│   │   ├── mod_phase_space.f90     位相空間積分 f
│   │   ├── mod_half_life.f90       半減期計算
│   │   └── beta_decay_main.f90     メインプログラム
│   ├── tests/
│   │   └── test_runner.f90  単体テスト
│   └── obj/                 コンパイル中間ファイル（自動生成）
└── PhysRevC.109.064319.pdf  参照論文
```

---

## 4. ビルド方法

### 4-1. `build.sh` を使う（Makefileなし・推奨）

```bash
cd fortran
bash build.sh
```

実行すると以下が順番に行われる：

1. `obj/` ディレクトリを作成
2. モジュールを依存順にコンパイル
3. 単体テスト（`test_runner`）をビルドして自動実行
4. 計算バイナリ（`beta_decay_calc`）をビルド

全テストが `PASS` になれば成功。

### 4-2. `make` を使う（make がインストールされている場合）

```bash
cd fortran
make all
```

| コマンド | 動作 |
|---|---|
| `make all` | テスト＋計算バイナリを一括ビルド |
| `make test` | テストのみビルド・実行 |
| `make calc` | 計算バイナリのみビルド |
| `make clean` | 生成ファイルをすべて削除してクリーンな状態に戻す |

### 4-3. コンパイルエラーが出たとき

`obj/` が壊れている可能性がある。一度クリーンしてからやり直す：

```bash
rm -rf obj beta_decay_calc test_runner
bash build.sh
```

---

## 5. 計算の実行

```bash
cd fortran
./beta_decay_calc
```

画面に結果が表示され、同時に `fortran/results.dat` が生成される。

---

## 6. 出力の読み方

### 画面出力

```
Nucleus         Z    A   Q[MeV]   T_calc[s]    T_exp[s]  Calc/Exp   logft
207Tl>207Pb    81  207    1.427     80467.0        286.2   281.157    6.85
```

| 列 | 内容 |
|---|---|
| Nucleus | 親核 > 娘核 |
| Z | 親核の陽子数 |
| A | 質量数 |
| Q[MeV] | β崩壊 Q 値（AME2020）|
| T_calc[s] | 計算半減期 [秒] |
| T_exp[s] | 実験半減期 [秒]（不明の場合 `-1`）|
| Calc/Exp | T_calc / T_exp の比（1に近いほど一致）|
| logft | log₁₀(f × T₁/₂) |

### results.dat

スペース区切りのテキストファイル。`#` で始まる行はコメント（ヘッダー）。  
Python や gnuplot で読み込んで描画できる。

---

## 7. 入力データの変更

`src/beta_decay_main.f90` の核種データ部分を編集する。

```fortran
! 207Tl (Z=81, N=126, A=207) → 207Pb
dat(1) = NuclideDatum('207Tl>207Pb ', 81, 207, 82, &
                       1.427d0, 0.0010d0, 0.0d0, 0.0d0, 286.2d0)
!                      Q[MeV]  B_GT      xi     zeta   T_exp[s]
```

| フィールド | 説明 | 出典 |
|---|---|---|
| Q_MeV | β崩壊 Q 値 [MeV] | AME2020 |
| B_GT | Gamow-Teller 遷移強度 | KSHELL の OBDM 出力 |
| xi | 一次禁止遷移の行列要素 | KSHELL の OBDM 出力 |
| zeta | 一次禁止遷移の行列要素 | KSHELL の OBDM 出力 |
| T_exp_s | 実験半減期 [秒]（不明なら `-1.0d0`）| NUBASE2020 |

編集後は再ビルドが必要：

```bash
bash build.sh
./beta_decay_calc
```

> **注意**: 現在の B_GT, xi, zeta はプレースホルダー値。  
> KSHELL の実際の計算結果に置き換えること。

---

## 8. 結果の描画

計算を実行して `results.dat` が生成された後：

```bash
cd HF_beta          # fortran/ の一つ上のディレクトリ
python3 plot_results.py
```

`results.pdf` と `results.png` が生成される。

| グラフ | 内容 |
|---|---|
| 左 | 各核種の計算値 vs 実験値の棒グラフ（対数スケール）|
| 右 | 計算値 vs 実験値の散布図（対角線に乗るほど計算が実験と一致）|

---

## 9. モジュール一覧

| ファイル | 役割 | 主な定数・関数 |
|---|---|---|
| `mod_constants.f90` | 物理定数の定義 | `D_CONST=6289 s`, `G_A=1.2756`, `ME_MEV`, `ALPHA` |
| `mod_wave_function.f90` | 波動関数・OBDM の格納 | `NuclearState` 型, OBDM の set/get |
| `mod_matrix_elements.f90` | 単粒子行列要素・B(GT) | `single_particle_me()`, `compute_bgt()` |
| `mod_quenching.f90` | クエンチング適用 | `effective_ga(q)`, `Q_STANDARD=0.74` |
| `mod_shape_factor.f90` | シェイプファクター C(W) | `make_allowed_sf()`, `make_ff_sf()` |
| `mod_phase_space.f90` | 位相空間積分 f（シンプソン法 200点）| `phase_space_integral()`, `fermi_function()` |
| `mod_half_life.f90` | 半減期・log ft の計算 | `Transition` 型, `total_half_life()` |

---

## 10. よくあるエラーと対処法

| エラーメッセージ | 原因 | 対処 |
|---|---|---|
| `make: command not found` | make が未インストール | `bash build.sh` を使う |
| `No rule to make target 'all'` | `fortran/` 以外のディレクトリで実行した | `cd fortran` してから実行 |
| `No such file or directory: obj/...` | `obj/` が存在しない | `bash build.sh` が自動作成するので `build.sh` から実行 |
| `results.dat` が存在しない | `./beta_decay_calc` を実行していない | 計算を先に実行する |
| Python 描画でエラー | matplotlib/numpy が未インストール | `pip install matplotlib numpy` |

---

## 11. beta_decay_main.f90 コード解説

### 11-1. 全体の流れ

```
① 使用モジュールの宣言
② 核種データ構造（NuclideDatum型）の定義
③ 8核種分のデータを配列に設定
④ クエンチング係数を計算
⑤ 核種ごとにループして半減期を計算
⑥ 画面と results.dat に結果を出力
```

---

### 11-2. 使用モジュールの宣言（1〜19行目）

```fortran
use mod_constants,    only: D_CONST
use mod_quenching,    only: effective_ga, Q_STANDARD
use mod_shape_factor, only: ShapeFactor, make_allowed_sf, make_ff_sf
use mod_phase_space,  only: phase_space_integral
use mod_half_life,    only: partial_half_life_gt, partial_half_life_ff, &
                             Transition, total_half_life
```

各モジュールから必要なものだけを `only:` で指定して取り込んでいる。

| 取り込んでいるもの | 内容 |
|---|---|
| `D_CONST` | 弱崩壊定数 D = 6289 s |
| `effective_ga` | g_A^eff = q × g_A を返す関数 |
| `Q_STANDARD` | 標準クエンチング因子 q = 0.74 |
| `ShapeFactor` | シェイプファクターを格納する型 |
| `make_allowed_sf` | 許容遷移（GT）のシェイプファクターを作る関数 |
| `make_ff_sf` | 一次禁止遷移（FF）のシェイプファクターを作る関数 |
| `phase_space_integral` | 位相空間積分 f を数値計算する関数 |
| `Transition` | 1つの遷移を表す型 |
| `total_half_life` | 複数遷移を合算して全半減期を返す関数 |

---

### 11-3. 核種データ構造（25〜35行目）

```fortran
type :: NuclideDatum
  character(len=12) :: label    ! 核種名（例: "207Tl>207Pb"）
  integer  :: Z_par             ! 親核の陽子数
  integer  :: A                 ! 質量数
  integer  :: Z_dau             ! 娘核の陽子数
  real(8)  :: Q_MeV            ! β崩壊 Q 値 [MeV]
  real(8)  :: B_GT             ! GT 遷移強度（クエンチング前）
  real(8)  :: xi               ! FF 行列要素 xi（クエンチング前）
  real(8)  :: zeta             ! FF 行列要素 zeta（クエンチング前）
  real(8)  :: T_exp_s          ! 実験半減期 [秒]（不明なら -1）
end type NuclideDatum
```

1つの核種に必要な全データをまとめた「型」。  
`dat(1)` 〜 `dat(8)` の配列として8核種分を保持する。

---

### 11-4. 核種データの設定（65〜103行目）

```fortran
dat(1) = NuclideDatum('207Tl>207Pb ', 81, 207, 82, &
                       1.427d0, 0.0010d0, 0.0d0, 0.0d0, 286.2d0)
```

引数の順番は型定義と同じ：  
`label, Z_par, A, Z_dau, Q_MeV, B_GT, xi, zeta, T_exp_s`

| 核種 | 中性子数 | 遷移種別 |
|---|---|---|
| dat(1)〜dat(4) | N=126 | GT のみ（xi=zeta=0） |
| dat(5)〜dat(6) | N=125 | GT のみ（xi=zeta=0） |
| dat(7)〜dat(8) | N=127 | GT + FF（論文との比較用）|

---

### 11-5. クエンチングの適用（108〜109行目）

```fortran
q_factor = Q_STANDARD      ! 0.74
ga_eff   = effective_ga(q_factor)   ! = 0.74 × 1.2756 = 0.9440
```

β崩壊の行列要素は核内媒質中でそのまま使えず、クエンチング因子 q を乗じて補正する。  
この補正後の値 `ga_eff` を以降の計算で使う。

---

### 11-6. 核種ごとの計算ループ（135〜222行目）

各核種について以下を順番に実行する。

#### GT 遷移の処理（140〜154行目）

```fortran
if (dat(i)%B_GT > 0.0d0) then
  sf_unit = make_allowed_sf(1.0d0)
  f0_gt   = phase_space_integral(dat(i)%Z_dau, dat(i)%Q_MeV, sf_unit)
  ...
end if
```

B_GT > 0 のときのみ実行。  
`make_allowed_sf(1.0d0)` で C(W)=1 のシェイプファクターを作り、  
`phase_space_integral` で位相空間積分 f を計算する。  
（B_GT は `total_half_life` 内で掛け合わせる）

#### FF 遷移の処理（157〜170行目）

```fortran
if (abs(dat(i)%xi) > 0.0d0 .or. abs(dat(i)%zeta) > 0.0d0) then
  sf_ff = make_ff_sf(dat(i)%xi, dat(i)%zeta)
  f0_ff = phase_space_integral(dat(i)%Z_dau, dat(i)%Q_MeV, sf_ff)
  ...
end if
```

xi または zeta が 0 でないときのみ実行（N=127 核種のみ該当）。  
FF 遷移のシェイプファクター C(W) = (ξ − ζW)² + ζ²(W²−1)/3 を使って f を計算する。

#### 全半減期の計算（173〜177行目）

```fortran
t_tot_s = total_half_life(tr(1:n_tr), n_tr)
```

GT 遷移・FF 遷移の崩壊率を足し合わせて全半減期を求める：

```
1/T_total = 1/T_GT + 1/T_FF
```

#### log ft の計算（181〜185行目）

```fortran
lft_gt = log10(D_CONST / (ga_eff**2 * dat(i)%B_GT))
```

log ft は遷移の「禁止度」を表す指標。  
許容 GT 遷移では log ft ≈ 4〜6 程度になる。

---

### 11-7. 出力処理（194〜220行目）

実験値の有無・計算値の有効性に応じて3パターンで分岐して出力する。

| 条件 | 出力内容 |
|---|---|
| 計算値あり・実験値あり | T_calc, T_exp, Calc/Exp, logft を全列表示 |
| 計算値あり・実験値なし | T_exp 列に `unknown` と表示 |
| 計算値なし（B_GT=0など）| T_calc 列に `(invalid)` と表示 |

同じ内容を画面（`write(*,...)`) と `results.dat`（`write(FUNIT,...)`) の両方に書き出す。  
`huge(1.0d0)` は Fortran の「実質無限大」で、計算が成立しない場合のフラグとして使っている。
