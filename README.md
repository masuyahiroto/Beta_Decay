# Beta Decay — N=126, 125 核種の β⁻ 崩壊半減期計算

Phys. Rev. C 109, 064319 (2024) の計算手順に基づき、N=126・125 近傍核種の
β⁻ 崩壊半減期を Fortran 90 で実装したコード群。

---

## 物理的背景

β⁻ 崩壊の半減期は以下の式で与えられる。

```
t_{1/2} = D / (g_A^eff² × f × B_GT)
```

| 記号 | 意味 |
|------|------|
| `D = 6289 s` | 弱崩壊定数 |
| `g_A^eff = q × g_A` | 有効軸性結合定数（クエンチング適用後）|
| `q ≈ 0.74` | クエンチング因子 |
| `f` | 位相空間積分 |
| `B_GT` | Gamow-Teller 遷移強度（シェル模型 KSHELL から取得）|

複数の遷移がある場合は崩壊率を足し合わせる。

```
1/t_{1/2} = Σ_i  1/t_i
```

### 対象核種

- **N=126** : ²⁰⁷Tl, ²⁰⁶Hg, ²⁰⁵Au, ²⁰⁴Pt
- **N=125** : ²⁰⁶Tl, ²⁰⁵Hg
- **N=127** (比較例) : ²⁰⁸Tl, ²⁰⁷Hg（GT + 一次禁止遷移）

---

## ディレクトリ構成

```
HF_beta/
├── fortran/
│   ├── src/                 # Fortran 90 モジュール群
│   │   ├── mod_constants.f90
│   │   ├── mod_quenching.f90
│   │   ├── mod_shape_factor.f90
│   │   ├── mod_phase_space.f90
│   │   ├── mod_half_life.f90
│   │   ├── mod_wave_function.f90
│   │   ├── mod_matrix_elements.f90
│   │   └── beta_decay_main.f90  # メインドライバ
│   ├── tests/
│   │   └── test_runner.f90
│   └── Makefile
├── docs/
│   └── Note.md
├── notes/
│   └── summary.tex
└── PhysRevC.109.064319.pdf  # 参照論文
```

---

## 計算の流れ

```
KSHELL 出力 (B_GT, xi, zeta)
        ↓
クエンチング適用  g_A^eff = 0.74 × g_A
        ↓
シェイプファクター C(W) の構成
  ・許容遷移 (GT) : C(W) = B_GT
  ・一次禁止遷移  : C(W) = (ξ - ζW)² + ζ²(W²-1)/3
        ↓
位相空間積分  f = ∫ C(W) F(Z,W) p W (W₀-W)² dW
        ↓
半減期  t_{1/2} = D / (g_A^eff² × f)
        ↓
実験値 (NUBASE2020) との比較・log ft 出力
```

---

## 使い方

### 必要なもの

- `gfortran`（GCC 9 以降推奨）
- `make`

インストール確認：

```bash
gfortran --version
make --version
```

---

### ディレクトリへ移動（必須）

**`make` はこのディレクトリで実行しないと動かない。** まず移動する：

```bash
cd fortran
```

移動できているか確認（`fortran` と表示されればOK）：

```bash
pwd   # 末尾が .../HF_beta/fortran になっていること
```

---

### ビルド

#### 全ターゲットを一括ビルド（推奨）

```bash
make all
```

内部では以下の順で処理される：

1. `obj/` ディレクトリを作成
2. `src/` 内のモジュール（`mod_*.f90`）を依存順にコンパイルし `obj/*.o` を生成
3. テストバイナリ `test_runner` をビルドしてそのまま実行
4. 計算バイナリ `beta_decay_calc` をビルド

#### 個別ターゲット

| コマンド | 動作 |
|----------|------|
| `make test` | モジュール＋テストをビルドし、`./test_runner` を自動実行 |
| `make calc` | モジュール＋`beta_decay_calc` のみビルド（テストは実行しない）|
| `make clean` | `obj/`・`test_runner`・`beta_decay_calc` を削除してクリーンな状態に戻す |

再ビルドしたいときは `make clean && make all`。

---

### 実行

#### 単体テスト

```bash
./test_runner
```

各モジュールの数値チェックを行い、PASS / FAIL を表示する。

#### 半減期計算

```bash
./beta_decay_calc
```

全対象核種（N=126, 125, 127）の計算結果を標準出力に表示する。

**出力例**

```
==============================================================================
  β⁻ 崩壊半減期計算  (N=126, 125 核種)
  Phys. Rev. C 109, 064319 (2024) の計算手順
  クエンチング q =  0.74   g_A^eff =  0.9440
==============================================================================
Nucleus         Z    A   Q[MeV]   T_calc[s]    T_exp[s]  Calc/Exp   logft
------------------------------------------------------------------------------
207Tl>207Pb    81  207    1.427       ...          286.2     ...      ...
...
```

結果をファイルに保存したい場合はリダイレクトを使う：

```bash
./beta_decay_calc > result.txt
```

---

## 入力データの更新（KSHELL 連携）

`beta_decay_main.f90` 内の `B_GT`, `xi`, `zeta` は現在プレースホルダー値。
KSHELL の OBDM 出力から以下の手順で更新する。

```
B_GT = |M_GT|² / (2J_i + 1)
xi   = Σ_{pn}  <p || r C¹ || n>  × ρ_{pn}       (rank 0,2 成分)
zeta = Σ_{pn}  <p || r C¹ × σ || n>  × ρ_{pn}   (rank 1 成分)
```

---

## 参考文献

- Phys. Rev. C **109**, 064319 (2024) — 本コードの計算手順の根拠
- Behrens & Bühring (1982) — 位相空間積分・シェイプファクターの定義
- Suhonen (2007) *From Nucleons to Nucleus* — β崩壊理論の教科書
- AME2020 — Q 値データ
- NUBASE2020 — 実験半減期データ
