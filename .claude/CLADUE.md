# CLAUDE.md — AI作業ルール

## プロジェクト概要
LaTeXノートを管理・PDF化するリポジトリ。

---

## ビルドコマンド

```bash
# 特定ノートをPDF化
cd notes && pdflatex <filename>.tex && mv <filename>.pdf ../pdf/

# 全ノートを一括ビルド
make all

# クリーンビルド（補助ファイル削除）
make clean && make all
```

---

## LaTeXルール
- エンコーディング: UTF-8（`\usepackage[utf8]{inputenc}`）
- 日本語使用時: `\usepackage{luatexja}` + LuaLaTeXでコンパイル
- 図は必ず `figures/` に置き、相対パスで参照
- 参考文献は BibTeX（`references.bib`）で管理

---

## AI作業後の必須手順

AIが何らかの作業（ファイル編集・作成・ビルドなど）を完了したら、**必ず以下の順序で実行すること**。

### 1. PROGRESS.md を更新する

`PROGRESS.md` を開き、以下を記録する：

```
## [YYYY-MM-DD] 作業タイトル

### ✅ やったこと
- （具体的に箇条書き）

### 🔜 次にやること
- （具体的に箇条書き）

### ⚠️ 注意・メモ
- （あれば）
```

### 2. git commit & push する

```bash
git add -A
git commit -m "chore: <作業内容を一行で要約>"
git push origin main
```

コミットメッセージの例：
- `feat: linear-algebra ノート追加`
- `fix: figures パスの修正`
- `chore: PROGRESS.md 更新`

---

## ディレクトリ構成

```
.
├── CLAUDE.md        # このファイル（AIへの指示）
├── PROGRESS.md      # 作業ログ（何をしたか・次に何をするか）
├── Makefile
├── notes/           # .tex ソースファイル
├── figures/         # 図・画像
├── pdf/             # 生成済みPDF
└── references.bib   # 参考文献
```

---

## コミットメッセージ規則

| プレフィックス | 用途 |
|---|---|
| `feat:` | 新しいノート・機能の追加 |
| `fix:` | バグ・誤りの修正 |
| `chore:` | ログ更新・設定変更など |
| `build:` | ビルド設定・Makefile変更 |
| `docs:` | ドキュメント更新 |
