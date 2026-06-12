#!/bin/bash
# notes/ 内の .tex をコンパイルする
# PDF は notes/ に、補助ファイルは build/ に出力される

set -e
LATEXMK="/mnt/c/texlive/2025/bin/windows/latexmk.exe"
NOTES_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$NOTES_DIR"

mkdir -p build

if [ "$#" -eq 0 ]; then
    # 引数なし → notes/ の全 .tex をコンパイル
    for tex in *.tex; do
        echo "=== Compiling $tex ==="
        "$LATEXMK" "$tex"
    done
else
    # 引数あり → 指定したファイルのみ
    for tex in "$@"; do
        echo "=== Compiling $tex ==="
        "$LATEXMK" "$tex"
    done
fi

echo ""
echo "完了: PDF は notes/ に、補助ファイルは build/ に出力されました"
