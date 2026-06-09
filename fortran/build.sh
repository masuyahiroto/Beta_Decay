#!/bin/bash
# Makefile不要でFortranをビルドするスクリプト
set -e

FC=gfortran
FLAGS="-O2 -Wall -Wextra"
SRC=src
OBJ=obj

mkdir -p "$OBJ"

compile() {
    echo "  compiling $1.f90"
    $FC $FLAGS -c "$SRC/$1.f90" -o "$OBJ/$1.o" -J"$OBJ" -I"$OBJ"
}

# モジュールを依存順にコンパイル
compile mod_constants
compile mod_wave_function
compile mod_matrix_elements
compile mod_quenching
compile mod_shape_factor
compile mod_phase_space
compile mod_half_life

MODS="$OBJ/mod_constants.o $OBJ/mod_wave_function.o $OBJ/mod_matrix_elements.o \
      $OBJ/mod_quenching.o $OBJ/mod_shape_factor.o $OBJ/mod_phase_space.o \
      $OBJ/mod_half_life.o"

# テストバイナリ
echo "  compiling tests/test_runner.f90"
$FC $FLAGS -c tests/test_runner.f90 -o "$OBJ/test_runner.o" -I"$OBJ" -J"$OBJ"
$FC $FLAGS -o test_runner $MODS "$OBJ/test_runner.o"
echo "--- running tests ---"
./test_runner

# 計算バイナリ
compile beta_decay_main
$FC $FLAGS -o beta_decay_calc $MODS "$OBJ/beta_decay_main.o"

echo ""
echo "Build done."
echo "  Run:  ./beta_decay_calc"
echo "  Plot: cd .. && python3 plot_results.py"
