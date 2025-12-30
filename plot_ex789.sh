#!/bin/bash
# ==========================================
# Plot EX7 / EX8 / EX9 residual curves only
# Usage:
#   ./plot_ex789.sh bin
# ==========================================

if [ $# -lt 1 ]; then
  echo "Usage: $0 <data_dir>"
  exit 1
fi

DATADIR=$1
OUTDIR=plots_ex789

mkdir -p "$OUTDIR" || {
  echo "Error: cannot create $OUTDIR (permission denied?)"
  exit 1
}

for i in 0 1 2; do
  if [ ! -f "$DATADIR/RESVEC_$i.dat" ]; then
    echo "Missing $DATADIR/RESVEC_$i.dat"
    exit 1
  fi
done

gnuplot << EOF
set term pngcairo size 900,600 noenhanced
set grid
set logscale y
set xlabel "Iteration"
set ylabel "Residual norm"

# EX7
set output "$OUTDIR/ex7_richardson.png"
set title "EX7 - Richardson"
plot "$DATADIR/RESVEC_0.dat" using 0:1 with lines lw 2 title "Richardson"

# EX8
set output "$OUTDIR/ex8_jacobi.png"
set title "EX8 - Richardson + Jacobi"
plot "$DATADIR/RESVEC_1.dat" using 0:1 with lines lw 2 title "Jacobi"

# EX9
set output "$OUTDIR/ex9_gs.png"
set title "EX9 - Richardson + Gauss-Seidel"
plot "$DATADIR/RESVEC_2.dat" using 0:1 with lines lw 2 title "Gauss-Seidel"

# Comparison
set output "$OUTDIR/ex789_compare.png"
set title "Comparison EX7 / EX8 / EX9"
plot \
  "$DATADIR/RESVEC_0.dat" using 0:1 with lines lw 2 title "Richardson", \
  "$DATADIR/RESVEC_1.dat" using 0:1 with lines lw 2 title "Jacobi", \
  "$DATADIR/RESVEC_2.dat" using 0:1 with lines lw 2 title "Gauss-Seidel"

EOF

echo "Plots generated in $OUTDIR/"
