if (!exists("spectrum_data")) spectrum_data = "results/phc30_np12_v1_100_cdw/data/spectrum.dat"
if (!exists("structure_data")) structure_data = "results/phc30_np12_v1_100_cdw/data/structure_factor_ground_states.dat"
if (!exists("figure_dir")) figure_dir = "reports/fqahc_results/figures/optical_response"

set datafile commentschars "#"
set term pngcairo size 1800,720 enhanced font "sans,14"
set output figure_dir . "/phc30_manybody_spectrum.png"
set multiplot layout 1,2 title "Tilted 30-site ED, strong-coupling period-three CDW (V1=100)"
set border linewidth 1.2
set grid ytics lc rgb "#d8d8d8" lw 1
set key opaque box top left
set xlabel "momentum sector k"
set ylabel "E-E0"
set xrange [-0.5:14.5]
set xtics 0,1,14
set yrange [-0.04:2.70]
set title "eight saved levels per sector"
plot spectrum_data using (($2 <= 0.0001245630) ? $1 : 1/0):2 \
        with points pt 7 ps 1.35 lc rgb "#c43b3b" \
        title "CDW ground-state triplet", \
     spectrum_data using (($2 > 0.0001245630) ? $1 : 1/0):2 \
        with points pt 7 ps 0.75 lc rgb "#55758f" \
        title "saved excitations"

set ylabel ""
set yrange [-0.002:0.050]
set title "low-energy zoom"
set key top right
set arrow 1 from graph 0, first 0.0001245629 to graph 1, first 0.0001245629 \
    nohead dt 2 lw 1.5 lc rgb "#c43b3b"
set arrow 2 from graph 0, first 0.0213937117 to graph 1, first 0.0213937117 \
    nohead dt 3 lw 1.5 lc rgb "#2b6ca3"
set label 1 "triplet width = 0.0001245629" at graph 0.04,0.92 front
set label 2 "gap to fourth = 0.0212691488" at graph 0.04,0.83 front
plot spectrum_data using (($2 <= 0.0001245630) ? $1 : 1/0):2 \
        with points pt 7 ps 1.35 lc rgb "#c43b3b" \
        title "k=0,5,10", \
     spectrum_data using (($2 > 0.0001245630) ? $1 : 1/0):2 \
        with points pt 7 ps 0.75 lc rgb "#55758f" \
        title "low-energy excitations"
unset multiplot
unset arrow 1
unset arrow 2
unset label 1
unset label 2
unset output

set term pngcairo size 1260,780 enhanced font "sans,15"
set output figure_dir . "/phc30_cdw_structure_factor.png"
set title "30-site CDW ground-state static structure factor"
set xlabel "momentum q"
set ylabel "N(q), unit-cell normalized"
set xrange [-0.3:14.3]
set xtics 0,1,14
set yrange [0:*]
set grid ytics lc rgb "#d8d8d8" lw 1
set key opaque box top right
plot structure_data using (strcol(1) eq "0" ? $2 : 1/0):3 \
        with linespoints pt 5 ps 0.9 lw 1.2 lc rgb "#2b6ca3" title "ground state k=0", \
     structure_data using (strcol(1) eq "5" ? $2 : 1/0):3 \
        with linespoints pt 7 ps 0.9 lw 1.2 lc rgb "#c43b3b" title "ground state k=5", \
     structure_data using (strcol(1) eq "10" ? $2 : 1/0):3 \
        with linespoints pt 9 ps 0.9 lw 1.2 lc rgb "#388e5c" title "ground state k=10", \
     structure_data using (strcol(1) eq "average" ? $2 : 1/0):3 \
        with linespoints pt 6 ps 1.15 lw 2.5 lc rgb "#202020" title "three-state average"
unset output
