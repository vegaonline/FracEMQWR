set xlab 'b'
set ylab 'm'
set key left
se xra[:12]
se gri
set title "Determination of b for different m with alpha as parameter"
pl  "data/mbValues_0.0.dat" u 2:3 t "alpha = 1.0e-5" w p, "data/mbValues_0.2.dat" u 2:3 t "alpha = 0.2" w p, "data/mbValues_0.4.dat" u 2:3 t "alpha = 0.4" w p, "data/mbValues_0.6.dat" u 2:3 t "alpha = 0.6" w p, "data/mbValues_0.8.dat" u 2:3 t "alpha = 0.8" w p

