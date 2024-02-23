#set terminal type 'qt'
set terminal pngcairo
set output 'Position.png'
set xlabel 'X'
set ylabel 'Y'
set title "Kepler's 2 Body Problem"

pl './Position-A.txt' u 2:3 w lp title "Paticle A",'./Position-B.txt' u 2:3 w lp title "Paticle B"
