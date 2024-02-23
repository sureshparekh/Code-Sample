#set terminal type 'qt'
set terminal pngcairo
set output 'Velocity.png'
set xlabel 'Vx'
set ylabel "Vy"
set title "Kepler's 2 Body Problem"

pl './Velocity-A.txt' u 2:3 w lp title "Paticle A",'./Velocity-B.txt' u 2:3 w lp title "Paticle B"
