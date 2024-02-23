#set terminal type 'qt'
set xlabel 'nmc'
set ylabel 'Total Potential Energy'
pl './LJa95.txt' w l title'a=0.95', './LJa112.txt' w l title 'a=1.12'
