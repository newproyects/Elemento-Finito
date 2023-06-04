set macro


nodes = "u 1:2 t 'Nodos' w p pt 7 ps 2 lc rgb 'blue'"

elems = "u 1:2 t 'Elementos' w l lc rgb 'gray'"

bcs_tras = "u 1:2 t 'Fijeza Traslacional' w l lw 6 lc rgb 'black'"

bcs_rot = "u 1:2 t 'Fijeza Rotacional' w l lw 3 lc rgb 'red'"

bc_vals_tras = "u 1:2 t 'Fuerzas' w l lc rgb 'black'"

bc_vals_rot = "u 1:2 t 'Momentos' w l lc rgb 'red'"

elems_deform = "u 1:2 t 'Elementos deformados' w l lc rgb 'black'"

# Interruptor de la leyenda
set key off

plot "draw_structure.dat" index 0 @nodes, "" index 1 @elems, "" index 2 @bcs_tras, "" index 3 @bcs_rot, "" index 4 @bc_vals_tras, "" index 5 @bc_vals_rot, "" index 6 @elems_deform
