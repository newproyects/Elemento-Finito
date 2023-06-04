disp('Escoja un archivo .mat')

load_data;

draw_elems = [1 1 1 1 1 1];
scale_factor = 1;
shrink_factor = 0.005;
draw_elem_nums = 1;
draw_nodes = 1;
draw_node_nums = 1;
draw_bcs = 1;
draw_bc_vals = 1;
nod_disps = [];
view_max_st = 0;
draw_def_elems = [1 1 1 1 1 1];
draw_structure;

disp('Pulse cualquier tecla para continuar')
pause

for i = 1:1


[ nod_disps, nodal_int_forces, elem_props] = perifea2d_linear( elem_con, elem_props, coords, bcs, bc_vals);


end


draw_structure;

