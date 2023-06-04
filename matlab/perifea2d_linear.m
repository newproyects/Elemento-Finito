function [ nodal_displacements, nodal_int_forces, elem_props] = perifea2d_linear( elem_con, elem_props, coords, bcs, bc_vals)
            

% This function determines the link or the membrane with maximum strain or
% stretch in a linear elastic, peridynamic, 2D problem
%
% Inputs:
%        elem_con: a cell array of element nodal connectivities, for example: {[1 3] [1 2 3] [2 3 4 5]};
%
%        elem_props: a cell array of element  properties: {[1 100 1] [2 29000 0.3 1 100 200 150 30] [1 100 1.]};
%               for truss: [elem_type_code=1, E, A]
%               for membrane: [elem_type_code=3, num_nodes, E, nu, t, plane_code]
%        E: Young's modulus
%        nu: Poisson's Ratio
%        A: cross-sectional area of the element
%        t: element thickness
%        Rotation_Angle: Angle of rotation about local x-axis, degrees
%        plane_code: 1: plane stress; 2: plane stain; 3: axisymmetric
%
%       coords: a (num_nodes rows x 2 columns) array of element global (X, Y) coordinates
%
%       bcs: [num_nodes x 3] array of fixity codes (x, y, theta)
%           0: free (applied load/moment)
%           1: fixed (applied displacement/rotation)
%
%       bc_vals: [num_nodes x 3] array of dof values (x, y, theta)
%            (either a specified load/moment or a specified displacement/rotation, 
%             depending on the corresponding BC code.) 
%
% Outputs:
%        nodal_displacements: a (num_nodes x 3) array containing nodal displacements and rotations (Dx Dy Rz)
%        external_forces: a (num_nodes x 3) array containing the external nodal forces and moments (Fx Fy Mz)
tic

ilink = 0;
elems_load_factor = 0;
worst_elem = 0;
num_nodes = size(coords, 1);
num_elems = size(elem_con, 2);
free_dof = zeros(num_nodes, 3);
num_free_dof = 0;

%
% Number the free degrees of freedom and fill the forcing vector
%
for node = 1 : num_nodes
    for node_dof = 1 : 3
        if(bcs(node, node_dof) == 0)
            num_free_dof = num_free_dof + 1;
            free_dof(node, node_dof) = num_free_dof;
            F_free(num_free_dof, 1) = bc_vals(node, node_dof);
        end
    end
end

%
% Number the fixed degrees of freedom and fill the forcing vector
%
fixed_dof = zeros(num_nodes, 3);
num_fixed_dof = 0;
for node = 1 : num_nodes
    for node_dof = 1 : 3
        if(bcs(node, node_dof) ~= 0)
            num_fixed_dof = num_fixed_dof + 1;
            fixed_dof(node, node_dof) = num_fixed_dof + num_free_dof;
            D_fixed(num_fixed_dof, 1) = bc_vals(node, node_dof);
        end
    end
end
%
% Assemble the global stiffness matrix
%
K_ff = zeros(num_free_dof, num_free_dof);
K_fs = zeros(num_free_dof, num_fixed_dof);
for ielem = 1 : num_elems
    ielem_type = elem_props{ielem}(1);
    ielem_num_nodes = elem_props{ielem}(2);
    ielem_coords(1:ielem_num_nodes, 1:2) = coords(elem_con{ielem}(1:ielem_num_nodes), 1:2);
    ielem_coords = ielem_coords';
    ielem_coords = ielem_coords(1:ielem_num_nodes*2)';
    if (ielem_num_nodes == 1)
        ielem_coords = coords(elem_con{ielem}(ielem_num_nodes),1:2)';
    end
    
    for jelem = 1 : num_elems
        jelem_type = elem_props{jelem}(1);
        jelem_num_nodes = elem_props{jelem}(2);
        jelem_coords(1:jelem_num_nodes, 1:2) = coords(elem_con{jelem}(1:jelem_num_nodes), 1:2);
        jelem_coords = jelem_coords';
        jelem_coords = jelem_coords(1:jelem_num_nodes*2)';
        if (jelem_num_nodes == 1)
            jelem_coords = coords(elem_con{jelem}(jelem_num_nodes),1:2)';
        end
        
        
        if(ielem ~= jelem)
            
            %if  ((ielem_type == 6)&(jelem_type == 6))    
           %[k_ijelem, link] = k_ij_2D_link(ielem_coords, jelem_coords,...
           %    elem_props{ielem}, elem_props{jelem}, link, ielem, jelem);
           % else
           [k_ijelem] = k_ij_2D(ielem_coords, jelem_coords, elem_props{ielem}, elem_props{jelem});             
           % end
             
            for el_row_node = 1 : ielem_num_nodes + jelem_num_nodes
                if(el_row_node <= ielem_num_nodes)
                    elem_row = ielem;
                    el_row_nod = el_row_node;
                else
                    elem_row = jelem;
                    el_row_nod = el_row_node - ielem_num_nodes;
                end
                for el_row_node_dof = 1 : 3
                    el_row = 3*(el_row_node - 1) + el_row_node_dof;
                    glob_row = free_dof(elem_con{elem_row}(el_row_nod), el_row_node_dof);
                    if(glob_row ~= 0)
                        for el_col_node = 1 : ielem_num_nodes + jelem_num_nodes
                            if(el_col_node <= ielem_num_nodes)
                                elem_col = ielem;
                                el_col_nod = el_col_node;
                            else
                                elem_col = jelem;
                                el_col_nod = el_col_node - ielem_num_nodes;
                            end
                            for el_col_node_dof = 1 : 3
                                el_col = 3*(el_col_node - 1) + el_col_node_dof;
                                glob_col = free_dof(elem_con{elem_col}(el_col_nod), el_col_node_dof);
                                if(glob_col ~= 0)
                                    K_ff(glob_row, glob_col ) = K_ff(glob_row, glob_col) ...
                                        + k_ijelem(el_row, el_col);
                                   
                                    
                                else
                                    glob_col = fixed_dof(elem_con{elem_col}(el_col_nod), el_col_node_dof);
                                    K_fs(glob_row, glob_col - num_free_dof ) = K_fs(glob_row, glob_col - num_free_dof) ...
                                        + k_ijelem(el_row, el_col);
                                end
                            end
                        end
                    end
                end
            end
           
          
        else
            k_elem = k_ii_2D(ielem_coords, elem_props{ielem});
                      
            for el_row_node = 1 : ielem_num_nodes
                for el_row_node_dof = 1 : 3
                    el_row = 3*(el_row_node - 1) + el_row_node_dof;
                    glob_row = free_dof(elem_con{ielem}(el_row_node), el_row_node_dof);
                    if(glob_row ~= 0)
                        for el_col_node = 1 : ielem_num_nodes
                            for el_col_node_dof = 1 : 3
                                el_col = 3*(el_col_node - 1) + el_col_node_dof;
                                glob_col = free_dof(elem_con{ielem}(el_col_node), el_col_node_dof);
                                if(glob_col ~= 0)
                                    K_ff(glob_row, glob_col ) = K_ff(glob_row, glob_col) ...
                                        + k_elem(el_row, el_col);
                                    
                                else
                                    glob_col = fixed_dof(elem_con{ielem}(el_col_node), el_col_node_dof);
                                    K_fs(glob_row, glob_col - num_free_dof ) = K_fs(glob_row, glob_col - num_free_dof) ...
                                        + k_elem(el_row, el_col);
                                end
                            end
                        end
                    end
                end
            end           
        end        
    end
end 

num_free_dof;
if(sprank(K_ff) == 0)
    disp_free = [];
elseif(sprank(K_ff) == num_free_dof)
    disp_free = K_ff\(F_free - K_fs*D_fixed);
else
    'the structure is unstable'
    nodal_displacements = [];
    external_forces = [];
    return
end

Qu = K_fs'*disp_free;
Qk = K_ff*disp_free;

nodal_disp2D = zeros(num_nodes, 3);
for node = 1 : num_nodes
    for node_dof = 1 : 3
        if(bcs(node, node_dof) == 1)
            nodal_disp2D(node, node_dof) = bc_vals(node, node_dof);
        else
            nodal_disp2D(node, node_dof) = disp_free(free_dof(node, node_dof));
        end
    end
end
nodal_displacements = zeros(num_nodes, 6);
nodal_displacements(:, 1:2) = nodal_disp2D(:, 1:2);
nodal_displacements(:, 6) = nodal_disp2D(:, 3);

external_forces = zeros(num_nodes, 3);

nodal_ext_forces2D = zeros(num_nodes, 3);
for node = 1 : num_nodes
    for node_dof = 1 : 3
        if(bcs(node, node_dof) == 1)
            nodal_ext_forces2D(node, node_dof) = bc_vals(node, node_dof);
        else
            nodal_ext_forces2D(node, node_dof) = F_free(free_dof(node, node_dof));
        end
    end
end
nodal_ext_forces = zeros(num_nodes, 6);
nodal_ext_forces(:, 1:2) = nodal_ext_forces2D(:, 1:2);
nodal_ext_forces(:, 6) = nodal_ext_forces2D(:, 3);

nodal_int_forces2D = zeros(num_nodes, 3);
for node = 1 : num_nodes
    for node_dof = 1 : 3
        if(bcs(node, node_dof) == 0)
            nodal_int_forces2D(node, node_dof) = bc_vals(node, node_dof);
        else
            nodal_int_forces2D(node, node_dof) = Qu(fixed_dof(node, node_dof) - num_free_dof);
        end
    end
end
nodal_int_forces = zeros(num_nodes, 6);
nodal_int_forces(:, 1:2) = nodal_int_forces2D(:, 1:2);
nodal_int_forces(:, 6) = nodal_int_forces2D(:, 3);


toc

max_damage_factor = 0;

for elem = 1:num_elems
    elem_type = elem_props{elem}(1);
    elem_num_nodes = size(elem_con{elem}, 2);
    elem_coords(1:elem_num_nodes, 1:2) = coords(elem_con{elem}(1:elem_num_nodes), 1:2);    
    elem_bcs(1:elem_num_nodes, 1:6) = bcs(elem_con{elem}(1:elem_num_nodes), 1:6);
    elem_disps = nodal_displacements(elem_con{elem}(1:elem_num_nodes), 1:3);
    elem_forces = nodal_int_forces(elem_con{elem}(1:elem_num_nodes),1:3);
    strain = stretch(elem_coords, elem_props{elem}, elem_disps);
    
    strain_mem(elem) = strain;

    if (elem_type == 3)
        strain_crit = elem_props{elem}(9);

        damage_factor = strain/strain_crit;
        

        if(damage_factor >= max_damage_factor)
            max_damage_factor = damage_factor;
            worst_elem = elem;
        end
    end
end

elems_load_factor = 1/max_damage_factor;

strain_energy = 0;
px = 0;
py = 0;
for inode = 1:num_nodes
    strain_energy = strain_energy + (1/2)*nodal_displacements(inode, :)*nodal_int_forces(inode, :)';
    %nodal_displacements(inode, :)
    %nodal_int_forces(inode, :)
    py = py + nodal_ext_forces(inode, 2);
    px = px + nodal_ext_forces(inode, 1);
end

strain_energy;



nom_disp = strain_energy/py;



return
