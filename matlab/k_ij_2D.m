function [el_stiff_C1] = k_ij_2D(icoords, jcoords, ielem_props, jelem_props)
%
% Compute the contribution to the stiffness matrix of two elements: j wrt
% i, using the peridynamic model
%

num_nodes_i = ielem_props(2);
num_nodes_j = jelem_props(2);

c_i = ielem_props(3);
d_i = ielem_props(4);
t_i = ielem_props(5);
n_i = ielem_props(6);
delta_i = ielem_props(7);
As_i = ielem_props(5);
As_j = jelem_props(5);
if (num_nodes_i == 1)
    area_i = ielem_props(10);   
end


c_j = jelem_props(3);
d_j = jelem_props(4);
t_j = jelem_props(5);
n_j = jelem_props(6);
delta_j = jelem_props(7);

if (num_nodes_j == 1)
    area_j = jelem_props(10);
end


c = (c_i + c_j)/2;
d = (d_i + d_j)/2;
mat_horizon = min([delta_i delta_j]);



if (num_nodes_i == 1) & (num_nodes_j == 1)
    el_stiff_C0 = zeros(6, 6);
elseif (num_nodes_i ~= 1)&(num_nodes_j == 1)
    el_stiff_C0 = zeros(2*num_nodes_i + 3, 2*num_nodes_i + 3);
elseif (num_nodes_i == 1)&(num_nodes_j ~= 1)
    el_stiff_C0 = zeros(2*num_nodes_j + 3, 2*num_nodes_j + 3);
else
    el_stiff_C0 = zeros(2*num_nodes_i+2*num_nodes_j, 2*num_nodes_i+2*num_nodes_j);
end

el_dist = elem_distance(icoords, jcoords);


if(el_dist > mat_horizon)
    el_stiff_C1 = zeros(3*num_nodes_i+3*num_nodes_j, 3*num_nodes_i+3*num_nodes_j);    
   return 
end


for i1 = 1  : n_i
    for i2 = 1 : n_i  
        
        if(num_nodes_i == 3)
            del_i = 1/n_i;
            xi_i1 = del_i*i1 - del_i/2;
            xi_i2 = del_i*i2 - del_i/2;
            if( i2 <= (n_i + 1 - i1))
                xi_i1 = xi_i1 - del_i/6;
                xi_i2 = xi_i2 - del_i/6;
            else
                xi_i1 = (1 - xi_i1) + del_i/6;
                xi_i2 = (1 - xi_i2) + del_i/6;
            end
        elseif(num_nodes_i == 4)
            del_i = 2/n_i;
            xi_i1 = -1 + del_i*i1 - del_i/2;
            xi_i2 = -1 + del_i*i2 - del_i/2;
        elseif(num_nodes_i == 2)
            del_i = 2/n_i;
            xi_i1 = -1 + del_i*i1 - del_i/2;
            xi_i2 = 0;
        elseif(num_nodes_i == 1)
            del_i = 2/n_i;
            xi_i1 = 0;
            xi_i2 = 0;
            
        else    
            return
        end

        [coor_i, N_i, det_J_i] = peri_quants([xi_i1 xi_i2], icoords);
        
        if (num_nodes_i == 1)
            dvol_i = area_i*t_i;
        elseif (num_nodes_i == 2)
            dvol_i = det_J_i*del_i*As_i;
        else
            dvol_i = det_J_i*del_i*del_i*t_i;
        end
        
       
        for j1 = 1  : n_j
            for j2 = 1  : n_j
                if(num_nodes_j == 3)
                    del_j = 1/n_j;
                    xi_j1 = del_j*j1 - del_j/2;
                    xi_j2 = del_j*j2 - del_j/2;
                    if( j2 <= (n_j + 1 - j1))
                        xi_j1 = xi_j1 - del_j/6;
                        xi_j2 = xi_j2 - del_j/6;
                    else
                        xi_j1 = (1 - xi_j1) + del_j/6;
                        xi_j2 = (1 - xi_j2) + del_j/6;
                    end
                elseif(num_nodes_j == 4)
                    del_j = 2/n_j;
                    xi_j1 = -1 + del_j*j1 - del_j/2;
                    xi_j2 = -1 + del_j*j2 - del_j/2;
                elseif(num_nodes_j == 2)
                    del_j = 2/n_j;
                    xi_j1 = -1 + del_j*j1 - del_j/2;
                    xi_j2 = 0;
                elseif(num_nodes_j == 1)
                    del_j = 1;
                    xi_j1 = 0;
                    xi_j2 = 0;                    
                else                    
                    return
                end 
                
                
                [coor_j, N_j, det_J_j] = peri_quants([xi_j1 xi_j2], jcoords);
                
                
                if (num_nodes_j == 1)         
                dvol_j = area_j*t_j;
                elseif (num_nodes_j == 2)
                dvol_j = det_J_j*del_j*As_j;
                else
                dvol_j = det_J_j*del_j*del_j*t_j;          
                end
                  
                                             
                
                del_coor = coor_j - coor_i;
                L = norm(del_coor);                            
                
               
                if((L <= mat_horizon) & (L > 0))                   
                    
                    
                    if ((num_nodes_i == 1)&(num_nodes_j == 1))
                       N_ij = [eye(3, 3) zeros(3, 3);zeros(3, 3) eye(3, 3)];                       
                    elseif((num_nodes_i ~= 1)&(num_nodes_j == 1))
                       N_ij = [N_i zeros(3, 3);zeros(3, 2*num_nodes_i) eye(3, 3)];
                    elseif((num_nodes_i == 1)&(num_nodes_j ~= 1))
                       N_ij = [eye(3, 3) zeros(3, 2*num_nodes_j);zeros(3, 3) N_j];  
                    else    
                       N_ij = [N_i zeros(3, 2*num_nodes_j);zeros(3, 2*num_nodes_i) N_j];
                    end      

                            k_ij_hat = [c/L  0         0        -c/L 0         0;
                                0    12*d/L^3  6*d/L^2  0    -12*d/L^3 6*d/L^2;
                                0    6*d/L^2   4*d/L    0    -6*d/L^2  2*d/L;
                                -c/L 0         0        c/L  0         0;
                                0    -12*d/L^3 -6*d/L^2 0    12*d/L^3  -6*d/L^2;
                                0    6*d/L^2   2*d/L    0    -6*d/L^2  4*d/L];
                    ct = del_coor(1)/L;
                    st = del_coor(2)/L;
                    T = [ct  st 0 0   0  0;
                        -st ct 0 0   0  0;
                        0   0  1 0   0  0;
                        0   0  0 ct  st 0;
                        0   0  0 -st ct 0;
                        0   0  0 0   0  1];
                    k_ij = T'*k_ij_hat*T;
                                    
                 
                    el_stiff_C0 = el_stiff_C0 + 0.5*N_ij'*k_ij*N_ij*dvol_i*dvol_j;
                    
                end
            end
        end
    end
end
%
% Expand matrix to C1 continuity (3*num_nodes DOF total)
%

el_stiff_C1 = zeros(3*num_nodes_i+3*num_nodes_j, 3*num_nodes_i+3*num_nodes_j);
for el_row_node = 1 : num_nodes_i + num_nodes_j
    k_row = (el_row_node-1)*3 + 1;
    kc0_row = (el_row_node-1)*2 + 1;
    for el_col_node = 1 : num_nodes_i + num_nodes_j
        k_col = (el_col_node-1)*3 + 1;
        kc0_col = (el_col_node-1)*2 + 1;
        el_stiff_C1( k_row : k_row + 1, k_col : k_col + 1) = ...
            el_stiff_C0(kc0_row : kc0_row + 1, kc0_col : kc0_col+1);
    end
end

if ((num_nodes_i == 1)&(num_nodes_j == 1))    
    el_stiff_C1 = el_stiff_C0;
end

if ((num_nodes_i ~= 1)&(num_nodes_j == 1))
    for el_row_node = 1 : num_nodes_i + num_nodes_j
        k_row = (el_row_node-1)*3 + 1;
        kc0_row = (el_row_node-1)*2 + 1;
        el_stiff_C1( k_row : k_row + 1, 3*num_nodes_i + 3*num_nodes_j) = ...
            el_stiff_C0(kc0_row : kc0_row + 1, 2*num_nodes_i + 3);
    end
    for el_col_node = 1 : num_nodes_i + num_nodes_j
        k_col = (el_col_node-1)*3 + 1;
        kc0_col = (el_col_node-1)*2 + 1;
        el_stiff_C1( 3*num_nodes_i + 3*num_nodes_j, k_col : k_col + 1) = ...
            el_stiff_C0( 2*num_nodes_i + 3, kc0_col : kc0_col+1);
    end
     el_stiff_C1( 3*num_nodes_i + 3*num_nodes_j,3*num_nodes_i + 3*num_nodes_j ) = ...
            el_stiff_C0( 2*num_nodes_i + 3,2*num_nodes_i + 3 );
    
end

if ((num_nodes_i == 1)&(num_nodes_j ~= 1))
    el_stiff_C1 = zeros(3*num_nodes_i+3*num_nodes_j, 3*num_nodes_i+3*num_nodes_j);

    for el_row_node = 2 : num_nodes_i + num_nodes_j
        k_row = (el_row_node-1)*3 + 1;
        kc0_row = (el_row_node-1)*2 + 2;
        for el_col_node = 2 : num_nodes_i + num_nodes_j
            k_col = (el_col_node-1)*3 + 1;
            kc0_col = (el_col_node-1)*2 + 2;
            el_stiff_C1( k_row : k_row + 1, k_col : k_col + 1) = ...
                el_stiff_C0(kc0_row : kc0_row + 1, kc0_col : kc0_col+1);
        end
    end
    el_stiff_C1(1 : 3, 1 : 3) = el_stiff_C0(1 : 3, 1 : 3);

    for el_col_node = 2 : num_nodes_i + num_nodes_j
        k_col = (el_col_node-1)*3 + 1;
        kc0_col = (el_col_node-1)*2 + 2;
        el_stiff_C1( 1 : 3, k_col : k_col + 1) = ...
            el_stiff_C0(1 : 3, kc0_col : kc0_col+1);
    end

    for el_row_node = 2 : num_nodes_i + num_nodes_j
        k_row = (el_row_node-1)*3 + 1;
        kc0_row = (el_row_node-1)*2 + 2;
        el_stiff_C1( k_row : k_row + 1, 1 : 3) = ...
            el_stiff_C0(kc0_row : kc0_row + 1, 1 : 3);
    end


end

return