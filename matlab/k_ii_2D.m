function [el_stiff_C1] = k_ii_2D(coords, elem_props)
%
% Compute the contribution to the stiffness matrix of two elements: j wrt
% i, using the peridynamic model
%
num_nodes = elem_props(2);
As = elem_props(5);
c = elem_props(3);
d = elem_props(4);
t = elem_props(5);
n = elem_props(6);
mat_horizon = elem_props(7);
el_stiff_C0 = zeros(2*num_nodes, 2*num_nodes);

del = 1/n;
number = 0;
nc = 0;
for i1 = 1  : n
    for i2 = 1 : n
        if(num_nodes == 3)
            xi_i1 = del*i1 - del/2;
            xi_i2 = del*i2 - del/2;
            if( i2 <= (n + 1 - i1))
                xi_i1 = xi_i1 - del/6;
                xi_i2 = xi_i2 - del/6;
            else
                xi_i1 = (1 - xi_i1) + del/6;
                xi_i2 = (1 - xi_i2) + del/6;
            end
        elseif(num_nodes == 4)
            del = 2/n;
            xi_i1 = -1 + del*i1 - del/2;
            xi_i2 = -1 + del*i2 - del/2;
        elseif(num_nodes == 2)
            del = 2/n;
            xi_i1 = -1 + del*i1 - del/2;
            xi_i2 = 0;
        elseif(num_nodes == 1)
            del = 2/n;
            xi_i1 = 0;
            xi_i2 = 0;       
        else    
            return
        end
        [coor_i, N_i, det_J_i] = peri_quants([xi_i1 xi_i2], coords);
        
        if (num_nodes == 2)
            dvol_i = det_J_i*del*As;        
        else
            dvol_i = det_J_i*del*del*t;
        end  

        for j1 = 1  : n
            for j2 = 1  : n
                if(num_nodes == 3)
                    xi_j1 = del*j1 - del/2;
                    xi_j2 = del*j2 - del/2;
                    if( j2 <= (n + 1 - j1))
                        xi_j1 = xi_j1 - del/6;
                        xi_j2 = xi_j2 - del/6;
                    else
                        xi_j1 = (1 - xi_j1) + del/6;
                        xi_j2 = (1 - xi_j2) + del/6;
                    end
                elseif(num_nodes == 4)
                    del = 2/n;
                    xi_j1 = -1 + del*j1 - del/2;
                    xi_j2 = -1 + del*j2 - del/2;
                elseif(num_nodes == 2)
                    del = 2/n;
                    xi_j1 = -1 + del*j1 - del/2;
                    xi_j2 = 0;
                elseif(num_nodes == 1)
                    del = 2/n;
                    xi_j1 = 0;
                    xi_j2 = 0;           
                else
                    return
                end
                [coor_j, N_j, det_J_j] = peri_quants([xi_j1 xi_j2], coords);               
               
                if (num_nodes == 2)
                    dvol_j = det_J_j*del*As;
                else
                    dvol_j = det_J_j*del*del*t;
                end                                 
                
                del_coor = coor_j - coor_i;                
                
                L = (del_coor*del_coor')^0.5;
                 
                 
                if((L <= mat_horizon) & (L > 0))
                    
                    N_ij = [N_i; N_j];
                                               
                                                            
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
%el_stiff_C0

%
% Expand matrix to C1 continuity (3*num_nodes DOF total)
%

el_stiff_C1 = zeros(3*num_nodes, 3*num_nodes);
for el_row_node = 1 : num_nodes
    k_row = (el_row_node-1)*3 + 1;
    kc0_row = (el_row_node-1)*2 + 1;
    for el_col_node = 1 : num_nodes
        k_col = (el_col_node-1)*3 + 1;
        kc0_col = (el_col_node-1)*2 + 1;
        el_stiff_C1( k_row : k_row + 1, k_col : k_col + 1) = ...
            el_stiff_C0(kc0_row : kc0_row + 1, kc0_col : kc0_col+1);
    end
end


return