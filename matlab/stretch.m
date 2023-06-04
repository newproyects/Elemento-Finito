function [strain, center_x, center_y, u_max, v_max, u_min, v_min] = stretch(coords, elem_props, displacements)

num_nodes = elem_props(2);


if (num_nodes == 3)
    x(1:3) = coords(1:2:5);
    y(1:3) = coords(2:2:6);
    u_i(1:2:5) = displacements(1:num_nodes,1);
    u_i(2:2:6) = displacements(1:num_nodes,2);

    J = [(x(1)-x(3)) (y(1)-y(3));
        (x(2)-x(3)) (y(2)-y(3));];

    det_J = det(J);
    inv_J = inv(J);

    B = [inv_J(1, 1) 0 inv_J(1, 2) 0 -(inv_J(1,1)+inv_J(1,2)) 0;
        0 inv_J(2, 1) 0 inv_J(2, 2) 0 -(inv_J(2,1)+inv_J(2,2))];

    B(3, 1) = inv_J(2, 1);
    B(3, 2) = inv_J(1, 1);
    B(3, 3) = inv_J(2, 2);
    B(3, 4) = inv_J(1, 2);
    B(3, 5) = (-(inv_J(2, 1) + inv_J(2, 2)));
    B(3, 6) = (-(inv_J(1, 1) + inv_J(1, 2)));


    e = B*u_i';

       
    theta_p(1) = (0.5)*atan(e(3)/(e(1) - e(2)));
    theta_p(2) = (0.5)*atan(e(3)/(e(1) - e(2))) - pi/2;
    
    p_strain(1) = (e(1) + e(2))/2 + 0.5*(e(1) - e(2))*cos(2*theta_p(1)) + (e(3)/2)*sin(2*theta_p(1));
    p_strain(2) = (e(1) + e(2))/2 + 0.5*(e(1) - e(2))*cos(2*theta_p(2)) + (e(3)/2)*sin(2*theta_p(2));


    max_strain = max(p_strain);
    min_strain = min(p_strain);
    
    strain = max_strain;

elseif (num_nodes == 1)
    u_i(1) = displacements(1,1);
    u_i(2) = displacements(1,2);
    center_x = 0;
    center_y = 0;
    u_max = 0;
    v_max = 0;
    u_min = 0;
    v_min = 0;

    strain = u_i(2);
elseif (num_nodes == 4)

    %Assume xi(1)=xi(2)=0
    
    u_i(1:2:7) = displacements(1:num_nodes,1);
    u_i(2:2:8) = displacements(1:num_nodes,2);

    xi = [0 0];
    plane_code = 1;
    B = membrane_B_matrix(coords, plane_code, xi);
    %{
    N_derivs = zeros(4, 2);
    N_derivs(1, 1) = -0.25*(1 - xi(2));
    N_derivs(2, 1) =  0.25*(1 - xi(2));
    N_derivs(3, 1) =  0.25*(1 + xi(2));
    N_derivs(4, 1) = -0.25*(1 + xi(2));
    N_derivs(1, 2) = -0.25*(1 - xi(1));
    N_derivs(2, 2) = -0.25*(1 + xi(1));
    N_derivs(3, 2) =  0.25*(1 + xi(1));
    N_derivs(4, 2) =  0.25*(1 - xi(1));
    
        
    B = zeros(3, 8);
    J = zeros(2, 2);
    for row = 1 : 2
        for col = 1 : 2
            for node = 1 : num_nodes
                J(row, col) = J(row, col)...
                    + N_derivs(node, row)*coords(2*(node-1) + col);
            end
        end
    end
    det_J = det(J);
    inv_J = inv(J);

    x(1:4) = coords(1:2:7);
    y(1:4) = coords(2:2:8);    
   
    u_i(1:2:7) = displacements(1:num_nodes,1);
    u_i(2:2:8) = displacements(1:num_nodes,2);
    

    dy_dxi2 = N_derivs(1, 2)*y(1) + N_derivs(2, 2)*y(2) + N_derivs(3, 2)*y(3) + N_derivs(4, 2)*y(4);
    dy_dxi1 = N_derivs(1, 1)*y(1) + N_derivs(2, 1)*y(2) + N_derivs(3, 1)*y(3) + N_derivs(4, 1)*y(4);
    dx_dxi2 = N_derivs(1, 2)*x(1) + N_derivs(2, 2)*x(2) + N_derivs(3, 2)*x(3) + N_derivs(4, 2)*x(4);
    dx_dxi1 = N_derivs(1, 1)*x(1) + N_derivs(2, 1)*x(2) + N_derivs(3, 1)*x(3) + N_derivs(4, 1)*x(4);     

    B(1, 1) = (1/det_J)*(dy_dxi2*N_derivs(1, 1) - dy_dxi1*N_derivs(1, 2));
    B(1, 3) = (1/det_J)*(dy_dxi2*N_derivs(2, 1) - dy_dxi1*N_derivs(2, 2));
    B(1, 5) = (1/det_J)*(dy_dxi2*N_derivs(3, 1) - dy_dxi1*N_derivs(3, 2));
    B(1, 7) = (1/det_J)*(dy_dxi2*N_derivs(4, 1) - dy_dxi1*N_derivs(4, 2));
    B(2, 2) = (1/det_J)*(dx_dxi1*N_derivs(1, 2) - dx_dxi2*N_derivs(1, 1));
    B(2, 4) = (1/det_J)*(dx_dxi1*N_derivs(2, 2) - dx_dxi2*N_derivs(2, 1));
    B(2, 6) = (1/det_J)*(dx_dxi1*N_derivs(3, 2) - dx_dxi2*N_derivs(3, 1));
    B(2, 8) = (1/det_J)*(dx_dxi1*N_derivs(4, 2) - dx_dxi2*N_derivs(4, 1));
    B(3, 1) = B(2, 2);
    B(3, 2) = B(1, 1);
    B(3, 3) = B(2, 4);
    B(3, 4) = B(1, 3);
    B(3, 5) = B(2, 6);
    B(3, 6) = B(1, 5);
    B(3, 7) = B(2, 8);
    B(3, 8) = B(1, 7);   
    %}
    
    x(1:4) = coords(1:4,1);
    y(1:4) = coords(1:4,2);
    
    

    e = B*u_i';
    
    eig([e(1) e(3); e(3) e(2)]);
        
    theta_p(1) = (0.5)*atan(e(3)/(e(1) - e(2)));
    theta_p(2) = (0.5)*atan(e(3)/(e(1) - e(2))) - pi/2;
    
    p_strain(1) = (e(1) + e(2))/2 + 0.5*(e(1) - e(2))*cos(2*theta_p(1)) + (e(3)/2)*sin(2*theta_p(1));
    p_strain(2) = (e(1) + e(2))/2 + 0.5*(e(1) - e(2))*cos(2*theta_p(2)) + (e(3)/2)*sin(2*theta_p(2));

    max_strain = max(p_strain);
    min_strain = min(p_strain);
    
    if (max_strain == p_strain(1))
        theta_max = theta_p(1);
        theta_min = theta_p(2);
    else
        theta_max = theta_p(2);
        theta_min = theta_p(1);
    end
    
    center_x = (x(2) + x(1))/2;
    center_y = (y(4) + y(2))/2;
    
    u_max = max_strain*cos(theta_max);
    v_max = max_strain*sin(theta_max);
    
    u_min = min_strain*cos(theta_min);
    v_min = min_strain*sin(theta_min); 
    
       
    strain = max_strain;
   
    
else

strain = 0;
center_x = 0;
center_y = 0;
u_max = 0;
v_max = 0;
u_min = 0;
v_min = 0;


end

return