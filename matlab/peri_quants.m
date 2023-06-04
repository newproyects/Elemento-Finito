function [xypos, N, det_J] = peri_quants(xi, elcoords)

num_nodes = size(elcoords, 1)/2;


if(num_nodes == 3)
    N1 = xi(1);
    N2 = xi(2);
    N3 = 1 - xi(1) - xi(2);
    N = [N1 0  N2 0  N3 0;
        0  N1 0  N2 0  N3;
        0  0  0  0  0  0];
    xypos(1) = N(1, 1:2:5)*elcoords(1:2:5);
    xypos(2) = N(2, 2:2:6)*elcoords(2:2:6);

    N_derivs = zeros(3, 2);
    N_derivs(1, 1) = 1;
    N_derivs(2, 1) =  0;
    N_derivs(3, 1) =  -1;
    N_derivs(1, 2) = 0;
    N_derivs(2, 2) = 1;
    N_derivs(3, 2) =  -1;

elseif(num_nodes == 4)
    ms = 1-xi(1);
    ps = 1+xi(1);
    mt = 1-xi(2);
    pt = 1+xi(2);
    N1 = ms*mt/4;
    N2 = ps*mt/4;
    N3 = ps*pt/4;
    N4 = ms*pt/4;
    N = [N1 0  N2 0  N3 0  N4 0;
        0  N1 0  N2 0  N3 0  N4;
        0  0  0  0  0  0  0  0];
    
    xypos(1) = N(1, 1:2:7)*elcoords(1:2:7);
    xypos(2) = N(2, 2:2:8)*elcoords(2:2:8);        
        
    N_derivs = zeros(4, 2);
    
    N_derivs(1, 1) = -0.25*(1 - xi(2));
    N_derivs(2, 1) =  0.25*(1 - xi(2));
    N_derivs(3, 1) =  0.25*(1 + xi(2));
    N_derivs(4, 1) = -0.25*(1 + xi(2));
    N_derivs(1, 2) = -0.25*(1 - xi(1));
    N_derivs(2, 2) = -0.25*(1 + xi(1));
    N_derivs(3, 2) =  0.25*(1 + xi(1));
    N_derivs(4, 2) =  0.25*(1 - xi(1));
    
    
elseif(num_nodes == 2)
    N1 = (1 - xi(1))/2;
    N2 = (1 + xi(1))/2;
    
    N = [N1 0 N2 0;
         0  N1 0 N2;
         0  0  0  0];
    N_derivs(1, 1) = 0.5;
    N_derivs(2, 1) = -0.5;
    N_derivs(1, 2) =  0;
    N_derivs(2, 2) =  0;
    
    xypos(1) = N(1, 1:2:3)*elcoords(1:2:3);
    xypos(2) = N(2, 2:2:4)*elcoords(2:2:4);    
    
    
elseif(num_nodes == 1)
    
    N = [1 0;
         0 1;
         0 0];
     
    xypos(1) = N(1, 1)*elcoords(1);
    xypos(2) = N(2, 2)*elcoords(2);   


else    
    xypos = [];
    N = [];
    det_J = [];
    return
end

if (num_nodes > 2)
J = zeros(2, 2);
for row = 1 : 2
    for col = 1 : 2
        for node = 1 : num_nodes
            J(row, col) = J(row, col)...
                + N_derivs(node, row)*elcoords(2*(node-1) + col);
        end
    end
end
det_J = det(J);
inv_J = inv(J);
elseif (num_nodes == 2)
    
    diff_x = -0.5*elcoords(1) + 0.5*elcoords(3);
    diff_y = -0.5*elcoords(2) + 0.5*elcoords(4);
    
    det_J = norm([diff_x diff_y]);  

else    
    det_J = 1;
end


if(num_nodes == 3)
    N(3, 1) = -inv_J(2, 1)/2;
    N(3, 2) = inv_J(1, 1)/2;
    N(3, 3) = -inv_J(2, 2)/2;
    N(3, 4) = inv_J(1, 2)/2;
    N(3, 5) = (inv_J(2, 1) + inv_J(2, 2))/2;
    N(3, 6) = (-(inv_J(1, 1)+ inv_J(1, 2)))/2;
elseif(num_nodes == 4)
    
    x(1:4) = elcoords(1:2:7);
    y(1:4) = elcoords(2:2:8);
    
    dy_dxi2 = N_derivs(1, 2)*y(1) + N_derivs(2, 2)*y(2) + N_derivs(3, 2)*y(3) + N_derivs(4, 2)*y(4);
    dy_dxi1 = N_derivs(1, 1)*y(1) + N_derivs(2, 1)*y(2) + N_derivs(3, 1)*y(3) + N_derivs(4, 1)*y(4);
    dx_dxi2 = N_derivs(1, 2)*x(1) + N_derivs(2, 2)*x(2) + N_derivs(3, 2)*x(3) + N_derivs(4, 2)*x(4);
    dx_dxi1 = N_derivs(1, 1)*x(1) + N_derivs(2, 1)*x(2) + N_derivs(3, 1)*x(3) + N_derivs(4, 1)*x(4);     
    
    N(3, 1) = -(1/2)*(1/det_J)*(dx_dxi1*N_derivs(1, 2) - dx_dxi2*N_derivs(1, 1));
    N(3, 2) =  (1/2)*(1/det_J)*(dy_dxi2*N_derivs(1, 1) - dy_dxi1*N_derivs(1, 2));
    N(3, 3) = -(1/2)*(1/det_J)*(dx_dxi1*N_derivs(2, 2) - dx_dxi2*N_derivs(2, 1));
    N(3, 4) =  (1/2)*(1/det_J)*(dy_dxi2*N_derivs(2, 1) - dy_dxi1*N_derivs(2, 2));
    N(3, 5) = -(1/2)*(1/det_J)*(dx_dxi1*N_derivs(3, 2) - dx_dxi2*N_derivs(3, 1));
    N(3, 6) =  (1/2)*(1/det_J)*(dy_dxi2*N_derivs(3, 1) - dy_dxi1*N_derivs(3, 2));
    N(3, 7) = -(1/2)*(1/det_J)*(dx_dxi1*N_derivs(4, 2) - dx_dxi2*N_derivs(4, 1));
    N(3, 8) =  (1/2)*(1/det_J)*(dy_dxi2*N_derivs(4, 1) - dy_dxi1*N_derivs(4, 2));

end

return