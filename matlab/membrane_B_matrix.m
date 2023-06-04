function [ B, det_J, radius] = membrane_B_matrix(coords, plane_code, xi)
% membrane_B_matrix: returns the B (strain to nodal displacement) matrix
% and the determinant of the Jacobian matrix, det_J,
% evaluated at natural coordinate, xi,
% within a membrane element with nodes at locations denoted by coords
%
num_nodes = size(coords, 1);
J = zeros(2, 2);
[N, shapedir] = shape_functions(xi, 'membrane', num_nodes);
for row = 1 : 2
    for col = 1 : 2
        for node = 1 : num_nodes
            J(row, col) = J(row, col)...
                + shapedir(node, row)*coords(node, col);
        end
    end
end

det_J = det(J);
J_inv = inv(J);
if(plane_code == 3)
    B = zeros(4, num_nodes*2);
    radius = 0;
    for(node = 1 : num_nodes)
        radius = radius + N(node)*coords(node, 1);
    end
else
    radius = 0;
    B = zeros(3, num_nodes*2);
end

for(node = 1 : num_nodes)
    B(1, (2*node) - 1) = shapedir(node, 1)*J_inv(1,1)...
        + shapedir(node, 2)*J_inv(1,2);
    B(2, (2*node)) = shapedir(node, 1)*J_inv(2,1)...
        + shapedir(node, 2)*J_inv(2,2);
    B(3, (2*node)) = B(1, (2*node) - 1);
    B(3, (2*node) - 1) =B(2, (2*node));
    if(plane_code == 3)
        B(4, (2*node) - 1) = N(node)/radius;
    end
end

return

    
