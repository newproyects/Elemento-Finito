function el_dist = elem_distance(icoords, jcoords)
%
% Assume, for now, the distance is between element nodes
%

dist_ij = sqrt((icoords(1)-jcoords(1))^2 + (icoords(2) - jcoords(2))^2);
    
el_dist = dist_ij;

return