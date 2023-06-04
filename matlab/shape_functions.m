function [N, N_derivs] = shape_functions(xi, elem_type, num_nodes)
%shape_functions: shape functions and shape function derivatives
% for all Co element types
% inputs: xi: vector of 1, 2, or 3 natural coordinates xi, eta, gamma
% outputs: N: vector of num_nodes shape functions evaluated at location xi
%          N_derivs: each row contains the derivative with respect to csi
%          and eta of the corresponding shape function, for example:
%  |N1,xi   N1,eta|
%  |N2,xi   N2,eta|
%  |N3,xi   N3,eta|
%  |N4,xi   N4,eta|
%
num_dimensions = size(xi, 2);
N = zeros(1, num_nodes);
N_derivs = zeros(num_nodes, num_dimensions);
if(elem_type(1:5) == 'membr')
    s = xi(1);
    t = xi(2);
    if(num_nodes == 3)
        N = zeros(1, 3);
        N(1, 1) = s;
        N(1, 2) = t;
        N(1, 3) = 1 - s - t;

        N_derivs = zeros(3, 2);
        N_derivs(1, 1) = 1;
        N_derivs(2, 1) =  0;
        N_derivs(3, 1) =  -1;
        N_derivs(1, 2) = 0;
        N_derivs(2, 2) = 1;
        N_derivs(3, 2) =  -1;
        
    elseif(num_nodes == 4)
        N(1, 1) = (1-s)*(1-t)/4;
        N(1, 2) = (1+s)*(1-t)/4;
        N(1, 3) = (1+s)*(1+t)/4;
        N(1, 4) = (1-s)*(1+t)/4;

        N_derivs(1, 1) = -0.25*(1 - t);
        N_derivs(2, 1) =  0.25*(1 - t);
        N_derivs(3, 1) =  0.25*(1 + t);
        N_derivs(4, 1) = -0.25*(1 + t);
        N_derivs(1, 2) = -0.25*(1 - s);
        N_derivs(2, 2) = -0.25*(1 + s);
        N_derivs(3, 2) =  0.25*(1 + s);
        N_derivs(4, 2) =  0.25*(1 - s);
        
    elseif(num_nodes == 6)
        u = 1 - s - t;
        N(1, 1) = 2*u*(u - 0.5);
        N(1, 2) = 4*s*u;
        N(1, 3) = 2*s*(s - 0.5);
        N(1, 4) = 4*s*t;
        N(1, 5) = 2*t*(t - 0.5);
        N(1, 6) = 4*t*u;
        
        N_derivs = [  -3+4*s+4*t   -3+4*s+4*t;
                      4-8*s-4*t    -4*s;
                      4*s-1        0;
                      4*t          4*s;
                      0            4*t-1;
                      -4*t         4-4*s-8*t];

    elseif(num_nodes == 8)
        N(1, 1) = (1 - s)*(1 - t)*(-s - t - 1)/4;
        N(1, 2) = (1 - t)*(1 + s)*(1 - s)/2;
        N(1, 3) = (1 + s)*(1 - t)*(s - t - 1)/4;
        N(1, 4) = (1 + s)*(1 + t)*(1 - t)/2;
        N(1, 5) = (1 + s)*(1 + t)*(s + t - 1)/4;
        N(1, 6) = (1 + t)*(1 + s)*(1 - s)/2;
        N(1, 7) = (1 - s)*(1 + t)*(-s + t - 1)/4;
        N(1, 8) = (1 - s)*(1 + t)*(1 - t)/2;
        
        N_derivs = [ -1/4*(-1+t)*(2*s+t)         -1/4*(-1+s)*(s+2*t);
                      -s+t*s                      -1/2+1/2*s^2;
                      -1/4*(-1+t)*(2*s-t)         -1/4*(1+s)*(s-2*t);
                      1/2-1/2*t^2                 -t-t*s;
                      1/4*(1+t)*(2*s+t)           1/4*(1+s)*(s+2*t);
                      -s-t*s                     1/2-1/2*s^2;
                      1/4*(1+t)*(2*s-t)           1/4*(-1+s)*(s-2*t);
                      -1/2+1/2*t^2                -t+s*t];
    else

    end
elseif(elem_type(1:5) == 'solid')
    s = xi(1);
    t = xi(2);
    u = xi(3);
    if(num_nodes == 4)
        N(1, 1) = s;
        N(1, 2) = t;
        N(1, 3) = u;
        N(1, 4) = 1 - s - t - u;
        N_derivs(1, 1) = 1;
        N_derivs(2, 1) = 0;
        N_derivs(3, 1) = 0;
        N_derivs(4, 1) = -1;
        N_derivs(1, 2) = 0;
        N_derivs(2, 2) = 1;
        N_derivs(3, 2) = 0;
        N_derivs(4, 2) = -1;
        N_derivs(1, 3) = 0;
        N_derivs(2, 3) = 0;
        N_derivs(3, 3) = 1;
        N_derivs(4, 3) = -1;
    elseif(num_nodes == 10)
        return
    elseif(num_nodes == 8)
        N(1, 1) = 0.125*(1-s)*(1-t)*(1-u);
        N(1, 2) = 0.125*(1+s)*(1-t)*(1-u);
        N(1, 3) = 0.125*(1+s)*(1+t)*(1-u);
        N(1, 4) = 0.125*(1-s)*(1+t)*(1-u);
        N(1, 5) = 0.125*(1-s)*(1-t)*(1+u);
        N(1, 6) = 0.125*(1+s)*(1-t)*(1+u);
        N(1, 7) = 0.125*(1+s)*(1+t)*(1+u);
        N(1, 8) = 0.125*(1-s)*(1+t)*(1+u);
        N_derivs(1, 1) = -0.125*(1-t)*(1-u);
        N_derivs(2, 1) = 0.125*(1-t)*(1-u);
        N_derivs(3, 1) = 0.125*(1+t)*(1-u);
        N_derivs(4, 1) = -0.125*(1+t)*(1-u);
        N_derivs(5, 1) = -0.125*(1-t)*(1+u);
        N_derivs(6, 1) = 0.125*(1-t)*(1+u);
        N_derivs(7, 1) = 0.125*(1+t)*(1+u);
        N_derivs(8, 1) = -0.125*(1+t)*(1+u);
        
        N_derivs(1, 2) = -0.125*(1-s)*(1-u);
        N_derivs(2, 2) = -0.125*(1+s)*(1-u);
        N_derivs(3, 2) = 0.125*(1+s)*(1-u);
        N_derivs(4, 2) = 0.125*(1-s)*(1-u);
        N_derivs(5, 2) = -0.125*(1-s)*(1+u);
        N_derivs(6, 2) = -0.125*(1+s)*(1+u);
        N_derivs(7, 2) = 0.125*(1+s)*(1+u);
        N_derivs(8, 2) = 0.125*(1-s)*(1+u);
        
        N_derivs(1, 3) = -0.125*(1-s)*(1-t);
        N_derivs(2, 3) = -0.125*(1+s)*(1-t);
        N_derivs(3, 3) = -0.125*(1+s)*(1+t);
        N_derivs(4, 3) = -0.125*(1-s)*(1+t);
        N_derivs(5, 3) = 0.125*(1-s)*(1-t);
        N_derivs(6, 3) = 0.125*(1+s)*(1-t);
        N_derivs(7, 3) = 0.125*(1+s)*(1+t);
        N_derivs(8, 3) = 0.125*(1-s)*(1+t);
        
    elseif(num_nodes == 20)
        return
    else

    end
end

return