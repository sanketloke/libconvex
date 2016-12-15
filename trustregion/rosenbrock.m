function [f,g,h] = rosenbrock()
    f=@banana;
    g=@bananad;
    h=@bananah;
end



function z = banana(xy)

x = xy(1);
y = xy(2);

z = 1 + 100*(x*x-y)^2 + (x-1)^2;

end

% bananad.m  -  Rosenbrock banana function gradient for M604 HW02
%
% Input parameter xy must be a column vector with two components.
% Output is a row vector of length 2.

function z = bananad(xy)

z = zeros(1,2);

x = xy(1);
y = xy(2);

z = [400*x^3-400*x*y+2*x-2,  200*(y-x^2)];
z=z';
end

% bananah.m  -  Rosenbrock banana function Hessian for M604 HW02
%
% Input parameter xy must be a column vector with two components.
% Output is a symmetric  2 x 2  matrix.

function z = bananah(xy)

z = zeros(2,2);

x = xy(1);
y = xy(2);

z12 = -400*x;

z = [1200*x^2-400*y+2, z12; ...
     z12, 200];
end