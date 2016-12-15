% Rosenbrock function
% output:
% f: function
% gradf: gradient
% hessf: hessian
% x: symbolic_array
function [f,gradf,hessf,x] = rosenbrock()
    % Calculate objective f
    syms x1 x2 a real
    x=[x1 x2];
    f = 100*(x2 - x1^2)^2 + (1-x1)^2;
    gradf = jacobian(f,x).'; % column gradf
    hessf = jacobian(gradf,x);
end