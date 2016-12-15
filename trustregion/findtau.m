function tau=findtau(a,b,trust_radius)
    syms x real;
    x = [x]; % column vector of symbolic variables
    f = a'*a + (x-1)*(a'*b - a'*a) + (x-1)*(b'*a-a'*a)+(x-1)^2*norm(b-a)==trust_radius;
    solx=solve(f,x);
    tau=solx;
end