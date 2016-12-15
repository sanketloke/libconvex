clear;
f='nocedalfuncCh3';
u.p=zeros(12,1)';
u.f=feval(f,u.p,1);
u.g=feval(f,u.p,2);
u.h=feval(f,u.p,3);

bfgsparams = struct('maxit',10000,'toler',1.0e-4);

fu=@(fun,p) feval(fun, p, 1);
gu=@(fun,p) feval(fun, p, 2);


 %darray contains the converging value
 %iterarray contains the no of iterations required to convergant.

[inform]=BFGS(f,fu,gu,u,bfgsparams);
darray(1)=inform.f
iterarray(1)=size(inform.fvals,2);


[inform]=DFP(f,fu,gu,u,bfgsparams);
darray(2)=inform.f
iterarray(2)=size(inform.fvals,2);

sr1params = struct('maxit',10000,'toler',1.0e-4);

[f,g,h] = hw3func();
[optimal_val,optimal_function_minimum,fvals] =SR1TrustRegion(u.p',f,g,h);
darray(3)=optimal_val;
iterarray(3)=size(fvals,1);


darray
iterarray
