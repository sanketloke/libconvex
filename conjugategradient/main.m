clear;
f='nocedalfuncCh3';
u.p=50*ones(12,1)';
u.f=feval(f,u.p,1);
u.g=feval(f,u.p,2);
u.h=feval(f,u.p,3);

bfgsparams = struct('maxit',10000,'toler',1.0e-4);

fu=@(fun,p) feval(fun, p, 1);
gu=@(fun,p) feval(fun, p, 2);


%%% darray contains the converging value
%%% iterarray contains the no of iterations required to convergant. 

[inform]=FletcherReevesCG(f,fu,gu,u,bfgsparams);
darray(1)=inform.f
iterarray(1)=size(inform.fvals,2);

[inform]=PolakReibereCG(f,fu,gu,u,bfgsparams);
darray(2)=inform.f;
iterarray(2)=size(inform.fvals,2);

[inform]=FRPRCG(f,fu,gu,u,bfgsparams);
darray(3)=inform.f;
iterarray(3)=size(inform.fvals,2);

[inform]=HestenesStiefelCG(f,fu,gu,u,bfgsparams);
darray(4)=inform.f;
iterarray(4)=size(inform.fvals,2);

[inform]=ACG1(f,fu,gu,u,bfgsparams);
darray(5)=inform.f;
iterarray(5)=size(inform.fvals,2);

[inform]=ACG2(f,fu,gu,u,bfgsparams);
darray(6)=inform.f;
iterarray(6)=size(inform.fvals,2);