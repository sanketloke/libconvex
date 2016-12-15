function [inform] = gaussnewton(fun,fu,gu,guJ, x, cgtparams)
    f=@(x) fu(fun,x);
    g=@(x) gu(fun,x);
    J=@(x) guJ(fun,x);
    inform=struct;
    %  Initialize parameter structure for StepSize function call.
    params = struct('ftol', 1e-4, 'gtol', 0.9, 'xtol', 1e-6, 'stpmin', 0, ...
                    'stpmax', 1e20, 'maxfev', 10000);
    %  Populate local caching of cgtparams parameters.
    toler = cgtparams.toler;  % Set gradient tolerance.
    maxit = cgtparams.maxit;  % Set maximum number of allowed iterations.
    xk=x.p;

    pk=-( J(xk)' * J(xk)) * (J(xk))' * f(xk)';
    tvals(1)=0
    k=0;
    for i= 1:maxit
        tic;
        u.p=xk;
        u.f=f(xk);
        fvals(i)=(f(xk));
        u.g=g(xk);
        pk=-( J(xk)' * J(xk)) * g(xk);
        %%%%%%%%%%%%%%%%%%%%% Compute alphak
        [alfa, x] = StepSizeSW(fun,u, pk, 1, params);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        xk=xk+alfa*pk';

        t=toc;
        tvals(i)=t;
        if i>1 & (fvals(i)>=fvals(i-1))
            break;
        end
        %  Check for termination condition: norm of gradient less than toler.
        if norm(g(xk)) < toler
          break;
        end
    end
    if norm(g(xk)) < toler
      inform.status=1;
      inform.fvals=fvals;
      inform.x=xk;
      inform.tvals=tvals;
      inform.f=f(xk);
    else
      inform.status=0;
      inform.fvals=fvals;
      inform.x=xk;
      inform.tvals=tvals;
      inform.f=f(xk);
    end
end
