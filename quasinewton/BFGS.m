function [inform] = FletcherReevesCG(fun,fu,gu, x, bfgsparams)
    f=@(x) fu(fun,x);
    g=@(x) gu(fun,x);
    inform=struct;
    %  Initialize parameter structure for StepSize function call.
    params = struct('ftol', 1e-4, 'gtol', 0.9, 'xtol', 1e-6, 'stpmin', 0, ...
                    'stpmax', 1e20, 'maxfev', 10000);
    %  Populate local caching of cgtparams parameters.
    toler = bfgsparams.toler;  % Set gradient tolerance.
    maxit = bfgsparams.maxit;  % Set maximum number of allowed iterations.
    xk=x.p;

    k=0;
    I= eye(size(x.p,1));
    H= I;
    for i = 1:maxit
        %i
        pk=-H*g(xk);

        u.p=xk;
        u.f=f(xk);
        fvals(i)=u.f;
        u.g=g(xk);
        %%%%%%%%%%%%%%%%%%%%% Compute alphak
        [alfa, x] = StepSizeSW(fun,u, pk, 1, params);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        xk1=xk+alfa*pk';

        sk=xk1 - xk;
        yk= (g(xk1) - g(xk))';
        rho = 1/ (sk*yk');


        H= ( I- rho* sk * yk')*H*( I- rho*  yk * sk') + rho * sk * sk';

        xk=xk1;

        %pk=pk1;
        if i>1 & (fvals(i)>fvals(i-1))
            break;
        end
        %  Check for termination condition: norm of gradient less than toler.
        if norm(g(xk1)) < toler
          break;
        end
    end
    if norm(g(xk)) < toler
      inform.status=1;
      inform.fvals=fvals;
      inform.x=xk;
      inform.f=f(xk);
    else
      inform.status=0;
      inform.fvals=fvals;
      inform.x=xk;
      inform.f=f(xk);
    end
end
