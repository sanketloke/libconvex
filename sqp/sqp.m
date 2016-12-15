function [inform] = sqp( x,l, cgtparams)

    inform=struct;
    %  Initialize parameter structure for StepSize function call.
    params = struct('ftol', 1e-4, 'gtol', 0.9, 'xtol', 1e-6, 'stpmin', 0, ...
                    'stpmax', 1e20, 'maxfev', 10000);
    %  Populate local caching of cgtparams parameters
    toler = cgtparams.toler;  % Set gradient tolerance.
    maxit = cgtparams.maxit;  % Set maximum number of allowed iterations.
    for i = 1:maxit
        df=hw9func(x,l,2)';
        Lkk= hw9func(x,l,3);
        Ak= hw9func(x,l,5)' ;
        ck=hw9func(x,l,4)';
        fvals(i)=(hw9func(x,l,1));
        cvals(i)=sum(ck);
        %keyboard
        KKTMatrix=[ df; ck];

        LMatrix=[[Lkk,-Ak'];[Ak,zeros(size(Ak,1))]];
        if any(eig(LMatrix))<0
          break
        end
        pvec= LMatrix\ (-KKTMatrix);
        pk=pvec(1:5)';
        u.f=fvals(i);
        u.p=x;
        u.g=hw9func(x,l,2);
        x=x+pk;
        l=pvec(6:end);
    end
    xk=x;
    inform.status=1;
    inform.fvals=fvals;
    inform.x=xk;
    inform.l=l;
    inform.cvals=cvals;
    inform.ck=hw9func(x,l,4);
    inform.f=hw9func(x,l,1);
end
