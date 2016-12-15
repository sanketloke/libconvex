function [inform]= levemar(point,fu,gu,guJ,fun)

    %%%%%%%%%%% Function assignment
  %   F = @(x) func(x);
	% dF = @(x) grad_func(x);
  %   d2F= @(x) hess_func(x);

    f=@(x) fu(fun,x);
    g=@(x) gu(fun,x);
    J=@(x) guJ(fun,x);
	%Initialize hyperparameters
    trust_radius_max = 2;    % step size (learning rate)
    n = 0.4;  % tolerance of error

    % progress tracking
    fvals = [];       % store F(x) values across iterations
	progress = @(iter,x) fprintf('iter = %3d: x = %-32s, F(x) = %f\n', ...
    iter, mat2str(x,6), f(x));
	% Iterate
	iter = 1;         % iterations counter
	x = point;    % initial guess
    fvals(iter) = f(x); % Collection of function vals
    tvals(iter)=0; % Collection of time taken vals
    trust_radius=trust_radius_max/10;


	while (norm(g(x))>=0.01) & iter<=100
        tic; % Counter Hits
        d2F=(J(x)'*J(x));  % B to be utilized later in the model _decrease
        B=d2F;

        F=f(x);
        dF=g(x);
        %%%% function returns the solution to the quadratic optimization
        %%%% problem under constraints of the trust-radius
        p=(approximate_solution(F,dF,d2F,x,trust_radius))';

        %Capturing rho in the trust-region
        model_decrease=model(F,dF,B,zeros(size(x,2),1))-model(F,dF,B,p');
        true_decrease=f(x)-f(x+p);
        rho=true_decrease/model_decrease;

        % Heuristical operation of controlling trust-region
        if rho<1/4
            trust_radius=0.4*trust_radius;
        else
           if rho>3/4 &  norm(p)==trust_radius
               trust_radius=min(2*trust_radius,trust_radius_max);
           end
        end
        if rho>0.2
            x = x + p;  % gradient descent
        end


        %Bookkeeping
        etime=toc;
        tvals(iter)=etime;
	       fvals(iter) = f(x);     % evaluate objective function
        iter=iter+1;
        try
            progress(iter, x);      % show progress
        catch
            x;
        end
    end
    optimal_val=f(x);
    optimal_function_minimum=x;
	% Plot
	% plot(1:(iter-1), fvals(1:iter-1), 'LineWidth',2); grid on;
	% title('Logarithmic Error'); xlabel('Iteration'); ylabel('Error');
  %   saveas(figure1,strcat('newton',int2str(int32(rand(1)*100)),'.jpg'),'jpg')  % here you save the figure
  %   %subplot(2,1,2);
	% %plot(1:(iter-1), tvals(1:iter-1), 'LineWidth',2); grid on;
	% title('Time Taken Per Step'); xlabel('Iteration'); ylabel('Time taken');
  %   %saveas(figure1,strcat('',int2str(int32(rand(1)*100)),'.jpg'),'jpg')  % here you save the figure
	% Evaluate final solution of system of equations G(x)=0
	disp('G(x) = '); disp(f(x));
  inform= struct;
  inform.fvals=fvals(1:iter-1);
  inform.status=1;
  inform.x=x;
  inform.tvals=tvals(1:iter-1);
end

function p=approximate_solution(F,dF,d2F,x,trust_radius)

    if norm(inv(d2F)*dF) <=trust_radius

        %%%%%%%%%%Calculating Full Length Step direction if less than the
        %%%%%%%%%%trust-radius
        pFS=-inv(d2F)*dF;
        p=pFS;
    else
        pFS=-inv(d2F)*dF;




        %%%%%%%%Calculate Cauchy point
        pks=-trust_radius*dF/norm(dF);
        temp=dF'*d2F*dF;
        % if temp<=0
        %     tau=1;
        % else
        %     tau=min(1,(norm(dF^3)/(trust_radius*temp)));
        % end
        %pcauchy=tau*pks;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Gradient descent of the quadratic model
        pku= - ((dF'*dF/(temp))) * dF;


        %Finding tau using equation solver
        tau=findtau(pku,pFS,trust_radius);
        %Tuning p
        if (min(tau)>1 & min(tau)<2)
            p=pku+(double(min(tau))-1)*(pFS-pku);
        elseif (max(tau)>1 & max(tau)<2)
            p=pku+(double(max(tau))-1)*(pFS-pku);
        else
            tau=trust_radius/norm(pku);
            p=tau*pku;
        end
    end
end

function val= model(F,dF,d2F,p)

    val= F+p'*dF+0.5*p'*d2F*p;
end
