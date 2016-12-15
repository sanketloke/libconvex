function [optimal_val,optimal_function_minimum]= newtons_method_dogleg(point,func,grad_func,hess_func)
    
    %%%%%%%%%%% Function assignment
    F = @(x) func(x);
	dF = @(x) grad_func(x);
    d2F= @(x) hess_func(x);

	%Initialize hyperparameters
    trust_radius_max = 2;    % step size (learning rate)
    n = 0.4;  % tolerance of error

    % progress tracking 
    fvals = [];       % store F(x) values across iterations
	progress = @(iter,x) fprintf('iter = %3d: x = %-32s, F(x) = %f\n', ...
    iter, mat2str(x,6), F(x));
	% Iterate
	iter = 1;         % iterations counter
	x = point;    % initial guess
    fvals(iter) = F(x); % Collection of function vals
    tvals(iter)=0; % Collection of time taken vals
    trust_radius=trust_radius_max/10;
    
    
	while (norm(dF(x))>=0.01) & iter<=100
        tic; % Counter Hits
        B=d2F(x);  % B to be utilized later in the model _decrease
        % Making sure B is positive definite 
        if min(eig(B))<0
            B=B-1.1*min(eig(B))*eye(size(B));
        end
        
        
        %%%% function returns the solution to the quadratic optimization
        %%%% problem under constraints of the trust-radius
        p=(approximate_solution(F,dF,d2F,x,trust_radius));
        
        %Capturing rho in the trust-region
        model(F(x),dF(x),d2F(x),zeros(size(x,1),1));
        model(F(x),dF(x),d2F(x),p);
        model_decrease=model(F(x),dF(x),B,zeros(size(x,1),1))-model(F(x),dF(x),B,p);
        true_decrease=F(x)-F(x+p);
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
	    fvals(iter) = log(F(x));     % evaluate objective function
        iter=iter+1;
        try
            progress(iter, x);      % show progress
        catch
            x;
        end
    end
    optimal_val=F(x);
    optimal_function_minimum=x;
	% Plot
    figure
    figure1=figure;
    subplot(2,1,1);
	plot(1:(iter-1), fvals(1:iter-1), 'LineWidth',2); grid on;
	title('Logarithmic Error'); xlabel('Iteration'); ylabel('Error');
    saveas(figure1,strcat('newton',int2str(int32(rand(1)*100)),'.jpg'),'jpg')  % here you save the figure
    subplot(2,1,2);
	plot(1:(iter-1), tvals(1:iter-1), 'LineWidth',2); grid on;
	title('Time Taken Per Step'); xlabel('Iteration'); ylabel('Time taken');
    saveas(figure1,strcat('newton',int2str(int32(rand(1)*100)),'.jpg'),'jpg')  % here you save the figure
	% Evaluate final solution of system of equations G(x)=0
	disp('G(x) = '); disp(F(x));

	% Output:
	%
	% iter =   1: x = [0;0;0]                         , F(x) = 58.456136
	% iter =   2: x = [0.0075;0.002;-0.20944]         , F(x) = 23.306394
	% iter =   3: x = [0.015005;0.0015482;-0.335103]  , F(x) = 10.617030
	% ...
	% iter = 187: x = [0.683335;0.0388258;-0.52231]   , F(x) = 0.101161
	% iter = 188: x = [0.684666;0.0389831;-0.522302]  , F(x) = 0.099372
	%
	% (converged in 188 iterations after exceeding termination tolerance for F(x))

end

function p=approximate_solution(F,dF,d2F,x,trust_radius)
    if norm(inv(d2F(x))*dF(x)) <=trust_radius 
        
        %%%%%%%%%%Calculating Full Length Step direction if less than the
        %%%%%%%%%%trust-radius
        pFS=-inv(d2F(x))*dF(x);
        p=pFS;
    else
        pFS=-inv(d2F(x))*dF(x);
        
        
        
        
        %%%%%%%%Calculate Cauchy point
        pks=-trust_radius*dF(x)/norm(dF(x));
        temp=dF(x)'*d2F(x)*dF(x);
        if temp<=0
            tau=1;
        else
            tau=min(1,(norm(dF(x))^3)/(trust_radius*temp));
        end
        pcauchy=tau*pks;
        
        
        
        %%%%%%%%%%%%%%%% true Decrease obtained by newtons method
        true_decrease_cauchy=F(x)-F(x+pcauchy);
        %%%%%%%%%%%%%%%% model Decrease obtained by newtons method
        model_decrease_cauchy=model(F(x),dF(x),d2F(x),zeros(size(x,1),1))-model(F(x),dF(x),d2F(x),pcauchy);
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Gradient descent of the quadratic model
        pku= - ((dF(x)'*dF(x))/(temp)) * dF(x);
        
        
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
        
        
        %%%% Logs
        %%%%%%%%%%%%%%%% true Decrease obtained by newtons method
        true_decrease=F(x)-F(x+p);
        
        %%%%%%%%%%%%%%%% model Decrease obtained by newtons method
        model_decrease=model(F(x),dF(x),d2F(x),zeros(size(x,1),1))-model(F(x),dF(x),d2F(x),p);
    end
end

function val= model(F,dF,B,p)
    if min(eig(B))<0
        B=B-1.1*min(eig(B))*eye(size(B));
    end
    val= F+p'*dF+0.5*p'*B*p;
end