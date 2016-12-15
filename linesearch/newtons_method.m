function [optimal_val,optimal_function_minimum] = newtons_method(func,jacobian,hessian,sym_x,point)
    
	F = @(x) evaluate_f(func,sym_x,x');

	% Gradient of F (partial derivatives)
	dF = @(x) evaluate_f(jacobian,sym_x,x');

    d2F= @(x) evaluate_f(hessian,sym_x,x');

	%Initialize hyperparameters
    alphamax=1;
    GAMMA = 0.001;    % step size (learning rate)
    MAX_ITER = 1000;  % maximum number of iterations
    FUNC_TOL = 0.1;   % termination tolerance for F(x)

    % progress tracking 
    fvals = [];       % store F(x) values across iterations
	progress = @(iter,x,calculations) fprintf('iter = %3d: x = %-32s, calculations=%f, F(x) = %f\n', ...
    iter, mat2str(x,6),calculations, F(x));
	% Iterate
	iter = 1;         % iterations counter
	x = point;    % initial guess
    fvals(iter) = F(x);
    
	while iter < MAX_ITER && fvals(end) > FUNC_TOL
        global calculations;
        calculations=0;
	    iter = iter + 1;
        
        %%%%%%%%%%%%%%%%%%%
        %%%% Backtracking line search
        %%%%
        rho=0.7;
        alpha=alphamax;
        c=0.7;
        while ~(F(x-alpha*dF(x))<=F(x)+c*alpha*dot(dF(x),dF(x)))
            F(x-alpha*dF(x));
            F(x)+c*alpha*dot(dF(x),dF(x));
            alpha;
            alpha=rho*alpha;
        end
        %%%%%%%%%%%%%%%%
        %Updating GAMMA to new step
        GAMMA=alpha;
        %%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%
	    x = x - GAMMA*inv(d2F(x)) * dF(x);  % gradient descent
	    fvals(iter) = F(x);     % evaluate objective function
	    progress(iter, x,calculations);      % show progress
%         if (fvals(iter)-fvals(iter-1))> 0
%             break
%         end
    end
    optimal_val=F(x);
    optimal_function_minimum=x;
	% Plot
    figure1=figure;
	plot(1:iter, fvals, 'LineWidth',2); grid on;
	title('Objective Function'); xlabel('Iteration'); ylabel('F(x)');
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


