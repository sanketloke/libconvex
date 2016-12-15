function [optimal_val,optimal_function_minimum] = gradientdescent(func,jacobian,sym_x,point)
    
	F = @(x) evaluate_f(func,sym_x,x');

	% Gradient of F (partial derivatives)
	dF = @(x) normr(evaluate_f(jacobian,sym_x,x'));


	%Initialize hyperparameters
    alphamax=0.01;
    GAMMA = 0.001;    % step size (learning rate)
    MAX_ITER = 1000;  % maximum number of iterations
    FUNC_TOL = 0.3;   % termination tolerance for F(x)

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
        %%%% Wolfe condition line search
        %%%%
        alpha0=0;
        syms a;
        temp=dF(x);
        [phi_f,der_phi_x,a] = getphi(func,sym_x,x(1)+a*temp(1),x(2)+a*temp(2),a);
        PHI = @(x) evaluate_f(phi_f,a,1);
        DER_PHI= @(x) normr(evaluate_f(der_phi_x,a,1));
        c1=10^-4;
        c2=0.9;
        falpha=[];
        i=1;
        alpha(i)=(alpha0 + alphamax)/2;
        alphastr=1;
        while 1==1
            if ( PHI(alpha(i))> PHI(0)+c1*alpha*DER_PHI(0) | (i>1 & PHI(alpha(i))>=PHI(alpha(i-1))))
                try
                    alphastr=zoom_opti(alpha(i-1),alpha(i),PHI,DER_PHI,c1,c2);
                catch
                    d=1;
                end
                    
                break;
            end
            if DER_PHI(alpha(i))<=-c2*DER_PHI(0)
                alphastr=alpha(i);
            end
            if DER_PHI(alpha(i))>=0 && i>1
                alphastr=zoom_opti(alpha(i),alpha(i-1),PHI,DER_PHI,c1,c2);
                break;
            end
            alpha(i+1)=(alpha(i)+alphamax)/2;
            i=i+1;
        end
        %%%%%%%%%%%%%%%%
        %Updating GAMMA to new step
        GAMMA=alphastr;
        %%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%
	    x = x - GAMMA * dF(x);  % gradient descent
	    fvals(iter) = F(x);     % evaluate objective function
	    progress(iter, x,calculations);      % show progress
        dF(x);
        if (fvals(iter)-fvals(iter-1))> 0
            break;
        end
        %if (fvals(iter)-fvals(iter-1))> 0
        %    alphamax=alphamax/2;
        %end
    end
    
    optimal_val=F(x);
    optimal_function_minimum=x;
	% Plot
    figure1=figure;
	plot(1:iter, fvals, 'LineWidth',2); grid on;
	title('Objective Function'); xlabel('Iteration'); ylabel('F(x)');
    saveas(figure1,strcat('gradientdescentwolfe',int2str(int32(rand(1)*100)),'.jpg'),'jpg')  % here you save the figure
	% Evaluate final solution of system of equations G(x)=0
	disp('G(x) = '); disp(F(x))
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

function alphastr= zoom_opti(alpha1,alpha2,PHI,DER_PHI,c1,c2)
    alpha_zoom=[];
    j=1;
    alpha_zoom(j)=(alpha1+alpha2)/2;
    while 1==1 
        j;
        alpha_zoom(j)=alpha1+(alpha2-alpha1)/10000;
       if PHI(alpha_zoom(j))>PHI(0)+c1*alpha_zoom(j)*DER_PHI(0) || PHI(alpha_zoom(j))>=PHI(alpha1)
           alpha2=alpha_zoom(j);
       else
           if abs(DER_PHI(alpha_zoom(j)))<=-c2*DER_PHI(0)
               alphastr=alpha_zoom(j);
               break;
           end
           if DER_PHI(alpha_zoom(j))*(alpha2-alpha1)>=0
               alpha2=alpha1;
           end
           alpha1=alpha_zoom(j);
       end
       if length(alpha_zoom)>2
           p=alpha_zoom(end);
           q=alpha_zoom(end-1);
           if abs(p-q)<0.001
               break;
           end
       end
       if rem(j,100)==0
           alpha_zoom(j);
       end
       j=j+1;
    end
    alphastr=alpha_zoom(j);
end
