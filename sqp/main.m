leastsquaresparams = struct('maxit',1000,'toler',1.0e-4);

x=[-1.8, 1.7, 1.9, -0.8 -0.8];
disp('SQP Initial Value of Input Vector')
disp(x)
disp('SQP Initial Value of Lagrange Multipliers')
l=[0,0,0]';
disp(l)
keyboard;
[inform]=sqp( x,l, leastsquaresparams);
%plot(inform.fvals)
disp('SQP Method converges to:')
disp(inform.x)
disp('Function Value')
disp(inform.f)
disp('Constraint Vector')
disp(inform.ck)
disp('SQP Method Lagrange Multipliers Vector')
disp(inform.l)
figure
figure1=figure;
plot(inform.fvals(1:200))
title('SQP Progress: Function Value'); xlabel('Iteration'); ylabel('f(x)');
saveas(figure1,strcat('SQPFunction',int2str(int32(rand(1)*100)),'.jpg'),'jpg')  % here you save the figure

figure
figure1=figure;
plot(inform.cvals(1:200))
title('SQP Progress: Sum of Constraints'); xlabel('Iteration'); ylabel('f(x)');
saveas(figure1,strcat('SQPConstraints',int2str(int32(rand(1)*100)),'.jpg'),'jpg')  % here you save the figure
