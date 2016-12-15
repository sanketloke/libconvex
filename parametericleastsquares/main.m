clear;


% gaussNewtonTemperatureHigh=
% gaussNewtonTemperatureLow=
% levenMarTemperatureHigh=
% levenMarTemperatureLow=
% gaussNewtonPrecipitation=
% levenMarPrecipitation=


random=rand(4,1)';

f='func_precipitation';
u.p = random;
u.f=sum(abs(feval(f,u.p,1)));
u.g=feval(f,u.p,2);

leastsquaresparams = struct('maxit',10000,'toler',1.0e-4);

fu=@(fun,p) feval(fun, p, 1);
gu=@(fun,p) feval(fun, p, 2);
guJ=@(fun,p) feval(fun, p, 3);


[gaussNewtonPrecipitation]=gaussnewton(f,fu,gu,guJ,u,leastsquaresparams);


f='func_precipitation';
[levenMarPrecipitation]= levenmar(u.p,fu,gu,guJ,f)

iter= length(gaussNewtonPrecipitation.fvals)+1
%Plot
figure;
figure1=figure;
subplot(2,1,1);
plot(1:(iter-1), log(gaussNewtonPrecipitation.fvals), 'LineWidth',2); grid on;
title('Log Error of GaussNewton on Precipitation'); xlabel('Iteration'); ylabel('Error');
subplot(2,1,2);
plot(1:(iter-1), gaussNewtonPrecipitation.tvals, 'LineWidth',2); grid on;
title('Time Taken Per Step'); xlabel('Iteration'); ylabel('Time taken');
pause(1)
saveas(figure1,strcat('gaussNewtonPrecipitation',int2str(int32(random(1)*100)),'.jpg'),'jpg')  % here you save the figure


iter= length(levenMarPrecipitation.fvals)+1
%Plot
figure;
figure1=figure;
subplot(2,1,1);
plot(1:(iter-1), log(levenMarPrecipitation.fvals), 'LineWidth',2); grid on;
title('Log Error of Levenberg-Marquardt on Precipitation'); xlabel('Iteration'); ylabel('Error');
subplot(2,1,2);
plot(1:(iter-1), levenMarPrecipitation.tvals, 'LineWidth',2); grid on;
title('Time Taken Per Step'); xlabel('Iteration'); ylabel('Time taken');
pause(1)
saveas(figure1,strcat('levenMarPrecipitation',int2str(int32(random(1)*100)),'.jpg'),'jpg')  % here you save the figure






f='func_temperature_high';
u.p = random;
u.f=sum(abs(feval(f,u.p,1)));
u.g=feval(f,u.p,2);

leastsquaresparams = struct('maxit',10000,'toler',1.0e-4);

fu=@(fun,p) feval(fun, p, 1);
gu=@(fun,p) feval(fun, p, 2);
guJ=@(fun,p) feval(fun, p, 3);


[gaussNewtonTemperatureHigh]=gaussnewton(f,fu,gu,guJ,u,leastsquaresparams);


f='func_temperature_high';
[levenMarTemperatureHigh]= levenmar(u.p,fu,gu,guJ,f)





iter= length(gaussNewtonTemperatureHigh.fvals)+1

%Plot
figure;
figure1=figure;
subplot(2,1,1);
plot(1:(iter-1), log(gaussNewtonTemperatureHigh.fvals), 'LineWidth',2); grid on;
title('Log Error of GaussNewton on High Temperature'); xlabel('Iteration'); ylabel('Error');
subplot(2,1,2);
plot(1:(iter-1), gaussNewtonTemperatureHigh.tvals, 'LineWidth',2); grid on;
title('Time Taken Per Step'); xlabel('Iteration'); ylabel('Time taken');
pause(1)
saveas(figure1,strcat('gaussNewtonTemperatureHigh',int2str(int32(random(1)*100)),'.jpg'),'jpg')  % here you save the figure


iter= length(levenMarTemperatureHigh.fvals)+1
%Plot
figure;
figure1=figure;
subplot(2,1,1);
plot(1:(iter-1), log(levenMarTemperatureHigh.fvals), 'LineWidth',2); grid on;
title('Log Error of Levenberg-Marquardt on High Temperature'); xlabel('Iteration'); ylabel('Error');
subplot(2,1,2);
plot(1:(iter-1), levenMarTemperatureHigh.tvals, 'LineWidth',2); grid on;
title('Time Taken Per Step'); xlabel('Iteration'); ylabel('Time taken');
pause(1)
saveas(figure1,strcat('levenMarTemperatureHigh',int2str(int32(random(1)*100)),'.jpg'),'jpg')  % here you save the figure

f='func_temperature_low';
u.p = random;
u.f=sum(abs(feval(f,u.p,1)));
u.g=feval(f,u.p,2);

leastsquaresparams = struct('maxit',10000,'toler',1.0e-4);

fu=@(fun,p) feval(fun, p, 1);
gu=@(fun,p) feval(fun, p, 2);
guJ=@(fun,p) feval(fun, p, 3);


[gaussNewtonTemperatureLow]=gaussnewton(f,fu,gu,guJ,u,leastsquaresparams);


f='func_temperature_low';
[levenMarTemperatureLow]= levenmar(u.p,fu,gu,guJ,f)




iter= length(gaussNewtonTemperatureLow.fvals)+1

%Plot
figure;
figure1=figure;
subplot(2,1,1);
plot(1:(iter-1), log(gaussNewtonTemperatureLow.fvals), 'LineWidth',2); grid on;
title('Log Error of GaussNewton on Low Temperature'); xlabel('Iteration'); ylabel('Error');
subplot(2,1,2);
plot(1:(iter-1), gaussNewtonTemperatureLow.tvals, 'LineWidth',2); grid on;
title('Time Taken Per Step'); xlabel('Iteration'); ylabel('Time taken');
pause(1)
saveas(figure1,strcat('gaussNewtonTemperatureLow',int2str(int32(random(1)*100)),'.jpg'),'jpg')  % here you save the figure


iter= length(levenMarTemperatureLow.fvals)+1
%Plot
figure;
figure1=figure;
subplot(2,1,1);
plot(1:(iter-1), log(levenMarTemperatureLow.fvals), 'LineWidth',2); grid on;
title('Log Error of Levenberg-Marquardt on Low Temperature'); xlabel('Iteration'); ylabel('Error');
subplot(2,1,2);
plot(1:(iter-1), levenMarTemperatureLow.tvals, 'LineWidth',2); grid on;
title('Time Taken Per Step'); xlabel('Iteration'); ylabel('Time taken');
pause(1)
saveas(figure1,strcat('levenMarTemperatureLow',int2str(int32(random(1)*100)),'.jpg'),'jpg')  % here you save the figure


disp('Initial Value')
disp(random)
disp('Converged Values for Gauss Newton on Precipitation Data');
disp(gaussNewtonPrecipitation.x)
disp(gaussNewtonPrecipitation.fvals(end))

disp('Converged Values for Gauss Newton on High Temperature Data');
disp(gaussNewtonTemperatureHigh.x);
disp(gaussNewtonTemperatureHigh.fvals(end));

disp('Converged Values for Gauss Newton on Low Temperature Data');;
disp(gaussNewtonTemperatureLow.x);
disp(gaussNewtonTemperatureLow.fvals(end));


disp('Converged Values for Levenberg-Marquardt on Precipitation Data');
disp(levenMarPrecipitation.x)
disp(levenMarPrecipitation.fvals(end))

disp('Converged Values for Levenberg-Marquardt on  High Temperature Data');
disp(levenMarTemperatureHigh.x);
disp(levenMarTemperatureHigh.fvals(end));

disp('Converged Values for Levenberg-Marquardt on Low Temperature Data');
disp(levenMarTemperatureLow.x);
disp(levenMarTemperatureLow.fvals(end));
