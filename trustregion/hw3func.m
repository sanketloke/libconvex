function [f,g,h]=hw3func()
    f=@func;
    g=@grad_func;
    h=@hess_func;
end


function func_value = func(x)
    func_value=0;
    for i=1:(size(x,1)/2)
        func_value =func_value + ( (1-x(2*i-1))^2+ 10*(x(2*i) - x(2*i-1))^2);
    end
end


function grad_func_value = grad_func(x)
    grad_func_value=zeros(size(x,1),1);
    for i=1:(size(x,1)/2)
        grad_func_value(2*i-1)=2*(1-x(2*i-1))-40*(x(2*i)-x(2*i-1)^2)*x(2*i-1);
        grad_func_value(2*i)=20*(x(2*i)-x(2*i-1)^2);
    end
end

function hess_func_value = hess_func(x)
    hess_func_value=zeros(size(x,1),size(x,1));
    for i=1:(size(x,1)/2)
        hess_func_value(2*i-1,2*i-1)=-2-40*x(2*i)+3*40*x(2*i-1)^2;%-2-40*(x(2*i))+40*3*x(2*i-1)^2;
        hess_func_value(2*i,2*i)=20;
        hess_func_value(2*i-1,2*i)=-40*x(2*i-1);
        hess_func_value(2*i,2*i-1)=-40*x(2*i-1);
    end
    
end