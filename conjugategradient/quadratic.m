function val=quadratic(x,order)
    if order==1
        val=func(x);
    elseif order==2
        val=grad_func(x);
    else
        val=hess_func(x);
    end
end

function func_value = func(x)
    func_value=x(1)^2+x(2)^2+x(3)^2+x(4)^2;
end
function grad_func_value = grad_func(x)
   grad_func_value=[ 2*(x(1)) 2*x(2) 2*x(3) 2*x(4)]';
end

function hess_func_value = hess_func(x)
    hess_func_value=[2,0,0,0;0,2,0,0;0,0,2,0;0,0,0,2];
end
