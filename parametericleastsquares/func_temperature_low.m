  function val=func_temperature_low(x,order)
    x=x';
    if order==1
        val=sum(abs(func(x)));
    elseif order==2
        val=grad_func(x)' * (func(x))';
    else
        val=grad_func(x);
    end
end

function out = func(x)
  t=[1,2,3,4,5,6,7,8,9,10,11,12];
  r= [21 24 30 38 47 56 60 59 51 39 31 24];
  out=zeros(1,12);
  for u=1:length(t)
    out(u)=x(1)*sin(x(2)*t(u) + x(3))+x(4)-r(u);
  end
end

function J = grad_func(x)
  J=zeros(12,4);
  t=[1,2,3,4,5,6,7,8,9,10,11,12];
  r= [21 24 30 38 47 56 60 59 51 39 31 24];
  for i = 1:length(J)
    J(i,1)= sin(x(2)*t(i)+x(3));
    J(i,2)= x(1) * cos(x(2) * t(i) + x(3)) * t(i);
    J(i,3)= x(1) * cos(x(2) * t(i) + x(3));
    J(i,4)= 1;
  end
end
