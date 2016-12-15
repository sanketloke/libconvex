function val=hw9func(x,l,order)
    x=x';
    if order==1 % Calculate function value
        val=f(x);
    elseif order==2 % Calculate function value
        val=df(x);
    elseif order==3
        val=DLxx(x,l);
    elseif order==4
        val=c(x);
    else
        val=A(x);
    end
end


function out = f(x)
  out= exp(x(1)*x(2)*x(3)*x(4)*x(5)) - 0.5 * (x(1)^3 + x(2)^3+ 1)^2;
end

function out = df(x)
  out=zeros(1,5);
  t=x;
  x=t(1);
  y=t(2);
  z=t(3);
  w=t(4);
  u=t(5);
  out(1)= u*w*y*z* exp(x*y*w*z*u) -3 *x^2*(x^3  + y^3 + 1);
  out(2)= u*w*x*z* exp(x*y*w*z*u) -3 *y^2*(x^3  + y^3 + 1);
  out(3)= u*w*x*y* exp(x*y*w*z*u);
  out(4)= u*z*x*y* exp(x*y*w*z*u);
  out(5)= w*z*x*y* exp(x*y*w*z*u);
  %out=zeros(1,5);
  x=t;
  % out(1)=x(2)*x(3)*x(4)*x(5)*exp(x(1)*x(2)*x(3)*x(4)*x(5))  - 3 *x(1) *(x(1)^2 + x(2)^2+ 1);
  % out(2)=x(1)*x(3)*x(4)*x(5)*exp(x(1)*x(2)*x(3)*x(4)*x(5))  - 3 *x(2) *(x(1)^2 + x(2)^2+ 1);
  % out(3)=x(1)*x(2)*x(4)*x(5)*exp(x(1)*x(2)*x(3)*x(4)*x(5));
  % out(4)=x(1)*x(3)*x(2)*x(5)*exp(x(1)*x(2)*x(3)*x(4)*x(5));
  % out(5)=x(1)*x(3)*x(4)*x(2)*exp(x(1)*x(2)*x(3)*x(4)*x(5));
  % keyboard;
end


function out= DLxx(x,l)
  t=x;
  x=t(1);
  y=t(2);
  z=t(3);
  w=t(4);
  u=t(5);
  out=zeros(5,5);
  % out(1,1)=w^2*x^2*y^2*z^2*exp(x*y*z*w*u) - 2*l(1);
  % out(1,2)=l(2)+u*w*x^2*y^2*z^2*exp(x*y*z*w*u) + x*y*z*exp(x*y*z*w*u);
  % out(1,3)=u*w^2*x*y^2*z^2*exp(x*y*z*w*u) + w*y*z*exp(x*y*z*w*u);
  % out(1,4)=u*w^2*x^2*y*z^2*exp(x*y*z*w*u) + w*x*z*exp(x*y*z*w*u);
  % out(1,5)=u*w^2*x^2*y^2*z*exp(x*y*z*w*u) + w*x*y*exp(x*y*z*w*u);
  %
  % out(2,1)=l(2) + u*w*x^2*y^2*z^2*exp(x*y*z*w*u)+x*y*z*exp(x*y*z*w*u);
  % out(2,2)=u^2*x^2*y^2*z^2*exp(x*y*z*w*u) -2*l(1);
  % out(2,3)=u^2*w*x*y^2*z^2*exp(x*y*z*w*u) + u*y*z*exp(x*y*z*w*u);
  % out(2,4)=u^2*w*x^2*y*z^2*exp(x*y*z*w*u) + u*x*z*exp(x*y*z*w*u);
  % out(2,5)=u^2*w*x^2*y^2*z*exp(x*y*z*w*u) + u*x*y*exp(x*y*z*w*u);
  %
  % out(3,1)=u*w^2*x*y^2*z^2*exp(x*y*z*w*u) + w*y*z*exp(x*y*z*w*u);
  % out(3,2)=u^2*w*x*y^2*z^2*exp(x*y*z*w*u) + u*y*z*exp(x*y*z*w*u);
  % out(3,3)=-2*l(1) - 6*l(3) *x +u^2*w^2*y^2*z^2*exp(x*y*z*w*u) - 6 *x* (x^3 + y^3 +1 ) -9*x^4;
  % out(3,4)=u^2*w^2*x*y*z^2*exp(x*y*z*w*u) + u*w*z*exp(x*y*z*w*u) -9 *x^2*y^2;
  % out(3,5)=u^2 * w^2 *x * y^2* z*exp(x*y*z*w*u) + u*w*y*exp(x*y*z*w*u);
  %
  % out(4,1)=u*w^2*x^2*y*z^2*exp(x*y*z*w*u) + w*x*z*exp(x*y*z*w*u);
  % out(4,2)=u^2*w*x^2*y*z^2*exp(x*y*z*w*u) + u*x*z*exp(x*y*z*w*u);
  % out(4,3)=u^2*w^2*x*y*z^2*exp(x*y*z*w*u) + u*w*z*exp(x*y*z*w*u) -9 *x^2*y^2;
  % out(4,4)=-2*l(1) - 6*l(3) *y +u^2*w^2*x^2*z^2*exp(x*y*z*w*u) - 6*y*(x^3 + y^3 +1 ) -9*y^4;
  % out(4,5)=-l(2) +u^2*w^2*x^2*y*z *exp(x*y*z*w*u) + u*w*x*exp(x*y*z*w*u);
  %
  %
  % out(5,1)=u*w^2*x^2*y^2*z*exp(x*y*z*w*u) + w*x*y*exp(x*y*z*w*u);
  % out(5,2)=u^2*w*x^2*y^2*z*exp(x*y*z*w*u) + u*x*y*exp(x*y*z*w*u);
  % out(5,3)=u^2 * w^2 *x * y^2* z*exp(x*y*z*w*u) + u*w*y*exp(x*y*z*w*u);
  % out(5,4)=-l(2) +u^2*w^2*x^2*y*z *exp(x*y*z*w*u) + u*w*x*exp(x*y*z*w*u);
  % out(5,5)=u^2*w^2*x^2*y^2 *exp(x*y*z*w*u) -2*l(1);
  x=t;
  ep = exp(prod(x));
  f = ep - 0.5*(x(1)^3+x(2)^3+1)^2;

  H=[ x(2)^2*x(3)^2*x(4)^2*x(5)^2*ep-9*x(1)^4-6*(x(1)^3+x(2)^3+1)*x(1), x(3)*x(4)*x(5)*ep+x(2)*x(3)^2*x(4)^2*x(5)^2*x(1)*ep-9*x(2)^2*x(1)^2, x(2)*x(4)*x(5)*ep+x(2)^2*x(3)*x(4)^2*x(5)^2*x(1)*ep, x(2)*x(3)*x(5)*ep+x(2)^2*x(3)^2*x(4)*x(5)^2*x(1)*ep, x(2)*x(3)*x(4)*ep+x(2)^2*x(3)^2*x(4)^2*x(5)*x(1)*ep;
  x(3)*x(4)*x(5)*ep+x(2)*x(3)^2*x(4)^2*x(5)^2*x(1)*ep-9*x(2)^2*x(1)^2, x(1)^2*x(3)^2*x(4)^2*x(5)^2*ep-9*x(2)^4-6*(x(1)^3+x(2)^3+1)*x(2), x(1)*x(4)*x(5)*ep+x(1)^2*x(3)*x(4)^2*x(5)^2*x(2)*ep, x(1)*x(3)*x(5)*ep+x(1)^2*x(3)^2*x(4)*x(5)^2*x(2)*ep, x(1)*x(3)*x(4)*ep+x(1)^2*x(3)^2*x(4)^2*x(5)*x(2)*ep;
  x(2)*x(4)*x(5)*ep+x(2)^2*x(3)*x(4)^2*x(5)^2*x(1)*ep,x(1)*x(4)*x(5)*ep+x(1)^2*x(3)*x(4)^2*x(5)^2*x(2)*ep,x(1)^2*x(2)^2*x(4)^2*x(5)^2*ep,x(1)*x(2)*x(5)*ep+x(1)^2*x(2)^2*x(4)*x(5)^2*x(3)*ep,x(1)*x(2)*x(4)*ep+x(1)^2*x(2)^2*x(4)^2*x(5)*x(3)*ep;
  x(2)*x(3)*x(5)*ep+x(2)^2*x(3)^2*x(4)*x(5)^2*x(1)*ep,x(1)*x(3)*x(5)*ep+x(1)^2*x(3)^2*x(4)*x(5)^2*x(2)*ep,x(1)*x(2)*x(5)*ep+x(1)^2*x(2)^2*x(4)*x(5)^2*x(3)*ep,x(1)^2*x(2)^2*x(3)^2*x(5)^2*ep,x(1)*x(2)*x(3)*ep+x(1)^2*x(2)^2*x(3)^2*x(5)*x(4)*ep;
  x(2)*x(3)*x(4)*ep+x(2)^2*x(3)^2*x(4)^2*x(5)*x(1)*ep,x(1)*x(3)*x(4)*ep+x(1)^2*x(3)^2*x(4)^2*x(5)*x(2)*ep,x(1)*x(2)*x(4)*ep+x(1)^2*x(2)^2*x(4)^2*x(5)*x(3)*ep,x(1)*x(2)*x(3)*ep+x(1)^2*x(2)^2*x(3)^2*x(5)*x(4)*ep, x(1)^2*x(2)^2*x(3)^2*x(4)^2*ep];

  c1= diag([2*x(1) 2*x(2) 2*x(3) 2*x(4) 2*x(5)]);
  c2= zeros(5,5);
  c2(2,3) =1;
  c2(3,2)=1;
  c2(4,5)=-1;
  c2(5,4)=-1;

  c3=zeros(5,5);
  c3(1,1)=6*x(1);
  c3(2,2) = 6*x(2);
  out=zeros(5,5);
  out=H - l(1) * c1 -l(2)* c2 -l(3)*c3;





end

function out= c(x)
  out=zeros(1,3);
  out(1)= (x(1)^2 + x(2)^2 +x(3)^2 +x(4)^2 +x(5)^2 -10 );
  out(2)= x(2) * x(3) - 5*x(4)*x(5);
  out(3)= x(1)^3 + x(2)^3 + 1;
end

function out= A(x)
  out=zeros(5,3);
  out(1,1)=2*x(1);
  out(2,1)=2*x(2);
  out(3,1)=2*x(3);
  out(4,1)=2*x(4);
  out(5,1)=2*x(5);

  out(1,2)=0;
  out(2,2)=x(3);
  out(3,2)=x(2);
  out(4,2)=-5*x(5);
  out(5,2)=-5*x(4);

  out(1,3)=3*x(1)^2;
  out(2,3)=3*x(2)^2;
  out(3,3)=0;
  out(4,3)=0;
  out(5,3)=0;
end
