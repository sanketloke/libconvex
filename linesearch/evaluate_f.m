% returns evaluation of function on value
%
function value= evaluate_f(func,sym_x,point)
    global calculations;
    calculations=calculations+1;
    value=double(subs(func,sym_x,point));
end