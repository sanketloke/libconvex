function [phi_f,der_phi_x,a] = getphi(f,sym_x,sym1,sym2,a)
    new_f=subs(subs(f,sym_x(1),sym1),sym_x(2),sym2);
    phi_f=new_f;
    der_phi_x=diff(phi_f);
    a=a;
end