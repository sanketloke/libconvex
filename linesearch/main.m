%point=input('Input a starting position');
global calculations;
calculations=0;
[func,jacobian,hessf,sym_x] = rosenbrock();
display('%%%%%%%%%%%%%%%%%%%%%%%%% Gradient Descent Backtracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display('%%%%%%%%%%%%%%%%%%%%%%%%% Gradient Descent Backtracking :: Start point: [-1.2,1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
[optimal_val,optimal_function_minimum]=gradientdescent(func,jacobian,sym_x,[-1.2,1]');

display('%%%%%%%%%%%%%%%%%%%%%%%%% Gradient Descent Backtracking :: Start point: [1.2,1.2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
[optimal_val,optimal_function_minimum] =gradientdescent(func,jacobian,sym_x,[1.2,1.2]');

display('%%%%%%%%%%%%%%%%%%%%%%%%% Gradient Descent Wolfe Conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display('%%%%%%%%%%%%%%%%%%%%%%%%% Gradient Descent Wolfe Conditions Start point: [-1.2,1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
[optimal_val,optimal_function_minimum] =gradientdescent_wolfe(func,jacobian,sym_x,[-1.2,1]');

display('%%%%%%%%%%%%%%%%%%%%%%%%% Gradient Descent Wolfe Conditions Start point: [1.2,1.2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
[optimal_val,optimal_function_minimum] =gradientdescent_wolfe(func,jacobian,sym_x,[1.2,1.2]');


display('%%%%%%%%%%%%%%%%%%%%%%%%% Newtons method Backtracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display('%%%%%%%%%%%%%%%%%%%%%%%%% Newtons method Backtrackin:Start point: [-1.2,1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
[optimal_val,optimal_function_minimum] =newtons_method(func,jacobian,hessf,sym_x,[-1.2,1]');


display('%%%%%%%%%%%%%%%%%%%%%%%%% Newtons method Backtracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display('%%%%%%%%%%%%%%%%%%%%%%%%% Newtons method Backtrackin:Start point: [1.2,1.2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
[optimal_val,optimal_function_minimum] =newtons_method(func,jacobian,hessf,sym_x,[1.2,1.2]');



display('%%%%%%%%%%%%%%%%%%%%%%%%% Newtons method Wolfe %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display('%%%%%%%%%%%%%%%%%%%%%%%%% Newtons method Wolfe Start point: [-1.2,1] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
[optimal_val,optimal_function_minimum] =newtons_method_wolfe(func,jacobian,hessf,sym_x,[-1.2,1]');

display('%%%%%%%%%%%%%%%%%%%%%%%%% Newtons method Wolfe Start point: [1.2,1.2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
[optimal_val,optimal_function_minimum] =newtons_method_wolfe(func,jacobian,hessf,sym_x,[1.2,1.2]');