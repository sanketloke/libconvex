[f,g,h]= hw3func();
x1=[0.0413;0.6432;0.6756;0.9006;0.4818;0.7941;0.8919;0.9035;0.7755;0.9288;0.8659;0.5903]
x2=[0.7967;0.0322;0.0121;0.5089;0.6617;0.9302;0.1893;0.9486;0.0772;0.3338;0.5090;0.1753]
x3=[0.9800;0.8121;0.3678;0.8734;0.8551;0.9032;0.3760;0.2975;0.0772;-0.1378;0.8106;0.5418]

display('%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display('%%%%%%%%%%%%%%%%%%%%%%%%% Newtons method Dogleg %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
[optimal_val,optimal_function_minimum]=newtons_method_dogleg(x1,f,g,h);
 
display('%%%%%%%%%%%%%%%%%%%%%%%%% Newtons method Subspace Minimization  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
[optimal_val,optimal_function_minimum] =newtons_method_subspace(x1,f,g,h);

display('%%%%%%%%%%%%%%%%%%%%%%%%% Newtons method Iterative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
[optimal_val,optimal_function_minimum] =newtons_method_iterative(x1,f,g,h);




display('%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

display('%%%%%%%%%%%%%%%%%%%%%%%%% Newtons method Dogleg %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
[optimal_val,optimal_function_minimum]=newtons_method_dogleg(x2,f,g,h);
 
display('%%%%%%%%%%%%%%%%%%%%%%%%% Newtons method Subspace Minimization  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
[optimal_val,optimal_function_minimum] =newtons_method_subspace(x2,f,g,h);

display('%%%%%%%%%%%%%%%%%%%%%%%%% Newtons method Iterative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
[optimal_val,optimal_function_minimum] =newtons_method_iterative(x2,f,g,h);


display('%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')



display('%%%%%%%%%%%%%%%%%%%%%%%%% Newtons method Dogleg %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
[optimal_val,optimal_function_minimum]=newtons_method_dogleg(x3,f,g,h);
 
display('%%%%%%%%%%%%%%%%%%%%%%%%% Newtons method Subspace Minimization  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
[optimal_val,optimal_function_minimum] =newtons_method_subspace(x3,f,g,h);

display('%%%%%%%%%%%%%%%%%%%%%%%%% Newtons method Iterative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
[optimal_val,optimal_function_minimum] =newtons_method_iterative(x3,f,g,h);
