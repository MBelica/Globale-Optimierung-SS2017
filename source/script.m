%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Globale Optimierung - Sonderuebung II
%
% Aufgabe S.2.3
% b): Testen Sie Ihren Algorithmus an dem Optimierungsproblem
%
% P: \min_{x \in \R^2} x_2 s.t. x_1^2 + x_2^2 \leq 1, e^{x_1} - x_2 \leq 0
% 	   x_1 - x_2 \leq 0,  x_1 \in [-1, -\frac{1}{2}],  x_2 \in [-1, 1] $$
%          fuer $\epsilon \coloneqq 10^{-7}$.
%%
% Nadine Kostorz (1972005)
% Martin Belica  (1775706)

epsilon = 10^(-7); 

g   = @(x) [x(1).^2+x(2).^2-1 exp(x(1))-x(2)];
ddg = @(x) [2*x(1) 2*x(2); exp(x(1)) -1]; 

clear linear_model; 
linear_model.A          = sparse([1 -1]);
linear_model.obj        = [0 1];
linear_model.modelsense = 'min';
linear_model.rhs        = [0];
linear_model.sense      = ['<'];
linear_model.lb         = [-1; -1]; 
linear_model.ub         = [-0.5; 0.9]; 

x_lsg = kelley_opt(linear_model, g, ddg, epsilon)