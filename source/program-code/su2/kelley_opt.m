%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Globale Optimierung - Sonderuebung II
%
% Aufgabe S.2.3
% b): Implementieren Sie das Schnittebenenverfahren von Kelley (vgl. Algorithmus 2.1 im Skript), das eine Approximation eines globalen Minimalpunktes 
%      eines konvexen Optimierungsproblems
%	 P:  \min_{x \in \R^n} c^T x \text{ s.t. } g_i(x) \leq 0, ~ i \in I, ~Ax \leq b 
%      mit der nicht-leeren kompakten Mengen $M	^0 \coloneqq \{ x \in \R^m | Ax \leq b \}$ berechnet. 
%      Die in jeder Iteration auftretenden linearen Optimierungsprobleme sollen mit dem Solver Gurobi geloest werden. 
%      Ihr Algorithmus erhaelt folgende Inputs:
%          - Das Gurobi-Modell zu Hilfsproblems $LP^0: \quad \min_{x \in \R^n} c^T x \text{ s.t. } x \in M^0$.
%          - Die Restriktionsfunktionen $g_1(x), \dotsc, g_p(x)$.
%          - Die Gradienten (bzw. Ableitungen) der Restriktionsfuntionen.
%          - Eine Toleranz $\epsilon$.
%%
% Nadine Kostorz (1972005)
% Martin Belica  (1775706)

function x_lsg = kelley_opt(linear_model, g, jacobi_g, eps)

% lineares Modell
lp     = linear_model; 
result = gurobi(lp); 

% Anfangsparameter
j          = 0; 
x          = result.x; 
val        = g(x);
[g_max, k] = max(val); 

while g_max > eps 
    jacobi_val = jacobi_g(x); 
    
    b       = -val(k) - dot(jacobi_val(k, :),-x); 
    fmatrix = full(lp.A);
    fmatrix(end+1:end+size(jacobi_val(k, :), 1), :) = jacobi_val(k, :); 
    j = j+1; 
  
    % Model-Update
    lp.A            = sparse(fmatrix); 
    lp.rhs(end+1)   = b; 
    lp.sense(end+1) = ['<'];
    
    result = gurobi(lp);
    
    x = result.x; 
    val = g(result.x); 
    [g_max, k] = max(val); 
end

x_lsg = x;

end