function J = jacs(~,Y, params, xspan, P, Vn, Vc, Ac)

% This function calculates the Jacobian of the functions from 
% ftn_dosage with respect to the State Variables

%-------------------------------
% Number of nuclei in present NC
M = length(Y)/2;

% ----------------
% Unpacking params
params0 = num2cell(params);
[lambdaU, lambdaW, KeqD, KeqDC, gamma, beta, phi, kappa ] = ...
        params0{:};
    
%---------------------------
% Unpack Structure Variables
u = Y(1:M);
w = Y(M+1:2*M);

%--------------------------------------------
% Make matrix versions of the State Variables
U   = spdiags(u,0,M,M);
W   = spdiags(w,0,M,M);    

%----------------
% Other Parameters
I    = speye(M);
v    = 1;
c_d  = KeqD* Vn + Vc;
c_dc = KeqDC*Vn + Vc;

%----------------------------------------------
% Derivative of the Michaelis-Menten function
% ftn = g*w./(kappa + dl0*w)
x   = xspan;
g   = exp(-x.^2/2/phi^2);
G   = g.*kappa./(kappa + w).^2;
G   = spdiags(G,0,M,M);

% =========================================================================
% Jacobian elements
% =========================================================================

%---------------------------------------------------
% First row block, derivatives of f1, the eqn for u.
% f1 = (-Vn*(gamma*KeqD*u) + lambdaU*Ac*P*u + Vc*(beta*g.*w./(kappa + w) -
% gamma*u.*v))/c_d

f1_u = (-Vn*gamma*KeqD*I + lambdaU*Ac*P)/c_d;
f1_w = beta*G/c_d;

J1 = [f1_u f1_w];

%-----------------------------------------------------
% Second row block, derivatives of f2, the eqn for uc.
% f2 = (Vn*(gamma*KeqD*u) + lambdaW*Ac*P*w - Vc*(beta*g.*w./(kappa + w) -
% gamma*u.*v))/c_dc

f2_u = (Vn*gamma*KeqD*I - Vc*gamma*v*I)/c_dc;
f2_w = (lambdaW*Ac*P - beta*G)/c_dc;

J2 = [f2_u f2_w];

J = [J1; J2];

%{

%}

end


