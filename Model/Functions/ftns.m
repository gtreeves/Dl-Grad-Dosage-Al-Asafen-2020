function dYdt = ftns(~,Y,params, xspan, P, Vn, Vc, Ac)

% This is the function-handle where all the differential equations reside
% for the dorsal dosage project.

%---------------------------
% Unpack Structure Variables
M = length(Y)/2;
u = Y(1:M);
w = Y(M+1:2*M);

%------------------
% Unpacking params
params0 = num2cell(params);
[lambdaU, lambdaW, KeqD, KeqDC, gamma, beta, phi, kappa ] = ...
        params0{:};
    
%-------------------------
% Defining some parameters
x    = xspan;
v    = 1;
c_d  = (KeqD * Vn) + Vc;
c_dc = (KeqDC* Vn) + Vc;
g    = exp(-x.^2/2/phi^2);
Bn   = Vn*(gamma*KeqD*u);                           % dl/Cact binding eqn in nucleus 
Bc   = Vc*( beta*g.*w./(kappa + w) - gamma*u.*v);   % dl/Cact binding eqn in cytoplasm
   
%--------------------------------------------------------------------------
% diffyQ's
f1 = (-Bn + lambdaU*Ac*P*u + Bc)/c_d;
f2 = ( Bn + lambdaW*Ac*P*w - Bc)/c_dc;

dYdt = [f1; f2];

end



