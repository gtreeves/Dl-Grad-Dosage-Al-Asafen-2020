function [u_space, w_space, u, w] = ftn_run(fhandle, tspan, options, params, xspan, P, Vn, Vc, Ac, m)
 
dl0 = [0.5, 1, 2];

for i = 1:3
     
    % Initial conditions
    u0      = zeros(m,1);
    w0      = dl0(i)*ones(m,1);
    Y0      = [u0 w0];    
    
    % Run
    soln    = ode15s(fhandle, tspan, Y0, options, params, xspan, P, Vn, Vc, Ac, dl0(i));
    Y       = deval(soln, tspan);
    u       = Y(1:m,:);
    w       = Y(m+1:2*m,:);
    
    % Store
    u_space{i}  = u(1:end,end);
    w_space{i}  = w(1:end,end);

end

