function [error_sna, error_sogv, error_sogd, hard_constraint,u_x0, u_x1, u_max] ...
    = ftn_errors(u_space, xspan)
% This function calculates the errors by minimizing the SSE in the x border
% of gene expression
    
load('Mat/Borders.mat');
fun = @SSE;
options = optimset;

u_space = cell2mat(u_space);

if max(max(u_space)) == 0
   error_sna  = inf;
   error_sogd = inf;
   error_sogv = inf;
else
    [p,  error_sna]  = fminbnd(fun, 0, 2, options, xsna,  sigma_xsna,  u_space, xspan);
    [q, error_sogv]  = fminbnd(fun, 0, 2, options, xsogv, sigma_xsogv, u_space, xspan);
    [r, error_sogd]  = fminbnd(fun, 0, 2, options, xsogd, sigma_xsogd, u_space, xspan);
end



% ------------------------------------
% This is the value of concentration for dl0 = 2, at x = 0, t = end
hard_constraint = u_space(1,3);

% -----------------------------------------------------
% Concentrations at 2 ends and max values and position!
u_x0 = [u_space(1,1) u_space(1,2) u_space(1,3)];
u_x1 = [u_space(end,1) u_space(end,2) u_space(end,3)];
[value,cell]  = max(u_space);
u_max_pos = [xspan(cell(1)), xspan(cell(2)), xspan(cell(3))];
u_max_val = value;
u_max = [u_max_val u_max_pos];
end
