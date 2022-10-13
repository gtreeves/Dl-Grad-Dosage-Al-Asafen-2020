function [error] = SSE(theta, xgene, sigma, u_space, xspan)

    a1 = interp1(u_space(:,1), xspan, theta, 'spline');
    a2 = interp1(u_space(:,2), xspan, theta, 'spline');
    a3 = interp1(u_space(:,3), xspan, theta, 'spline');
    
    e1 = (xgene(1) - a1)/sigma(1);
    e2 = (xgene(2) - a2)/sigma(2);
    e3 = (xgene(3) - a3)/sigma(3);

    error = sqrt(e1^2 + e2^2 + e3^2);

end