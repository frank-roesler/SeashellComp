function val = u_D(r,theta,alpha,r_outer,R)

% Sets up an extension to the annulus [r_outer, R] of the boundary data exp(i*alpha*theta)

    val = zeros(length(r),length(alpha));
    [ind, y] = find(r>r_outer);
    val(ind,:) = abs(r(ind)-r_outer)/(R-r_outer).*exp(1i*alpha.*theta(ind));
    val = val/sqrt(2*pi*R);
end