%% Analytical Solution for membrane covered electrode using infinite series
% and bessel functions of the first kind. Wakil Sarfaraz 12/08/2014

function [I] = Analytical(N)
ze = 0.5;
zm = 1;
epsilon1 = 0.9214;
epsilon2 = 0.9214;
L = zeros(N+1, N+1);
F = zeros(N+1,1);
F(1) = sqrt(2/pi);



for m = 0:N
    for n = 0:N
        L(m+1,n+1) = (4*n+1)*integral(@(k) (2./k).*(epsilon1*exp(-2*k*ze)-epsilon2*exp(-2*k*zm))...
            ./(1-epsilon1*epsilon2 *exp(-2*k*(zm-ze))-epsilon1*exp(-2*k*ze)+epsilon2*exp(-2*k*zm))...
            .*besselj(2*m+1/2,k).*besselj(2*n+1/2,k), 0 , Inf);

    end
end

A = L+ eye(N+1,N+1);
a = A\F;

I = sqrt(pi/2)*a(1);
