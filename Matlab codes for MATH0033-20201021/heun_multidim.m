function [t, u] = heun_multidim(f, tspan, y0, N)
% HEUN_MULTIDIM Solve systems of differential equations using the improved
% Euler method, often called Heun's method.
% [T, U] = heun_multidim(F, TSPAN, Y0, N) with TSPAN = [T0 TFINAL]
% integrates the system of differential equations y' = f(t,y) from
% time T0 to TFINAL with initial condition Y0 using the Heun
% method on an equispaced grid of N intervals.
% F should be a function handle to F(T,Y).
% The output T is a row vector of time values tn, n=0,...,N and the output 
% U is a PxN matrix, where P is the dimension of the system of ODEs. So 
% the ith row of U represents the time evolution of the jth unknown in the 
% solution, and the jth column of U represents a snapshot of all of the 
% unknowns at time tj.
% The input Y0 should be a column vector of length P, and the function F
% should, given a scalar t and column vector v, return a column vector
% corresponding to f(t,v).

h = (tspan(2)-tspan(1))/N;
t = linspace(tspan(1),tspan(2),N+1);
u = zeros(length(y0),length(t));
u(:,1) = y0;
for n = 1:N
    u(:,n+1)=u(:,n)+h/2*(f(t(n),u(:,n))+f(t(n+1),u(:,n)+h*f(t(n),u(:,n))));
end
