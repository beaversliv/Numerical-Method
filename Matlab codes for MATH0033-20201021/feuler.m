function [t, u] = feuler(f, tspan, y0, N)
% FEULER Solve differential equations using the forward Euler method.
% [T, U] = feuler(F, TSPAN, Y0, N) with TSPAN = [T0 TFINAL]
% integrates the differential equation y' = f(t,y) from
% time T0 to TFINAL with initial condition Y0 using the forward Euler
% method on an equispaced grid of N intervals.
% Input F should be a function handle to a function F(T,Y).

h = (tspan(2)-tspan(1))/N;
t = linspace(tspan(1),tspan(2),N+1);
u = zeros(size(t));
u(1) = y0;
for n = 1:N
    u(n+1)=u(n)+h*f(t(n),u(n));
end
