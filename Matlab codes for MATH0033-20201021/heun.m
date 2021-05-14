function [t, u] = heun(f, tspan, y0, N)
% HEUN Solve differential equations using the improved Euler method,
% Heun's method.
% [T, U] = heun(F, TSPAN, Y0, N) with TSPAN = [T0 TFINAL]
% integrates the differential equation y' = f(t,y) from
% time T0 to TFINAL with initial condition Y0 using the Heun
% method on an equispaced grid of N intervals.
% F should be a function handle to F(T,Y).

h = (tspan(2)-tspan(1))/N;
t = linspace(tspan(1),tspan(2),N+1);
u = zeros(size(t));
u(1) = y0;
for n = 1:N
    u(n+1)=u(n)+h/2*(f(t(n),u(n))+f(t(n+1),u(n)+h*f(t(n),u(n))));
end
