function [t, u] = beuler(f, dfdy, tspan, y0, N, varargin)
% BEULER Solve differential equations using the backward Euler method.
% [T, U] = beuler(F, DFDY, TSPAN, Y0, N) with TSPAN = [T0 TFINAL]
% integrates the differential equation y' = f(t,y) from
% time T0 to TFINAL with initial condition Y0 using the backward Euler
% method on an equispaced grid of N intervals.
% F and DFDY should be function handles to F(T,Y) and DFDY(T,Y).
% The nonlinear solve at each iteration is carried out using Newton's
% method.

h = (tspan(2)-tspan(1))/N;
t = linspace(tspan(1),tspan(2),N+1);
u = zeros(size(t));
u(1) = y0;
for n = 1:N
    fnewt=@(y)y-u(n)-h*f(t(n+1),y);
    dfnewt=@(y)1-h*dfdy(t(n+1),y);
    u(n+1)=newton(fnewt,dfnewt,u(n),1e-5,100);
end
