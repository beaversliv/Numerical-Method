function [t, u] = cn(f, dfdy, tspan, y0, N, varargin)
% CN Solve differential equations using the Crank Nicolson method.
% [T, U] = CN(F, DFDY, TSPAN, Y0, N) with TSPAN = [T0 TFINAL]
% integrates the differential equation y' = f(t,y) from
% time T0 to TFINAL with initial condition Y0 using the Crank Nicolson
% method on an equispaced grid of N intervals.
% F and DFDY should be function handles to F(T,Y) and DFDY(T,Y).
% The nonlinear solve at each iteration is carried out using Newton's
% method.

h = (tspan(2)-tspan(1))/N;
t = linspace(tspan(1),tspan(2),N+1);
u = zeros(size(t));
u(1) = y0;
for n = 1:N
    fnewt=@(y)y-u(n)-h/2*(f(t(n),u(n))+f(t(n+1),y));
    dfnewt=@(y)1-h/2*dfdy(t(n+1),y);
    u(n+1)=newton(fnewt,dfnewt,u(n),1e-5,100);
end
