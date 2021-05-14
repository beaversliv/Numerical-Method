function [zero,res,niter,iters]=chord(f,a,b,x0,tol,nmax,varargin)
% CHORD Chord metod.
%   ZERO=CHORD(FUN,A,B,X0,TOL,NMAX) tries to find the zero ZERO of the 
%   continuous function FUN in the interval [A,B] using the chord method 
%   starting from X0, A < X0 < B. FUN is an inline function which accept 
%   real scalar input x and returns a real scalar value. 
%   If the search fails an errore message is displayed.
%
%   [ZERO,RES,NITER]= CHORD(FUN,...) returns the value of the residual in ZERO
%   and the iteration number at which ZERO was computed.

q=(f(b)-f(a))/(b-a);
x = x0;
niter = 0;
iters = x;
incr = tol+1;
while incr >= tol && niter < nmax
    niter = niter + 1;
    fx = f(x);    
    incr = - fx/q;
    x = x + incr;
    iters=[iters; x];
    incr = abs(incr);    
end
if niter == nmax && incr >= tol
    fprintf(['chord stopped without converging to the desired tolerance',...
        ' because the maximum number of iterations was reached\n']);    
end
zero = x;
res = f(x);