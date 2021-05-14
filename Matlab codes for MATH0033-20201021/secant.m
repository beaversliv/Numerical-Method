function [zero,res,niter,iters]=secant(f,x0,x1,tol,nmax)
%SECANT Find function zeros.
%   ZERO=SECANT(FUN,DFUN,X0,TOL,NMAX) tries to find the zero ZERO of the
%   continuous and differentiable function FUN starting from X0 and X1 using the
%   secant method. FUN should be a function handle accepting
%   real scalar input x and returning a real scalar value. If the search fails
%   an error message is displayed.
%
%   [ZERO,RES,NITER]= SECANT(FUN,...) returns the value of the residual in
%   ZERO, the iteration number at which ZERO was computed, and the iteration history ITERS.

xm = x0;
x = x1;
niter = 1;
iters = [xm;x];
incr = tol+1;
fxm=f(xm);
while incr >= tol && niter < nmax
    niter = niter + 1;
    fx = f(x);    
    incr = - fx*(x-xm)/(fx-fxm);
    xm = x;
    fxm=fx;
    x = x + incr;
    iters=[iters; x];
    incr = abs(incr);    
end
if niter == nmax && incr >= tol
    fprintf(['newton stopped without converging to the desired tolerance',...
        ' because the maximum number of iterations was reached\n']);    
end
zero = x;
res = f(x);