function [zero,res,niter,iters]=newton_multidim(f,Jf,x0,tol,nmax)
%NEWTON Find function zeros.
%   ZERO=NEWTON(FUN,DFUN,X0,TOL,NMAX) tries to find the zero ZERO of the
%   continuous and differentiable function FUN nearest to X0 using the Newton
%   method. FUN and its Jacobian JFUN should be function handles accepting
%   real column vector input x and returning a real column vector value and real matrix 
%   value respectively. If the search fails an error message is displayed.
%
%   [ZERO,RES,NITER]= NEWTON(FUN,...) returns the value of the residual in
%   ZERO, the iteration number at which ZERO was computed, and the iteration history ITERS.

x = x0;
niter = 0;
iters = x';
incr = tol+1;
while incr >= tol && niter < nmax
    niter = niter + 1;
    fx = f(x);
    Jfx = Jf(x);
    incr = - Jfx\fx;
    x = x + incr;
    iters=[iters; x'];
    incr = norm(incr);    
end
if niter == nmax && incr >= tol
    fprintf(['newton stopped without converging to the desired tolerance',...
        ' because the maximum number of iterations was reached\n']);    
end
zero = x;
res = f(x);