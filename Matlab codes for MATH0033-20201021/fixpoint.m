function [p,res,niter,iters]=fixpoint(phi,x0,tol,nmax)
% FIXPOINT Fixed point iteration.
%   P=FIXPOINT(PHI,X0,TOL,NMAX) tries to find the fixed point P of the
%   function PHI using a fixed point iteration starting from X0.
%   PHI should be a function handle which accept real scalar input x
%   and returns a real scalar value. If the search fails an error message is displayed.
%
%   [P,RES,NITER]= FIXPOINT(PHI,...) returns the value of the residual in ZERO,
% the iteration number at which ZERO was computed, and the iteration history ITERS.

x = x0;
niter = 0;
iters = x;
incr = tol+1;
while incr >= tol && niter < nmax
    niter = niter + 1;
    phix = phi(x);
    iters=[iters; phix];
    incr = abs(phix - x);
    x = phix;
end
if niter == nmax && incr >= tol
    fprintf(['fixpoint stopped without converging to the desired tolerance',...
        ' because the maximum number of iterations was reached\n']);
end
p = x;
res = phi(x) - x;