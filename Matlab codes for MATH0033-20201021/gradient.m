function [x, niter, relresiter, xiter] = gradient(A, b, x0, nmax, tol, P)
%
% [X, ITER] = GRADIENT(A, B, X0, N_MAX, TOL, P)
%
%  Preconditioned gradient method (Richardson dynamic method).
%  It attemps to solve the linear system A*X = B for X.
%  P is a N-by-N preconditioner. TOL specifies the
%  tolerance of the method. The method stops when the relative residual ||r^k||/||b|| falls
%  below TOL. NMAX specifies the maximum
%  number of iterations
n=size(A,1);
if nargin < 6
    P = eye(n);
end
nb = norm(b);
if nb == 0
    nb = 1;
end
niter = 0;
x=x0;
xiter = [];
relresiter = [];
r = b-A*x;
relres = norm(r)/nb;
while (relres > tol) && (niter <  nmax)
    niter = niter + 1;
    p  = P\r;
    rho = (r'*p);
    q = A*p;
    alpha = rho / (p'*q);
    x = x + alpha*p;
    r = r - alpha*q;
    relres = norm(r)/nb;
    relresiter = [relresiter, relres];
    xiter = [xiter,x];
end

if (niter >= nmax)
    disp(sprintf('gradient reached the maximum iteration number without converging.'));
else
    disp(sprintf('gradient converged in %i iterations.',niter));
end