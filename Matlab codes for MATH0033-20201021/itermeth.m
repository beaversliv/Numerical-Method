function [x, niter, relresiter, xiter] = itermeth(A,b,x0,nmax,tol,P)
%ITERMETH    Stationary iterative method
%   X = ITERMETH(A,B,X0,NMAX,TOL,P) attempts to solve the system of
%   linear equations A*X=B for X.
%   [X, ITER] = ITERMETH(A,B,X0,NMAX,TOL,P) additionally returns the number
%   of iterations ITER.
%   [X, ITER, RES] = ITERMETH(A,B,X0,NMAX,TOL,P) additionally returns the
%   vector RES of residuals at each iteration.
%   The N-by-N coefficient matrix A must be not singular
%   and the right hand side column vector B must have length N.
%   If P='J' the Jacobi method is used.
%   If P='G' the Gauss-Seidel method is used. 
%   If P=alpha is a scalar then the stationary Richardson method is used.
%   Otherwise, P is a N-by-N matrix that plays the role of a preconditioner. 
%   If no P is supplied then the basic stationary method is used, i.e. P=I.
%   TOL specifies the tolerance of the method. The method stops when the 
%   relative residual ||r^k||/||b|| falls below TOL.
%   NMAX specifies the maximum number of iterations.
n=size(A,1);
if nargin == 6
    if ischar(P)==1
        if P=='J'
            P = diag(diag(A));
        elseif P == 'G'
            P = tril(A);
        end
    elseif isscalar(P)
        P = (1/P)*eye(n);
    end
else
    P = eye(n);
end
N = P-A;
nb = norm(b);
if nb == 0
    nb = 1;
end
niter = 0;
x = x0;
xiter = [];
relresiter = [];
r = b-A*x;
relres = norm(r)/norm(b);
while (relres > tol) && (niter <  nmax)
    niter = niter + 1;
    x = P\(N*x + b);
    r = b - A*x;
    relres = norm(r)/nb;
    relresiter = [relresiter, relres]; 
    xiter = [xiter,x];
end

if (niter >= nmax)
    disp(sprintf('itermeth reached the maximum iteration number without converging.'));
else
    disp(sprintf('itermeth converged in %i iterations.',niter));
end