function [zero,res,niter,iters]=bisection(f,a,b,tol,nmax)
%BISECTION Find function zeros.
%   ZERO=BISECTION(F,A,B,TOL,NMAX) tries to find a zero ZERO of the continuous
%   function F in the interval [A,B] using the bisection method. F should
%   be a function handle accepting real scalar input x and return a real
%   scalar value. If the search fails an error message is displayed.
%
%   [ZERO,RES,NITER,ITERS]= BISECTION(FUN,...) returns the value of the residual
%   at ZERO, the iteration number NITER at which ZERO was computed, and the iteration
%   history ITERS. The flag NITER=-1 means the zero was one of the original
%   endpoints supplied (either a or b).

fa=f(a);
fb=f(b);
if fa*fb>0
    error(' The sign of FUN at the extrema of the interval must be different');
elseif fa == 0
    zero = a; res = 0; niter = -1; iters=[];    
    return
elseif fb == 0
    zero = b; res = 0; niter = -1; iters=[];
    return
end
x = 0.5*(a+b);          % Initial guess 
niter = 0; iters=x;     % Initialise the counter niter and the vector of iterates iters
errbnd = 0.5*(b-a);     % Initialise the error bound
while errbnd >= tol && niter < nmax
    fx=f(x);
    if fx==0
        errbnd=0;       % If f(x)=0 this will force the iteration to stop
    else
        if fa*fx < 0
            b = x;
        else
            a = x;
            fa=fx;
        end
        x = 0.5*(a+b);              % New iterate        
        errbnd = 0.5*errbnd;
        niter = niter + 1;
        iters=[iters; x];           % Add new iterate to the vector iters
    end
end
if niter == nmax && errbnd >= tol
    fprintf(['bisection stopped without converging to the desired tolerance',...
        ' because the maximum number of iterations was reached\n']);
end
zero = x;
res = f(x);