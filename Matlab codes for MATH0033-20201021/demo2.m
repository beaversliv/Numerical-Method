% D. Hewett 
% MATH0033 Numerical Methods 
% Demo 2 - linear systems
%% Setup
close all, clear all, clc
format long, format compact
tol=1e-8;
nmax=100;
fs=16;
set(0,'defaulttextfontsize',fs);
set(0,'defaultaxesfontsize',fs);
%% Example 4.8.5
% Setup
A=[2,1;1,3]
b=[1;0];
x_ex=[3/5;-1/5];
eigA=eig(A)
cond(A)
x0=[0;1/2];     
%% Basic stationary method
rhoBSM=max(abs(eig(eye(2)-A)))
[x, niter, relresiter, xiter] = itermeth(A,b,x0,nmax,tol);
normres=norm(relresiter(end))
%% Stationary Richardson method
lmin=min(eigA);
lmax=max(eigA);
alphastar=2/(lmin+lmax)    % Optimal choice of alpha
% You can check analytically that alphastar=2/5 for this example,
% and that the resulting convergence constant
% rho(B_alphastar)=1/sqrt(5) = 0.4472...
[x, niter, relresiter, xiter] = itermeth(A,b,x0,nmax,tol,alphastar);
err_SR=zeros(1,niter);
for n=1:niter, err_SR(n)=norm(xiter(:,n)-x_ex); end
%% Jacobi method
[x, niter, relresiter, xiter] = itermeth(A,b,x0,nmax,tol,'J');
err_J=zeros(1,niter);
for n=1:niter, err_J(n)=norm(xiter(:,n)-x_ex); end
%% Gauss-Seidel method
[x, niter, relresiter, xiter] = itermeth(A,b,x0,nmax,tol,'G');
err_GS=zeros(1,niter);
for n=1:niter, err_GS(n)=norm(xiter(:,n)-x_ex); end
relresGS=relresiter;
%% Gradient method
x0=[0;1/2];     % Gives gradient convergence rate close to theoretical estimate
%x0=[0;0];       % Gives gradient convergence rate faster than theoretical estimate 
%x0=[1;1/2];     % Gives gradient convergence rate much faster than theoretical estimate
[x, niter, relresiter, xiter] = gradient(A,b,x0,nmax,tol);
err_G=zeros(1,niter);
for n=1:niter, err_G(n)=norm(xiter(:,n)-x_ex); end
figure
N=501;
x=linspace(0,1,N); y=linspace(-.5,.5,N);
[xx,yy]=meshgrid(x,y);
Phifun=@(x).5*x'*A*x-x'*b;     % x needs to be a column vector for this to work
Phi=zeros(N);
for i=1:N
    for j=1:N
        X=[xx(i,j);yy(i,j)];
        Phi(i,j)=Phifun(X);
    end
end
surf(xx,yy,Phi,'Linestyle','none')
xlim([x(1) x(end)])
ylim([y(1) y(end)])
view(2)
axis square
hold on
Phiiter=zeros(niter,1);
for i=1:niter
    X=[xiter(1,i);xiter(2,i)];
    Phiiter(i)=Phifun(X);
end
plot3([x0(1) xiter(1,:)],[x0(2) xiter(2,:)],[Phifun(x0); Phiiter],'o-r')
%% Conjugate gradient method
% The conjugate gradient method is already included as an inbuilt Matlab
% command, called "pcg" (for "preconditioned conjugate gradient"). We are
% going to use the unpreconditioned method, hence the "[],[]" where one
% would normally supply preconditioners. We know that for this
% 2-dimensional problem the cg algorithm should converge in 2 steps. This
% is indeed the case.
xcg1 = pcg(A,b,tol,1,[],[],x0)     % 1 step of cg
xcg2 = pcg(A,b,tol,2,[],[],x0)     % 2 steps of cg.
% (Note that pcg doesn't give an option to output the whole iteration
% history, so I have computed x1 and x2 manually by setting nmax=1 and
% nmax=2 respectively.
plot3([x0(1) xcg1(1) xcg2(1)],[x0(2) xcg1(2) xcg2(2)],[Phifun(x0); Phifun(xcg1); Phifun(xcg2)],'o-g')
%% Plot convergence for each method
figure
semilogy(err_SR)
hold on
semilogy(err_J)
semilogy(err_GS)
semilogy(err_G)
xlabel('k')
ylabel('||x^k-x||_2')
legend('Stat Rich','Jacobi','Gauss-Seidel','Gradient')
%% Investigate convergence constants.
% Suppose that the error E^k=C^k*E^0 for some 0<C<1. Then
% log(E^k)=k*log(C)+log(E^0). So the slope of the semilogy plot
% tells us log(C). We can extract this from the data in a very 
% crude way as follows:
Calphastar=(cond(A)-1)/(cond(A)+1)
shift=4;
logC_SR=(log(err_SR(end))-log(err_SR(end-shift)))/shift;
C_SR=exp(logC_SR)
logC_J=(log(err_J(end))-log(err_J(end-shift)))/shift;
C_J=exp(logC_J)
logC_GS=(log(err_GS(end))-log(err_GS(end-shift)))/shift;
C_GS=exp(logC_GS)
logC_G=(log(err_G(end))-log(err_G(end-shift)))/shift;
C_G=exp(logC_G)
% Note that one can also use other tools e.g. polyfit for linear
% regression.
%% Reliability of relative residual as stopping criterion
% First investigate this issue for the same matrix as above, which is
% well-conditioned.
A=[2,1;1,3]
cond(A)
b=[0.5;1];
x_ex=A\b;
x0=[0;0];
[x, niter, relresiter, xiter] = itermeth(A,b,x0,nmax,tol,'G');
err_GS=zeros(1,niter);
for n=1:niter, err_GS(n)=norm(xiter(:,n)-x_ex); end
relresGS=relresiter;
figure
semilogy(err_GS/norm(x_ex))
hold on
semilogy(relresGS)
ylim([1e-8,1e2])
xlabel('k')
legend('relerr','relres')
title('Well-conditioned example')
ratio_wellcond=(err_GS(end)/norm(x_ex))/relresGS(end)
% Now change A to something with worse conditioning. A from
% Example 5.10.1
A=[5,7;7,10]
nmax=1000;
%b=[1;1];
theta=0.9569;
b=[cos(theta);sin(theta)];
x_ex=A\b;
eigA=eig(A)
cond(A)
[x, niter, relresiter, xiter] = itermeth(A,b,x0,nmax,tol,'G');
err=zeros(1,niter);
err_GS=zeros(1,niter);
for n=1:niter, err_GS(n)=norm(xiter(:,n)-x_ex); end
relresGS=relresiter;
figure
semilogy(err_GS/norm(x_ex))
hold on
semilogy(relresGS)
ylim([1e-8,1e2])
xlabel('k')
legend('relerr','relres')
title('Ill-conditioned example')
ratio_illcond=(err_GS(end)/norm(x_ex))/relresGS(end)
%% Optimising over b
% In the previous section we saw that if the matrix is
% ill-conditioned then for certain right-hand sides b the
% relative residual in GS can be much smaller than the relative
% error, meaning that the relative residual is a bad choice for a
% stopping criterion. (Cf. Fig 5.9 in notes.) But how did I know
% that
%        theta=0.9569; b=[cos(theta);sin(theta)];
% gives a particularly large ratio between relative residual and
% relative error, for this particular matrix A? Well, I evaluated
% this ratio for lots of different thetas and chose the one
% giving the largest ratio. I first started looking over the
% whole interval (0,pi), and then focussed in on (0.9,1).
% Question for you: why is it sufficient to just consider b of
% the form
%        b=[cos(theta);sin(theta)];
% i.e. with norm 1?
A=[5,7;7,10]
%A=[2,1;1,3]
eigA=eig(A);
x0=[0;1];
nmax=10;
t=linspace(0,pi,1000);
%t=linspace(0.9,1,10000);
err=zeros(size(t));
relres=zeros(size(t));
for i=1:length(t)
    b=[cos(t(i));sin(t(i))];
    x_ex=A\b;
    [x, niter, relresiter, xiter] = itermeth(A,b,x0,nmax,tol,'G');
    err(i)=norm(xiter(:,end)-x_ex)/norm(x_ex);
    relres(i)=relresiter(end);
end
figure
plot(t,err./relres)
xlabel('t')
ylabel('relerr/relres')
[M,m]=max(err./relres)
tm=t(m)
% Note: this is a very inefficient way of finding the maximum
% because we compute the value of the objective function at lots
% of points we don't need. Better to use an optimization routine
% like fminbnd. But you would first need to write a function file
% that computed relerr/relres for any given value of t; you can
% then feed the function handle to fminbnd.