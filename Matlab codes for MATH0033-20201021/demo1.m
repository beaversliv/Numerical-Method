% D. Hewett 
% MATH0033 Numerical Methods 
% Demo 1 - nonlinear equations
%% Setup
close all, clear all, clc
format long, format compact     % Display all digits and condense output in the command window
tol=1e-4;                           % Tolerance for stopping criteria
nmax=20;                           % Maximum number of iterations (in case stopping criterion is not achieved)
fs=16;                              % Setting font size in plots (the default is a bit small)
set(groot,'defaulttextfontsize',fs);
set(groot,'defaultaxesfontsize',fs);
set(groot,'defaultLineLineWidth',2) % Set line width in plots (the default is a bit thin)
set(groot,'defaultContourLineWidth',2) % Set line width in contour plots (the default is a bit thin)
set(0,'DefaultLegendAutoUpdate','off') % Stops Matlab adding extra legend entries automatically
%% Bisection
% Example 3.2.3
f=@(x)sin(2*x)-1+x;                 % Define the function f
x=linspace(-3,3,101);               % Define a set of x values for plotting
figure                              % Create a new figure
plot(x,f(x))          % Plot f
grid on
xlabel('x')
ylabel('f(x)')
hold on                             % Hold allows you to overlay things on top of the current plot
[zero,res,niter,itersb]=bisection(f,-1,1,tol,nmax)    % Run the bisection method with a=-1, b=1, and tol and nmax as specified above
%scatter(zero,res)                  % Plot the approximate root outputted by bisection (i.e. the last iterate)
for i=1:10
    scatter(itersb(i,1),f(itersb(i,1)),'r')   % Plot the first 10 iterates from bisection
    pause                                     % pausing after each one
end
% Exercise - play with nmax, tol, a, b to see how this affects things
%% Fixed point iterations
% Example 3.3.8
phi1=@(x)1-sin(2*x);                % Define the two fixed point functions from Example 3.3.8
phi2=@(x).5*asin(1-x);
figure
x=linspace(0,2,101);
plot(x,phi1(x),'-k')
hold on
plot(x,phi2(x),'--b')
plot(x,x,'-.r')
grid on
xlabel('x')
legend('\phi_1(x)','\phi_2(x)','x','Location','NorthWest')
x0=1;
[fixp,res,niter,itersfp1]=fixpoint(phi1,x0,tol,nmax)   % Run the fixed point iteration for phi1 starting from x0
[fixp,res,niter,itersfp2]=fixpoint(phi2,x0,tol,nmax)   % Run the fixed point iteration for phi2 starting from x0
for i=1:niter+1
    % Comment out one of the following two lines to plot one or other of
    % the sets of iterates:
    
    scatter(itersfp1(i,1),itersfp1(i,1),100,'r')     % Not convergent
    %scatter(itersfp2(i,1),itersfp2(i,1),100,'r')      % Convergent
    pause
end
%% Newton method
% Example 3.6.1
df = @(x)1/2-cos(x)                                % Define the derivative df/dx (Newton needs this!)
x0=0.7;
[zero,res,niter,itersn]=newton(f,df,x0,tol,nmax)    % Run Newton 
%% Chord method
% Example 3.6.1
[zero,res,niter,itersc]=chord(f,-1,1,x0,tol,nmax)   % Run the chord method
%% Secant method
x1=chord(f,-1,1,x0,tol,1)                          % For secant we need two initial values x0 and x1. Generate x1 using one step of the chord method
[zero,res,niter,iterss]=secant(f,x0,x1,tol,nmax)   % Run the secant method
%% Comparison of convergence
zeroex=newton(f,df,x0,1e-15,100) % Generate an "exact" solution. Actually a highly accurate Newton approximation, since we can't find the root analytically.
errb=abs(itersb-zeroex);           % Compute the error in each of the bisection iterates (note that itersb is an array, and hence so too is errb).
errfp1=abs(itersfp1-zeroex);       % Compute the error in the convergent fixed point method with phi1
errfp2=abs(itersfp2-zeroex);       % Compute the error in the convergent fixed point method with phi2
errn=abs(itersn-zeroex);           % Compute the error in the Newton iterates 
errc=abs(itersc-zeroex);           % Compute the error in the chord iterates
errs=abs(iterss-zeroex);           % Compute the error in the secant iterates
figure
semilogy(errb,'-k')  % Plot bisection error against iteration number on a semilogy scale, i.e. x axis is normal but y axis is logarithmic
hold on
semilogy(errfp1,'--c') % Plot fixed point error for phi1
semilogy(errfp2,'--b') % Plot fixed point error for phi2
semilogy(errn,':m')   % Plot Newton error
semilogy(errc,'-.r')  % Plot fixed point error
semilogy(errs,'-g')   % Plot secant error
xlabel('iteration')
ylabel('error')
legend('bisection','fixed point \phi_1','fixed point \phi_2','Newton','chord','secant','Location','SouthEast')
ylim([1e-10,10])  % Set y axis limits
%% Multi-dim Newton
% Example 3.7.9
f1 = @(x1,x2)x1.^2-2*x1.*x2-2;    % Define the first component of f (a scalar-valued function of two variables)
f2 = @(x1,x2)x1+x2.^2+1;          % Define the first component of f (a scalar-valued function of two variables) 
[x1,x2] = meshgrid(-6:.01:2,-4:.01:4);
figure
contour(x1,x2,f1(x1,x2),[0,0],'k'); hold on;
contour(x1,x2,f2(x1,x2),[0,0],'b--')
legend('f_1(x_1,x_2)=0','f_2(x_1,x_2)=0','Location','South')
f=@(x)[f1(x(1),x(2));f2(x(1),x(2))];            % Define the (vector-valued) function f
Jf = @(x)[2*x(1)-2*x(2),-2*x(1);1,2*x(2)];       % Define Jacobian matrix
[zero,res,niter,itersn1]=newton_multidim(f,Jf,[-1;0],tol,nmax)
for i=1:niter+1
    scatter(itersn1(i,1),itersn1(i,2),100,'r')
    pause
end
[zero,res,niter,itersn2]=newton_multidim(f,Jf,[-5;-3],tol,nmax)
for i=1:niter+1
    scatter(itersn2(i,1),itersn2(i,2),100,'b')
    pause
end