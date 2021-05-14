% D. Hewett MATH0033 Numerical Methods 
% Demo 3 - differential equations
%% Setup
close all, clear all, clc
format long, format compact
tol=1e-8;
nmax=100;
fs=16;
set(0,'defaulttextfontsize',fs);
%% 
% The next three sections define three different cases. Run one of these 
% sections and then run either the "Plotting solutions" or "Error
% comparison" sections.
%% Case1 - Model problem for absolute stability
f=@(t,y)-y;     % lambda = -1
titlestring = 'Case 1: y^\prime(t) =-y';
dfdy=@(t,y)-1;
y0=1;
y_ex=@(t)exp(-t);
% To investigate convergence, comment out line 27 and take N larger in line 
% 22 to see better approximations:
N=80; tmax=5; h=tmax/N;
% To investigate stability to rounding errors and the link to
% absolute stability, comment out line 22 and experiment with line 27, to
% see how big can you take h before you get blow-up in forward Euler? 
% The critical value from theory is h=2 for this example.
%h=1; N=10; tmax=N*h;
%% Case2 - Example 5.3.1
f=@(t,y)-t.*y.^2;
titlestring = 'Case 2: y^\prime(t) =- ty^2';
dfdy=@(t,y)-2*t.*y;
y0=2;
y_ex=@(t)2./(1+t.^2);
% To investigate convergence comment out line 41 and take N larger in line 
% 36 to see better approximations:
N=10; tmax=7; h=tmax/N;
% To investigate stability to rounding errors and the link to
% absolute stability, comment out line 36 and experiment with line 41, to
% see how big can you take h before you get blow-up in forward Euler? 
% The critical value seems to be h=0.7 for this example.
%h=.5; N=10; tmax=N*h;
%% Case3 - Example 5.6.1
c=-3;
f=@(t,y)(cos(t)+c)*y;
titlestring = ['Case 3: y^\prime(t)=(' num2str(c) '+cos(t)) y'];
dfdy=@(t,y)cos(t)+c;
y0=1;
y_ex=@(t)exp(c*t+sin(t));
% To investigate convergence comment out line 57 and take N larger in line 
% 51 to see better approximations:
N=40; tmax=12; h=tmax/N;
% To investigate stability to rounding errors, comment out line 51 and 
% experiment with line 57. For c=-3 the theoretical critical value is h=0.5. 
% But numerically the critical value seems to be h=0.7. This is not a 
% contradiction - the theoretical value provided a sufficient but not 
% necessary condition for stability.
%h=.2; N=50; tmax=N*h;
%% Plotting solutions
[t_fe, u_fe] = feuler(f,[0,tmax],y0,N);
[t_be, u_be] = beuler(f,dfdy,[0,tmax],y0,N);
[t_cn, u_cn] = cn(f,dfdy,[0,tmax],y0,N);
[t_he, u_he] = heun(f,[0,tmax],y0,N);
t=tmax*(0:0.01:1);     % t values for plotting exact solution
y=y_ex(t);
figure
%plot(t_be, u_be','go-') plot(t,y,'k',t_fe,u_fe,'ro-',t_be,
%u_be','go-',t_cn,u_cn,'bo-')
plot(t,y,'k',t_fe,u_fe,'ro-',t_be, u_be','go-',t_cn,u_cn,'bo-',t_he,u_he,'mo-')
%legend('exact','f Euler','b Euler','CN')
legend('exact','f Euler','b Euler','CN','Heun')
title(titlestring)
%% Error comparison
% Here I am computing the error at the final time tmax.
% Alternatively you could compute the max error over the whole
% interval, by e.g.
%      err_fe(i) = max(abs(y_ex(t_fe) - u_fe));
% or the root-mean-square error by e.g.
%      err_fe(i) = norm(y_ex(t_fe) - u_fe));
NN=10;
Nvec=10*2.^(1:NN);
err_fe=zeros(NN,1);
err_be=zeros(NN,1);
err_cn=zeros(NN,1);
err_he=zeros(NN,1);
ytmax=y_ex(tmax);
for i=1:NN
    N=Nvec(i);
    [t_fe, u_fe] = feuler(f,[0,tmax],y0,N);
    [t_be, u_be] = beuler(f,dfdy,[0,tmax],y0,N);
    [t_cn, u_cn] = cn(f,dfdy,[0,tmax],y0,N);
    [t_he, u_he] = heun(f,[0,tmax],y0,N);
    err_fe(i) = abs(ytmax - u_fe(end));
    err_be(i) = abs(ytmax - u_be(end));
    err_cn(i) = abs(ytmax - u_cn(end));
    err_he(i) = abs(ytmax - u_he(end));
end
hvec=1./Nvec;
figure
loglog(hvec,err_fe,'r',hvec,err_be,'g',hvec,err_cn,'b',hvec,err_he,'m')
grid on
xlabel('h')
ylabel('error')
legend('f Euler - O(h)','b Euler - O(h)','CN - O(h^2)','Heun - O(h^2)','Location','SouthEast')
title(titlestring)
%% Systems
% Let's solve the system from Example 6.8.1, modelling the
% predator-prey interactions between a group of foxes and a group
% of rabbits. We're going to use Heun's method. First we'll
% implement the method inline in our script. Then we'll move the
% implementation to a function file heun_multidim.m. The two
% approaches give exactly the same result, but the latter is
% slightly more elegant.
f=@(t,y)[0.08*y(1) - 0.004*y(1)*y(2);-0.06*y(2) + 0.002*y(1)*y(2)];
y0=[40;20]; tmax=120;
M=20;
t = linspace(0,tmax,M+1);
h = tmax/M;
u = zeros(2,M+1);
u(:,1) = y0;
for n = 1:M
    % Heun
    u(:,n+1)=u(:,n)+h/2*(f(t(n),u(:,n))+f(t(n+1),u(:,n)+h*f(t(n),u(:,n))));
    % FE (in case you want to play with that)
%    u(:,n+1)=u(:,n)+h*(f(t(n),u(:,n)));
end
figure
plot(t,u(1,:),'bx-', t,u(2,:),'rx-')
legend('rabbits','foxes')
xlabel('time (months)')
ylabel('population')
title('Created without heun-multidim.m')
%hold on
[t_he, u_he] = heun_multidim(f, [0,tmax], y0, M);
figure
plot(t_he,u_he(1,:),'bx-', t_he,u_he(2,:),'rx-')
legend('rabbits','foxes')
xlabel('time (months)')
ylabel('population')
title('Created using heun-multidim.m')
%% Some hints/tips for computational homework 3
% For homework 3 you have to solve a IBVP for a PDE using the
% method of lines, which reduces the problem to the solution of a
% system of IVPs. You are asked to solve this system using FE, BE
% and CN. To do this, you have two choices: either implement the
% methods directly in your code, or write function files to do
% this. (These two choices are illustrated above for heun)
%
% The first approach is the simpler in the sense that everything is
% in a single script. So if you are a beginner, this is probably the way to
% go.

% The second is more elegant, but more advanced, as it would require you 
% to spend some time writing separate function files. 
% In particular, you would need to modify feuler.m, beuler.m and cn.m (which
% currently only work for scalar IVPs) to work for systems of
% IVPs, following the example in heun_multidim.m. But I only
% recommend doing this if you are confident in using function
% files, and if you can work out what to supply as the Jacobian
% for the implicit methods. (Note that in these implicit methods
% you will need to replace newton.m with newton_multidim.m.) Note
% that since the RHS of the system f is linear, this Newton
% business is actually unnecessary for this problem - the update
% step for the implicit methods just requires the solution of a
% LINEAR system of equations. So if you wanted, you could write
% special files that only work for linear f.
%% Practice with nested loops: 
N=5; M=4;
U=zeros(M+1,N-1);
for m=1:M+1
    for n=1:N-1
        U(m,n)=sqrt(m^2+3*n^2+5);
    end
end
% Practice with array indexing: Reference the mth row (a row
% vector) by A(m,:)
v=rand(1,N-1);
w=rand(M+1,1);
for m=1:M+1
    U(m,:)=m*v;
end
% Reference the nth column (a column vector) by A(:,n)
for n=1:N-1
    U(:,n)=n*w;
end