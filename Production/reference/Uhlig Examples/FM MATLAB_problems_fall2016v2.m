%% ECON 200, Fall 2016
%Flavien Moreau, based on Robert Kurtzman's code

%Auxiliary files: hitzero.m , hitzero1.m , 

clear;
clc;
%% Problem 1: Convergence
% a) analytic solution
quad=@(x) x^2 - x - 7
x1=(1+sqrt(1+32))/2;
x2=(1-sqrt(1+32))/2;
[quad(x1) quad(x2)]             %check that they are roots

% b) while loop
% positive solution
tolq = 10^-20;
iterq=1;
zhigh = 20;
zlow = -20;
diff = 1;
z1=1;
while abs(diff) > tolq && iterq < 1000;
    diff=quad(z1);
    if diff> 0;
    zhigh = z1;
    z1 = .5*z1 + .5* zlow;
    else
    zlow = z1;
    z1 = .5*z1 + .5*zhigh;
    end;
    iterq = iterq + 1;
end;
z1
% negative solution
iterq=1;
zhigh = 20;
zlow = -20;
diff = 1;
z2=-2;
while abs(diff) > tolq && iterq < 1000;
    diff=quad(z2);
    if diff< 0;
    zhigh = z2;
    z2 = .5*z2 + .5* zlow;
    else
    zlow = z2;
    z2 = .5*z2 + .5*zhigh;
    end;
    iterq = iterq + 1;
end;
z2

%c) with fsolve
fsolve(@(x) quad(x), -3)
fsolve(@(x) quad(x), 3)

%% Problem 2 Finite Horizon Life-cycle 
%Parameters
delta=.98;
r=.05;
W1=100;
T=40;
W = zeros(1,T+1);
x1 = zeros(1,T+1);
x2 = zeros(1,T+1);
p2=ones(T+1,1)*1.5;
a=.6;
b=.4;

%(Part d)
x11 = 1;            % some initial guess
W(1)=W1;

x1=fsolve(@(z) hitzero(z,W,x1,x2,r,delta,T,p2,a,b),x11);    %find the right initial condition to use up all wealth
[F,W,x1,x2]=hitzero(x1,W,x1,x2,r,delta,T,p2,a,b);

plot(1:length(x1),x1,1:length(x1),x2,1:length(x1),W); legend('x','x2','W');
[x1' x2' W']

% (Part e) Riley Problem - Chapter 6
W = zeros(1,T+1);
x = zeros(1,T+1);
x1 = 2;
W(1)=W1;
x1=fsolve(@(z) hitzero1(z,W,x,r,delta,T),x1);
[F,W,x]=hitzero1(x1,W,x,r,delta,T);

%% Problem 3 (Fast version)
% Initialize
delta=.08;
alpha=1/3;
beta=.96;
tol=10^-10;

% Define a grid
max_k = delta ^ (1 / (alpha - 1));
n = 100; 
k = linspace(0,max_k,n);

tic                     %ask Matlab to start a timer
G=ones(n,1)*(k.^alpha+(1-delta).*k)-k'*ones(1,n);
T=log(max(G,10^-10000));

V3=zeros(n,1);
V=V3;
Vdif=1;
iter=1;
while Vdif > tol && iter < 10000;

    [V3a,I]=max(T+beta*V*ones(1,n)); 
    V3=max(T+beta*V*ones(1,n))';

    Vdif = max(abs(V3-V));

    V=V3;

    iter = iter + 1;

end;
toc %stops Matlab's timer

figure()
plot(k,V)

%% ECON 200 - Problem 3 (Slow version)

Gamma = @(x) x ^ alpha + (1 - delta) * x;
% Solve VFI
tic
for i=1:length(k);
    for j=1:length(k);
        R(i,j) = -Inf;
        if(k(j) <= Gamma(k(i)))
        R(i,j)=log(k(i)^alpha+(1-delta)*k(i)-k(j));
        end;
    end;
end;


V2=zeros(length(k),1);
V=V2;
Vdif=1;
iter=1;
while Vdif > tol && iter < 10000;

    for j=1:length(k);
        [V2(j),kj]=max(R(j,:)+beta*V');
    end;

    Vdif = max(abs(V2-V));
    V=V2;
    iter = iter + 1;

end;
toc