clc
clear all
close all
load('main')
N=length(dd); % number of agents
T=10;

p=(1-2*rand(N,1)); % N random points in [-1,1]

% interaction kernel
a=bb;

summ=@(x) 0;
for j=1:N
    summ=@(x) summ(x)+a(j).*(x(j)-x);
end
f=@(t,x) (1/N)*summ(x);

% solution of the ODE with rk4
[t,X] = rk4(f,[0,T],p);

figure(1)
plot(t,X)
%title('Evolution with a(x)')
%xlabel('x')
%ylabel('y')
axis([0 T -1 1])
