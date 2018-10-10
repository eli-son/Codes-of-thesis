clc
clear all
close all

%% INTERACTION KERNEL

N=100; % number of cells/agents
T=250;
bound_a=-1;
bound_b=1;

dx=(bound_b-bound_a)/N;

x=linspace(bound_a,bound_b,N);

vmax=4;
% CFL condition
dt=dx/vmax;
m=ceil(T/dt); %time steps
dt=T/m;

xd=0.4;
% penalization of the interaction kernel
gamma=0.1;
% gamma=1;

C=zeros(N^2,N^2);
b=zeros(N^2,1);
for l=1:N % rows
    for m1=1:N % columns
        C(((l-1)*N+1):(l*N),((l-1)*N+m1))=(x(:)-x(l))*(x(m1)-x(l));
        b(((l-1)*N+1):(l*N))=(N/(N*gamma+dt))*(x(:)-x(l))*(xd-x(l));
    end 
end

cvx_begin
    cvx_precision high
    variable a1(N^2,1)

    minimize ( (norm(a1,2)) )
    subject to
        C*a1 == b
        abs(a1) <= 2
cvx_end

a=reshape(a1,N,N);
a=a';

y=x;
[X,Y]=meshgrid(x,y);
figure(1);
mesh(X,Y,a)
print('-depsc2','aij_l2')
% title('Plot of interaction kernel')

save('i1l2','N','x','a','bound_a','bound_b','T','m','dt')