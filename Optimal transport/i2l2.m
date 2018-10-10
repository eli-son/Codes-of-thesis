clc
clear all
close all

load('i1l2')
xd=0.4;

T=250;
NN=100;
% T=7;
% NN=80;
%a=ones(NN);
dx=(bound_b-bound_a)/NN;
xx=bound_a+(dx/2):dx:bound_b-(dx/2); % cells' centers
xx_g=[xx(1)-dx,xx,xx(end)+dx]; % vector with ghosts cells
x_g=[x(1)-dx,x,x(end)+dx];

vmax=2;
% CFL condition
dt=1*dx/(2*vmax);
m=ceil(T/dt); %time steps
dt=T/m;

% initial condition
f=@(x) ((4)/3).*(x<0).*(x>-0.75)+(((1-x)/2)).*(x>=0).*(x<0);

f0=f(xx);
f_c=f(xx); % function evaluated on the centers of cells

f_c=[f_c(1),f_c,f_c(end)];

a=[a(1,:);a;a(end,:)];
speed=zeros(1,m);
U = zeros(m,N);
U(1,:)=f_c(2:N+1);
funz=f_c;
for k=1:m
    intf=interp1(xx,f_c(2:end-1),x,'pchip');
    if sum(intf)==Inf || sum(isnan(intf))~=0
        f_c=funz;
        break
    end
    s=zeros(N+2,1);
    for i=1:N+2
        for j=1:N
            s(i)=s(i)+(dx*a(i,j)*(x(j)-x_g(i))*intf(j));
        end
    end
    ints=interp1(x_g,s',xx_g,'pchip');

    for i=2:NN+1
        if ints(i)>0
            fp=f_c(i);
        else
            fp=f_c(i+1); 
        end
        if ints(i-1)>0
            fm=f_c(i-1);
        else
            fm=f_c(i);
        end
        f_c(i)=f_c(i)-(dt/dx)*(ints(i)*fp-ints(i-1)*fm);
    end
    f_c(1)=f_c(2); f_c(end)=f_c(end-1);
    funz=f_c;
    speed(k)=dx*sum(((xx-xd).^2).*f_c(2:end-1));

    U(k,:) = f_c(2:N+1);
end

figure(1)
plot(xx,f0,'b:',xx,f_c(2:end-1),'r--','LineWidth',1.5)
% title(['t = ' num2str(dt*k)])
hold on
plot([0.4,0.4],[0,100],'m','LineWidth',0.6)
axis([-1 1 0 max(f_c)*(1.1)])
hold off
hh=legend('$u_0(x)$','$u(T,x)$','$\overline x$');
set(hh,'FontSize',19,'interpreter','latex','Location','NorthWest')
print('-depsc2','PDE3_')

figure(2)
plot(linspace(0,T,m),speed,'b','LineWidth',1)
print('-depsc2','rcost3')
save('i2l2')

figure(3)
t = (1:m)*dt;
[XX,TT]=meshgrid(x,t);
contourf(XX,TT,U,100,'LineStyle','none')
colorbar
hold on
plot(xd*ones(1,m),linspace(0,150,m),'r--','LineWidth',2)
print('-depsc2','density3')