clc
clear all
close all

load('i1')
xd=0.4;

%% SOLUTION OF PDE
T=400;
NN=200;

dx=(bound_b-bound_a)/NN;
xx=bound_a+(dx/2):dx:bound_b-(dx/2); % cells' centers
xx_g=[xx(1)-dx,xx,xx(end)+dx]; % vector with ghosts cells
x_g=[x(1)-dx,x,x(end)+dx];

vmax=0.5;
% CFL condition
dt=0.1*dx/(2*vmax);

m=ceil(T/dt); %time steps
dt=T/m;

% initial condition
% f=@(x) exp(-((x).^2)/0.01);
% f=@(x) 1/NN(count)+0*x;
f=@(x) ((x+1)/2).*(x<0)+(((1-x)/2)).*(x>=0);

f0=f(xx);
f_c=f(xx); % function evaluated on the centers of cells

f_c=[f_c(1),f_c,f_c(end)];

a=[a(1,:);a;a(end,:)];
speed=zeros(1,m);
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
end
figure(1)
plot(xx,f0,'b',xx,f_c(2:end-1),'r','LineWidth',1)
hold on
plot([0.4,0.4],[0,2.5],'m','LineWidth',0.6)
% title(['N= ',num2str(NN(count)),'  T = ',num2str(T)])
% legend('f(0)','f(T)')
% axis([-1 1 0 1])
print('-depsc2','PDE_bb12')

figure(2)
plot(linspace(0,T,m),speed,'b','LineWidth',1)
print('-depsc2','rcost3')
save('i2')