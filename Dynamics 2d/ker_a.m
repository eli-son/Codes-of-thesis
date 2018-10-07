clc
clear all
close all

N=50; % number of particles
p=(rand(2,N)-1)+rand(2,N); p=p'; % N random points in the square [-1,1]x[-1,1]
x0=p(:,1); % x-coord
y0=p(:,2); % y-coord

gamma=0.45; % >=0 parameter of the interaction kernel
% gamma=0.1 : the distances between the particles go to 0 quickly
% gamma=0.9 : the distances between the particles decreases but slowly
% a=@(x) 1/((1+sum(x.^2))^gamma); % interaction kernel a
a=@(x) 1/((1+norm(x,2)^2)^gamma); % interaction kernel a

alpha=0.5; % parameter of the repulsion  
gamma=5; % parameter of the attraction
eps=0.01; % to exclude the singularities


% plot of the initial positions
figure(1)
plot(x0,y0,'b*')
title([num2str(N) ' random points in the square [-1,1]x[-1,1]'])
xlabel('x')
ylabel('y')
axis([-1 1 -1 1])
print('-depsc2','attr1')

%% solution of the ODEs with rk4

T=10; % time horizon
xM = mean(x0);
yM = mean(y0);

% definition of the ODE function
summ=@(x) 0;
for j= 1:N
    q=p-p(j,:); q1=q.^2; % N x 2
    q2=sum(q1,2); q3=q2.^(1/2); % N x 1
    summ=@(x) summ(x)+a(q3)*[(x(j)-x(1:N));(x(N+j)-x(N+1:2*N))];
end
f=@(t,x) (1/N)*summ(x);

% solution of the ODE with rk4
[t,x_y] = rk4(f,[0,T],[x0;y0]);
% % solution of the ODE with ode45
% [tout,X_Yode45] = ode45(f,[0,T],[x0' y0']);

% X_Y : matrix where we save the solutions
% note that the particles converge to the mean of their initial
% positions

% plot of the evolution of the positions withh rk4
px=x_y(1:N,:);
py=x_y(N+1:2*N,:);

figure(2)
for i=1:N
    plot(px(i,:),py(i,:),'b.')
    title('Evolution')
    axis([-1 1 -1 1])
    hold on
end
print('-depsc2','attr2')
% plot(px,py,'b.')
% title('Evolution')
% xlabel('x')
% ylabel('y')
% axis([-1 1 -1 1])



figure(3)
plot(px(:,end),py(:,end),'b*')
hold on
plot(xM,yM,'r*')

% plot(X(:,end),Y(:,end),'b*')
title('Final configuration')
xlabel('x')
ylabel('y')
axis([-1 1 -1 1])
print('-depsc2','attr3')

save('ker_a','N','p','x_y','a')

% % plot of the evolution of the positions with ode45
% pxode45=X_Yode45(:,1:N);
% pyode45=X_Yode45(:,N+1:2*N);
% figure(4)
% plot(pxode45,pyode45,'b.')
% title('Evolution ode45')
% xlabel('x')
% ylabel('y')
% axis([-1 1 -1 1])
% 
% figure(5)
% plot(pxode45(end,:),pyode45(end,:),'b*')
% hold on
% plot(xM,yM,'ro')
% % plot(X(:,end),Y(:,end),'b*')
% title('Final configuration ode45')
% xlabel('x')
% ylabel('y')
% axis([-1 1 -1 1])