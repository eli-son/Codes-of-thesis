clc
clear all
close all

N=50; % number of particles
p=(rand(2,N)-1)+rand(2,N); p=p'; % N random points in the square [-1,1]x[-1,1]
x0=p(:,1); % x-coord
y0=p(:,2); % y-coord

% plot of the initial positions
figure(1)
plot(x0,y0,'b*')
title([num2str(N) ' random points in the square [-1,1]x[-1,1]'])
xlabel('x')
ylabel('y')
axis([-1 1 -1 1])
print('-depsc2','rep1')

% a=2.5 g=15,  a=1.5 g=7 o 5, a=0.5 g=5
alpha=1.5; % parameter of the repulsion  
gamma=5; % parameter of the attraction
eps=0.01; % to exclude the singularities

w=@(x) -(((x+eps).^(alpha))/alpha)+(((x+eps).^(gamma))/gamma);
dw=@(x) -(((x+eps).^(alpha-2)))+((x+eps).^(gamma-2));

x_y=[x0;y0];

%% initial values of the energy
de= 0; e= 0;
dt=0.1;
count=1;
for j=1:N
    qx=x_y(1:N,1)-x_y(j,1); qy=x_y(N+1:2*N,1)-x_y(N+j,1);
    q1x=qx.^2; q1y=qy.^2;
    q2=q1x+q1y; q3=q2.^(1/2); % q3 is the vector of the distances |X-Xj|
    q5=[qx;qy];
    de= de+dw([q3;q3]).*q5;
    e= e+sum(w(q3));
end
x_y(:,count+1)= x_y(:,count)-dt*(1/(2*N))*de;
e1=(1/(2*(N^2)))*e;
e0=e1; % initial values

count=count+1;
m=100; % number of time steps
while (count<m) && (e1<=e0)
    de= 0; e= 0;    
    for j=1:N
        qx=x_y(1:N,count)-x_y(j,count); qy=x_y(N+1:2*N,count)-x_y(N+j,count);
        q1x=qx.^2; q1y=qy.^2;
        q2=q1x+q1y; q3=q2.^(1/2);
        q5=[qx;qy];
        de= de+dw([q3;q3]).*q5; % definition of the ODE function
        e= e+sum(w(q3));
    end

    % GRADIENT METHOD
    x_y(:,count+1)= x_y(:,count)-dt*(1/(2*N))*de;
    
    e0=e1;
    e1=(1/(2*(N^2)))*e;
    c=count;
    % if we are not minimizing the interaction energy, we reduce the time step
    while (e1>e0) && (dt>0)
        dt = dt/2;
        de= 0; e= 0;    
        for j=1:N
            qx=x_y(1:N,c)-x_y(j,c); qy=x_y(N+1:2*N,c)-x_y(N+j,c);
            q1x=qx.^2; q1y=qy.^2;
            q2=q1x+q1y; q3=q2.^(1/2);
            q5=[qx;qy];
            de= de+dw([q3;q3]).*q5; % definition of the ODE function
            e= e+sum(w(q3));
        end
        % GRADIENT METHOD
        x_y(:,c+1)= x_y(:,c)-dt*(1/(2*N))*de;
        e1=(1/(2*(N^2)))*e;
        c=c+1;
    end
    
    count=count+1;
end

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
print('-depsc2','rep2')

figure(3)
plot(px(:,end),py(:,end),'b*')
title(['Final configuration, alpha=' num2str(alpha) ', gamma=' num2str(gamma)])
xlabel('x')
ylabel('y')
axis([-1 1 -1 1])
print('-depsc2','rep3')
%% plot of the graphs of W and dW

% maximum initial distance between two particles
dist=zeros(N);
for j=1:N
    q=p-(ones(N,1)*p(j,:)); q1=q.^2; q2=sum(q1,2); q3=q2.^(1/2);
    dist(:,j)=q3; % matrix of distances between the particles (symmetric)
end
R=max(max(dist)); % maximum distance between two particles

DN=50; % discretization of the interval [0,2R]
supp=linspace(0,2*R,DN);

% the interaction kernel that we used
vect_dw=zeros(DN,1); vect_w=zeros(DN,1);
for s=1:DN
    vect_dw(s)=dw(supp(s));
    vect_w(s)=w(supp(s));
end
figure(4)
plot(supp,vect_w,'b',supp,vect_dw,'r','LineWidth',1)
title('Interaction kernel')
legend('W','dW')
print('-depsc2','rep4')
save('ker_w','N','p','x_y','dw','m','dt')