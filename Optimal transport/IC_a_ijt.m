clc
clear all
close all

N=80; % number of agents
T=100;

p=(1-2*rand(N,1)); % N random points in [-1,1]
% wanted consensus
xd=0.4;

m=500;
dt = T/(m-1);
t=linspace(0,T,m);
% penalization of the interaction kernel
gamma=0.01;

X_a=zeros(N,m);
X_a(:,1)=p;
for k=1:m-1
    average=sum(X_a(:,k))/N;
    C=zeros(N^2,N);
    b=zeros(N^2,1);

    for p1=1:N % construction of the system to compute a2
        for m1=1:N   
            C(((p1-1)*N+1):(p1*N),m1)=(dt/N)*(X_a(:,k)-X_a(p1,k)).*(X_a(m1,k)-X_a(p1,k));
            C(((p1-1)*N+m1),m1)=C(((p1-1)*N+m1),m1)+gamma;
        end
        b(((p1-1)*N+1):(p1*N),1)=(xd-X_a(p1,k)).*(X_a(:,k)-X_a(p1,k));
    end
    
    a2=zeros(N,N);
    for i=1:N
        cc=C(((i-1)*N+1):(i*N),:);
        bb=b(((i-1)*N+1):(i*N));
        a2(i,:)=cc\bb;
    end

    % dynamics
    s1=0;
    for i=1:N
        s1=s1+a2(:,i).*(X_a(i,k)-X_a(:,k));
    end
    % evolution
    X_a(:,k+1)=X_a(:,k)+((dt)/(N))*s1;
end

figure(1)
plot(t,X_a)
% title('Evolution')
% xlabel('time')
% ylabel('dynamic')
axis([0 T -1 1])

disp(['Initial average: ',num2str(sum(p)/N)])
disp(['Final average: ',num2str(sum(X_a(:,end))/N)])

save('ICaij')