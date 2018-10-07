clc
clear all
% close all

%load('not_cons')
%load('CVX_1D')
%load('m_dt')
load('main')

figure;
plot(dd,bb,'g')

figure;
plot(dd,muhat,'b')

% X_1=zeros(N,m); X_1(:,1)=p;
% 
% for k=1:m-1
%     phi_1=0;
%     
%     for lambda=1:DN
%         
%         for j=1:N
%             
%             % reproduction with a1
%             dist_1=abs(X_1(:,k)-X_1(j,k));
%             dist11=dist_1<((lambda)*deltat);
%             dist21=dist_1>((lambda-2)*deltat);
%             
%             dist_1=(dist11.*dist21)*a1(lambda);
%             
%             phi_1=phi_1+dist_1.*(X_1(j,k)-X_1(:,k));
%         end
%         
%     end
%     
%     X_1(:,k+1)=X_1(:,k)+dt*(1/N)*phi_1;
% end
% 
% 
% t_a=linspace(0,T,m);
% 
% figure(3)
% plot(t_a,X_1,'b')
% %title('Evolution with a1')
% %xlabel('x')
% %ylabel('y')
% axis([0 T -1 1])