clc
close all

clear all
load('MPC_aij1')
xd=0.4;
speed=zeros(1,length(t));
for k=1:length(t)
    speed(k)=sum((WW(:,k)-xd).^2);
end
figure(1)
semilogy(t,speed,'b','LineWidth',1.5)
hold on

clear all
load('MPC_aij_3')
xd=0.4;
speed=zeros(1,length(t));
for k=1:length(t)
    speed(k)=sum((WW(:,k)-xd).^2);
end
figure(1)
semilogy(t,speed,'m','LineWidth',1.5)
hold on

clear all
xd=0.4;
load('MPC5_aij')
speed=zeros(1,length(t));
for k=1:length(t)
    speed(k)=sum((WW(:,k)-xd).^2);
end
semilogy(t,speed,'r','LineWidth',1.5)
print('-depsc2','MPCaij_rcost')