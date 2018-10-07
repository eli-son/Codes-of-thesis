clc
clear all
load('i2')
figure(1)
semilogx(linspace(0,T,m),speed,'b','LineWidth',1)
hold on

clear all
load('i2l2')
semilogx(linspace(0,T,m),speed,'r','LineWidth',1)
hold on

clear all
load('i2_ICaij')
semilogx(linspace(0,T,m),speed,'m','LineWidth',1)
print('-depsc2','allcost')