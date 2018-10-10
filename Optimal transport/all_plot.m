clc
clear all
close all

figure(1)
load('i2_ICaij2')
semilogy(linspace(0,T,m),speed,'m--','LineWidth',2)
hold on

clear all
load('i2l2')
semilogy(linspace(0,T,m),speed,'b:','LineWidth',2)
hold on

clear all
load('i2')
semilogy(linspace(0,T,m),speed,'r','LineWidth',2.5)
hh=legend('$V(t)$ control problem','$V(t)$ with $a\in\L^{2}$','$V(t)$ with $a\in\L^{\infty}$');
set(hh,'FontSize',19,'interpreter','latex','Location','SouthWest')
axis([0 T 1e-5 1])
print('-depsc2','allcost')