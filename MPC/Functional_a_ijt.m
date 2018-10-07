function [val]= Functional_a_ijt(U)
% INPUT
% U : variable of the minimization
% OUTPUT
% val : value of the functional for the minumized U

global winit ti tf Nsteps dt wd gamma N i

[W] = DynEvol_a_ijt(winit,U,i,ti,tf,Nsteps,@dynamics_a_ijt); % evolution of the system

val = (dt/(2*N))*(sum(sum((W-wd).^2))+(gamma/N)*(sum(sum(U.^2)))); % functional J(u,w)