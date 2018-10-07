function [val]= Functional_aj(U)
% INPUT
% U : variable of the minimization
% OUTPUT
% val : value of the functional for the minumized U

global winit ti tf Nsteps dt wd gamma N i

[W] = DynEvol_aj(winit,U,i,ti,tf,Nsteps,@dynamics_aj); % evolution of the system

val = (dt/(2*N))*(sum(sum((W-wd).^2))+gamma*(sum(sum(U.^2)))); % functional J(u,w)