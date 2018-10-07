function [w1] = dynamics_a_ijt(w0,u,n)
% INPUT
% w0 : initial position of the particles
% u : control of the dynamic
% n : variable of the current time
% OUTPUT
% w1 : dynamic of the particles

global dt

w     = w0; % initial
N     = length(w0); 

summ= 0;
for j=1:N
    summ= summ+u(j,:)'.*(w(j)-w);
end

w1 = (dt/N).*summ;

end