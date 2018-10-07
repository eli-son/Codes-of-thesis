function [x]=DynEvol_a_ijt(x0,u,n,t0,T,Nt,fun)
% INPUT
% x0 : the position of the particles at the time in which we are working
% u : the control for which we have to minimize the functional
% n : the time in which we are working
% t0 : the initial time of the prediction
% T : the final time of the prediction
% Nt : number of steps for the current prediction horizon
% fun : is the dynamic that we have to compute
% OUTPUT
% x : matrix with the computed dynamic
global N

sw = 2; % we are selecting the runge kutta
h  = (T-t0)/Nt; % time steps

x(:,1) = x0;
for k=1:Nt
    ni = n+k-1;
    switch sw % switch variable
        
        case 1 % Euler Explicit
            
            e0   = x(:,k);
            cont = u(:,k);
            x(:,k+1) = x(:,k)+h*feval(fun,e0,cont,ni);
            
        case 2  % Runge - Kutta O(4)
            
            e0=x(:,k);
            cont=u(:,(k-1)*N+1:k*N);
            % cont=u(:,k);
            p1=feval(fun,e0,cont,ni);
            p2=feval(fun,e0+(h/2).*p1,cont,ni);
            p3=feval(fun,e0+h/2*p2,cont,ni);
            p4=feval(fun,e0+h*p3,cont,ni);
            x(:,k+1)=e0+h*(1/6*p1+1/3*p2+1/3*p3+1/6*p4);
    end
end