function [tout,yout] = rk4(odefun,tspan,y0,options)
% INPUT:
% odefun: function of the ODE that we have to solve
% tspan: [0,T] time interval
% y0: initial conditions
% dt: time step
% OUTPUT:
% tout: vector of the times (column vector)
% yout: matrix of the solutions at every time : in the rows there are the
%       partcles and in the columns, the evolutions in time of each
%       particle
% global m dt

%m=100;
%dt = tspan(2)/(m-1);
% dt=0.05;
dt=0.1;
m = ceil(tspan(2)/dt); % number of time steps

% Runge-Kutta of order 4
a(2,1)=1/2;
a(3,2)=1/2;
a(4,3)=1;
% a(4,4)=0;

c=[0 1/2 1/2 1];
b=[1/6 1/3 1/3 1/6];

n = 1; % time position
tout(n) = tspan(1); % sets the initial time

yout(:,n) = y0(:); % initial data in the first column

% the xi are the stages (we have 4 stages)
xi = zeros(length(y0),4);

% we work until we reach the time horizon T
for iter=1:m-1
    % stage1: set the initial data
    xi(:,1) = yout(:,n);
    % stage2
    xi(:,2) = yout(:,n)+dt*a(2,1)*odefun(tout(n)+c(2)*dt,xi(:,1)); % forse c(2)
    % stage3
    xi(:,3) = yout(:,n)+dt*(a(3,1)*odefun(tout(n)+c(2)*dt,xi(:,1))+a(3,2)*odefun(tout(n)+c(3)*dt,xi(:,2)));
    % stage4
    xi(:,4) = yout(:,n)+dt*(a(4,1)*odefun(tout(n)+c(2)*dt,xi(:,1))+a(4,2)*odefun(tout(n)+c(3)*dt,xi(:,2))+a(4,3)*odefun(tout(n)+c(4)*dt,xi(:,3)));
    
    % update of the solution at the n-th step
    yout(:,n+1) = yout(:,n)+dt*(b(1)*odefun(tout(n)+c(1)*dt,xi(:,1))+b(2)*odefun(tout(n)+c(2)*dt,xi(:,2))+b(3)*odefun(tout(n)+c(3)*dt,xi(:,3))+b(4)*odefun(tout(n)+c(4)*dt,xi(:,4)));
    
    tout(n+1) = tout(n)+dt; % update of the time
    n = n+1; % update of the time position
end
save('m_dt','m','dt')