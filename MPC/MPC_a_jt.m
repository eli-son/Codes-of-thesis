clc
close all
clear all
warning off

global  N gamma wd T dt winit ti tf Nsteps Ntot i
N = 20;    % total number vertices
gamma = 0.01;
wd = 0.4;

% Time evolution:

T = 100;
t0 = 0;

Ntot =500;
dt = T/(Ntot-1);

% MPC steps:
Np = 1;  % time horizon

NH = Np*dt; % length of the prediction horizon
ti = t0;
tf = ti+NH; % [ti,tf] time frame on which the prediction is computed
% Nsteps = 1; % number of time steps on the prediction horizon
Nsteps = round((tf-ti)/dt); % numb of time steps on the prediction horizon
U = [];% storage variable for control on [ti, tf]
uinit = 0*randn(N,Nsteps);% initialization of the control
u = uinit;

M = 5;  % Control bounds [-50, 50]
lb = -M*ones(1,Nsteps);% lower bound of the control
ub = M*ones(1,Nsteps); % upper bound

% winitt = 1-2*rand(N,1); % initial data w^0 in [-1,0]
winitt = [-0.8675;-0.6219;0.0309;-0.5135;0.1659;-0.9436;-0.9759;-0.7283;0.2222;0.0905;0.5066;-0.5688;-0.7657;-0.8274;-0.1166;-0.1977;0.7022;-0.7994;0.0992;0.5887];
WW(:,1) = winitt; % storage matrix for the solution of the optimization
winit = winitt;

% evolution of the control dynamics
for i=1:Ntot-1
    options = optimset('Display','final-detailed','Algorithm','SQP');
    options.MaxFunEvals = 100000; options.MaxIter = 100;
    
    [u,Fval,exitflag] = fmincon(@Functional_aj,uinit,[],[],[],[],lb,ub,[],options);
    
    [W1] = DynEvol_aj(winit,u.*ones(N,Nsteps),i,ti,ti+dt,Nsteps,@dynamics_aj); % evolution of the system where we apply the control u
    
    winit = W1(:,2); % we save the solution at the time step ti+dt (control horizon [ti, ti+dt])
    WW = [WW winit]; % store the data inside of the matrix WW
    
    %% store of the control in U
    U = [U u(:,1)];
    uinit = u;
    ti = ti + dt; % update
    tf = ti + NH;
    
end

t=linspace(t0,T,Ntot);
figure
plot(t,WW,'LineWidth',1.5)
% title('Evolution')
% xlabel('time')
% ylabel('dynamic')
plot(t,WW,'LineWidth',1.5)
axis([0 T -1 1])
print('-depsc2','MPC5_ai')
disp(['Initial average: ',num2str(sum(winitt)/N)])
disp(['Final average: ',num2str(sum(WW(:,end))/N)])

save('MPC_ai1')
% for j=1:N
%     
%     plot(WW(j,:),'Color',[0.1 0. 0.55],'LineWidth',1);
%     title('Evolution')
%     xlabel('time')
%     ylabel('dynamic')
%     hold on
%     drawnow
%     % axis([0 Ntot -1 1])
% 
% end