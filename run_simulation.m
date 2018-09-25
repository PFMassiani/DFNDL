clc;
clear;
close all;

params = params2();

%% Getting the basic operators
matrices = matrices_linearized_t_indep(params);

NTOT = (matrices.N_s + 1)* (matrices.N_elyte - 1) + matrices.N_elyte - 1;
Nem2 = matrices.N_elyte - 1;
Nsm2 = matrices.N_s - 1;
Ns = matrices.N_s + 1;

%% Initial conditions
init = initial_conditions(matrices,params);
y0 = init.state;

%% Computing the discretized system in state-space form
% [D,exo_bdry,nonlinear] = flatten_matrices(matrices,params);
% dynamic = @(t,x) (full(D(t)) * x  + exo_bdry(t) + nonlinear(t,x));

%% Running the simulation
tf = 0.1;
tspan = [0 tf];
Mass = blkdiag(speye(NTOT - Nem2), speye(Nem2,Nem2));
options = odeset('NormControl','on','MassSingular','yes','Mass',Mass,'AbsTol',1e-9,'RelTol',1e-7);
[t,state] = ode15s(@(t,x)(compute_dynamic(t,x,matrices,params)),tspan,y0,options);
[t,x,r,dl,elyte,us,cs] = postproc(t,state,params,matrices);
params.K = 0;

%% Plots

figure(1);
for k = 1:length(t)
    subplot(2,1,1);
    title(sprintf('time = %0.7f',t(k)));
    plot(x,elyte(k,:),'o-');grid on
    ylabel('CE [mol/m^3]');
%     ylim([990 1010]);
    subplot(2,1,2);
    plot(x,100*(elyte(k,:) - elyte(1,:))./elyte(1,:),'o-'); grid on;
    ylabel('CE relative variation (%)')
	ylim([-100 100]);
    pause(0.01);
end
figure(2)
for k = 1:length(t)
    subplot(2,1,1);
    title(sprintf('time = %0.7f',t(k)));
    plot(x,dl(k,:),'o-');grid on
    ylabel('PHI [V]');
    subplot(2,1,2);
    plot(x,100*(dl(k,:) - dl(1,:))./dl(1,:),'o-'); grid on;
    ylabel('PHI relative variation (%)')
% 	ylim([-1 1]);
    pause(0.01);
end
figure(3)
for k = 1:length(t)
    subplot(2,1,1);
    title(sprintf('time = %0.7f',t(k)));
    plot(reshape(cs(k,:,:),Nem2*Ns,1),'o-'); grid on;
    ylabel('CS [mol/m^3]');
%     ylim([0.5*params.csmax 1.5*params.csmax]);
    subplot(2,1,2);
    plot(100*reshape((cs(k,:,:) - cs(1,:,:))./cs(1,:,:),Nem2*Ns,1),'o-'); grid on;
    ylabel('CS relative variation (%)')
    pause(0.001);
end