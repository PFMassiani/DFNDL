function [ p ] = params2
%PARAMS2 Summary of this function goes here
%   Detailed explanation goes here

% Values taken from params_LCO.m on Github - FastDFN
p.n.L = 100e-6;
p.s.L = 25e-6;
p.p.L = 100e-6;
p.n.R_s = 10e-6;
p.p.R_s = 10e-6;
p.n.epsilon_e = 0.3;
p.s.epsilon_e = 1;
p.p.epsilon_e = 0.3;
p.n.epsilon_s = 0.6;
p.p.epsilon_s = 0.5;
p.n.a_s = 3*p.n.epsilon_s / p.R_s;
p.p.a_s = 3*p.p.epsilon_s / p.p.R_s;
p.n.D_s = 3.9e-14;
p.p.D_s = 1e-13;
p.n.sigma = 100;
p.p.sigma = 10;
p.t_0 = 0.4;
p.Faraday = 96487;
p.R = 8.314472;
p.alpha = 0.5; 
p.T = 298.15;
p.n.csmax =  3.6e3 * 372 * 1800 / p.Faraday;
p.p.csmax =  3.6e3 *  274 * 5010 / p.Faraday;
p.volt_max = 4.7; %4.1113; %4.7;
p.volt_min = 3.105; %2.6;
p.n.k = 1e-5;
p.p.k = 3e-7;


% Values taken from the simulations sent by Saehong
p.n.kappa = (1.5 + 0.1) / 2; % Mean value of the conductivity of electrolyte [S/m]
p.p.kappa = (1.5 + 0.1) / 2;
p.n.D_e = (3.5 + 2.4)*10^(-10) / 2; % Mean value of the diffusivity [m^2/s]
p.p.D_e = (3.5 + 2.4)*10^(-10) / 2; 

% Value taken from "State-of-Charge Estimation with a
% Doyle-Fuller-Newman Li-ion Battery Model", R. Drummond
p.n.aC = 22*1e3; % [F/m^2]
p.p.aC = 22*1e3;


%Values to find
p.n.dlnfdlnce = 0;
p.p.dlnfdlnce = 0;

%Computing OCP characteristics
% p.OCP_0 = 4;
% p.beta = 0.1; % "Identifiability and Parameter Estimation of the Single Particle Lithium-Ion Battery Model", David Howey
% p.csmax = 30000; % [mol/m^3], "Identifiability and Parameter Estimation of the Single Particle Lithium-Ion Battery Model", David Howey
% p.A = p.epsilon_s * p.L / (4*pi*p.R_s^3 / 3);
% p.Qth = p.epsilon_s * p.L * p.csmax * p.Faraday * p.A;
% p.OCP_slope_wrt_stoechio = -p.epsilon_s * p.L * p.csmax * p.Faraday * p.A * p.beta;% "Identifiability and Parameter Estimation of the Single Particle Lithium-Ion Battery Model", David Howey
% p.OCP_slope = p.OCP_slope_wrt_stoechio / (p.R_s * p.csmax);
%p.OCP_slope = -1e4;

% Parameters appearing in the system
p.n.V_0 = (p.Faraday / p.n.aC) * (p.n.k * (p.n.R_s ^2 / p.n.D_s) * (p.n.a_s*(p.n.csmax^2)*(1 - p.t_0))^p.alpha)^(1/(1 - p.alpha));
p.p.V_0 = (p.Faraday / p.p.aC) * (p.p.k * (p.p.R_s ^2 / p.p.D_s) * (p.n.a_s*(p.p.csmax^2)*(1 - p.t_0))^p.alpha)^(1/(1 - p.alpha));
p.n.i_0_ref = p.n.aC * p.n.D_s * p.n.V_0/ (p.n.a_s * p.n.R_s^2);
p.p.i_0_ref = p.p.aC * p.p.D_s * p.p.V_0/ (p.p.a_s * p.p.R_s^2);
p.n.theta_c = (p.n.R_s^2 * p.n.sigma * p.n.kappa ) / (p.n.L^2*p.n.aC * (p.n.sigma + p.n.kappa) * p.n.D_s);
p.p.theta_c = (p.p.R_s^2 * p.p.sigma * p.p.kappa ) / (p.p.L^2*p.p.aC * (p.p.sigma + p.p.kappa) * p.p.D_s);
p.n.theta_d = (p.n.D_e * p.n.R_s^2) / (p.n.D_s * p.n.L^2 * p.n.epsilon_e);
p.p.theta_d = (p.p.D_e * p.p.R_s^2) / (p.p.D_s * p.p.L^2 * p.p.epsilon_e);
p.n.E = p.Faraday * p.n.V_0 / (p.R * p.T);
p.p.E = p.Faraday * p.p.V_0 / (p.R * p.T);
%p.E = p.E *10^(-floor(log10(p.E)));
p.n.K = 2 * (1 - p.t_0) * (1 + p.n.dlnfdlnce) / p.n.E;
p.p.K = 2 * (1 - p.t_0) * (1 + p.p.dlnfdlnce) / p.p.E;
p.n.mu = p.n.R_s^2*p.n.i_0_ref/(p.n.D_s*p.Faraday);
p.p.mu = p.p.R_s^2*p.p.i_0_ref/(p.p.D_s*p.Faraday);
% p.i_0_ref = p.k * p.Faraday * (p.Faraday / ((1 - p.t_0) * p.V_0 * p.aC)) / p.csmax; % Expression to check

%% Electrodes comparison ratios
p.rho_aCV_0 = sqrt(p.p.aC*p.p.V_0 / (p.n.aC*p.n.V_0));
p.rho_L = sqrt(p.p.L / p.n.L);
p.rho_kappa = sqrt(p.p.kappa / p.n.kappa);
p.rho_V_0 = sqrt(p.p.V_0 / p.n.V_0);

%% Electrodes-separator comparison ratios
p.p.lambda = p.s.L / p.p.L;
p.n.lambda = p.s.L / p.n.L;
p.p.delta = p.s.D_e / p.p.D_e;
p.n.delta = p.s.D_2 / p.n.D_e;

%% Simulation parameters
p.n.dscrtzn.N_e = 10; % Number of nodes in the electrolyte
p.p.dscrtzn.N_e = 10; 
p.n.dscrtzn.N_s = 5;% Number of nodes in the solid phase
p.n.dscrtzn.N_s = 5;
p.s.dscrtzn.N_e = 5;
end

