function [ x_dot ] = state_t_deriv( t,x,matrices,masks,params )
%STATE_T_DERIV Summary of this function goes here
%   Detailed explanation goes here

x_dot = zeros(size(x));

%% Getting the boundary values

[ elyte0p,elyte1p,~,~,elyte1n,elyte0n,dl0p,dl1p,dl1n,dl0n ] = bdry_values(t,x,masks,matrices);

%% Electrolyte concentration
% Cathode

x_dot(masks.p.elyte) = params.p.theta_d * D2_BC(x(masks.p.elyte),elyte1p,elyte0p,matrices.p.adim.D2_noBC)...
                                        + params.p.theta_c * D2_BC(x(masks.p.dl),dl1p,dl0p,matrices.p.adim.D2_noBC)...
                                        + params.p.K * params.p.theta_c * logterm(x(masks.p.elyte),elyte1p,elyte0p,matrices.p.adim.D2_noBC);
% Anode
x_dot(masks.n.elyte) = params.n.theta_d * D2_BC(x(masks.n.elyte),elyte1n,elyte0n,matrices.n.adim.D2_noBC)...
                                        + params.n.theta_c * D2_BC(x(masks.n.dl),dl1n,dl0n,matrices.n.adim.D2_noBC)...
                                        + params.n.K * params.n.theta_c * logterm(x(masks.n.elyte),elyte1n,elyte0n,matrices.n.adim.D2_noBC);

% Separator
x_dot(masks.s.elyte) = params.s.theta_d * D2_BC(x(masks.s.elyte),elyte1p,elyte0p,matrices.s.adim.D2_noBC)...
                                        + params.s.K * params.s.theta_c * logterm(x(masks.s.elyte),elyte1p,elyte0p,matrices.s.adim.D2_noBC);

%% Double layer
% Cathode
x_dot(masks.p.dl) = params.p.theta_c * D2_BC(x(masks.p.dl),dl1p,dl0p,matrices.p.adim.D2_noBC)...
                                  + params.p.K * params.p.theta_c * logterm(x(masks.p.elyte),elyte1p,elyte0p,matrices.p.adim.D2_noBC)...
                                  - sinhterm(x(masks.p.elyte),x(masks.p.dl),x(masks.p.ussurf),params.p.alpha,params.p.E,params.p.V_0,params.p.ocpname);
                              
%Anode
x_dot(masks.n.dl) = params.n.theta_c * D2_BC(x(masks.n.dl),dl1p,dl0p,matrices.n.adim.D2_noBC)...
                                  + params.n.K * params.n.theta_c * logterm(x(masks.n.elyte),elyte1p,elyte0p,matrices.n.adim.D2_noBC)...
                                  - sinhterm(x(masks.n.elyte),x(masks.n.dl),x(masks.n.ussurf),params.n.alpha,params.n.E,params.n.V_0,params.n.ocpname);
                              
%% Particles concentration
%Cathode
for i = 1:params.p.dscrtzn.N_s - 1
    x_dot(masks.p.us(i)) = D2_BC(x(masks.p.us(i)),x(masks.p.ussurf(i)),zeros(size(x(masks.p.ussurf(i)))),matrices.p.adim.sphase.D2_noBC);
end

%Anode
for i = 1:params.n.dscrtzn.N_s - 1
    x_dot(masks.n.us(i)) = D2_BC(x(masks.n.us(i)),x(masks.n.ussurf(i)),zeros(size(x(masks.n.ussurf(i)))),matrices.n.adim.sphase.D2_noBC);
end

%% Algebraic equation
% Cathode
x_dot(mask.p.ussurf) = algebraicterm(x(masks.p.elyte),x(masks.p.dl),x(masks.p.us),x(masks.p.ussurf),...
                                                                        matrices.p.adim.sphase.D1_noBC(1,:),params.p.mu,params.alpha,params.p.E,params.p.V_0,params.p.ocpname);

% Anode
x_dot(mask.s.ussurf) = algebraicterm(x(masks.s.elyte),x(masks.s.dl),x(masks.s.us),x(masks.s.ussurf),...
                                                                        matrices.s.adim.sphase.D1_noBC(1,:),params.s.mu,params.alpha,params.s.E,params.s.V_0,params.s.ocpname);

end

function [logterm] = logterm(elyte,elyte1,elyte0,D2_noBC)
    logterm = D2_BC(ln(elyte),ln(elyte1),ln(elyte0),D2_noBC);
end

function [d2_bc] = D2_BC(x,x1,x0,D2_noBC)
    d2_bc = D2_noBC(2:end-1,2:end-1) * x...
                                + D2_noBC(2:end-1,1) * x1...
                                +  D2_noBC(2:end-1,end) * x0;
end

function [sinhterm] = sinhterm(elyte,dl,ussurf,alpha,E,V_0,ocpname)
    sinhterm = 2 .* (elyte .* (1 - ussurf) .* ussurf) .^ alpha ...
                            .* sinh(alpha * E * (dl - feval(ocpname,ussurf) / V_0));
end

function [algebraicterm] = algebraicterm(elyte,dl,us,ussurf,D1,mu,alpha,E,V_0,ocpname)
    algebraicterm = ussurf + mu / (D1(1) - 1) * sinhterm(elyte,dl,ussurf,alpha,E,V_0,ocpname)...
                                    + D1(2:end-1) / (D1(1) - 1) * us;
end