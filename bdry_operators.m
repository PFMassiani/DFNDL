function [ operators ] = bdry_operators( D1p,D0p,D1s,Dm1s,D1n,D0n,params )
% Summary :
%       Given the following boundary conditions:
%                D1p DLp = ... * i_adim + params.p.K * D1p ln(ELYTEp)
%                D1n DLn = ... * i_adim + params.p.K * D1p ln(ELYTEn)
%                D0p DLp = ... * i_adim
%                D0n DLn = ... * i_adim
%                D1p ELYTEp = params.p.delta * Dm1s ELYTEs
%                D1n ELYTEn = params.n.delta * D1s ELYTEs
%                D0p ELYTEp = 0
%                D0n ELYTEn = 0
%       Returns the operators st.
%               [ELYTEp(1);ELYTEn(1);ELYTEp(end);ELYTEn(end)] = bdry_electrodes * [ELYTEp(2:end-1);ELYTEn(2:end-1)]
%                                                                                   + bdry_separator * ELYTEs(2:end-1)
%               [DLp(1);DLn(1);DLp(end);DLn(end)] = bdry_dl * [DLp(2:end-1);DLn(2:end-1)]
%                                                                                    + bdry_i_adim * i_adim
%                                                                                    + bdry_ln_internal * [ln(ELYTEp(2:end-1));ln(ELYTEn(2:end-1))]
%                                                                                    + bdry_ln_1 * [ln(ELYTEp(1));ln(ELYTEn(1))];
%                                                                                    + bdry_ln_0 * [ln(ELYTEp(end));ln(ELYTEn(end))];
%       This script assumes the continuity of the electrolyte concentration
%       at the boundaries.
%       The implicit assumptions are : 
%               - the boundary cathode/separator is at index 1 in the cathode(x=1),
%                           index end in the separator(x=-1)
%               - the boundary anode/separator is at index 1 in the anode(x=1),
%                           index end in the separator(x=1).

M = [D1p(1) - params.p.delta/params.p.lambda / params.rho_aCV_0 * Dm1s(end) ,   -params.p.delta/params.p.lambda/params.rho_aCV_0 * Dm1s(1)                 ,   D1p(end)      ,    0;
        -params.n.delta/params.n.lambda * params.rho_aCV_0 * D1s(end)                    ,   D1n(1) - params.n.delta / params.n.lambda * params.rho_aCV_0 * D1s(1)  ,   0                    ,   D1n(end);
        D0p(1)                                                                                                                                 ,   0                                                                                                                                     ,   D0p(end)     ,   0;
        0                                                                                                                                            ,   D0n(1)                                                                                                                          ,   0                    ,   D0n(end)];

A_electrodes = -[blkdiag(D1p(2:end-1),D1n(2:end-1));
         blkdiag(D0p(2:end-1),D0n(2:end-1))];

A_separator = kron([1;0],[params.p.delta/params.p.lambda / params.rho_aCV_0 * Dm1s(2:end-1);
                            params.n.delta/params.n.lambda * params.rho_aCV_0 * D1s(2:end-1)]);
                        
operators.elyte.bdry_electrodes = M\A_electrodes;
operators.elyte.bdry_separator = M\A_separator;

clear M A_electrodes A_separator;

M = [D1p(end)   ,   0                   ,   D1p(1)  ,   0;
        0                   ,   D1n(end)    ,   0            ,   D1n(1);
        D0p(end)     ,   0                  ,   D0P(1)  ,   0;
        0                    ,   D0n(end)   ,   0             ,   D0n(1)];

A_dl = -[blkdiag(D1p(2:end-1),D1n(2:end-1));
               blkdiag(D0p(2:end-1),D0n(2:end-1))];

A_i_adim = [params.rho_L / params.rho_kappa / params.rho_V_0;
                       params.rho_kappa * params.rho_V_0 / params.rho_L;
                       - params.rho_L * sqrt(params.p.kappa * params.p.kappa) / params.p.sigma / params.rho_V_0;
                       - params.rho_V_0 * sqrt(params.p.kappa * params.n.kappa) / params.n.sigma / params.rho_L ];

A_ln_internal = kron([1;0], ...
                                        blkdiag(params.p.K * D1p(2:end-1),...
                                                                                        params.n.K * D1n(2:end-1)));
A_ln_1 = [params.p.K * D1p(1)   ,   0;
                  0   ,   params.n.K * D1n(1);
                  0    ,    0;
                  0    ,    0];
A_ln_0 = [params.p.K * D1p(end) ,   0;
                  0 ,   params.n.K * D1n(end);
                  0 ,   0;
                  0 ,   0];
              
operators.dl.bdry_dl = M\A_dl;
operators.dl.bdry_i_adim = M\A_i_adim;
operators.dl.bdry_ln_internal = M\A_ln_internal;
operators.dl.bdry_ln_1 = M\A_ln_1;
operators.dl.bdry_ln_0 = M\A_ln_0;
end

