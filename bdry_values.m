function [ elyte0p,elyte1p,elytem1s,elyte1s,elyte1n,elyte0n,dl0p,dl1p,dl1n,dl0n ] = bdry_values( t,state,masks,matrices )
%BDRY_VALYES Summary of this function goes here
%   Detailed explanation goes here

ELYTE_bdrys = matrices.bdry_operators.elyte.bdry_electrodes * [state(masks.p.elyte);state(masks.n.elyte)] ...
                            + matrices.elyte.bdry_operators.bdry_separator * state(masks.s.elyte);

elyte1p = ELYTE_bdrys(1);
elyte1n = ELYTE_bdrys(2);
elyte0p = ELYTE_bdrys(3);
elyte0n = ELYTE_bdrys(4);
elytem1s = ELYTE_bdrys(1);
elyte1s = ELYTE_bdrys(2);
                        
DL_bdrys = bdry_dl * [state(masks.p.dl);state(masks.n.dl)]...
                             + bdry_i_adim * adim_current(t)...
                             + bdry_ln_internal * [ln(state(masks.p.elyte));ln(state(masks.n.elyte))]...
                             + bdry_ln_1 * [ln(elyte1p);ln(elyte1n)]...
                             + bdry_ln_0 * [ln(elyte0p);ln(elyte0n)];
                         
dl1p = DL_bdrys(1);
dl1n = DL_bdrys(2);
dl0p = DL_bdrys(3);
dl0n = DL_bdrys(4);
end

