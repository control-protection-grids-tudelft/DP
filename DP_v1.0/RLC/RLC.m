close all;
clear all;

N = 5;
Ts = 20e-6;
f = 50;
Ws = 2*pi*f;
Vnom_prim = 240e3;

R_src= 0.1;
L_src=0.00222;

R_lo = 50;
L_lo = 0.002;
 

%%
% load rlc.mat;
load Srlc.mat;
run('DQsymRLC')
sim('DQsymRLC')
