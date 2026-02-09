close all;
clear all;

N = 1;
Ts = 20e-6;
Ts_Power = 20e-6;
f = 50;
Ws = 2*pi*f;
Vnom_prim = 24e3;

%%
load rlc.mat;
load Srlc.mat;
run('DQsymRLC')
sim('DQsymRLC')
