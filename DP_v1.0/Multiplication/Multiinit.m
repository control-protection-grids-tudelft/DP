clc ; clear all; 


N= 5;
N2= 2*N;
Ts= 2e-5;

f= 50 ; 

load Xdcpnz_c.mat
load Ydcpnz_c.mat
load Ydcpnz_c2.mat

open Multiplication.slx
sim Multiplication.slx