%% Open-loop MMC: abc, HSS, DP and DQsym models
% Builds the harmonic-domain models of the open-loop MMC benchmark (Section 3)
% and produces the model-based results of the paper:
%
%   Fig. 2   sparsity of the linearised abc model
%   Fig. 3   sparsity of the HSS model
%   Fig. 4   sparsity of the DP model
%   Fig. 5   sparsity of the DQsym solver (pnz coordinates)
%   Fig. 9   eigenvalues of DP, HSS and DQsym
%   Table 3  state-space dimensions and sparsity
%
% Figs. 6-8 and Table 4 need the DQsym Simulink model and are produced by
% mmc_dp_hss_vs_simulink_cos.m, which uses the workspace of this script.
%
% Required: Symbolic Math Toolbox, Optimization Toolbox (fsolve).
% Figures are exported as vector PDFs to ./figs at their printed size.

clear; clc; close all;

%% 0. Settings

K  = 4;                         % harmonic truncation order
ph = 'abc';

pub.export = true;
pub.dir    = 'figs';
pub.fs     = 8;                 % font size [pt]
pubStyle(pub.fs);

% colours of the state groups i_u, i_l, v_Cu, v_Cl, i_n (Okabe-Ito)
cGrp = [  0 114 178;
        213  94   0;
          0 158 115;
        230 159   0;
        204 121 167]/255;

%% 1. Benchmark parameters (Table 2)

fg     = 50;   wg_v = 2*pi*fg;   Tp = 1/fg;
Ra     = 1.07;              % arm resistance [Ohm]
La     = 50e-3;             % arm inductance [H]
Csm    = 10e-3;             % submodule capacitance [F]
Nsm    = 400;               % submodules per arm
Vdc    = 640e3;             % pole-to-pole dc voltage [V]
Rst    = 10;                % star-point resistance [Ohm]
Lst    = 5;                 % star-point inductance [H]
Vac_pk = sqrt(2/3)*333e3;   % ac phase voltage, peak [V]
m      = 0.85;              % modulation amplitude
delta  = deg2rad(-2);       % modulation angle

%% 2. DP model: symbolic parameters and phasors
% Arm parameters are kept per arm and phase. Each signal x is a (K+1)x2 array
% whose row k+1 is [Re<x>_k, Im<x>_k], with Im<x>_0 = 0.

wg    = sym('w_g',    'positive');
Lstar = sym('L_star', 'positive');
Rstar = sym('R_star', 'positive');

Rau = sym(zeros(1,3)); Lau = Rau; Ral = Rau; Lal = Rau;
Cu  = Rau; Cl = Rau; Nu = Rau; Nl = Rau;
for p = 1:3
    Rau(p) = sym(['R_au_' ph(p)], 'positive');
    Lau(p) = sym(['L_au_' ph(p)], 'positive');
    Ral(p) = sym(['R_al_' ph(p)], 'positive');
    Lal(p) = sym(['L_al_' ph(p)], 'positive');
    Cu(p)  = sym(['C_SMu_' ph(p)], 'positive');
    Cl(p)  = sym(['C_SMl_' ph(p)], 'positive');
    Nu(p)  = sym(['N_u_' ph(p)], 'positive');
    Nl(p)  = sym(['N_l_' ph(p)], 'positive');
end

IU = cell(1,3); IL = IU; VCU = IU; VCL = IU; IN = IU;
NU = IU; NL = IU; VAC = IU;
for p = 1:3
    IU{p}  = dpvar(['iu_'  ph(p)], K);
    IL{p}  = dpvar(['il_'  ph(p)], K);
    VCU{p} = dpvar(['vCu_' ph(p)], K);
    VCL{p} = dpvar(['vCl_' ph(p)], K);
    IN{p}  = dpvar(['in_'  ph(p)], K);
    NU{p}  = dpvar(['nu_'  ph(p)], K);
    NL{p}  = dpvar(['nl_'  ph(p)], K);
    VAC{p} = dpvar(['vac_' ph(p)], K);
end
VUDC = dpvar('vudc', K);
VLDC = dpvar('vldc', K);
VN   = sym(zeros(K+1,2));                % v_n = 0

%% 3. DP equations (13), h = 0..K
% Convolutions are truncated to |l| <= K and |h-l| <= K.

dIU = cell(1,3); dIL = dIU; dVCU = dIU; dVCL = dIU; dIN = dIU;
for p = 1:3
    dIU{p} = sym(zeros(K+1,2)); dIL{p} = dIU{p};
    dVCU{p} = dIU{p}; dVCL{p} = dIU{p}; dIN{p} = dIU{p};
    for k = 0:K
        r = k + 1;
        dIU{p}(r,:) = ( dp(VUDC,k) - dp(VN,k) - dp(VAC{p},k) ...
                        - Rau(p)*dp(IU{p},k) - dpconv(NU{p},VCU{p},k) ) / Lau(p) ...
                      + mjkw(IU{p},k,wg);
        dIL{p}(r,:) = ( dp(VN,k) - dp(VLDC,k) + dp(VAC{p},k) ...
                        - Ral(p)*dp(IL{p},k) - dpconv(NL{p},VCL{p},k) ) / Lal(p) ...
                      + mjkw(IL{p},k,wg);
        dVCU{p}(r,:) = Nu(p)/Cu(p) * dpconv(NU{p},IU{p},k) + mjkw(VCU{p},k,wg);
        dVCL{p}(r,:) = Nl(p)/Cl(p) * dpconv(NL{p},IL{p},k) + mjkw(VCL{p},k,wg);
        sumIN = dp(IN{1},k) + dp(IN{2},k) + dp(IN{3},k);
        dIN{p}(r,:) = ( dp(VAC{p},k) - Rstar*sumIN ) / Lstar + mjkw(IN{p},k,wg);
    end
end

imag0 = [];
for p = 1:3
    imag0 = [imag0; dIU{p}(1,2); dIL{p}(1,2); dVCU{p}(1,2); dVCL{p}(1,2); dIN{p}(1,2)]; %#ok<AGROW>
end
assert(all(isAlways(expand(imag0) == 0)), 'The h = 0 phasors are not real.');

%% 4. Real-valued state and input vectors (15)

x = sym([]); f = sym([]);
for p = 1:3, [x,f] = addBlock(x,f,{IU{p},IL{p}},  {dIU{p},dIL{p}},   0);   end
for p = 1:3, [x,f] = addBlock(x,f,{IU{p},IL{p}},  {dIU{p},dIL{p}},   1:K); end
for p = 1:3, [x,f] = addBlock(x,f,{VCU{p},VCL{p}},{dVCU{p},dVCL{p}}, 0);   end
for p = 1:3, [x,f] = addBlock(x,f,{VCU{p},VCL{p}},{dVCU{p},dVCL{p}}, 1:K); end
for p = 1:3, [x,f] = addBlock(x,f,{IN{p}},        {dIN{p}},          0);   end
for p = 1:3, [x,f] = addBlock(x,f,{IN{p}},        {dIN{p}},          1:K); end

u = sym([]);
for p = 1:3, u = addBlock(u,sym([]),{VAC{p},NU{p},NL{p}},{VAC{p},NU{p},NL{p}},0);   end
for p = 1:3, u = addBlock(u,sym([]),{VAC{p},NU{p},NL{p}},{VAC{p},NU{p},NL{p}},1:K); end
u = addBlock(u,sym([]),{VUDC,VLDC},{VUDC,VLDC},0);
u = addBlock(u,sym([]),{VUDC},{VUDC},1:K);
u = addBlock(u,sym([]),{VLDC},{VLDC},1:K);

nx = numel(x);
stateNames = string(x);
inputNames = string(u);

%% 5. Linearisation of the DP model

A_sym = jacobian(f, x);
B_sym = jacobian(f, u);

parSyms = [wg, Lstar, Rstar, Rau, Lau, Ral, Lal, Cu, Cl, Nu, Nl];
parVals = [wg_v, Lst, Rst, Ra*ones(1,3), La*ones(1,3), Ra*ones(1,3), La*ones(1,3), ...
           Csm*ones(1,3), Csm*ones(1,3), Nsm*ones(1,3), Nsm*ones(1,3)];

f_fun = matlabFunction(subs(f,     parSyms, parVals), 'Vars', {x, u});
A_fun = matlabFunction(subs(A_sym, parSyms, parVals), 'Vars', {x, u});
B_fun = matlabFunction(subs(B_sym, parSyms, parVals), 'Vars', {x, u});

%% 6. Operating point: ac source and insertion indices (6)

nVAC = cell(1,3); nNU = nVAC; nNL = nVAC;
for p = 1:3
    th = -(p-1)*2*pi/3;
    nVAC{p} = zeros(K+1,2); nVAC{p}(2,:) = Vac_pk/2*[cos(th), sin(th)];
    nNU{p}  = zeros(K+1,2); nNU{p}(1,1) = 0.5; nNU{p}(2,:) = -m/4*[cos(th+delta), sin(th+delta)];
    nNL{p}  = zeros(K+1,2); nNL{p}(1,1) = 0.5; nNL{p}(2,:) = +m/4*[cos(th+delta), sin(th+delta)];
end
nVUDC = zeros(K+1,2); nVUDC(1,1) = +Vdc/2;
nVLDC = zeros(K+1,2); nVLDC(1,1) = -Vdc/2;

uSym = []; uVal = [];
for p = 1:3
    uSym = [uSym; dpflat(VAC{p}); dpflat(NU{p}); dpflat(NL{p})];      %#ok<AGROW>
    uVal = [uVal; dpflat(nVAC{p}); dpflat(nNU{p}); dpflat(nNL{p})];   %#ok<AGROW>
end
uSym = [uSym; dpflat(VUDC); dpflat(VLDC)];
uVal = [uVal; dpflat(nVUDC); dpflat(nVLDC)];
u0 = double(subs(u, uSym, uVal));

%% 7. DP steady state and linearised DP model

xinit = zeros(nx,1);
xinit(startsWith(stateNames,"vC") & endsWith(stateNames,"_R0")) = Vdc;

opts = optimoptions('fsolve', 'Algorithm','levenberg-marquardt', ...
    'SpecifyObjectiveGradient',true, 'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12, 'OptimalityTolerance',1e-12, ...
    'MaxIterations',500, 'Display','off');
obj = @(xx) deal(f_fun(xx,u0), A_fun(xx,u0));
x0 = fsolve(obj, xinit, opts);

A0 = A_fun(x0, u0);
B0 = B_fun(x0, u0);
lamDP = eig(A0);

%% 8. abc model (4): periodic steady state

prm.Rau = Ra*ones(3,1); prm.Lau = La*ones(3,1); prm.Ral = Ra*ones(3,1); prm.Lal = La*ones(3,1);
prm.Cu = Csm*ones(3,1); prm.Cl = Csm*ones(3,1); prm.Nu = Nsm*ones(3,1); prm.Nl = Nsm*ones(3,1);
prm.Lstar = Lst; prm.Rstar = Rst;
Uin.vac = nVAC; Uin.nu = nNU; Uin.nl = nNL; Uin.vudc = nVUDC; Uin.vldc = nVLDC;

xTDnames = strings(15,1); c = 0;
for s = ["iu","il","vCu","vCl","in"]
    for p = 1:3, c = c + 1; xTDnames(c) = s + "_" + ph(p); end
end
uTDnames = strings(11,1); c = 0;
for s = ["vac","nu","nl"]
    for p = 1:3, c = c + 1; uTDnames(c) = s + "_" + ph(p); end
end
uTDnames(10) = "vudc"; uTDnames(11) = "vldc";

z0 = zeros(15,1);
for i = 1:15
    z0(i) = dp2time(dpget(xTDnames(i), x0, stateNames, K), 0, wg_v);
end

Tsim    = 1.0;
odeOpts = odeset('RelTol',1e-9, 'AbsTol',1e-6, 'MaxStep',1e-4);
solTD   = ode15s(@(t,z) mmcTD(t,z,prm,Uin,wg_v), [0 Tsim], z0, odeOpts);

M  = 512;
tg = Tsim - Tp + (0:M-1)*Tp/M;           % last period of the run
xg_TD = deval(solTD, tg);
ug = zeros(11, M);
for p = 1:3
    ug(p,:)   = dp2time(nVAC{p}, tg, wg_v);
    ug(3+p,:) = dp2time(nNU{p},  tg, wg_v);
    ug(6+p,:) = dp2time(nNL{p},  tg, wg_v);
end
ug(10,:) = dp2time(nVUDC, tg, wg_v);
ug(11,:) = dp2time(nVLDC, tg, wg_v);

%% 9. HSS model (9)-(11)
% A(t) and B(t) of (9) are evaluated on the abc periodic trajectory, expanded in
% Fourier series up to |h| <= 2K and arranged in block-Toeplitz form.

iuS  = sym('iu',[3 1],'real');  ilS  = sym('il',[3 1],'real');
vCuS = sym('vCu',[3 1],'real'); vClS = sym('vCl',[3 1],'real');
inS  = sym('istar',[3 1],'real');
vacS = sym('vac',[3 1],'real'); nuS = sym('nu',[3 1],'real'); nlS = sym('nl',[3 1],'real');
vudcS = sym('vudc','real');     vldcS = sym('vldc','real');

xTD = [iuS; ilS; vCuS; vClS; inS];
uTD = [vacS; nuS; nlS; vudcS; vldcS];
fTD = [ (vudcS - vacS - prm.Rau.*iuS - nuS.*vCuS) ./ prm.Lau;
        (-vldcS + vacS - prm.Ral.*ilS - nlS.*vClS) ./ prm.Lal;
        prm.Nu ./ prm.Cu .* nuS .* iuS;
        prm.Nl ./ prm.Cl .* nlS .* ilS;
        (vacS - prm.Rstar*sum(inS)) / prm.Lstar ];

A_TD = jacobian(fTD, xTD);
B_TD = jacobian(fTD, uTD);
Atd_fun = matlabFunction(A_TD, 'Vars', {xTD, uTD});
Btd_fun = matlabFunction(B_TD, 'Vars', {xTD, uTD});

[Ah, Bh] = hssCoeffs(Atd_fun, Btd_fun, xg_TD, ug, tg, wg_v, 2*K);
[A_hss, B_hss, hssXnames, hssUnames] = hssBuild(Ah, Bh, K, wg_v, xTDnames, uTDnames);
lamHSS = eig(A_hss);

%% 10. DP inputs in HSS coordinates (used by the Simulink comparison)

[Tu, hssUnamesDP] = dp2hssT(inputNames, K);
[okU, pU] = ismember(hssUnames, hssUnamesDP);
assert(all(okU), 'Input names differ between DP and HSS.');

%% 11. DQsym model (16)-(22)
% Solver (16)-(20): arm and neutral currents in sequence coordinates, per
% harmonic h: d/dt X_h = (A_dq - j h w0 I) X_h + B_dq U_h, with
%   X = [i_u,pnz; i_l,pnz; i_n,pnz],
%   U = [v^ac_pnz; v^dc_u,pnz; v^dc_l,pnz; v_u,pnz; v_l,pnz].
% Closed loop (22): the capacitor voltages and the products m.*v_C, m.*i are
% added around the solver in harmonic form, h = -K..K.

I3 = eye(3);  Z3 = zeros(3);
A_dq = blkdiag(-Ra/La*eye(6), -Rst/Lst*diag([0 0 3]));
B_dq = [ [-I3,  I3,  Z3, -I3,  Z3]/La;
         [ I3,  Z3, -I3,  Z3, -I3]/La;
         [ I3/Lst, zeros(3,12)] ];

nG = 2*K + 1;  ka = -K:K;  Ih = eye(nG);
Nh = diag(1j*wg_v*ka);
A_RL = kron(A_dq, Ih) - kron(eye(9), Nh);           % solver, harmonics -K..K

a1  = exp(1j*2*pi/3);
Tsq = [1 a1 a1^2; 1 a1^2 a1; 1 1 1]/3;               % abc -> pnz (17)
T6  = kron(blkdiag(Tsq, Tsq), Ih);                   % six arm signals
Mn  = blkdiag(toeplitzDP(nNU{1},K), toeplitzDP(nNU{2},K), toeplitzDP(nNU{3},K), ...
              toeplitzDP(nNL{1},K), toeplitzDP(nNL{2},K), toeplitzDP(nNL{3},K));
B_v  = kron(B_dq(:, 10:15), Ih);                     % arm voltages into the solver
C_RL = T6 \ [eye(6*nG), zeros(6*nG, 3*nG)];          % arm currents (abc) from the solver
A_C  = -kron(eye(6), Nh);                            % capacitor integrators
B_C  = Nsm/Csm*Mn;                                   % m .* i
C_C  = T6*Mn;                                        % v_m = m .* v_C, in pnz
A_cl = [A_RL, B_v*C_C; B_C*C_RL, A_C];

lamRL = eig(A_RL);
lamCL = eig(A_cl);

%% 12. Table 3: state-space dimensions and sparsity

nzTol = @(Mx) abs(Mx) > 1e-9*max(abs(Mx(:)));
Sdq_A = nzTol(A_dq - 1j*wg_v*eye(9));                % pattern for h ~= 0
Sdq_B = nzTol(B_dq);
S.A = {Sdq_A, kron(eye(K+1), Sdq_A), nzTol(A0), nzTol(A_hss)};
S.B = {Sdq_B, kron(eye(K+1), Sdq_B), nzTol(B0), nzTol(B_hss)};

nA = cellfun(@nnz, S.A).';   nB = cellfun(@nnz, S.B).';
nX = cellfun(@(Q) size(Q,1), S.A).';   nU = cellfun(@(Q) size(Q,2), S.B).';
table3 = table(nX, nU, nA, nB, 100*nA./nX.^2, 100*nB./(nX.*nU), nA + nB, ...
    'VariableNames', {'states','inputs','nnz_A','nnz_B','density_A','density_B','total_nnz'}, ...
    'RowNames', {'DQsym','DQsym (all K+1)','DP','HSS'})

%% 13. Figs. 2-5: sparsity patterns

nzA_TD = ~arrayfun(@(e) isequal(e, sym(0)), A_TD);
nzB_TD = ~arrayfun(@(e) isequal(e, sym(0)), B_TD);

pDx = signalOrder(stateNames, xTDnames, K, "dp");
pDu = signalOrder(inputNames, uTDnames, K, "dp");
pHx = signalOrder(hssXnames,  xTDnames, K, "hss");
pHu = signalOrder(hssUnames,  uTDnames, K, "hss");

labX   = {'$i_{u,abc}$', '$i_{l,abc}$', '$v_{Cu,abc}$', '$v_{Cl,abc}$', '$i_{n,abc}$'};
labU   = {'$v^{ac}_{abc}$', '$m_{u,abc}$', '$m_{l,abc}$', '$v^{dc}_{u,l}$'};
legAll = {'$i_u$', '$i_l$', '$v_{Cu}$', '$v_{Cl}$', '$i_n$'};

sp = struct('rows', 3*[1 1 1 1 1], 'cols', [3 3 3 2], 'w', 1, 'colour', 1:5, ...
            'labX', {labX}, 'labU', {labU}, 'leg', {legAll}, 'ms', 7, ...
            'ttl', {{'abc, $\mathbf{A}_{\mathrm{abc}}$', 'abc, $\mathbf{B}_{\mathrm{abc}}$'}});
fig = pubFig('Fig. 2 sparsity abc', 8.8, 6.4);
sparsityFigure(fig, nzA_TD, nzB_TD, sp, cGrp);
saveFig(fig, 'fig2_sparsity_abc', pub);

sp.rows = 3*nG*[1 1 1 1 1];  sp.cols = nG*[3 3 3 2];  sp.w = nG;  sp.ms = 3;
sp.ttl  = {'HSS, $\mathbf{A}_{\mathrm{HSS}}-\mathcal{N}$', 'HSS, $\mathbf{B}_{\mathrm{HSS}}$'};
fig = pubFig('Fig. 3 sparsity HSS', 8.8, 6.4);
sparsityFigure(fig, S.A{4}(pHx,pHx), S.B{4}(pHx,pHu), sp, cGrp);
saveFig(fig, 'fig3_sparsity_hss', pub);

sp.ttl  = {'DP, $\mathbf{A}_{\mathrm{DP}}$', 'DP, $\mathbf{B}_{\mathrm{DP}}$'};
fig = pubFig('Fig. 4 sparsity DP', 8.8, 6.4);
sparsityFigure(fig, S.A{3}(pDx,pDx), S.B{3}(pDx,pDu), sp, cGrp);
saveFig(fig, 'fig4_sparsity_dp', pub);

sp5 = struct('rows', [3 3 3], 'cols', [3 3 3 3 3], 'w', 1, 'colour', [1 2 5], ...
             'labX', {{'$i_{u,pnz}$', '$i_{l,pnz}$', '$i_{n,pnz}$'}}, ...
             'labU', {{'$v^{ac}_{pnz}$', '$v^{dc}_{u,pnz}$', '$v^{dc}_{l,pnz}$', ...
                       '$v_{u,pnz}$', '$v_{l,pnz}$'}}, ...
             'leg', {{'$i_u$', '$i_l$', '$i_n$'}}, 'ms', 7, ...
             'ttl', {{'DQsym, $\mathbf{A}_{\mathrm{DQsym}}$', 'DQsym, $\mathbf{B}_{\mathrm{DQsym}}$'}});
fig = pubFig('Fig. 5 sparsity DQsym', 8.8, 4.6);
sparsityFigure(fig, Sdq_A, Sdq_B, sp5, cGrp);
saveFig(fig, 'fig5_sparsity_dqsym', pub);

%% 14. Fig. 9: eigenvalues

sArm  = -Ra/La;                 % solver arm pole
sArm2 = sArm/2;                 % arm pole with the capacitors connected
sStar = -3*Rst/Lst;             % star branch, zero sequence

fs   = 8;   fsS = 7;
cDP  = [0.000 0.447 0.698];   cHS = [0.835 0.369 0.000];
cDQ  = [0.000 0.620 0.451];   cRL = [0.800 0.475 0.655];
cStr = [0.93 0.93 0.93];      cFam = [0.45 0.45 0.45];
xL   = [sArm - 2.5, 2];
fMax = ceil(1.1*max(abs(imag([lamCL; lamRL])))/(2*pi)/fg)*fg;
yL   = [-fMax, fMax];
fz   = @(l) imag(l)/(2*pi);

fig = pubFig('Fig. 9 eigenvalues', 8.8, 10.5);
tl  = tiledlayout(fig, 2, 1, 'TileSpacing','tight', 'Padding','compact');

% (a) DP, HSS and DQsym closed loop
ax1 = eigAxes(tl, xL, yL, fg, cStr, fs);
h = gobjects(1,4);
h(1) = plot(ax1, real(lamDP),  fz(lamDP),  'o', 'MarkerSize',5,   'Color',cDP, 'LineWidth',0.7);
h(2) = plot(ax1, real(lamHSS), fz(lamHSS), 'x', 'MarkerSize',5,   'Color',cHS, 'LineWidth',0.8);
h(3) = plot(ax1, real(lamCL),  fz(lamCL),  '+', 'MarkerSize',5.5, 'Color',cDQ, 'LineWidth',1.0);
h(4) = patch(ax1, nan, nan, cStr, 'EdgeColor','none');
set(ax1, 'XTickLabel', []);
text(ax1, 0.015, 0.96, '(a)', 'Units','normalized', 'FontSize',fs, 'VerticalAlignment','top');
legend(ax1, h, {'DP', 'HSS', 'DQsym, $A_{\mathrm{cl}}$', '$|f|<f_g/2$'}, ...
       'FontSize',fsS, 'Orientation','horizontal', 'Location','northoutside', ...
       'Box','off', 'NumColumns',4);

% (b) DQsym solver alone and closed loop
ax2 = eigAxes(tl, xL, yL, fg, cStr, fs);
fam = {sArm,  '$-\frac{R_{arm}}{L_{arm}}$';
       sArm2, '$-\frac{R_{arm}}{2L_{arm}}$';
       sStar, '$-\frac{3R_{star}}{L_{star}}$';
       0,     '$0$'};
for j = 1:size(fam,1)
    xline(ax2, fam{j,1}, '--', 'Color',cFam, 'LineWidth',0.5, 'HandleVisibility','off');
    text(ax2, fam{j,1}, yL(2), fam{j,2}, 'FontSize',fsS, 'Color',cFam, ...
         'HorizontalAlignment','center', 'VerticalAlignment','bottom');
end
g = gobjects(1,2);
g(1) = plot(ax2, real(lamRL), fz(lamRL), 's', 'MarkerSize',4.5, 'Color',cRL, 'LineWidth',0.7);
g(2) = plot(ax2, real(lamCL), fz(lamCL), '+', 'MarkerSize',5.5, 'Color',cDQ, 'LineWidth',1.0);

ya = 0.80*yL(2);  xa = [sArm + 0.4, sArm2 - 0.4];
plot(ax2, xa, [ya ya], 'k-', 'LineWidth',0.7, 'HandleVisibility','off');
plot(ax2, xa(2), ya, 'k>', 'MarkerFaceColor','k', 'MarkerSize',3.5, 'HandleVisibility','off');
text(ax2, mean(xa), ya, 'closing the $C$ loop', 'FontSize',fsS, ...
     'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
     'BackgroundColor','w', 'Margin',0.5);
text(ax2, 0.015, 0.96, '(b)', 'Units','normalized', 'FontSize',fs, 'VerticalAlignment','top');
xlabel(ax2, '$\sigma$ (s$^{-1}$)', 'FontSize',fs);
legend(ax2, g, {'solver alone, $A_{RL}$', 'closed loop, $A_{\mathrm{cl}}$'}, ...
       'FontSize',fsS, 'Location','southwest', 'Box','on');
uistack(findobj(ax2, 'Type','line', 'Marker','+'), 'top');
saveFig(fig, 'fig9_eigenvalues', pub);

%% Local functions: figures

function pubStyle(fs)
set(groot, ...
    'defaultTextInterpreter',               'latex', ...
    'defaultAxesTickLabelInterpreter',      'latex', ...
    'defaultLegendInterpreter',             'latex', ...
    'defaultColorbarTickLabelInterpreter',  'latex', ...
    'defaultAxesFontSize',                  fs, ...
    'defaultTextFontSize',                  fs, ...
    'defaultLegendFontSize',                fs, ...
    'defaultAxesLabelFontSizeMultiplier',   1, ...
    'defaultAxesTitleFontSizeMultiplier',   1, ...
    'defaultAxesTitleFontWeight',           'normal', ...
    'defaultAxesLineWidth',                 0.5, ...
    'defaultLineLineWidth',                 1.0, ...
    'defaultAxesBox',                       'on', ...
    'defaultAxesXGrid',                     'on', ...
    'defaultAxesYGrid',                     'on', ...
    'defaultAxesGridAlpha',                 0.15, ...
    'defaultAxesMinorGridAlpha',            0.08, ...
    'defaultAxesTickLength',                [0.015 0.015], ...
    'defaultLegendBox',                     'off', ...
    'defaultFigureColor',                   'w');
end

function fig = pubFig(name, w, h)
% figure at its printed size [cm]
fig = figure('Name', name, 'Color', 'w', 'Units', 'centimeters');
fig.Position = [2 2 w h];
end

function saveFig(fig, fname, pub)
if ~pub.export, return; end
if ~exist(pub.dir, 'dir'), mkdir(pub.dir); end
drawnow;
exportgraphics(fig, fullfile(pub.dir, [fname '.pdf']), ...
    'ContentType', 'vector', 'BackgroundColor', 'none');
end

function sparsityFigure(fig, SA, SB, sp, cGrp)
% A and B side by side, nonzeros coloured by the state group of their row
tl = tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');
nexttile(tl);
sparsityPanel(SA, sp.rows, sp.rows, sp.w, sp.colour, cGrp, sp.labX, sp.labX, ...
              sp.ttl{1}, 'states', true, sp.ms);
ax = nexttile(tl);
sparsityPanel(SB, sp.rows, sp.cols, sp.w, sp.colour, cGrp, sp.labX, sp.labU, ...
              sp.ttl{2}, 'inputs', false, sp.ms);
hL = gobjects(1, numel(sp.leg));
for g = 1:numel(sp.leg)
    hL(g) = plot(ax, nan, nan, '.', 'Color', cGrp(sp.colour(g),:), 'MarkerSize', 10);
end
lg = legend(hL, sp.leg, 'Orientation','horizontal', 'Box','off');
lg.Layout.Tile = 'south';
lg.ItemTokenSize = [8 8];
end

function sparsityPanel(S, rG, cG, w, cIdx, cGrp, rLab, cLab, ttl, xl, showY, ms)
% S: logical pattern; rG, cG: rows / columns per group; w: rows per signal
ax = gca;  hold(ax, 'on');
[nr, nc] = size(S);
rE = [0 cumsum(rG)] + 0.5;
cE = [0 cumsum(cG)] + 0.5;
tint = @(c) 1 - 0.15*(1 - c);
for g = 1:numel(rG)
    patch(ax, [0.5 nc+0.5 nc+0.5 0.5], rE([g g g+1 g+1]), tint(cGrp(cIdx(g),:)), ...
          'EdgeColor','none', 'HandleVisibility','off');
end
gMi = [0.85 0.85 0.85];  gMa = [0.45 0.45 0.45];
for e = setdiff((w:w:nc-1) + 0.5, cE)
    plot(ax, [e e], [0.5 nr+0.5], '-', 'Color', gMi, 'LineWidth', 0.3, 'HandleVisibility','off');
end
for e = setdiff((w:w:nr-1) + 0.5, rE)
    plot(ax, [0.5 nc+0.5], [e e], '-', 'Color', gMi, 'LineWidth', 0.3, 'HandleVisibility','off');
end
for e = cE(2:end-1)
    plot(ax, [e e], [0.5 nr+0.5], '-', 'Color', gMa, 'LineWidth', 0.5, 'HandleVisibility','off');
end
for e = rE(2:end-1)
    plot(ax, [0.5 nc+0.5], [e e], '-', 'Color', gMa, 'LineWidth', 0.5, 'HandleVisibility','off');
end
[ii, jj] = find(S);
grp = discretize(ii, rE);
for g = 1:numel(rG)
    k = grp == g;
    plot(ax, jj(k), ii(k), '.', 'Color', cGrp(cIdx(g),:), 'MarkerSize', ms, ...
         'HandleVisibility','off');
end
set(ax, 'YDir','reverse', 'XGrid','off', 'YGrid','off', 'Box','on', 'Layer','top', ...
        'TickLength', [0 0], 'FontSize', 6, 'XTickLabelRotation', 0, ...
        'XTick', (cE(1:end-1) + cE(2:end))/2, 'XTickLabel', cLab, ...
        'YTick', (rE(1:end-1) + rE(2:end))/2, 'YTickLabel', rLab);
axis(ax, 'image');  xlim(ax, [0.5 nc+0.5]);  ylim(ax, [0.5 nr+0.5]);
if showY
    ylabel(ax, 'states');
else
    set(ax, 'YTickLabel', []);
end
xlabel(ax, xl);
title(ax, ttl);
end

function ax = eigAxes(tl, xL, yL, fg, cStr, fs)
% eigenvalue axes with the strip |f| < fg/2 shaded, y ticks at multiples of fg
ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');
patch(ax, [xL(1) xL(2) xL(2) xL(1)], [-fg -fg fg fg]/2, cStr, ...
      'EdgeColor','none', 'HandleVisibility','off');
xlim(ax, xL);  ylim(ax, yL);
step = fg*max(1, round(diff(yL)/fg/8));
set(ax, 'YTick', -floor(yL(2)/step)*step:step:yL(2), 'FontSize', fs, ...
        'GridAlpha', 0.12, 'TickLength', [0.008 0.008], 'Layer','top');
ylabel(ax, '$f$ (Hz)');
end

%% Local functions: models

function X = dpvar(name, K)
X = sym(zeros(K+1,2));
X(1,1) = sym([name '_R0'], 'real');
for k = 1:K
    X(k+1,1) = sym(sprintf('%s_R%d', name, k), 'real');
    X(k+1,2) = sym(sprintf('%s_I%d', name, k), 'real');
end
end

function z = dp(X, k)
% <x>_k for k = -K..K, as [Re Im]
K = size(X,1) - 1;
if abs(k) > K
    z = sym([0 0]);
elseif k >= 0
    z = X(k+1,:);
else
    z = [X(-k+1,1), -X(-k+1,2)];
end
end

function z = cmul(a, b)
z = [a(1)*b(1) - a(2)*b(2), a(1)*b(2) + a(2)*b(1)];
end

function z = dpconv(A, B, k)
% <a b>_k = sum_l <a>_{k-l} <b>_l
K = size(A,1) - 1;
z = sym([0 0]);
for l = -K:K
    if abs(k-l) <= K
        z = z + cmul(dp(A,k-l), dp(B,l));
    end
end
end

function z = mjkw(X, k, w)
% -j k w <x>_k
z = [k*w*X(k+1,2), -k*w*X(k+1,1)];
end

function [x, f] = addBlock(x, f, Xs, Fs, ks)
for k = ks
    for s = 1:numel(Xs)
        if k == 0
            x = [x; Xs{s}(1,1)];        f = [f; Fs{s}(1,1)];        %#ok<AGROW>
        else
            x = [x; Xs{s}(k+1,:).'];    f = [f; Fs{s}(k+1,:).'];    %#ok<AGROW>
        end
    end
end
end

function v = dpflat(X)
v = X(1,1);
for k = 1:size(X,1)-1
    v = [v; X(k+1,1); X(k+1,2)]; %#ok<AGROW>
end
end

function X = dpget(base, vec, names, K)
base = char(base);
X = zeros(K+1,2);
X(1,1) = vec(names == [base '_R0']);
for k = 1:K
    X(k+1,1) = vec(names == sprintf('%s_R%d', base, k));
    X(k+1,2) = vec(names == sprintf('%s_I%d', base, k));
end
end

function y = dp2time(X, t, w)
y = X(1,1)*ones(size(t));
for k = 1:size(X,1)-1
    y = y + 2*(X(k+1,1)*cos(k*w*t) - X(k+1,2)*sin(k*w*t));
end
end

function dz = mmcTD(t, z, prm, Uin, w)
% abc arm-averaged model (4)
iu = z(1:3); il = z(4:6); vCu = z(7:9); vCl = z(10:12); in = z(13:15);
vac = zeros(3,1); nu = vac; nl = vac;
for p = 1:3
    vac(p) = dp2time(Uin.vac{p}, t, w);
    nu(p)  = dp2time(Uin.nu{p},  t, w);
    nl(p)  = dp2time(Uin.nl{p},  t, w);
end
vud = dp2time(Uin.vudc, t, w);
vld = dp2time(Uin.vldc, t, w);
dz = [ (vud - vac - prm.Rau.*iu - nu.*vCu) ./ prm.Lau;
       (-vld + vac - prm.Ral.*il - nl.*vCl) ./ prm.Lal;
       prm.Nu ./ prm.Cu .* nu .* iu;
       prm.Nl ./ prm.Cl .* nl .* il;
       (vac - prm.Rstar*sum(in)) / prm.Lstar ];
end

function [Ah, Bh] = hssCoeffs(Afun, Bfun, xg, ug, tg, w, H)
% Fourier coefficients h = -H..H of A(t), B(t) sampled over one period
M  = numel(tg);
A1 = Afun(xg(:,1), ug(:,1)); B1 = Bfun(xg(:,1), ug(:,1));
Ah = zeros([size(A1), 2*H+1]); Bh = zeros([size(B1), 2*H+1]);
for q = 1:M
    Aq = Afun(xg(:,q), ug(:,q));
    Bq = Bfun(xg(:,q), ug(:,q));
    for h = -H:H
        e = exp(-1j*h*w*tg(q)) / M;
        Ah(:,:,h+H+1) = Ah(:,:,h+H+1) + Aq*e;
        Bh(:,:,h+H+1) = Bh(:,:,h+H+1) + Bq*e;
    end
end
end

function [Ahss, Bhss, xNames, uNames] = hssBuild(Ah, Bh, K, w, xn, un)
% block-Toeplitz A - N and B (11), harmonic-outer ordering
H  = (size(Ah,3) - 1)/2;
nx = size(Ah,1); nu = size(Bh,2);
ka = -K:K; nH = numel(ka);
Acal = zeros(nx*nH); Bcal = zeros(nx*nH, nu*nH);
for a = 1:nH
    ra = (a-1)*nx + (1:nx);
    for b = 1:nH
        h = ka(a) - ka(b);
        if abs(h) <= H
            Acal(ra, (b-1)*nx + (1:nx)) = Ah(:,:,h+H+1);
            Bcal(ra, (b-1)*nu + (1:nu)) = Bh(:,:,h+H+1);
        end
    end
end
Ahss = Acal - kron(diag(1j*w*ka), eye(nx));
Bhss = Bcal;
xNames = strings(nx*nH,1); uNames = strings(nu*nH,1);
for a = 1:nH
    xNames((a-1)*nx + (1:nx)) = xn(:) + "_" + ka(a);
    uNames((a-1)*nu + (1:nu)) = un(:) + "_" + ka(a);
end
end

function [T, names] = dp2hssT(realNames, K)
% real DP coordinates [Re, Im] -> complex harmonics -K..K
realNames = realNames(:);
bases = unique(regexprep(realNames, '_[RI]\d+$', ''), 'stable');
n = numel(realNames);
T = zeros(n); names = strings(n,1);
r = 0;
for b = 1:numel(bases)
    for k = -K:K
        r = r + 1;
        names(r) = bases(b) + "_" + k;
        T(r, realNames == bases(b) + "_R" + abs(k)) = 1;
        if k ~= 0
            T(r, realNames == bases(b) + "_I" + abs(k)) = sign(k)*1j;
        end
    end
end
assert(r == n, 'Inconsistent DP naming.');
end

function p = signalOrder(names, sigs, K, form)
% permutation that groups the harmonic components of each signal
p = [];
for s = sigs(:).'
    if form == "dp"
        c = s + "_R0";
        for k = 1:K, c = [c, s + "_R" + k, s + "_I" + k]; end %#ok<AGROW>
    else
        c = s + "_" + (-K:K);
    end
    [ok, i] = ismember(c, names);
    assert(all(ok), 'Signal %s not found.', s);
    p = [p, i]; %#ok<AGROW>
end
end

function T = toeplitzDP(X, K)
% Toeplitz matrix (harmonics -K..K) of a periodic signal given by its DP array
c = zeros(1, 2*K+1);
c(K+1) = X(1,1);
for k = 1:min(K, size(X,1)-1)
    c(K+1+k) = X(k+1,1) + 1j*X(k+1,2);
    c(K+1-k) = conj(c(K+1+k));
end
ka = -K:K;
T = zeros(2*K+1);
for a = 1:2*K+1
    for b = 1:2*K+1
        h = ka(a) - ka(b);
        if abs(h) <= K, T(a,b) = c(h+K+1); end
    end
end
end
