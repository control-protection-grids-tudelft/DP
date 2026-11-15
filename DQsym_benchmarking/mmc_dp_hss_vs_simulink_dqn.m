%% Open-loop MMC: DQsym (Simulink) vs abc ODE, DP and HSS
% Simulates the DQsym model of the open-loop MMC and compares it with the abc
% ODE, DP and HSS models built by OL_TD_DP_HSSv2_cos.m:
%
%   Fig. 6   periodic steady state over one fundamental period
%   Fig. 7   response over three fundamental periods after a 1 % step in the
%            modulation amplitude, and error w.r.t. the abc ODE
%   Fig. 8   harmonic magnitudes of i_u,a and v_Cu,a in periodic steady state
%   Table 4  the same harmonic magnitudes
%
% Signals logged in Simulink (To Workspace, Timeseries, sample time Ts_log):
%   sl_iu  3x1  upper arm currents [ITa; ITb; ITc]
%   sl_il  3x1  lower arm currents [IBa; IBb; IBc]
%   sl_vc  6x1  capacitor voltages [u_a u_b u_c l_a l_b l_c]
% Current convention: i_u from dc+ to the ac node, i_l from the ac node to dc-.
clc; clear all

if ~exist('A_hss','var') || ~exist('solTD','var')
    run('OL_TD_DP_HSS_dqn.m');
end
load('SOL.mat', 'SOL');                          % DQsym_SS solver matrices

mdl    = 'OL_MMC_dqn_GIThub';
figDir = 'figs';

%% 1. Simulink model parameters

phMod  = pi/2 + delta + [0, -2*pi/3, 2*pi/3];    % Sine Wave block phases
zSS0   = devalPeriodic(solTD, Tsim, Tp, 0);      % initial conditions from the abc ODE
iu0    = zSS0(1:3);
il0    = zSS0(4:6);
vC0_sl = zSS0(7:12);

Ts_log   = Tp/M;          % logging sample time
Ts_dqn   = 2e-5;
Tss      = 0.6;           % run time before the step [s]
Tobs     = 1.5;           % run time after the step [s]
tStep    = Tss;
stepRel  = 0.01;          % relative step in m
mStepRel = 0;
Vpeak    = Vac_pk;
N        = 1;
Vg_sim   = complex([[0 1i*Vpeak; 0 0; 0 0], zeros(3, N+2)]);
Low_arm_gain = Nsm/Csm;
kdcL     = 1;             % lower dc voltage after T_dc, fraction of nominal
T_dc     = 0.66;

% harmonic injection in the grid source (switched off: flag_hAC = 0)
f1 = fg;  su = -1;  fh = 150;  su_f = -1;
flag_hAC = 0;  t_entranceHF = 0.65;

sgn.iu = 1;  sgn.il = 1;  % measurement signs

%% 2. Simulink runs: base case and step in m

load_system(mdl);
in = Simulink.SimulationInput(mdl);
in = in.setModelParameter('StopTime', num2str(Tss + Tobs), ...
    'SolverType', 'Variable-step', 'Solver', 'ode23tb', ...
    'RelTol', '1e-6', 'MaxStep', num2str(Tp/400));
in = in.setVariable('tStep', tStep).setVariable('flag_hAC', 0).setVariable('Vpeak', Vac_pk);
outSL = sim([in.setVariable('mStepRel', 0), in.setVariable('mStepRel', stepRel)], ...
            'ShowProgress', 'off');
for r = 1:2
    if ~isempty(outSL(r).ErrorMessage), error(outSL(r).ErrorMessage); end
end

%% 3. Periodic steady state

tgS = Tss - Tp + (0:M-1)*Tp/M;                   % last period before the step

SL = slStates(outSL(1), tgS, sgn);
OD = armStates(devalPeriodic(solTD, Tsim, Tp, tgS));

zD = zeros(15, M);
for i = 1:15
    zD(i,:) = dp2time(dpget(xTDnames(i), x0, stateNames, K), tgS, wg_v);
end
DP = armStates(zD);

% HSS: with n(t) prescribed, (4) is linear in the states, so the harmonic steady
% state solves 0 = (A_HSS - N) X + B_HSS U over the inputs v^ac and v^dc. The
% minimum-norm solution sets the undamped star-branch dc modes to zero.
U_hss  = Tu*u0;  U_hss = U_hss(pU);
linCol = startsWith(hssUnames, ["vac_", "vudc_", "vldc_"]);
X_hss  = lsqminnorm(A_hss, -B_hss(:, linCol)*U_hss(linCol));
zH = zeros(15, M);
for i = 1:15
    zH(i,:) = hssTraj2time(X_hss, hssXnames, xTDnames(i), tgS, wg_v, K);
end
HS = armStates(zH);

% colours and axes style shared by Figs. 6-8
cODE = [213  94   0]/255;                        % vermillion
cDP  = [230 159   0]/255;                        % orange
cHSS = [204 121 167]/255;                        % reddish purple
cDQ  = [0.25 0.25 0.25];
styleAx = @(ax, fsT) set(ax, 'TickLabelInterpreter','latex', 'FontSize', fsT, 'Box','on', ...
    'LineWidth', 0.5, 'XGrid','on', 'YGrid','on', 'GridAlpha', 0.15, ...
    'Layer','top', 'TickLength', [0.01 0.01]);

%% 4. Fig. 6: periodic steady state, phase a

fs = 8;  lw = 1.0;
tms  = (tgS - tgS(1))*1e3;
sig6 = {'iu',  1,    '$i_{u,a}$ [A]';
        'il',  1,    '$i_{l,a}$ [A]';
        'vCu', 1e-3, '$v_{Cu,a}$ [kV]';
        'vCl', 1e-3, '$v_{Cl,a}$ [kV]'};
nS = size(sig6, 1);

fig6 = figure('Name','Fig. 6 periodic steady state', 'Color','w', ...
              'Units','centimeters', 'Position',[2 2 8.4 6.5]);
tl6  = tiledlayout(fig6, nS, 1, 'TileSpacing','tight', 'Padding','compact');
for s = 1:nS
    q = sig6{s,1};  c = sig6{s,2};
    ax = nexttile(tl6);  hold(ax, 'on');
    h = gobjects(1,4);
    h(2) = plot(ax, tms, c*SL.(q)(1,:), '-',  'Color', 'k',  'LineWidth', 1.2);
    h(1) = plot(ax, tms, c*OD.(q)(1,:), '--', 'Color', cODE, 'LineWidth', lw);
    h(3) = plot(ax, tms, c*DP.(q)(1,:), '-.', 'Color', cDP,  'LineWidth', lw);
    h(4) = plot(ax, tms, c*HS.(q)(1,:), ':',  'Color', cHSS, 'LineWidth', 1.3);
    ylabel(ax, sig6{s,3}, 'Interpreter','latex', 'FontSize', fs);
    xlim(ax, [0 Tp*1e3]);  styleAx(ax, fs);

    % y range with a small margin and three round ticks
    y  = c*[SL.(q)(1,:), OD.(q)(1,:), DP.(q)(1,:), HS.(q)(1,:)];
    yr = [min(y) max(y)];  pad = 0.08*diff(yr);
    ylim(ax, yr + [-pad pad]);
    st = 10^floor(log10(diff(yr)/2));
    tk = unique(round(linspace(yr(1), yr(2), 3)/st)*st);
    ax.YTick = tk(tk >= ax.YLim(1) & tk <= ax.YLim(2));

    if s < nS, ax.XTickLabel = []; end
    if s == 1, hLeg = h; end
end
lg6 = legend(hLeg, {'abc ODE', 'DQsym', 'DP', 'HSS'}, 'Interpreter','latex', ...
             'FontSize', fs, 'Orientation','horizontal', 'Box','off');
lg6.Layout.Tile = 'north';
lg6.ItemTokenSize = [14 8];
xlabel(tl6, '$T$ [ms]', 'Interpreter','latex', 'FontSize', fs);
exportFig(fig6, figDir, 'fig6_steady_state');

%% 5. Fig. 7: step in m, DQsym vs nonlinear abc vs linear DP and HSS
% The abc ODE and the linear models get the same step as DQsym, at the same time. Top: total waveforms over the first three
% periods after the step. Bottom: error w.r.t. the abc ODE, e = x - x_abc.

dtl = Tp/200;
Nt  = round(Tobs/dtl) + 1;
tl  = tStep + (0:Nt-1)*dtl;
S1  = slStates(outSL(2), tl, sgn);

% step in the DP and HSS input coordinates
dm   = stepRel*m;
z3   = zeros(K+1,2);
Uin3 = Uin;  duS = [];
for p = 1:3
    th   = -(p-1)*2*pi/3;
    dNUs = z3;  dNUs(2,:) = -dm/4*[cos(th+delta), sin(th+delta)];
    dNLs = -dNUs;
    Uin3.nu{p} = nNU{p} + dNUs;
    Uin3.nl{p} = nNL{p} + dNLs;
    duS = [duS; dpflat(z3); dpflat(dNUs); dpflat(dNLs)]; %#ok<AGROW>
end
duS = [duS; dpflat(z3); dpflat(z3)];
dUs = double(subs(u, uSym, duS));
dUs_hss = Tu*dUs;  dUs_hss = dUs_hss(pU);

% nonlinear abc model
zE   = devalPeriodic(solTD, Tsim, Tp, tStep);
solS = ode15s(@(t,z) mmcTD(t,z,prm,Uin3,wg_v), [tStep tStep+Tobs], zE, odeOpts);
zS   = deval(solS, tl);
zB   = devalPeriodic(solTD, Tsim, Tp, tl);

% linear models
dX_DP  = stepSim(A0,    B0*dUs,          dtl, Nt);
dX_HSS = stepSim(A_hss, B_hss*dUs_hss,   dtl, Nt);
zDb = zeros(12, Nt);
linDP = zeros(12, Nt);  linHS = zeros(12, Nt);
for i = 1:12
    zDb(i,:)   = dp2time(dpget(xTDnames(i), x0, stateNames, K), tl, wg_v);
    linDP(i,:) = dpTraj2time(dX_DP, stateNames, xTDnames(i), tl, wg_v, K);
    linHS(i,:) = hssTraj2time(dX_HSS, hssXnames, xTDnames(i), tl, wg_v, K);
end

% rows: i_u,a, i_l,a, v_Cu,a, v_Cl,a
ia = [1 4 7 10];
selQ = @(Q) [Q.iu(1,:); Q.il(1,:); Q.vCu(1,:); Q.vCl(1,:)];
tot.SL = selQ(S1);
tot.OD = zS(ia,:);
tot.DP = zDb(ia,:) + linDP(ia,:);
tot.HS = zB(ia,:)  + linHS(ia,:);

fs = 9;  fsT = 8;  fsL = 7;  lw = 0.8;
sc  = [1e-3 1e-3 1e-3 1e-3];                     % totals: kA, kA, kV, kV
scE = [1    1    1e-3 1e-3];                     % errors: A, A, kV, kV
lbT = {'$i_{u,a}$ [kA]', '$i_{l,a}$ [kA]', '$v_{Cu,a}$ [kV]', '$v_{Cl,a}$ [kV]'};
lbE = {'$e_{i_{u,a}}$ [A]', '$e_{i_{l,a}}$ [A]', '$e_{v_{Cu,a}}$ [kV]', '$e_{v_{Cl,a}}$ [kV]'};
tr  = (tl - tStep)*1e3;
nCyc = 3;
iw   = tr <= nCyc*Tp*1e3 + 1e-9;
trP  = tr(iw);
totP = structfun(@(X) X(:, iw), tot, 'UniformOutput', false);

fig7 = figure('Name','Fig. 7 step response', 'Color','w', ...
              'Units','centimeters', 'Position',[2 2 8.5 16.5]);
tl7  = tiledlayout(fig7, 2, 1, 'TileSpacing','compact', 'Padding','compact');
tTop = tiledlayout(tl7, 4, 1, 'TileSpacing','tight', 'Padding','none');  tTop.Layout.Tile = 1;
tBot = tiledlayout(tl7, 4, 1, 'TileSpacing','tight', 'Padding','none');  tBot.Layout.Tile = 2;

for s = 1:4
    ax = nexttile(tTop);  hold(ax, 'on');
    h = gobjects(1,4);
    h(2) = plot(ax, trP, sc(s)*totP.SL(s,:), '-',  'Color', 'k',  'LineWidth', 1.0);
    h(1) = plot(ax, trP, sc(s)*totP.OD(s,:), '--', 'Color', cODE, 'LineWidth', lw);
    h(3) = plot(ax, trP, sc(s)*totP.DP(s,:), '-.', 'Color', cDP,  'LineWidth', lw);
    h(4) = plot(ax, trP, sc(s)*totP.HS(s,:), ':',  'Color', cHSS, 'LineWidth', 1.1);
    ylabel(ax, lbT{s}, 'Interpreter','latex', 'FontSize', fs);
    xlim(ax, [0 nCyc*Tp*1e3]);  styleAx(ax, fsT);
    if s == 1, hLeg = h; end
    if s < 4, set(ax, 'XTickLabel', []); end
end
for s = 1:4
    ax = nexttile(tBot);  hold(ax, 'on');
    yline(ax, 0, '-', 'Color', cODE, 'LineWidth', 0.6);
    plot(ax, trP, scE(s)*(totP.DP(s,:) - totP.OD(s,:)), '-.', 'Color', cDP,  'LineWidth', lw);
    plot(ax, trP, scE(s)*(totP.HS(s,:) - totP.OD(s,:)), ':',  'Color', cHSS, 'LineWidth', 1.1);
    plot(ax, trP, scE(s)*(totP.SL(s,:) - totP.OD(s,:)), '-',  'Color', 'k',  'LineWidth', 1.0);
    ylabel(ax, lbE{s}, 'Interpreter','latex', 'FontSize', fs);
    xlim(ax, [0 nCyc*Tp*1e3]);  styleAx(ax, fsT);
    if s < 4, set(ax, 'XTickLabel', []); end
end
title(tBot, 'Error w.r.t.\ abc ODE, $e = x - x_{\mathrm{abc}}$', ...
      'Interpreter','latex', 'FontSize', fs);
xlabel(tBot, 'T [ms]', 'Interpreter','latex', 'FontSize', fs);
lg7 = legend(hLeg, {'abc', 'DQsym', 'DP', 'HSS'}, 'Interpreter','latex', 'FontSize', fsL, ...
             'Orientation','horizontal', 'Box','off');
lg7.ItemTokenSize = [14 8];
lg7.Layout.Tile = 'north';
exportFig(fig7, figDir, 'fig7_step_response');

%% 6. Fig. 8 and Table 4: harmonic content, h = 0..K
% Amplitude: dc value for h = 0, peak value 2|X_h| for h >= 1.

amp = @(y) [abs(y(1)); 2*abs(y(2:end))];
hAmp = @(y) amp(fourierK(y, tgS, wg_v, K));

fs  = 12;
cBar = [cODE; cDQ; cDP; cHSS];
sig8 = {'iu',  '$|i_{u,a}|$ [A]';
        'vCu', '$|v_{Cu,a}|$ [V]'};

fig8 = figure('Name','Fig. 8 harmonic content', 'Color','w', ...
              'Units','centimeters', 'Position',[2 2 18 11]);
tl8  = tiledlayout(fig8, 2, 1, 'TileSpacing','compact', 'Padding','compact');
for s = 1:size(sig8,1)
    q  = sig8{s,1};
    Ab = [hAmp(OD.(q)(1,:)), hAmp(SL.(q)(1,:)), hAmp(DP.(q)(1,:)), hAmp(HS.(q)(1,:))];
    ax = nexttile(tl8);
    hb = bar(ax, 0:K, Ab, 'grouped', 'BarWidth', 0.85);
    for b = 1:numel(hb)
        hb(b).FaceColor = cBar(b,:);
        hb(b).EdgeColor = 'none';
    end
    set(ax, 'YScale','log', 'XTick', 0:K, 'TickLabelInterpreter','latex', ...
            'FontSize', fs, 'Box','on', 'LineWidth', 0.5, ...
            'YGrid','on', 'XGrid','off', 'GridAlpha', 0.15, 'YMinorGrid','off', ...
            'Layer','top');
    ylim(ax, [max(Ab(:))*1e-6, max(Ab(:))*3]);
    xlim(ax, [-0.6, K + 0.6]);
    ylabel(ax, sig8{s,2}, 'Interpreter','latex', 'FontSize', fs);
    if s == 1, hLeg = hb; end
end
lg8 = legend(hLeg, {'abc ODE', 'DQsym', 'DP', 'HSS'}, 'Interpreter','latex', ...
             'FontSize', fs, 'Orientation','horizontal', 'Box','off');
lg8.Layout.Tile = 'north';
xlabel(tl8, 'harmonic order $h$', 'Interpreter','latex', 'FontSize', fs);
exportFig(fig8, figDir, 'fig8_harmonics');

table4 = array2table([hAmp(OD.iu(1,:)).';  hAmp(SL.iu(1,:)).';  hAmp(DP.iu(1,:)).'; ...
                      [hAmp(OD.vCu(1,:)).'; hAmp(SL.vCu(1,:)).'; hAmp(DP.vCu(1,:)).']/1e3], ...
    'VariableNames', compose('h%d', 0:K), ...
    'RowNames', {'|i_ua| [A] abc ODE', '|i_ua| [A] DQsym', '|i_ua| [A] DP/HSS', ...
                 '|v_Cua| [kV] abc ODE', '|v_Cua| [kV] DQsym', '|v_Cua| [kV] DP/HSS'})

%% Local functions

function exportFig(fig, folder, name)
if ~exist(folder, 'dir'), mkdir(folder); end
drawnow;
exportgraphics(fig, fullfile(folder, [name '.pdf']), 'ContentType','vector');
end

function Q = slStates(out, t, sgn)
% arm currents and capacitor voltages of a Simulink run, resampled on t
Q.iu = sgn.iu * slres(out, 'sl_iu', t);
Q.il = sgn.il * slres(out, 'sl_il', t);
vc   = slres(out, 'sl_vc', t);
Q.vCu = vc(1:3,:);
Q.vCl = vc(4:6,:);
end

function Q = armStates(z)
% arm currents and capacitor voltages from the 15 abc states
Q.iu = z(1:3,:);  Q.il = z(4:6,:);  Q.vCu = z(7:9,:);  Q.vCl = z(10:12,:);
end

function y = slres(out, name, t)
% logged Simulink signal (To Workspace or logsout), resampled on t
v = [];
if any(strcmp(out.who, name))
    v = out.get(name);
elseif any(strcmp(out.who, 'logsout'))
    try
        v = out.get('logsout').getElement(name);
    catch
        v = [];
    end
end
if isempty(v)
    error('Simulink signal "%s" not found: add a To Workspace block (Timeseries).', name);
end
if isa(v, 'Simulink.SimulationData.Signal'), v = v.Values; end
if isa(v, 'timeseries')
    ts = v.Time(:); d = v.Data;
elseif isstruct(v)
    ts = v.time(:); d = v.signals.values;
else
    error('Unsupported format for "%s".', name);
end
if ndims(d) == 3                                 % n x 1 x N  ->  N x n
    d = permute(d, [3 1 2]); d = reshape(d, size(d,1), []);
end
[ts, ia] = unique(ts, 'last');
y = interp1(ts, double(d(ia,:)), t(:), 'linear').';
if any(isnan(y(:)))
    error('Requested times for "%s" are outside the simulated range.', name);
end
end

function X = fourierK(x, t, w, H)
% complex Fourier coefficients X_0..X_H over one period (uniform t)
X = zeros(H+1, 1);
for k = 0:H
    X(k+1) = mean(x .* exp(-1j*k*w*t));
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

function v = dpflat(X)
v = X(1,1);
for k = 1:size(X,1)-1
    v = [v; X(k+1,1); X(k+1,2)]; %#ok<AGROW>
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

function z = devalPeriodic(sol, Tend, Tp, t)
% periodic steady state of the abc ODE at arbitrary times t
tb = Tend - Tp + mod(t - Tend, Tp);
[tbs, ord] = sort(tb);
zs = deval(sol, tbs);
z  = zeros(size(zs));
z(:, ord) = zs;
end

function X = stepSim(A, b, dt, Nt)
% step response of dx/dt = A x + b, exact discretisation
n   = size(A,1);
Phi = expm([A, b; zeros(1, n+1)]*dt);
X   = zeros(n, Nt);
z   = [zeros(n,1); 1];
for q = 2:Nt
    z = Phi*z;
    X(:,q) = z(1:n);
end
end

function y = dpTraj2time(Xr, names, base, t, w, K)
base = char(base);
y = Xr(names == [base '_R0'], :);
for k = 1:K
    R = Xr(names == sprintf('%s_R%d', base, k), :);
    I = Xr(names == sprintf('%s_I%d', base, k), :);
    y = y + 2*(R.*cos(k*w*t) - I.*sin(k*w*t));
end
end

function y = hssTraj2time(Xc, names, base, t, w, K)
y = zeros(1, numel(t));
for k = -K:K
    y = y + Xc(names == base + "_" + k, :) .* exp(1j*k*w*t);
end
y = real(y);
end
