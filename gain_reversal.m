clear; clc; close all;
run('dab_params.m');
root = fileparts(mfilename('fullpath'));
fdir = fullfile(root, 'figures');
if ~isfolder(fdir), mkdir(fdir); end

om = omega;
FILES = {'pointwise_fullTPS_NG70.mat'};
NAMES = {'ZVS-constrained optimal modulation, full TPS, NG70'};

H_GAIN    = 0.5;
ROWSTEP   = 3;
I_MAX     = 125;
P_RATED   = 50e3;
RAMP_RATE = 8e3;
N_TRAJ    = 600;

S  = load(FILES{1});
Ig = S.mod_cmd(:);
Vg = S.mod_V2(:);
D1 = S.mod_D1;  D2 = S.mod_D2;  D3 = S.mod_D3;
NI = numel(Ig);  NV = numel(Vg);
rr = S.IRMS_T(isfinite(S.IRMS_T));

fprintf('%s\n%s\n', NAMES{1}, FILES{1});
fprintf('table: %d x %d\n', NI, NV);
fprintf('stored RMS: median %.6f A | p90 %.6f A | max %.6f A\n', median(rr), pct(rr,90), max(rr));

G = [];  GI = [];  GV = [];
t0 = tic;
for j = 1:NV
    P2 = param;  P2.V2 = Vg(j);
    for i = 1:ROWSTEP:NI
        for f = [0 0.5]
            c  = Ig(i) + f*(Ig(min(i+1,NI)) - Ig(i));
            cp = min(max(c + H_GAIN, Ig(1)), Ig(end));
            cm = min(max(c - H_GAIN, Ig(1)), Ig(end));
            if cp-cm < 1e-9, continue, end
            ip = avg_exact(look2(Vg,Ig,D1,D2,D3,Vg(j),cp), Vg(j), P2, om);
            im = avg_exact(look2(Vg,Ig,D1,D2,D3,Vg(j),cm), Vg(j), P2, om);
            if isfinite(ip) && isfinite(im)
                G(end+1,1)  = (ip-im)/(cp-cm);
                GI(end+1,1) = c;
                GV(end+1,1) = Vg(j);
            end
        end
    end
end
gpos = G(G > 0);
negSurface = G < 0;
inEnv = GI.*GV <= P_RATED + 1e-6;

fprintf('\n2-D composed gain over surface:\n');
fprintf('  samples                    : %d\n', numel(G));
fprintf('  median gain, all samples   : %.6f\n', median(G));
fprintf('  negative gain              : %d / %d = %.3f%%\n', sum(negSurface), numel(G), 100*mean(negSurface));
fprintf('  minimum sampled gain       : %.6f\n', min(G));
if ~isempty(gpos)
    fprintf('  positive gain median       : %.6f\n', median(gpos));
    fprintf('  positive gain p10 / p90    : %.6f / %.6f\n', pct(gpos,10), pct(gpos,90));
end
fprintf('  runtime                    : %.2f s\n', toc(t0));
fprintf('  negative inside 50 kW      : %d   (outside, never commanded: %d)\n', ...
        sum(negSurface & inEnv), sum(negSurface & ~inEnv));
for q = find(negSurface).'
    fprintf('     V2 %5.0f V  I %7.2f A  gain %9.3f  %s\n', GV(q), GI(q), G(q), ...
            ternary_s(inEnv(q), 'inside', 'outside 50 kW'));
end

Vtr = linspace(min(Vg), max(Vg), N_TRAJ);
gtr = nan(size(Vtr));
Itr = nan(size(Vtr));
for q = 1:numel(Vtr)
    P2 = param;  P2.V2 = Vtr(q);
    c  = min(min(I_MAX, Ig(end)-H_GAIN), P_RATED/Vtr(q));
    Itr(q) = c;
    cp = min(max(c + H_GAIN, Ig(1)), Ig(end));
    cm = min(max(c - H_GAIN, Ig(1)), Ig(end));
    if cp-cm < 1e-9, continue, end
    ip = avg_exact(look2(Vg,Ig,D1,D2,D3,Vtr(q),cp), Vtr(q), P2, om);
    im = avg_exact(look2(Vg,Ig,D1,D2,D3,Vtr(q),cm), Vtr(q), P2, om);
    if isfinite(ip) && isfinite(im)
        gtr(q) = (ip-im)/(cp-cm);
    end
end
valid = isfinite(gtr);
neg = valid & (gtr < 0);
d  = diff([false neg false]);
s0 = find(d == 1);
s1 = find(d == -1) - 1;

fprintf('\ncharging-profile gain, Iref = min(122, 50 kW/V2), %g-%g V:\n', min(Vg), max(Vg));
fprintf('  valid profile samples      : %d / %d\n', sum(valid), numel(gtr));
fprintf('  median profile gain        : %.6f\n', median(gtr(valid)));
fprintf('  negative profile fraction  : %.3f%%\n', 100*mean(neg(valid)));
fprintf('  reversal stretches         : %d\n', numel(s0));
if isempty(s0)
    widths = [];  widest = 0;  trev = 0;  fcmax = inf;
    fprintf('  widest reversal            : NONE\n');
else
    widths = Vtr(min(s1+1,numel(Vtr))) - Vtr(max(s0-1,1));
    widest = max(widths);
    trev   = widest/RAMP_RATE;
    fcmax  = 1/(2*pi*trev);
    fprintf('  widest reversal            : %.3f V\n', widest);
    fprintf('  crossing time @ 8 V/ms     : %.3f ms\n', trev*1e3);
    fprintf('  corresponding fc bound     : %.3f Hz\n', fcmax);
    for r = 1:numel(s0)
        fprintf('    %2d) %.3f .. %.3f V  (%.3f V)\n', r, Vtr(s0(r)), Vtr(s1(r)), Vtr(s1(r)) - Vtr(s0(r)));
    end
end

RES = struct();
RES.name = NAMES{1};               RES.file = FILES{1};
RES.medianRMS = median(rr);        RES.p90RMS = pct(rr,90);          RES.maxRMS = max(rr);
RES.G = G;  RES.G_Icmd = GI;  RES.G_V2 = GV;
RES.surfaceMedianGain = median(G); RES.surfaceNegCount = sum(negSurface);
RES.surfaceNegInEnv = sum(negSurface & inEnv);
RES.surfaceNSamples = numel(G);    RES.surfaceNegPct = 100*mean(negSurface);
RES.surfaceMinGain = min(G);
RES.Vtraj = Vtr;  RES.Itraj = Itr;  RES.Gtraj = gtr;
RES.profileNegPct = 100*mean(neg(valid));
RES.nStretches = numel(s0);  RES.stretchStartIdx = s0;  RES.stretchEndIdx = s1;
RES.stretchWidthsV = widths; RES.widestV = widest;  RES.trev_s = trev;  RES.fcmax_Hz = fcmax;

ctrl_ok = RES.surfaceNegInEnv == 0 && RES.nStretches == 0;
fprintf('\nVERDICT (gain)              : %s\n', ternary_s(ctrl_ok, ...
        'CONTROLLABLE - no negative gain, no reversal', 'NOT CONTROLLABLE - negative gain / reversal present'));
save('gain_reversal.mat','RES','FILES','NAMES','H_GAIN','ROWSTEP','I_MAX','P_RATED','RAMP_RATE','N_TRAJ','-v7.3');

f1 = figure('Color','w','Position',[100 100 700 420]);
ax = axes(f1);  hold(ax,'on');
plot(ax, Vtr, gtr, '-', 'Color', [0 0 0.90], 'LineWidth', 1.4);
yline(ax, 0, 'k--', 'LineWidth', 1);
xline(ax, [500 700], 'k:', 'LineWidth', 1);
xlabel(ax, 'V_2 (V)');  ylabel(ax, 'd\itI\rm_{del} / d\itI\rm_{cmd}');
style(ax);
save_png(f1, fullfile(fdir, 'gain_profile.png'));

fprintf('saved: gain_reversal.mat\n');

function style(ax)
set(ax, 'FontName','Times New Roman', 'FontSize', 11, 'Box','on', 'TickDir','in', ...
        'LineWidth', 0.75, 'XColor','k', 'YColor','k');
ax.Toolbar.Visible = 'off';
end

function save_png(f, fn)
drawnow;
exportgraphics(f, fn, 'Resolution', 300);
fprintf('saved %s\n', fn);
end

function D = look2(Vg,Ig,T1,T2,T3,V2,c)
    vq = min(max(V2,min(Vg)),max(Vg));
    cq = min(max(c,min(Ig)),max(Ig));
    D = [interp2(Vg(:).',Ig(:),T1,vq,cq,'linear'), ...
         interp2(Vg(:).',Ig(:),T2,vq,cq,'linear'), ...
         interp2(Vg(:).',Ig(:),T3,vq,cq,'linear')];
end

function I = avg_exact(D,V2,P,om)
    [tk,S1,S2] = sched(jitter(D),om);
    x = fixed_point(tk,S1,S2,[P.V1;V2],P);
    if isempty(x), I = NaN; return, end
    Q = 0;  T = sum(tk);
    for k = 1:8
        if tk(k) <= 0, continue, end
        dt = tk(k)/6;
        [E,p] = imap(S1(k),S2(k),dt,[P.V1;V2],P);
        for q = 1:6
            xn = E*x + p;
            Q  = Q + 0.5*(x(3) + xn(3))*dt;
            x  = xn;
        end
    end
    I = -Q/T;
end

function D = jitter(D)
    D = D(:).';
    D(1) = min(max(D(1),3e-9),1-3e-9) - 1.3e-9;
    D(2) = min(max(D(2),3e-9),1-3e-9) - 0.7e-9;
    D(3) = D(3) + 0.9e-9;
end

function [tk,S1,S2] = sched(D,om)
    p1 = D(1)*pi;  p2 = D(2)*pi;  p3 = mod(D(3),2)*pi;
    e  = [0, mod(p1,2*pi), pi, mod(pi+p1,2*pi), mod(p3,2*pi), ...
          mod(p3+p2,2*pi), mod(pi+p3,2*pi), mod(pi+p3+p2,2*pi)];
    es = sort(e);  tk = zeros(1,8);  S1 = tk;  S2 = tk;
    for k = 1:8
        lo = es(k);  hi = es(mod(k,8)+1);
        if k == 8, hi = hi + 2*pi; end
        md = (lo+hi)/2;
        S1(k) = 0.5*(sign(sin(md)) - sign(sin(md-p1)));
        S2(k) = 0.5*(sign(sin(md-p3)) - sign(sin(md-p3-p2)));
        tk(k) = (hi-lo)/om;
    end
end

function x = fixed_point(tk,S1,S2,Vg,P)
    G = eye(5);  H = zeros(5,1);
    for k = 1:8
        [E,p] = imap(S1(k),S2(k),tk(k),Vg,P);
        G = E*G;  H = E*H + p;
    end
    if max(abs(eig(G))) > 1-1e-12, x = []; return, end
    x = (eye(5) - G)\H;
end

function [E,p] = imap(s1,s2,dt,Vg,P)
    fv = [0; Vg(1)/P.L1; Vg(2)/P.L2; 0; 0];
    A  = [-P.R/P.L, 0, 0, s1/P.L, -P.n*s2/P.L
           0, -P.R1/P.L1, 0, -1/P.L1, 0
           0, 0, -P.R2/P.L2, 0, -1/P.L2
          -s1/P.C1, 1/P.C1, 0, 0, 0
           P.n*s2/P.C2, 0, 1/P.C2, 0, 0];
    M = expm([A, fv; zeros(1,6)]*dt);
    E = M(1:5,1:5);
    p = M(1:5,6);
end

function q = pct(x,p)
    x = sort(x(isfinite(x)));
    if isempty(x), q = NaN; return, end
    if numel(x) == 1, q = x; return, end
    r = 1 + (numel(x)-1)*p/100;
    lo = floor(r);  hi = ceil(r);
    if lo == hi, q = x(lo); else, q = x(lo) + (r-lo)*(x(hi)-x(lo)); end
end

function s = ternary_s(c,a,b)
    if c, s = a; else, s = b; end
end