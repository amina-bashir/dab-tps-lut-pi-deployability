clear; clc; close all;
run('dab_params.m');
root = fileparts(mfilename('fullpath'));
fdir = fullfile(root, 'figures');
if ~isfolder(fdir), mkdir(fdir); end

om = omega;
Ts = 1/f_sw;
FILES = {'pointwise_fullTPS_NG70.mat'};
NAMES = {'ZVS-constrained optimal modulation, full TPS, NG70'};

H_GAIN    = 0.5;
ROWSTEP   = 3;
FC_LIST   = [300 200 150 100 75 50 30 20];
PM_MIN    = 45;
GM_MIN    = 6;
I_MAX     = 125;
P_RATED   = 50e3;
RAMP_RATE = 8e3;
U_LIM     = 20;
N_DELAY   = 1;

BANDS = [200 500; 500 700; 700 1000];
RAMP  = [true true true];
CHECK = [false true true];
NB = size(BANDS,1);

S  = load(FILES{1});
Ig = S.mod_cmd(:);
Vg = S.mod_V2(:);
D1 = S.mod_D1;  D2 = S.mod_D2;  D3 = S.mod_D3;
NI = numel(Ig);  NV = numel(Vg);
rv = S.IRMS_T(isfinite(S.IRMS_T));

fprintf('%s\n%s\n', NAMES{1}, FILES{1});
fprintf('stored RMS: median %.6f A | p90 %.6f A | max %.6f A\n', median(rv), pct(rv,90), max(rv));

G = [];  GP = [];  GV = [];
for j = 1:NV
    P2 = param;  P2.V2 = Vg(j);
    for i = 1:ROWSTEP:NI
        for f = [0 0.5]
            c  = Ig(i) + f*(Ig(min(i+1,NI)) - Ig(i));
            cp = min(max(c + H_GAIN, Ig(1)), Ig(end));
            cm = min(max(c - H_GAIN, Ig(1)), Ig(end));
            if cp-cm < 1e-9, continue, end
            ip = samp_exact(look2(Vg,Ig,D1,D2,D3,Vg(j),cp), Vg(j), P2, om);
            im = samp_exact(look2(Vg,Ig,D1,D2,D3,Vg(j),cm), Vg(j), P2, om);
            if isfinite(ip) && isfinite(im)
                G(end+1,1)  = (ip-im)/(cp-cm);
                GP(end+1,1) = c*Vg(j);
                GV(end+1,1) = Vg(j);
            end
        end
    end
end

Vtr = linspace(min(Vg), max(Vg), 600);
gtr = nan(size(Vtr));
for q = 1:numel(Vtr)
    P2 = param;  P2.V2 = Vtr(q);
    c  = min(min(I_MAX, Ig(end)-H_GAIN), P_RATED/Vtr(q));
    cp = min(max(c + H_GAIN, Ig(1)), Ig(end));
    cm = min(max(c - H_GAIN, Ig(1)), Ig(end));
    if cp-cm < 1e-9, continue, end
    ip = samp_exact(look2(Vg,Ig,D1,D2,D3,Vtr(q),cp), Vtr(q), P2, om);
    im = samp_exact(look2(Vg,Ig,D1,D2,D3,Vtr(q),cm), Vtr(q), P2, om);
    if isfinite(ip) && isfinite(im)
        gtr(q) = (ip-im)/(cp-cm);
    end
end

ii0 = max(3, round(NI/4));
RES = struct([]);
for b = 1:NB
    lo = BANDS(b,1);  hi = BANDS(b,2);
    fprintf('\n=== BAND %d: %g - %g V ===\n', b, lo, hi);

    sb = GV >= lo & GV <= hi;
    Gb = G(sb);  GPb = GP(sb);
    gpos = Gb(Gb > 0);
    gmed = median(gpos);  glo = pct(gpos,10);  ghi = pct(gpos,90);
    fprintf('surface gain: median(all) %.6f | negative %.3f%% (%d/%d) | negative inside 50 kW %d\n', ...
            median(Gb), 100*mean(Gb<0), sum(Gb<0), numel(Gb), sum(Gb<0 & GPb<=P_RATED+1e-6));
    fprintf('positive gain: median %.6f | p10 %.6f | p90 %.6f | spread p90/p10 %.3f\n', gmed, glo, ghi, ghi/glo);

    mb = Vtr >= lo & Vtr <= hi;
    Vb = Vtr(mb);  gb = gtr(mb);
    valid = isfinite(gb);
    neg = valid & (gb < 0);
    dd = diff([false neg false]);
    s0 = find(dd == 1);
    s1 = find(dd == -1) - 1;
    if isempty(s0)
        widest = 0;  fcmax = inf;
        fprintf('profile reversal: NONE\n');
    else
        widths = Vb(min(s1+1,numel(Vb))) - Vb(max(s0-1,1));
        widest = max(widths);
        trev   = widest/RAMP_RATE;
        fcmax  = 1/(2*pi*trev);
        fprintf('profile reversal: %.3f%% | %d stretches | widest %.3f V\n', 100*mean(neg(valid)), numel(s0), widest);
        fprintf('crossing %.3f ms @ %.1f V/ms -> fc_max %.3f Hz\n', trev*1e3, RAMP_RATE/1e3, fcmax);
    end

    [~,iv0] = min(abs(Vg - (lo+hi)/2));
    P0  = param;  P0.V2 = Vg(iv0);
    Vg0 = [Vi; P0.V2];
    Jc  = [gradient(D1(:,iv0),Ig), gradient(D2(:,iv0),Ig), gradient(D3(:,iv0),Ig)];
    Dk  = [D1(ii0,iv0); D2(ii0,iv0); mod(D3(ii0,iv0),2)];
    [xk,Gk] = fixed_point_G(Dk,Vg0,P0,om);
    eD = 1e-6;
    Bk = zeros(5,3);
    for j = 1:3
        Dp = Dk;  Dm = Dk;
        Dp(j) = Dp(j) + eD;
        Dm(j) = Dm(j) - eD;
        Bk(:,j) = (period_map(xk,Dp,Vg0,P0,om) - period_map(xk,Dm,Vg0,P0,om))/(2*eD);
    end
    Gd_raw  = ss(Gk, -Bk*Jc(ii0,:).', [0 0 1 0 0], 0, Ts) * tf(1, [1 zeros(1,N_DELAY)], Ts);
    k_local = dcgain(Gd_raw);
    Gd = Gd_raw*(gmed/k_local);

    fprintf('design point %g V, %.2f A | local Jacobian %.6f; scaled to %.6f\n', Vg(iv0), Ig(ii0), k_local, gmed);
    fprintf('%8s %12s %12s | %8s %8s | %8s %8s | %7s %s\n', 'fc','Kp','Ki','PMlo','GMlo','PMhi','GMhi','<=rev','usable');
    fprintf('%s\n', repmat('-',1,88));
    pick = [];
    for fc = FC_LIST
        C0  = pidtune(Gd, 'PI', 2*pi*fc);
        Cpi = pid(C0.Kp, C0.Ki, 0, 0, Ts);
        [pm1,gm1] = margins_local(Cpi*Gd*(glo/gmed));
        [pm2,gm2] = margins_local(Cpi*Gd*(ghi/gmed));
        good = pm1>=PM_MIN && gm1>=GM_MIN && pm2>=PM_MIN && gm2>=GM_MIN && fc<=fcmax;
        fprintf('%8.0f %12.5g %12.5g | %8.1f %8.1f | %8.1f %8.1f | %7s %s\n', ...
                fc, C0.Kp, C0.Ki, pm1, gm1, pm2, gm2, yn(fc<=fcmax), yn(good));
        if good && isempty(pick)
            pick = struct('fc',fc,'kp',C0.Kp,'ki',C0.Ki,'pmlo',pm1,'gmlo',gm1,'pmhi',pm2,'gmhi',gm2,'pass',true);
        end
    end
    if isempty(pick)
        fc  = FC_LIST(end);
        C0  = pidtune(Gd, 'PI', 2*pi*fc);
        Cpi = pid(C0.Kp, C0.Ki, 0, 0, Ts);
        [pm1,gm1] = margins_local(Cpi*Gd*(glo/gmed));
        [pm2,gm2] = margins_local(Cpi*Gd*(ghi/gmed));
        pick = struct('fc',fc,'kp',C0.Kp,'ki',C0.Ki,'pmlo',pm1,'gmlo',gm1,'pmhi',pm2,'gmhi',gm2,'pass',false);
        fprintf('no FC_LIST entry passed all constraints; reporting %g Hz\n', fc);
    end
    fprintf('SELECTED: fc %.0f Hz | Kp %.8g | Ki %.8g\n', pick.fc, pick.kp, pick.ki);

    if RAMP(b)
        I_EFF = min(I_MAX, Ig(end)-H_GAIN);
        Nsim  = ceil((hi-lo)/RAMP_RATE/Ts);
        [t,y_ff,r_t] = sweep(lo,RAMP_RATE,Nsim,0,0,Vg,Ig,D1,D2,D3,param,Vi,om,Ts,N_DELAY,U_LIM,I_EFF,P_RATED);
        [~,y_cl,~,u_cl,m_cl] = sweep(lo,RAMP_RATE,Nsim,pick.kp,pick.ki,Vg,Ig,D1,D2,D3,param,Vi,om,Ts,N_DELAY,U_LIM,I_EFF,P_RATED);
        eff = y_ff - r_t;
        ecl = y_cl - r_t;
        fprintf('ramp %g -> %g V @ %.1f V/ms | %d periods (%.1f ms)\n', lo, hi, RAMP_RATE/1e3, Nsim, Nsim*Ts*1e3);
        fprintf('open loop   : mean|e| %.6f A | worst|e| %.6f A | peak I %.6f A\n', mean(abs(eff)), max(abs(eff)), max(y_ff));
        fprintf('closed loop : mean|e| %.6f A | worst|e| %.6f A | peak I %.6f A\n', mean(abs(ecl)), max(abs(ecl)), max(y_cl));
        fprintf('closed loop : peak|u| %.6f A | clamp %d/%d | mode changes %d\n', ...
                max(abs(u_cl)), sum(abs(u_cl)>=U_LIM-1e-9), Nsim, sum(diff(m_cl)~=0));

        f = figure('Color','w','Position',[100 100 700 420]);
        ax = axes(f);  hold(ax,'on');  grid(ax,'on');
        h1 = plot(ax, t*1e3, y_ff, '-', 'Color', [0.85 0.10 0.10], 'LineWidth', 1.4);
        h2 = plot(ax, t*1e3, y_cl, '-', 'Color', [0 0.447 0.741], 'LineWidth', 1.1);
        h3 = plot(ax, t*1e3, r_t, 'k--', 'LineWidth', 1.4);
        uistack(h1,'top');  uistack(h3,'top');
        xlim(ax, [0 t(end)*1e3]);
        xlabel(ax, 't (ms)');  ylabel(ax, 'i_{2,s} (A)');
        legend(ax, [h1 h2 h3], {'open loop','PI','I_{ref}'}, 'Location','best', 'FontName','Times New Roman', 'FontSize', 10);
        set(ax, 'FontName','Times New Roman', 'FontSize', 11, 'Box','on', 'TickDir','in', ...
                'LineWidth', 0.75, 'XColor','k', 'YColor','k');
        ax.Toolbar.Visible = 'off';
        drawnow;
        fn = fullfile(fdir, sprintf('pi_ramp_%d_%d.png', lo, hi));
        exportgraphics(f, fn, 'Resolution', 300);
        fprintf('saved %s\n', fn);
    else
        eff = NaN;  ecl = NaN;  y_cl = NaN;  u_cl = NaN;  m_cl = NaN;  Nsim = 0;
        t = [];  r_t = [];  y_ff = [];
        fprintf('ramp skipped: CC part sits on the table top, no LUT headroom for the PI\n');
    end

    R = struct();
    R.band = [lo hi];
    R.Vdesign = Vg(iv0);  R.Idesign = Ig(ii0);
    R.gainNegPct   = 100*mean(Gb<0);
    R.gainNegCount = sum(Gb<0);
    R.gainNegInEnv = sum(Gb<0 & GPb<=P_RATED+1e-6);
    R.gainN = numel(Gb);
    R.gmed = gmed;  R.glo = glo;  R.ghi = ghi;
    R.profileNegPct = 100*mean(neg(valid));
    R.nStretches = numel(s0);
    R.widestV = widest;
    R.fcmax = fcmax;
    R.kLocal = k_local;
    R.fc = pick.fc;  R.kp = pick.kp;  R.ki = pick.ki;
    R.pmlo = pick.pmlo;  R.gmlo = pick.gmlo;
    R.pmhi = pick.pmhi;  R.gmhi = pick.gmhi;
    R.marginsPass = pick.pass;
    R.rampMeanFF  = mean(abs(eff));
    R.rampWorstFF = max(abs(eff));
    R.rampMeanCL  = mean(abs(ecl));
    R.rampWorstCL = max(abs(ecl));
    R.rampPeakI   = max(y_cl);
    R.rampPeakU   = max(abs(u_cl));
    R.rampClamp   = sum(abs(u_cl) >= U_LIM-1e-9) * RAMP(b);
    R.rampN       = Nsim;
    R.rampModeChanges = sum(diff(m_cl) ~= 0);
    R.ok = pick.pass && R.nStretches == 0 && (~CHECK(b) || (R.rampWorstCL <= U_LIM && R.rampClamp == 0));
    R.t = t;  R.ref = r_t;  R.yff = y_ff;  R.ycl = y_cl;  R.ucl = u_cl;
    RES = [RES, R];
end

fprintf('\n%-12s %8s %8s %8s %10s %8s %12s %12s %8s %8s %10s %7s  %s\n', 'band (V)', 'gmed', 'p10', 'p90', 'reversals', 'fc (Hz)', 'Kp', 'Ki', 'PMmin', 'GMmin', 'worst|e|', 'clamp', 'verdict');
fprintf('%s\n', repmat('-',1,137));
for b = 1:NB
    R = RES(b);
    if CHECK(b), es = sprintf('%10.3f %7d', R.rampWorstCL, R.rampClamp); else, es = sprintf('%10s %7s', '-', '-'); end
    fprintf('%4g - %-5g %8.4f %8.4f %8.4f %10d %8.0f %12.5g %12.5g %8.1f %8.1f %s  %s\n', ...
            R.band(1), R.band(2), R.gmed, R.glo, R.ghi, R.nStretches, R.fc, R.kp, R.ki, ...
            min(R.pmlo,R.pmhi), min(R.gmlo,R.gmhi), es, yn2(R.ok, 'OK', 'PI FAILS'));
end
save('controller_ramp.mat','RES','BANDS','FILES','NAMES','H_GAIN','ROWSTEP','FC_LIST','PM_MIN','GM_MIN', ...
     'I_MAX','P_RATED','RAMP_RATE','U_LIM','N_DELAY','-v7.3');

fprintf('saved: controller_ramp.mat\n');

function I = samp_exact(D,V2,P,om)
    [tk,S1,S2] = sched(jitter(D),om);
    x = fixed_point(tk,S1,S2,[P.V1;V2],P);
    if isempty(x), I = NaN; else, I = -x(3); end
end

function [pm,gm] = margins_local(L)
    S = allmargin(L);
    if isempty(S.PhaseMargin), pm = NaN; else, pm = min(S.PhaseMargin); end
    if isempty(S.GainMargin), gm = Inf; else, gm = 20*log10(min(S.GainMargin)); end
    if ~S.Stable, pm = -inf; gm = -inf; end
end

function y = yn(x)
    if x, y = 'yes'; else, y = 'no'; end
end

function s = yn2(c,a,b)
    if c, s = a; else, s = b; end
end

function D = look2(Vg,Ig,T1,T2,T3,V2,c)
    vq = min(max(V2,min(Vg)),max(Vg));
    cq = min(max(c,min(Ig)),max(Ig));
    D = [interp2(Vg(:).',Ig(:),T1,vq,cq,'linear'), ...
         interp2(Vg(:).',Ig(:),T2,vq,cq,'linear'), ...
         interp2(Vg(:).',Ig(:),T3,vq,cq,'linear')];
end

function [x,G] = fixed_point_G(D,Vg,P,om)
    [tk,S1,S2] = sched(jitter(D),om);
    G = eye(5);  H = zeros(5,1);
    for k = 1:8
        [E,p] = imap(S1(k),S2(k),tk(k),Vg,P);
        G = E*G;  H = E*H + p;
    end
    x = (eye(5) - G)\H;
end

function x = period_map(x,D,Vg,P,om)
    [tk,S1,S2] = sched(jitter(D),om);
    for k = 1:8
        [E,p] = imap(S1(k),S2(k),tk(k),Vg,P);
        x = E*x + p;
    end
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

function [t,y,r,u,md] = sweep(V20,RATE,N,kp,ki,V2v,cmd,T1,T2,T3,P0,V1,om,Ts,ND,UL,Imax,Prated)
    P = P0;
    t = (0:N-1)*Ts;
    v = V20 + RATE*t;
    r = min(Imax, Prated./v);
    y = zeros(1,N);  u = y;  md = y;
    P.V2 = v(1);
    D0 = look2(V2v,cmd,T1,T2,T3,v(1),r(1));
    x  = fixed_point_l([D0(1);D0(2);mod(D0(3),2)],[V1;v(1)],P,om);
    ui  = 0;
    buf = repmat(r(1),1,max(ND,1));
    for k = 1:N
        P.V2 = v(k);
        Vg = [V1; v(k)];
        e  = r(k) - (-x(3));
        uu = kp*e + ui;
        us = max(min(uu,UL),-UL);
        if abs(uu-us) < 1e-12 || sign(e) ~= sign(uu)
            ui = ui + ki*Ts*e;
        end
        ck = min(max(r(k)+us,min(cmd)),max(cmd));
        if ND > 0
            ca  = buf(end);
            buf = [ck, buf(1:end-1)];
        else
            ca = ck;
        end
        D = look2(V2v,cmd,T1,T2,T3,v(k),ca);
        x = period_map_l(x,[D(1);D(2);mod(D(3),2)],Vg,P,om);
        y(k)  = -x(3);
        u(k)  = us;
        md(k) = mode_of(D(1),D(2),D(3));
    end
end

function m = mode_of(D1,D2,D3)
    t = 1e-9;  m = 0;
    if D3 >= -t
        if     D1>=D2-t && D3<=D1-D2+t,                  m=1;
        elseif D2>=D1-t && D3>=1+D1-D2-t && D3<=1+t,     m=2;
        elseif D2<=1-D1+t && D3>=D1-t && D3<=1-D2+t,     m=3;
        elseif D2>=1-D3-t && D2<=1-D3+D1+t && D3>=D1-t,  m=4;
        elseif D2>=D1-D3-t && D2<=1-D3+t && D3<=D1+t,    m=5;
        elseif D1>=1-D2-t && D3>=1-D2-t && D3<=D1+t,     m=6;
        end
    end
end

function x = fixed_point_l(D,Vg,P,om)
    [tk,S1,S2] = sched(D,om);
    G = eye(5);  H = zeros(5,1);
    for k = 1:8
        [E,p] = imap(S1(k),S2(k),tk(k),Vg,P);
        G = E*G;  H = E*H + p;
    end
    x = (eye(5) - G)\H;
end

function x = period_map_l(x,D,Vg,P,om)
    [tk,S1,S2] = sched(D,om);
    for k = 1:8
        [E,p] = imap(S1(k),S2(k),tk(k),Vg,P);
        x = E*x + p;
    end
end