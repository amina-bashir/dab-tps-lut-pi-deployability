clear; clc; close all;
run('dab_params.m');
root = fileparts(mfilename('fullpath'));
fdir = fullfile(root, 'figures');
if ~isfolder(fdir), mkdir(fdir); end

V1   = param.V1;
Coss = 1e-9;
TJ   = 125;
P_RATED = 50e3;
NAM  = {'ZVS-constrained optimal modulation','SPS'};

ds  = load_cab320(root);
T13 = readtable(fullfile(root,'fig13_Esw_vs_Tj.csv'),'TextType','string');
f13 = @(nm) interp1(T13.x(T13.curve==nm), T13.y(T13.curve==nm), [25 TJ], 'linear','extrap');
a = f13("Eon");   kTon  = a(2)/a(1);
a = f13("Eoff");  kToff = a(2)/a(1);
T14 = readtable(fullfile(root,'fig14_Err_vs_Tj.csv'),'TextType','string');
a = interp1(T14.x(T14.curve=="Err_900V"), T14.y(T14.curve=="Err_900V"), [25 TJ], 'linear','extrap');
kTrr = a(2)/a(1);
ES.on  = @(i,v) max(Esw(ds.Eon900, ds.Eon1200, i, v)*kTon, ds.Eoss(v));
ES.off = @(i,v) Esw(ds.Eoff900, ds.Eoff1200, i, v)*kToff;
ES.rr  = @(i,v) Esw(ds.Err900,  ds.Err1200,  i, v)*kTrr;
RP = ds.Rpkg_hi + ds.Rpkg_lo;

K  = load('pointwise_fullTPS_NG70.mat');
Vk = double(K.mod_V2(:));  Ik = double(K.mod_cmd(:));
K1 = double(K.mod_D1);  K2 = double(K.mod_D2);  K3 = double(K.mod_D3);

V  = (min(Vk):10:max(Vk)).';
I  = min(Ik(end), P_RATED./V);
NV = numel(V);
ETA = nan(NV,2);  PC = ETA;  PS = ETA;  IR = ETA;  NH = ETA;
for q = 1:NV
    D = look2(Vk,Ik,K1,K2,K3,V(q),I(q));
    Dq = {[D(1), D(2), D(3)-(D(1)-D(2))/2], [1, 1, sps_shift(V(q)*I(q),V1,n,V(q),Ls,f_sw)]};
    for m = 1:2
        d = Dq{m};
        if ~all(isfinite(d)), continue, end
        [Irms,iE,ok] = tps_edges(d(1),d(2),d(3),V1,n,V(q),Ls,f_sw,Coss);
        I2 = n*Irms;
        Pc = (2*ds.Rds(Irms,TJ) + RP)*Irms^2 + (2*ds.Rds(I2,TJ) + RP)*I2^2 + Rs*Irms^2;
        E = 0;
        for k = 1:8
            if k <= 4, ik = abs(iE(k));  vb = V1;  else, ik = n*abs(iE(k));  vb = V(q); end
            if ok(k), E = E + ES.off(ik,vb); else, E = E + ES.on(ik,vb) + ES.rr(ik,vb); end
        end
        Ps = f_sw*E;  P = V(q)*I(q);
        PC(q,m) = Pc;  PS(q,m) = Ps;  IR(q,m) = Irms;  NH(q,m) = 8 - nnz(ok);
        ETA(q,m) = P/(P + Pc + Ps);
    end
end

fprintf('T_vj %d C, R_DS(on) from Fig. 2, energies Figs. 11-14 (CAB320M17XM3)\n\n', TJ);
q6 = find(abs(V-600) < 1e-9, 1);
fprintf('at 600 V, %.1f A (%.1f kW):\n', I(q6), V(q6)*I(q6)/1e3);
for m = 1:2
    fprintf('  %-36s I_RMS %6.1f A | conduction %7.1f W | switching %6.1f W (%d hard) | eta %.2f %%\n', ...
            NAM{m}, IR(q6,m), PC(q6,m), PS(q6,m), NH(q6,m), 100*ETA(q6,m));
end
fprintf('\nalong the profile:\n');
for m = 1:2
    fprintf('  %-36s eta %.2f - %.2f %% | mean loss %.0f W | hard edges at %d of %d points\n', NAM{m}, ...
            100*min(ETA(:,m)), 100*max(ETA(:,m)), mean(PC(:,m)+PS(:,m),'omitnan'), nnz(NH(:,m)>0), NV);
end
save('efficiency_estimate.mat','V','I','ETA','PC','PS','IR','NH','NAM','TJ');

f = figure('Color','w','Position',[100 100 700 420]);
ax = axes(f);  hold(ax,'on');
h1 = plot(ax, V, 100*ETA(:,1), '-', 'Color', [0 0 0.90], 'LineWidth', 1.4);
h2 = plot(ax, V, 100*ETA(:,2), '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.4);
legend(ax, [h1 h2], NAM, 'Location','best', 'FontName','Times New Roman', 'FontSize', 10);
xlim(ax, [min(V) max(V)]);
set(ax, 'XTick', 200:200:1000);
xlabel(ax, 'V_2 (V)');  ylabel(ax, '\eta (%)');
set(ax, 'FontName','Times New Roman', 'FontSize', 11, 'Box','on', 'TickDir','in', ...
        'LineWidth', 0.75, 'XColor','k', 'YColor','k');
ax.Toolbar.Visible = 'off';
drawnow;
exportgraphics(f, fullfile(fdir,'efficiency_estimate.png'), 'Resolution', 300);
fprintf('saved %s\nsaved: efficiency_estimate.mat\n', fullfile(fdir,'efficiency_estimate.png'));

function E = Esw(f900,f1200,i,v)
    i40 = max(i,40);
    e9 = f900(i40);  e12 = f1200(i40);
    k  = min(max(log(e12/e9)/log(1200/900), 0.5), 2);
    E  = e9*(v/900)^k * min(i/40,1);
end

function d3 = sps_shift(Ptar,V1,n,V2,Ls,fs)
    lo = 0;  hi = 0.5;
    if tps_power(1,1,hi,V1,n,V2,Ls,fs) < Ptar, d3 = NaN; return, end
    for it = 1:60
        md = (lo+hi)/2;
        if tps_power(1,1,md,V1,n,V2,Ls,fs) < Ptar, lo = md; else, hi = md; end
    end
    d3 = (lo+hi)/2;
end

function P = tps_power(D1,D2,D3c,V1,n,V2,Ls,fs)
    nV = n*V2;  s0 = D3c + (D1-D2)/2;  Th = 1/(2*fs);
    E  = mod([0, D1, 1, 1+D1, s0, s0+D2, s0+1, s0+1+D2],2);
    Es = sort(E);
    hseg = diff([Es, Es(1)+2]);  mid = Es + hseg/2;
    vp = vprim(mid,D1,V1);  vs = vsec(mid,s0,D2,nV);
    sl = (vp - vs)/Ls*Th;
    i0 = [0, cumsum(sl(1:7).*hseg(1:7))];
    i0 = i0 - sum(i0.*hseg + 0.5*sl.*hseg.^2)/2;
    P  = sum(vp.*(i0.*hseg + 0.5*sl.*hseg.^2))/2;
end

function [Irms,iE,ok] = tps_edges(D1,D2,D3c,V1,n,V2,Ls,fs,Coss)
    nV = n*V2;  s0 = D3c + (D1-D2)/2;  Th = 1/(2*fs);
    E  = mod([0, D1, 1, 1+D1, s0, s0+D2, s0+1, s0+1+D2],2);
    [Es,ord] = sort(E);
    hseg = diff([Es, Es(1)+2]);  mid = Es + hseg/2;
    vp = vprim(mid,D1,V1);  vs = vsec(mid,s0,D2,nV);
    sl = (vp - vs)/Ls*Th;
    i0 = [0, cumsum(sl(1:7).*hseg(1:7))];
    i0 = i0 - sum(i0.*hseg + 0.5*sl.*hseg.^2)/2;
    Irms = sqrt(max(sum(i0.^2.*hseg + i0.*sl.*hseg.^2 + sl.^2.*hseg.^3/3)/2,0));
    iE = zeros(1,8);  iE(ord) = i0;
    Zp = sqrt(Ls/(2*Coss));  Zs = Zp*n;
    va = [0 1 0 -1 0 1 0 -1];  vb = [1 0 -1 0 1 0 -1 0];  ep = 1e-9;
    ok = false(1,8);
    for k = 1:8
        x = E(k);  ik = iE(k);
        if k <= 4
            Va = va(k)*V1;  Vb = vb(k)*V1;
            vo = (vsec(x-ep,s0,D2,nV) + vsec(x+ep,s0,D2,nV))/2;
            need = sqrt(max((Vb-vo)^2 - (Va-vo)^2, 0))/Zp;
            if Vb > Va, ok(k) = ik <= -need + 1e-2; else, ok(k) = ik >= need - 1e-2; end
        else
            Va = va(k)*nV;  Vb = vb(k)*nV;
            vo = (vprim(x-ep,D1,V1) + vprim(x+ep,D1,V1))/2;
            need = sqrt(max((Vb-vo)^2 - (Va-vo)^2, 0))/Zs;
            if Vb > Va, ok(k) = ik >= need - 1e-2; else, ok(k) = ik <= -need + 1e-2; end
        end
    end
end

function v = vprim(x,D1,V1)
    x = mod(x,2);
    v = V1*((x < D1) - (x >= 1 & x < 1+D1));
end

function v = vsec(x,s0,D2,nV)
    y = mod(x - s0,2);
    v = nV*((y < D2) - (y >= 1 & y < 1+D2));
end

function D = look2(Vg,Ig,T1,T2,T3,V2,c)
    vq = min(max(V2,min(Vg)),max(Vg));
    cq = min(max(c,min(Ig)),max(Ig));
    D = [interp2(Vg(:).',Ig(:),T1,vq,cq,'linear'), ...
         interp2(Vg(:).',Ig(:),T2,vq,cq,'linear'), ...
         interp2(Vg(:).',Ig(:),T3,vq,cq,'linear')];
end