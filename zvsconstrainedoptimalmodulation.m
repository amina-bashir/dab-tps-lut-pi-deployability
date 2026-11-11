clear; clc;
run('dab_params.m');
V1   = param.V1;
Coss = 1e-9;
NG   = 70;
V2_grid = 200:50:1000;
NI = 40;
I_grid = linspace(5, 0.98*125, NI);
NV = numel(V2_grid);
OUT_FILE = 'pointwise_fullTPS_NG70.mat';

fprintf('ZVS-constrained optimal modulation table, NG = %d\n', NG);
fprintf('V1 %g V, n %.2f, Ls %.1f uH, fs %.0f kHz, Coss %.0f pF\n\n', V1, n, Ls*1e6, f_sw/1e3, Coss*1e12);

mod_D1 = nan(NI,NV);  mod_D2 = mod_D1;  mod_D3 = mod_D1;  IRMS_T = mod_D1;
SCORE_T = mod_D1;  MODE_T = mod_D1;  ZVS_OK = false(NI,NV);
t0 = tic;
for iv = 1:NV
    for ii = 1:NI
        H = tps_cloud(I_grid(ii), V2_grid(iv), V1, n, Ls, f_sw, Coss, NG);
        if isempty(H), error('No triplet at %.1f A, %d V.', I_grid(ii), V2_grid(iv)); end
        smax = max(H(:,5));
        Hs = H(H(:,5) == smax, :);
        [r,k] = min(Hs(:,4));
        mod_D1(ii,iv) = Hs(k,1);  mod_D2(ii,iv) = Hs(k,2);  mod_D3(ii,iv) = Hs(k,3);
        IRMS_T(ii,iv) = r;  SCORE_T(ii,iv) = smax;  MODE_T(ii,iv) = Hs(k,6);  ZVS_OK(ii,iv) = (smax == 8);
    end
    fprintf('  V2 = %4d V (%2d/%d)   %.0f s\n', V2_grid(iv), iv, NV, toc(t0));
end
mod_cmd = I_grid;
mod_V2  = V2_grid;
save(OUT_FILE,'mod_cmd','mod_V2','mod_D1','mod_D2','mod_D3','IRMS_T','ZVS_OK','SCORE_T','MODE_T','NG','Coss');

R = IRMS_T(:);
fprintf('\nnodes                : %d\n', numel(R));
fprintf('8/8 ZVS              : %d / %d\n', nnz(ZVS_OK), numel(ZVS_OK));
fprintf('RMS median / mean    : %.2f / %.2f A\n', median(R), mean(R));
fprintf('RMS p10 / p90        : %.2f / %.2f A\n', prctile(R,10), prctile(R,90));
fprintf('RMS min / max        : %.2f / %.2f A\n', min(R), max(R));
fprintf('modes used           : %s   (7 = outside the six classical modes: %d nodes)\n', ...
        mat2str(unique(MODE_T(:))'), nnz(MODE_T == 7));
fprintf('D3 range             : %+.3f .. %+.3f\n', min(mod_D3(:)), max(mod_D3(:)));
fprintf('saved: %s\n', OUT_FILE);

function H = tps_cloud(Icmd, V2, V1, n, Ls, f_sw, Coss, NG)
    Ptar = Icmd*V2;
    g = (1:NG)/NG;
    [A,B] = ndgrid(g,g);  A = A(:);  B = B(:);
    NS = 241;  s = linspace(0,1,NS);
    H = zeros(0,6);
    CH = 400;
    for c0 = 1:CH:numel(A)
        idx = c0:min(c0+CH-1,numel(A));
        a = A(idx);  b = B(idx);  np = numel(idx);
        aa = repmat(a,1,NS);  bb = repmat(b,1,NS);  ss = repmat(s,np,1);
        P = tps_eval(aa(:),bb(:),ss(:),V1,n,V2,Ls,f_sw);
        F = reshape(P,np,NS) - Ptar;
        sg = sign(F);
        [r,q] = find(sg(:,1:end-1).*sg(:,2:end) < 0);
        if isempty(r), continue, end
        lo = s(q).';  hi = s(q+1).';  flo = F(sub2ind(size(F),r,q));
        d1 = a(r);  d2 = b(r);
        for it = 1:45
            md = (lo+hi)/2;
            fm = tps_eval(d1,d2,md,V1,n,V2,Ls,f_sw) - Ptar;
            same = sign(fm) == sign(flo);
            lo(same) = md(same);  hi(~same) = md(~same);
        end
        d3c = (lo+hi)/2;
        [~,Ir,sc] = tps_eval(d1,d2,d3c,V1,n,V2,Ls,f_sw,Coss);
        D3 = d3c + (d1-d2)/2;
        H = [H; d1 d2 D3 Ir sc zeros(numel(d1),1)];
    end
    for k = 1:size(H,1)
        H(k,6) = mode_of(H(k,1),H(k,2),H(k,3));
    end
end

function [P, Irms, score] = tps_eval(D1,D2,D3c,V1,n,V2,Ls,fs,Coss)
    N  = numel(D1);  nV = n*V2;
    D1 = D1(:);  D2 = D2(:);  s0 = D3c(:) + (D1-D2)/2;
    Th = 1/(2*fs);
    E  = [zeros(N,1), D1, ones(N,1), 1+D1, s0, s0+D2, s0+1, s0+1+D2];
    E  = mod(E,2);
    [Es,ord] = sort(E,2);
    hseg = diff([Es, Es(:,1)+2],1,2);
    mid  = Es + hseg/2;
    vp = vprim(mid,D1,V1);  vs = vsec(mid,s0,D2,nV);
    sl = (vp - vs)/Ls*Th;
    i0 = [zeros(N,1), cumsum(sl(:,1:7).*hseg(:,1:7),2)];
    i0 = i0 - sum(i0.*hseg + 0.5*sl.*hseg.^2,2)/2;
    P  = sum(vp.*(i0.*hseg + 0.5*sl.*hseg.^2),2)/2;
    if nargout < 2, return, end
    Irms = sqrt(max(sum(i0.^2.*hseg + i0.*sl.*hseg.^2 + sl.^2.*hseg.^3/3,2)/2,0));
    if nargout < 3, return, end
    Zp = sqrt(Ls/(2*Coss));  Zs = Zp*n;
    iE = zeros(N,8);
    for k = 1:8
        iE(sub2ind([N 8],(1:N).',ord(:,k))) = i0(:,k);
    end
    br = [1 1 1 1 2 2 2 2];
    va = [0 1 0 -1 0 1 0 -1];
    vb = [1 0 -1 0 1 0 -1 0];
    ep = 1e-9;
    score = zeros(N,1);
    for k = 1:8
        x = E(:,k);  ik = iE(:,k);
        if br(k) == 1
            Va = va(k)*V1;  Vb = vb(k)*V1;
            vo = (vsec(x-ep,s0,D2,nV) + vsec(x+ep,s0,D2,nV))/2;
            need = sqrt(max((Vb-vo).^2 - (Va-vo).^2, 0))/Zp;
            if Vb > Va, ok = ik <= -need + 1e-9; else, ok = ik >= need - 1e-9; end
        else
            Va = va(k)*nV;  Vb = vb(k)*nV;
            vo = (vprim(x-ep,D1,V1) + vprim(x+ep,D1,V1))/2;
            need = sqrt(max((Vb-vo).^2 - (Va-vo).^2, 0))/Zs;
            if Vb > Va, ok = ik >= need - 1e-9; else, ok = ik <= -need + 1e-9; end
        end
        score = score + ok;
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

function m = mode_of(D1,D2,D3)
    t = 1e-6;  m = 7;
    if     D1>=D2-t && D3>=-t && D3<=D1-D2+t,                    m=1;
    elseif D2>=D1-t && D3>=1+D1-D2-t && D3<=1+t,                 m=2;
    elseif D2<=1-D1+t && D3>=D1-t && D3<=1-D2+t,                 m=3;
    elseif D2>=1-D3-t && D2<=1-D3+D1+t && D3>=D1-t && D3<=1+t,   m=4;
    elseif D2>=D1-D3-t && D2<=1-D3+t && D3>=-t && D3<=D1+t,      m=5;
    elseif D1>=1-D2-t && D3>=1-D2-t && D3<=D1+t,                 m=6;
    end
end
