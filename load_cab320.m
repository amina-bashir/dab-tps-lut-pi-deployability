function ds = load_cab320(folder)
if nargin < 1, folder = fileparts(mfilename('fullpath')); end
rd = @(f) readtable(fullfile(folder, [f '.csv']), 'TextType', 'string');

T = readtable(fullfile(folder, 'coss_merged_Qoss_Eoss.csv'));
ds.Coss = @(v) interp1(T.VDS_V, T.Coss_nF*1e-9, v, 'linear', 'extrap');
ds.Qoss = @(v) interp1(T.VDS_V, T.Qoss_uC*1e-6, v, 'linear', 'extrap');
ds.Eoss = @(v) interp1(T.VDS_V, T.Eoss_uJ*1e-6, v, 'linear', 'extrap');

T = rd('fig02_Rdson_norm_vs_ID');
temps = [-40 25 100 125 150 175];
ig = linspace(0, 640, 129);
R = zeros(numel(temps), numel(ig));
for k = 1:numel(temps)
    s = T(T.curve == sprintf('%dC', temps(k)), :);
    [x, iu] = unique(s.x);
    R(k,:) = interp1(x, s.y(iu), ig, 'linear', 'extrap');
end
F = griddedInterpolant({temps, ig}, R * 3.5e-3, 'linear', 'linear');
ds.Rds = @(i, tj) F(tj, abs(i));
ds.Rpkg_hi = 0.72e-3;
ds.Rpkg_lo = 0.63e-3;

for vset = [900 1200]
    T = rd(sprintf('fig%d_Esw_vs_I_%dV', 11 + (vset==1200), vset));
    for nm = ["Eon", "Eoff", "Err"]
        s = T(T.curve == nm, :);
        [x, iu] = unique(s.x);
        ds.(sprintf('%s%d', nm, vset)) = @(i) interp1(x, s.y(iu)*1e-3, abs(i), 'linear', 'extrap');
    end
end
end
