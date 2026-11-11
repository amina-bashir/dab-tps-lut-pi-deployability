clear; clc;
run('dab_params.m');

BAND = 2;

Rds = 1e-6;  Vth = 0.5;  G_mosfet = 1e-9;
Tsolver = 0.5e-6;
tdead = 0;

L = load('pointwise_fullTPS_NG70.mat','mod_cmd','mod_V2','mod_D1','mod_D2','mod_D3');
mod_cmd = double(L.mod_cmd(:));
mod_V2  = double(L.mod_V2(:));
mod_D1  = double(L.mod_D1);
mod_D2  = double(L.mod_D2);
mod_D3  = double(L.mod_D3);
i2_cmd_min = min(mod_cmd);
i2_cmd_max = max(mod_cmd);

C = load('controller_ramp.mat','RES','RAMP_RATE','U_LIM');
R = C.RES(BAND);
kp_cmd   = R.kp;
ki_cmd   = R.ki;
umax_cmd = C.U_LIM;

V_cap_init = Vi;

fprintf('PI of band : %g - %g V\n', R.band(1), R.band(2));
fprintf('PI         : fc %g Hz | Kp %.9g | Ki %.9g | clamp +/- %g A\n', R.fc, kp_cmd, ki_cmd, umax_cmd);
fprintf('LUT        : %d x %d | I %.2f .. %.2f A\n', numel(mod_cmd), numel(mod_V2), i2_cmd_min, i2_cmd_max);