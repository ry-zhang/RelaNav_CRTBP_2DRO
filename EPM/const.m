%%% 一些常用的常数

% 会合坐标系的单位
LU         = 389703e3;   %m, lenght unit
% LU=38439e3;
%TU         = 3.829808977771119e+05; %s, time unit
TU_day=27.284429211881122;
TU=TU_day*86400;
% T0_day=6.276897008112217*TU_day;
% T0=6.276897008112217*TU;
T0_day=0.5*TU_day;
T0=0.5*TU;
SU         = LU/TU; %m/s, velocity convertion

R_Earth    = 6378137.0/LU;
R_Moon     = 1738000.0/LU;
mu         = 0.012150584269940; 
% posEarth   = [-mu;0;0];
% posMoon    = [1-mu;0;0];
posEarth   = [-1;0;0];
posMoon    = [0;0;0];