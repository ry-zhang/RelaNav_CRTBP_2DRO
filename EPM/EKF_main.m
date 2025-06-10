% LiAISON--CRTBP--自主定轨
% 说明：
% 1.不考虑时间相关项；
% 2.考虑遮挡；
% 3.定轨过程是在以地月质心为原点的地月旋转系下完成的。

clc;clear;close all;

warning('off')

%% Simulation parameters
% 加载常数
const

% 初始误差协方差矩阵, sigma, meter, meter/s (diagonal elements)
% P0 = [1000 1000 1000 1 1 1 1000 1000 1000 1 1 1];
%  P0 = [1 1 1 0.001 0.001 0.001];

% P0 = [0 0 0 0 0 0];
% P0 = [1000 1000 1000 1 1 1];
% P0 = [1000 1000 1000 0.1 0.1 0.1];
% P0 = P0.^2;
% C0=[10^-5 10^-5 10^-5 10^-5 10^-5 10^-5]*0.2;
C0=[10^-5 10^-5 10^-5 10^-5 10^-5 10^-5]*0.2;
P0 = C0.^2;



% 初始状态偏差
Yerr = sqrt(P0);
% Yerr = [1000 1000 1000 1 1 1 1 1 1 0.001 0.001 0.001];
% Yerr = zeros(1,12);

%1:EML2south, 
%2:EML2south,
%3:EML2north, 
%4:EML2north 
%5:EML1south 
%6:EML1south,
%7:EML1north 
%8:EML1north
%9:Lunar Elliptic
%10:Lunar Polar-Circular
%11:NRHO EML2south,
%12:NRHO EML1south 
%13:NRHO EML2north, 
%14:NRHO EML1north 
%15:EML2 Lyapunov 
%16:EML1 Lyapunov
%17:DRO 2:1->1.3670569939965526E+1天
%18:DRO 3:1->9.1105936571618997E+0天
%19:DRO 4:1->6.8342565493518084E+0天
%20:DRO 5:1->5.4773934847797472E+0天
%21:LEO

% S/C-1
SC1 =22;%DRO

% S/C-2
SC2 = 21;%LEO

% 测量类型 (make it 1 in case )
rm = 2;          % 测距1，测距测角

% 仿真时长
simdur =13.64;    % [days]

% 仿真步长
simstep = 60;   % [s] 3min

% 测量误差的标准差
sigma_range = 1;            % 测距误差 SD, m

% 测量误差的期望
bias_range=0;

% 主函数

% [structRss, structState, distance] = CRTBP_OD_EKF_func(P0, Yerr, SC1, SC2, rm, simdur, simstep,sigma_range, bias_range);%CRTBP相对状态估计
% [structRss, structState, distance] = CRTBP_CD_EKF_func(P0, Yerr, SC1, SC2, rm, simdur, simstep, sigma_range, bias_range);%CRTBP模态估计
[structRss, structState, distance] =EPM_CD_EKF_func(P0, Yerr, SC1, SC2, rm, simdur, simstep, sigma_range, bias_range);%星历下模态估计
% [structRss, structState, distance] =CRTBP_OD_UKF_func(P0, Yerr, SC1, SC2, rm, simdur, simstep, sigma_range, bias_range);


% 收敛点
endpoint = 200; % 需要根据实际情况修改
disp('Simulation is done.')
disp('-------------------------------------------------------------------')
disp(['SC 1st mean RMS Position error              : ' num2str(mean(structRss.rsspos(end-endpoint,1))*LU) ' m'])
disp(['SC 1st mean RMS Position uncertainty (1sig) : ' num2str(mean(structRss.rssposunc(end-endpoint,1))*LU) ' m'])
disp(['SC 1st mean RMS Velocity error              : ' num2str(mean(structRss.rssvel(end-endpoint,1))*SU) ' m/s'])
disp(['SC 1st mean RMS Velocity uncertainty (1sig) : ' num2str(mean(structRss.rssvelunc(end-endpoint,1))*SU) ' m/s'])

% disp('---------------------------------------------------------')
% disp(['SC 2nd mean RMS Position error              : ' num2str(mean(structRss.rsspos(end-endpoint,2))*LU) ' m'])
% disp(['SC 2nd mean RMS Position uncertainty (1sig) : ' num2str(mean(structRss.rssposunc(end-endpoint,2))*LU) ' m'])
% disp(['SC 2nd mean RMS Velocity error              : ' num2str(mean(structRss.rssvel(end-endpoint,2))*SU) ' m/s'])
% disp(['SC 2nd mean RMS Velocity uncertainty (1sig) : ' num2str(mean(structRss.rssvelunc(end-endpoint,2))*SU) ' m/s'])
% disp('-------------------------------------------------------------------')

%% Plot
% if full state results are required = 1
fullerr = 1; % 绘制误差图
fullunc = 1; % 绘制不确定图

% RSS uncertainty plots
Plot_RSS_EKF_detailed

% figure()
% plot(distance*LU)

% figure()
% plot3(structState.trueState(:,1),structState.trueState(:,2),structState.trueState(:,3))
% hold on
% plot3(structState.estState(:,1+6),structState.estState(:,2+6),structState.estState(:,3+6))
