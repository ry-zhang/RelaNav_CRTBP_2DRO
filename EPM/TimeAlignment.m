%CRTBP时间对齐
function[r_MCRLVLH_rel_ref]=TimeAlignment(x0_DRO,x0_REL,x_MCR_target)
mu=0.012150584269940;
T0=0.5*2*pi;
dt =1*T0;
length_t = 2000;
t_sample = linspace(0,dt,length_t);
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);


sol=ode113(@(t,x)CRTBP_RelaDyn(t,x,mu),[0,dt],[x0_DRO;x0_REL],opts);
sol_sample = deval(sol,t_sample)';
abs_motion_M= sol_sample(:,1:6);
rr_abs_motion_M= sol_sample(:,1:3);
vv_abs_motion_M= sol_sample(:,4:6);
% rel_motion_L= sol_sample(7:12,:);
rr_MCRLVLH_rel_CRTBP=sol_sample(:,7:9);
vv_MCRLVLH_rel_CRTBP=sol_sample(:,10:12);

%% 拟合CRTBP下的参考周期轨道
% ----------------------------拟合自变量与因变量--------------------------
theta_all = mod(atan2(-abs_motion_M(:,2),abs_motion_M(:,1)),2*pi);%计算主星在 MCR 坐标系中的几何相位。
% -----------------------------------拟合（最小二乘）---------------------------------
ft = fittype( 'fourier8' );%使用 8 项的 Fourier 系列进行拟合。
opts_fit = fitoptions( 'Method', 'NonlinearLeastSquares','Display','Off' );
[fitresult_x, gof_x] = fit( theta_all, rr_MCRLVLH_rel_CRTBP(:,1), ft, opts_fit );
[fitresult_y, gof_y] = fit( theta_all, rr_MCRLVLH_rel_CRTBP(:,2), ft, opts_fit );
[fitresult_vx, gof_vx] = fit( theta_all, vv_MCRLVLH_rel_CRTBP(:,1), ft, opts_fit );
[fitresult_vy, gof_vy] = fit( theta_all, vv_MCRLVLH_rel_CRTBP(:,2), ft, opts_fit );

%----------------- ------------将CRTBP轨道转换成以theta为自变量的相对轨道--------------
theta_MCR_targe_eph = mod(atan2(-x_MCR_target(2),x_MCR_target(1)),2*pi);%主星在 MCR 坐标系中的方向角/相位。
fit_x =feval(fitresult_x,theta_MCR_targe_eph);%基于主星的轨道递推结果，计算参考相对轨道
fit_y =feval(fitresult_y,theta_MCR_targe_eph);
fit_vx =feval(fitresult_vx,theta_MCR_targe_eph);%基于主星的轨道递推结果，计算参考相对轨道
fit_vy =feval(fitresult_vy,theta_MCR_targe_eph);
r_MCRLVLH_rel_ref = [fit_x,fit_y,zeros(size(fit_x)),fit_vx,fit_vy,zeros(size(fit_vx))];%计算CRTBP下相对参考轨道（转换成eph相位）。


    
