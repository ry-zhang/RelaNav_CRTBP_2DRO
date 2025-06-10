%CRTBP时间对齐
function[XA,r_MCRLVLH_rel_ref,Sia]=TimeAlignmentSia(C,x0_DRO,x0_REL,x_MCR_target)
load('Meigva_diag.mat','Meigva_diag');
mu=0.012150584269940;
T0=0.5*2*pi;
dt =1*T0;
length_t = 2000;
t_sample = linspace(0,dt,length_t);
options= odeset('RelTol',1e-10,'AbsTol',1e-10);


% 平面拟周期
alpha1 = atan2(imag(Meigva_diag(2)),real(Meigva_diag(2)));
T1 = 2*pi*T0/alpha1;
sol1 = ode113(@(t,x)CRTBP_RelaDyn(t,x,mu),t_sample, [x0_DRO;(x0_REL(1,:))'], options);
sol1_sample = deval(sol1,t_sample)';
rel_motion_L_linear_1(:,1:6) = sol1_sample(:,7:12);

sol2= ode113(@(t,x)CRTBP_RelaDyn(t,x,mu),t_sample, [x0_DRO;(x0_REL(2,:))'], options);
sol2_sample = deval(sol2,t_sample)';
rel_motion_L_linear_2(:,1:6) = sol2_sample(:,7:12);

% e_1 e_2
e1_hat = cos(-alpha1*t_sample'/T0).*rel_motion_L_linear_1 + sin(-alpha1*t_sample'/T0).*rel_motion_L_linear_2;
e2_hat = -sin(-alpha1*t_sample'/T0).*rel_motion_L_linear_1 + cos(-alpha1*t_sample'/T0).*rel_motion_L_linear_2;

X1=cos(alpha1*t_sample'/T0).*e1_hat + sin(alpha1*t_sample'/T0).*e2_hat;
X2=-sin(alpha1*t_sample'/T0).*e1_hat + cos(alpha1*t_sample'/T0).*e2_hat;

% 周期轨道
sol3 = ode113(@(t,x)CRTBP_RelaDyn(t,x,mu),t_sample, [x0_DRO;(x0_REL(3,:))'], options);
sol3_sample = deval(sol3,t_sample)';
e3_hat(:,1:6) = sol3_sample(:,7:12);
X3=e3_hat;

% 发散轨道
sol4 = ode113(@(t,x)CRTBP_RelaDyn(t,x,mu),t_sample, [x0_DRO;(x0_REL(4,:))'], options);
sol4_sample = deval(sol4,t_sample)';
Y1=sol4_sample(1:6,end);
e4_hat(:,1:6) = sol4_sample(:,7:12);
X4=e4_hat+(t_sample'/T0).*e3_hat;

%法向拟周期
alpha2 = atan2(imag(Meigva_diag(6)),real(Meigva_diag(6)));
T2 = 2*pi*T0/alpha2;
sol5 = ode113(@(t,x)CRTBP_RelaDyn(t,x,mu),t_sample, [x0_DRO;(x0_REL(5,:))'], options);
sol5_sample = deval(sol5,t_sample)';
rel_motion_L_linear_5(:,1:6) = sol5_sample(:,7:12);

sol6= ode113(@(t,x)CRTBP_RelaDyn(t,x,mu),t_sample, [x0_DRO;(x0_REL(6,:))'], options);
sol6_sample = deval(sol6,t_sample)';
abs_motion_M= sol6_sample(:,1:6);
rel_motion_L_linear_6(:,1:6) = sol6_sample(:,7:12);

% e_1 e_2
e5_hat = cos(-alpha2*t_sample'/T0).*rel_motion_L_linear_5 + sin(-alpha2*t_sample'/T0).*rel_motion_L_linear_6;
e6_hat = -sin(-alpha2*t_sample'/T0).*rel_motion_L_linear_5 + cos(-alpha2*t_sample'/T0).*rel_motion_L_linear_6;

X5=cos(alpha2*t_sample'/T0).*e5_hat + sin(alpha2*t_sample'/T0).*e6_hat;
X6=-sin(alpha2*t_sample'/T0).*e5_hat + cos(alpha2*t_sample'/T0).*e6_hat;

rr_A_MCR_CRTBP=abs_motion_M(:,1:3);
vv_A_MCR_CRTBP=abs_motion_M(:,4:6);


X1_rr_MCRLVLH_rel_CRTBP=X1(:,1:3);
X1_vv_MCRLVLH_rel_CRTBP=X1(:,4:6);

X2_rr_MCRLVLH_rel_CRTBP=X2(:,1:3);
X2_vv_MCRLVLH_rel_CRTBP=X2(:,4:6);

X3_rr_MCRLVLH_rel_CRTBP=X3(:,1:3);
X3_vv_MCRLVLH_rel_CRTBP=X3(:,4:6);

X4_rr_MCRLVLH_rel_CRTBP=X4(:,1:3);
X4_vv_MCRLVLH_rel_CRTBP=X4(:,4:6);

X5_rr_MCRLVLH_rel_CRTBP=X5(:,1:3);
X5_vv_MCRLVLH_rel_CRTBP=X5(:,4:6);

X6_rr_MCRLVLH_rel_CRTBP=X6(:,1:3);
X6_vv_MCRLVLH_rel_CRTBP=X6(:,4:6);



%% 拟合CRTBP下的参考周期轨道
% ----------------------------拟合自变量与因变量--------------------------
theta_all = mod(atan2(-abs_motion_M(:,2),abs_motion_M(:,1)),2*pi);%计算主星在 MCR 坐标系中的几何相位。
% -----------------------------------拟合（最小二乘）---------------------------------
ft = fittype( 'fourier8' );%使用 8 项的 Fourier 系列进行拟合。
opts_fit = fitoptions( 'Method', 'NonlinearLeastSquares','Display','Off' );
[fitresult_xA, gof_xA] = fit( theta_all, rr_A_MCR_CRTBP(:,1), ft, opts_fit );
[fitresult_yA, gof_xA] = fit( theta_all, rr_A_MCR_CRTBP(:,2), ft, opts_fit );
[fitresult_zA, gof_zA] = fit( theta_all, rr_A_MCR_CRTBP(:,3), ft, opts_fit );
[fitresult_vxA, gof_vxA] = fit( theta_all, rr_A_MCR_CRTBP(:,1), ft, opts_fit );
[fitresult_vyA, gof_vyA] = fit( theta_all, rr_A_MCR_CRTBP(:,2), ft, opts_fit );
[fitresult_vzA, gof_vzA] = fit( theta_all, rr_A_MCR_CRTBP(:,3), ft, opts_fit );

[fitresult_x1, gof_x1] = fit( theta_all, X1_rr_MCRLVLH_rel_CRTBP(:,1), ft, opts_fit );
[fitresult_y1, gof_x1] = fit( theta_all, X1_rr_MCRLVLH_rel_CRTBP(:,2), ft, opts_fit );
[fitresult_z1, gof_z1] = fit( theta_all, X1_rr_MCRLVLH_rel_CRTBP(:,3), ft, opts_fit );
[fitresult_vx1, gof_vx1] = fit( theta_all, X1_vv_MCRLVLH_rel_CRTBP(:,1), ft, opts_fit );
[fitresult_vy1, gof_vy1] = fit( theta_all, X1_vv_MCRLVLH_rel_CRTBP(:,2), ft, opts_fit );
[fitresult_vz1, gof_vz1] = fit( theta_all, X1_vv_MCRLVLH_rel_CRTBP(:,3), ft, opts_fit );

[fitresult_x2, gof_x2] = fit( theta_all, X2_rr_MCRLVLH_rel_CRTBP(:,1), ft, opts_fit );
[fitresult_y2, gof_x2] = fit( theta_all, X2_rr_MCRLVLH_rel_CRTBP(:,2), ft, opts_fit );
[fitresult_z2, gof_z2] = fit( theta_all, X2_rr_MCRLVLH_rel_CRTBP(:,3), ft, opts_fit );
[fitresult_vx2, gof_vx2] = fit( theta_all, X2_vv_MCRLVLH_rel_CRTBP(:,1), ft, opts_fit );
[fitresult_vy2, gof_vy2] = fit( theta_all, X2_vv_MCRLVLH_rel_CRTBP(:,2), ft, opts_fit );
[fitresult_vz2, gof_vz2] = fit( theta_all, X2_vv_MCRLVLH_rel_CRTBP(:,3), ft, opts_fit );

[fitresult_x3, gof_x3] = fit( theta_all, X3_rr_MCRLVLH_rel_CRTBP(:,1), ft, opts_fit );
[fitresult_y3, gof_x3] = fit( theta_all, X3_rr_MCRLVLH_rel_CRTBP(:,2), ft, opts_fit );
[fitresult_z3, gof_z3] = fit( theta_all, X3_rr_MCRLVLH_rel_CRTBP(:,3), ft, opts_fit );
[fitresult_vx3, gof_vx3] = fit( theta_all, X3_vv_MCRLVLH_rel_CRTBP(:,1), ft, opts_fit );
[fitresult_vy3, gof_vy3] = fit( theta_all, X3_vv_MCRLVLH_rel_CRTBP(:,2), ft, opts_fit );
[fitresult_vz3, gof_vz3] = fit( theta_all, X3_vv_MCRLVLH_rel_CRTBP(:,3), ft, opts_fit );

[fitresult_x4, gof_x4] = fit( theta_all, X4_rr_MCRLVLH_rel_CRTBP(:,1), ft, opts_fit );
[fitresult_y4, gof_x4] = fit( theta_all, X4_rr_MCRLVLH_rel_CRTBP(:,2), ft, opts_fit );
[fitresult_z4, gof_z4] = fit( theta_all, X4_rr_MCRLVLH_rel_CRTBP(:,3), ft, opts_fit );
[fitresult_vx4, gof_vx4] = fit( theta_all, X4_vv_MCRLVLH_rel_CRTBP(:,1), ft, opts_fit );
[fitresult_vy4, gof_vy4] = fit( theta_all, X4_vv_MCRLVLH_rel_CRTBP(:,2), ft, opts_fit );
[fitresult_vz4, gof_vz4] = fit( theta_all, X4_vv_MCRLVLH_rel_CRTBP(:,3), ft, opts_fit );


[fitresult_x5, gof_x5] = fit( theta_all, X5_rr_MCRLVLH_rel_CRTBP(:,1), ft, opts_fit );
[fitresult_y5, gof_x5] = fit( theta_all, X5_rr_MCRLVLH_rel_CRTBP(:,2), ft, opts_fit );
[fitresult_z5, gof_z5] = fit( theta_all, X5_rr_MCRLVLH_rel_CRTBP(:,3), ft, opts_fit );
[fitresult_vx5, gof_vx5] = fit( theta_all, X5_vv_MCRLVLH_rel_CRTBP(:,1), ft, opts_fit );
[fitresult_vy5, gof_vy5] = fit( theta_all, X5_vv_MCRLVLH_rel_CRTBP(:,2), ft, opts_fit );
[fitresult_vz5, gof_vz5] = fit( theta_all, X5_vv_MCRLVLH_rel_CRTBP(:,3), ft, opts_fit );

[fitresult_x6, gof_x6] = fit( theta_all, X6_rr_MCRLVLH_rel_CRTBP(:,1), ft, opts_fit );
[fitresult_y6, gof_x6] = fit( theta_all, X6_rr_MCRLVLH_rel_CRTBP(:,2), ft, opts_fit );
[fitresult_z6, gof_z6] = fit( theta_all, X6_rr_MCRLVLH_rel_CRTBP(:,3), ft, opts_fit );
[fitresult_vx6, gof_vx6] = fit( theta_all, X6_vv_MCRLVLH_rel_CRTBP(:,1), ft, opts_fit );
[fitresult_vy6, gof_vy6] = fit( theta_all, X6_vv_MCRLVLH_rel_CRTBP(:,2), ft, opts_fit );
[fitresult_vz6, gof_vz6] = fit( theta_all, X6_vv_MCRLVLH_rel_CRTBP(:,3), ft, opts_fit );


%----------------- -----------将CRTBP轨道转换成以theta为自变量的相对轨道--------------
theta_MCR_targe_eph = mod(atan2(-x_MCR_target(2),x_MCR_target(1)),2*pi);%主星在 MCR 坐标系中的方向角/相位。

fit_xA=feval(fitresult_xA,theta_MCR_targe_eph);%基于主星的轨道递推结果，计算参考相对轨道
fit_yA =feval(fitresult_yA,theta_MCR_targe_eph);
fit_zA =feval(fitresult_zA,theta_MCR_targe_eph);
fit_vxA =feval(fitresult_vxA,theta_MCR_targe_eph);%基于主星的轨道递推结果，计算参考相对轨道
fit_vyA =feval(fitresult_vyA,theta_MCR_targe_eph);
fit_vzA =feval(fitresult_vzA,theta_MCR_targe_eph);
XA=[fit_xA,fit_yA,fit_zA,fit_vxA,fit_vyA,fit_vzA]';

fit_x1 =feval(fitresult_x1,theta_MCR_targe_eph);%基于主星的轨道递推结果，计算参考相对轨道
fit_y1 =feval(fitresult_y1,theta_MCR_targe_eph);
fit_z1 =feval(fitresult_z1,theta_MCR_targe_eph);
fit_vx1 =feval(fitresult_vx1,theta_MCR_targe_eph);%基于主星的轨道递推结果，计算参考相对轨道
fit_vy1 =feval(fitresult_vy1,theta_MCR_targe_eph);
fit_vz1 =feval(fitresult_vz1,theta_MCR_targe_eph);

X1=[fit_x1,fit_y1,fit_z1,fit_vx1,fit_vy1,fit_vz1]';

fit_x2 =feval(fitresult_x2,theta_MCR_targe_eph);%基于主星的轨道递推结果，计算参考相对轨道
fit_y2 =feval(fitresult_y2,theta_MCR_targe_eph);
fit_z2 =feval(fitresult_z2,theta_MCR_targe_eph);
fit_vx2 =feval(fitresult_vx2,theta_MCR_targe_eph);%基于主星的轨道递推结果，计算参考相对轨道
fit_vy2 =feval(fitresult_vy2,theta_MCR_targe_eph);
fit_vz2 =feval(fitresult_vz2,theta_MCR_targe_eph);

X2=[fit_x2,fit_y2,fit_z2,fit_vx2,fit_vy2,fit_vz2]';

fit_x3 =feval(fitresult_x3,theta_MCR_targe_eph);%基于主星的轨道递推结果，计算参考相对轨道
fit_y3 =feval(fitresult_y3,theta_MCR_targe_eph);
fit_z3 =feval(fitresult_z3,theta_MCR_targe_eph);
fit_vx3 =feval(fitresult_vx3,theta_MCR_targe_eph);%基于主星的轨道递推结果，计算参考相对轨道
fit_vy3 =feval(fitresult_vy3,theta_MCR_targe_eph);
fit_vz3 =feval(fitresult_vz3,theta_MCR_targe_eph);

X3=[fit_x3,fit_y3,fit_z3,fit_vx3,fit_vy3,fit_vz3]';

fit_x4 =feval(fitresult_x4,theta_MCR_targe_eph);%基于主星的轨道递推结果，计算参考相对轨道
fit_y4 =feval(fitresult_y4,theta_MCR_targe_eph);
fit_z4 =feval(fitresult_z4,theta_MCR_targe_eph);
fit_vx4 =feval(fitresult_vx4,theta_MCR_targe_eph);%基于主星的轨道递推结果，计算参考相对轨道
fit_vy4 =feval(fitresult_vy4,theta_MCR_targe_eph);
fit_vz4 =feval(fitresult_vz4,theta_MCR_targe_eph);
X4=[fit_x4,fit_y4,fit_z4,fit_vx4,fit_vy4,fit_vz4]';

fit_x5 =feval(fitresult_x5,theta_MCR_targe_eph);%基于主星的轨道递推结果，计算参考相对轨道
fit_y5 =feval(fitresult_y5,theta_MCR_targe_eph);
fit_z5 =feval(fitresult_z5,theta_MCR_targe_eph);
fit_vx5 =feval(fitresult_vx5,theta_MCR_targe_eph);%基于主星的轨道递推结果，计算参考相对轨道
fit_vy5 =feval(fitresult_vy5,theta_MCR_targe_eph);
fit_vz5 =feval(fitresult_vz5,theta_MCR_targe_eph);
X5=[fit_x5,fit_y5,fit_z5,fit_vx5,fit_vy5,fit_vz5]';

fit_x6 =feval(fitresult_x6,theta_MCR_targe_eph);%基于主星的轨道递推结果，计算参考相对轨道
fit_y6 =feval(fitresult_y6,theta_MCR_targe_eph);
fit_z6 =feval(fitresult_z6,theta_MCR_targe_eph);x0_DRO,x0_REL
fit_vx6 =feval(fitresult_vx6,theta_MCR_targe_eph);%基于主星的轨道递推结果，计算参考相对轨道
fit_vy6 =feval(fitresult_vy6,theta_MCR_targe_eph);
fit_vz6 =feval(fitresult_vz6,theta_MCR_targe_eph);
X6=[fit_x6,fit_y6,fit_z6,fit_vx6,fit_vy6,fit_vz6]'; 

r_MCRLVLH_rel_ref=C(1)*X1+C(2)*X2 +C(3)*X3 +C(4)*X4 +C(5)*X5+C(6)*X6;%计算CRTBP下相对参考轨道（转换成eph相位）。
Sia=[X1,X2,X3,X4,X5,X6];


    
