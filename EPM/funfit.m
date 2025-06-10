function [Funfit]=funfit(x0_DRO,x0_REL)
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

Funfit=struct();
Funfit.XA.fitresult_xA=fitresult_xA;
Funfit.XA.fitresult_yA=fitresult_yA;
Funfit.XA.fitresult_zA=fitresult_zA;
Funfit.XA.fitresult_vxA=fitresult_vxA;
Funfit.XA.fitresult_vyA=fitresult_vyA;
Funfit.XA.fitresult_vzA=fitresult_vzA;

Funfit.X1.fitresult_x1 = fitresult_x1;
Funfit.X1.fitresult_y1 = fitresult_y1;
Funfit.X1.fitresult_z1 = fitresult_z1;
Funfit.X1.fitresult_vx1 = fitresult_vx1;
Funfit.X1.fitresult_vy1 = fitresult_vy1;
Funfit.X1.fitresult_vz1 = fitresult_vz1;

Funfit.X2.fitresult_x2 = fitresult_x2;
Funfit.X2.fitresult_y2 = fitresult_y2;
Funfit.X2.fitresult_z2 = fitresult_z2;
Funfit.X2.fitresult_vx2 = fitresult_vx2;
Funfit.X2.fitresult_vy2 = fitresult_vy2;
Funfit.X2.fitresult_vz2 = fitresult_vz2;


Funfit.X3.fitresult_x3 = fitresult_x3;
Funfit.X3.fitresult_y3 = fitresult_y3;
Funfit.X3.fitresult_z3 = fitresult_z3;
Funfit.X3.fitresult_vx3 = fitresult_vx3;
Funfit.X3.fitresult_vy3 = fitresult_vy3;
Funfit.X3.fitresult_vz3 = fitresult_vz3;

Funfit.X4.fitresult_x4 = fitresult_x4;
Funfit.X4.fitresult_y4 = fitresult_y4;
Funfit.X4.fitresult_z4 = fitresult_z4;
Funfit.X4.fitresult_vx4 = fitresult_vx4;
Funfit.X4.fitresult_vy4 = fitresult_vy4;
Funfit.X4.fitresult_vz4 = fitresult_vz4;

Funfit.X5.fitresult_x5 = fitresult_x5;
Funfit.X5.fitresult_y5 = fitresult_y5;
Funfit.X5.fitresult_z5 = fitresult_z5;
Funfit.X5.fitresult_vx5 = fitresult_vx5;
Funfit.X5.fitresult_vy5 = fitresult_vy5;
Funfit.X5.fitresult_vz5 = fitresult_vz5;

Funfit.X6.fitresult_x6 = fitresult_x6;
Funfit.X6.fitresult_y6 = fitresult_y6;
Funfit.X6.fitresult_z6 = fitresult_z6;
Funfit.X6.fitresult_vx6 = fitresult_vx6;
Funfit.X6.fitresult_vy6 = fitresult_vy6;
Funfit.X6.fitresult_vz6 = fitresult_vz6;










