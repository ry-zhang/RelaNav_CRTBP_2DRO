%CRTBP时间对齐
function[XA,r_MCRLVLH_rel_ref,Sia]=TimeAlignmentSiaFit(C,x_MCR_target)
load('Funfit.mat','Funfit');


%% 拟合CRTBP下的参考周期轨道


fitresult_xA = Funfit.XA.fitresult_xA;
fitresult_yA = Funfit.XA.fitresult_yA;
fitresult_zA = Funfit.XA.fitresult_zA;
fitresult_vxA = Funfit.XA.fitresult_vxA;
fitresult_vyA = Funfit.XA.fitresult_vyA;
fitresult_vzA = Funfit.XA.fitresult_vzA;

fitresult_x1 = Funfit.X1.fitresult_x1;
fitresult_y1 = Funfit.X1.fitresult_y1;
fitresult_z1 = Funfit.X1.fitresult_z1;
fitresult_vx1 = Funfit.X1.fitresult_vx1;
fitresult_vy1 = Funfit.X1.fitresult_vy1;
fitresult_vz1 = Funfit.X1.fitresult_vz1;

fitresult_x2 = Funfit.X2.fitresult_x2;
fitresult_y2 = Funfit.X2.fitresult_y2;
fitresult_z2 = Funfit.X2.fitresult_z2;
fitresult_vx2 = Funfit.X2.fitresult_vx2;
fitresult_vy2 = Funfit.X2.fitresult_vy2;
fitresult_vz2 = Funfit.X2.fitresult_vz2;

fitresult_x3 = Funfit.X3.fitresult_x3;
fitresult_y3 = Funfit.X3.fitresult_y3;
fitresult_z3 = Funfit.X3.fitresult_z3;
fitresult_vx3 = Funfit.X3.fitresult_vx3;
fitresult_vy3 = Funfit.X3.fitresult_vy3;
fitresult_vz3 = Funfit.X3.fitresult_vz3;

fitresult_x4 = Funfit.X4.fitresult_x4;
fitresult_y4 = Funfit.X4.fitresult_y4;
fitresult_z4 = Funfit.X4.fitresult_z4;
fitresult_vx4 = Funfit.X4.fitresult_vx4;
fitresult_vy4 = Funfit.X4.fitresult_vy4;
fitresult_vz4 = Funfit.X4.fitresult_vz4;

fitresult_x5 = Funfit.X5.fitresult_x5;
fitresult_y5 = Funfit.X5.fitresult_y5;
fitresult_z5 = Funfit.X5.fitresult_z5;
fitresult_vx5 = Funfit.X5.fitresult_vx5;
fitresult_vy5 = Funfit.X5.fitresult_vy5;
fitresult_vz5 = Funfit.X5.fitresult_vz5;

fitresult_x6 = Funfit.X6.fitresult_x6;
fitresult_y6 = Funfit.X6.fitresult_y6;
fitresult_z6 = Funfit.X6.fitresult_z6;
fitresult_vx6 = Funfit.X6.fitresult_vx6;
fitresult_vy6 = Funfit.X6.fitresult_vy6;
fitresult_vz6 = Funfit.X6.fitresult_vz6;



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
fit_z6 =feval(fitresult_z6,theta_MCR_targe_eph);
fit_vx6 =feval(fitresult_vx6,theta_MCR_targe_eph);%基于主星的轨道递推结果，计算参考相对轨道
fit_vy6 =feval(fitresult_vy6,theta_MCR_targe_eph);
fit_vz6 =feval(fitresult_vz6,theta_MCR_targe_eph);
X6=[fit_x6,fit_y6,fit_z6,fit_vx6,fit_vy6,fit_vz6]'; 

r_MCRLVLH_rel_ref=C(1)*X1+C(2)*X2 +C(3)*X3 +C(4)*X4 +C(5)*X5+C(6)*X6;%计算CRTBP下相对参考轨道（转换成eph相位）。
Sia=[X1,X2,X3,X4,X5,X6];


    
