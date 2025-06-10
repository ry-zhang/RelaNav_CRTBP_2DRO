function[Y1,Y,Sia]=LT(C,YL_true,x0_REL,tspan,tf)
load('Meigva_diag.mat','Meigva_diag');
const
mu= 0.012150584269940; 
options=odeset('RelTol',1e-13,'AbsTol',1e-20);

% 平面拟周期
alpha1 = atan2(imag(Meigva_diag(2)),real(Meigva_diag(2)));
T1 = 2*pi*T0/alpha1;
sol1 = ode113(@(t,x)CRTBP_RelaDyn(t,x,mu),tspan, [YL_true;(x0_REL(1,:))'], options);
sol1_sample = deval(sol1,tspan);
Y1=sol1_sample(1:6,end);
rel_motion_L_linear_1(:,1) = sol1_sample(7:12,end);

sol2= ode113(@(t,x)CRTBP_RelaDyn(t,x,mu),tspan, [YL_true;(x0_REL(2,:))'], options);
sol2_sample = deval(sol2,tspan);
rel_motion_L_linear_2(:,1) = sol2_sample(7:12,end);

% e_1 e_2
e1_hat = cos(-alpha1*tf/T0).*rel_motion_L_linear_1 + sin(-alpha1*tf/T0).*rel_motion_L_linear_2;
e2_hat = -sin(-alpha1*tf/T0).*rel_motion_L_linear_1 + cos(-alpha1*tf/T0).*rel_motion_L_linear_2;

X1=cos(alpha1*tf/T0).*e1_hat + sin(alpha1*tf/T0).*e2_hat;
X2=-sin(alpha1*tf/T0).*e1_hat + cos(alpha1*tf/T0).*e2_hat;

% 周期轨道
sol3 = ode113(@(t,x)CRTBP_RelaDyn(t,x,mu),tspan, [YL_true;(x0_REL(3,:))'], options);
sol3_sample = deval(sol3,tspan);
e3_hat(:,1) = sol3_sample(7:12,end);
X3=e3_hat;


% 发散轨道
sol4 = ode113(@(t,x)CRTBP_RelaDyn(t,x,mu),tspan, [YL_true;(x0_REL(4,:))'], options);
sol4_sample = deval(sol4,tspan);
Y1=sol4_sample(1:6,end);
e4_hat(:,1) = sol4_sample(7:12,end);
X4=e4_hat+(tf/T0).*e3_hat;

%法向拟周期
alpha2 = atan2(imag(Meigva_diag(6)),real(Meigva_diag(6)));
T2 = 2*pi*T0/alpha2;
sol5 = ode113(@(t,x)CRTBP_RelaDyn(t,x,mu),tspan, [YL_true;(x0_REL(5,:))'], options);
sol5_sample = deval(sol5,tspan);
rel_motion_L_linear_5(:,1) = sol5_sample(7:12,end);

sol6= ode113(@(t,x)CRTBP_RelaDyn(t,x,mu),tspan, [YL_true;(x0_REL(6,:))'], options);
sol6_sample = deval(sol6,tspan);
rel_motion_L_linear_6(:,1) = sol6_sample(7:12,end);

% e_1 e_2
e5_hat = cos(-alpha2*tf/T0).*rel_motion_L_linear_5 + sin(-alpha2*tf/T0).*rel_motion_L_linear_6;
e6_hat = -sin(-alpha2*tf/T0).*rel_motion_L_linear_5 + cos(-alpha2*tf/T0).*rel_motion_L_linear_6;

X5=cos(alpha2*tf/T0).*e5_hat + sin(alpha2*tf/T0).*e6_hat;
X6=-sin(alpha2*tf/T0).*e5_hat + cos(alpha2*tf/T0).*e6_hat;

Y=C(1)*X1+C(2)*X2 +C(3)*X3 +C(4)*X4 +C(5)*X5+C(6)*X6;
Sia=[X1,X2,X3,X4,X5,X6];


