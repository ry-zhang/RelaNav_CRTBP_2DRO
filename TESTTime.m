% clc;
% clear;
close all;
const
step=60/TU;
n_obs=27.284429211881122*86400/TU/step;
options=odeset('RelTol',1e-16,'AbsTol',1e-22);

IC=[ 0.808936204314439-(1-mu),	0,	0,0, 	0.515632164677207,0]; % 2:1平面DRO
Y0_rel_3=[0,0.576120458608156,0,0.817364800546933,0,0].*10^-05;%相对状态初值
% N=20000/ Y0_rel_3(2)/LU;
Y0_true=Y0_rel_3';
YL_true = [IC(1,:)]';%dro主星初值 
Y=Y0_true;
estState=[];
trueState=[];
k=0;
t=0;
for i=1:n_obs

    k=k+1;
    % Previous step
    t_old = t;

    % Propagation to measurement epoch
    t = i*step;     % Time since epoch
    tspan = linspace(t_old, t, 2);


    % 主星轨道和相对轨道
   sol=ode113(@(t,x)CRTBP_RelaDyn(t,x,mu),tspan, [YL_true;Y0_true],options);
   sol_t=deval(sol,tspan);
    Y=sol_t(7:12,end); %相对轨道(lvlh)
    YL_true=sol_t(1:6,end);
    % 2nd S/C
    trueABState(k,:) = Y;
    trueAState(k,:) = YL_true;

end
figure()

plot3(trueABState(:,1)*LU,trueABState(:,2)*LU,trueABState(:,3)*LU,'b','LineWidth', 2);
hold on;
% % plot3(trueState(:,1)*LU,trueState(:,2)*LU,trueState(:,3)*LU,'r','LineWidth', 2);
% % hold on;

% plot3(YAState(:,1)*LU,YAState(:,2)*LU,YAState(:,3)*LU,'g');
% hold on;
plot3(posMoon(1)*LU, posMoon(2)*LU,posMoon(3)*LU,'o', 'MarkerSize', 8, 'MarkerFaceColor', 'y');hold on;
plot3(trueABState(1,1)*LU,trueABState(1,2)*LU,trueABState(1,3)*LU,'g^','LineWidth', 2);
hold on;
plot3(trueABState(end,1)*LU,trueABState(end,2)*LU,trueABState(end,3)*LU,'r*','LineWidth', 2);
hold on;
legend('估计相对轨迹','真实相对轨迹','月球位置','开始位置','结束位置','Location', 'northeast');
hold on;


