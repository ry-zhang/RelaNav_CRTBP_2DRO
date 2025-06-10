close all;

load('x_rel_L_force.mat','x_rel_L_force');

figure()
estState=x_rel_L_force;
% estState=structState.estState;
% trueState=structState.trueState;
% YAState=structState.yA;
% estState=trueABState;

plot3(estState(:,2)*LU,estState(:,3)*LU,estState(:,4)*LU,'b','LineWidth', 2);
hold on;
% plot3(trueState(:,1)*LU,trueState(:,2)*LU,trueState(:,3)*LU,'r','LineWidth', 2);
% hold on;

% plot3(YAState(:,1)*LU,YAState(:,2)*LU,YAState(:,3)*LU,'g');
% hold on;
plot3(posMoon(1)*LU, posMoon(2)*LU,posMoon(3)*LU,'o', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
% plot3(trueState(1,1)*LU,trueState(1,2)*LU,trueState(1,3)*LU,'g^','LineWidth', 2);
% hold on;
% plot3(trueState(end,1)*LU,trueState(end,2)*LU,trueState(end,3)*LU,'r*','LineWidth', 2);
% hold on;
legend('估计相对轨迹','真实相对轨迹','月球位置','开始位置','结束位置','Location', 'northeast');
hold on;

% 显示网格
grid on;

xlabel('x[m]', 'fontsize',14,'interpreter','latex')
ylabel('y[m]', 'fontsize',14,'interpreter','latex')
zlabel('z[m]', 'fontsize',14,'interpreter','latex')
title('月心LVLH下相对轨迹 ', 'fontsize',16,'interpreter','latex')
hold on;