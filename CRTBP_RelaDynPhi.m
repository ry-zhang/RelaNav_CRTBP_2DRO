function ydot=CRTBP_RelaDynPhi(t,y,mu)
%% 主星/标称轨道在月球会合坐标系M下的动力学（包含状态转移矩阵）
% 原点位于月球旋转坐标系M，x轴从月球指向地球，y轴地球绕月球旋转的方向
r_C = y(1:3); % 主星/标称轨道在月球旋转坐标系M下的位置速度
v_C = y(4:6);
r_C_dot = v_C;
I = eye(3);
% 圆型限制性三体下
 % 到地球及月球距离取决于会合坐标系的原点是放在质心还是放在月球等
r_E_C = r_C-[-1;0;0]; r_E_C_norm  = norm(r_E_C); r_E_C3 = r_E_C_norm^3;
r_M_C_norm  = norm(r_C); r_M_C3 = r_M_C_norm^3;

Ar_C = [1-mu/r_M_C3-(1-mu)/r_E_C3, 0, 0;
        0, 1-mu/r_M_C3-(1-mu)/r_E_C3, 0;
        0, 0, -mu/r_M_C3-(1-mu)/r_E_C3];
Av_C = [0, 2, 0;
        -2, 0, 0;
        0, 0, 0];
B = [(1-mu)*(1-1/r_E_C3); 0; 0];
v_C_dot = Ar_C*r_C + Av_C*v_C + B;

M1 = 1/r_M_C_norm^5 * (-3*r_C*r_C'+I*r_M_C_norm^2);
M2 = 1/r_E_C_norm^5 * (-3*r_E_C*r_E_C'+I*r_E_C_norm^2);

omega_mi = [0;0;1];
v_C_ddot = - 2*cross(omega_mi,v_C_dot) - cross(omega_mi,cross(omega_mi,v_C)) ...
           - ( mu*M1 + (1-mu)*M2)*v_C;


%% 副星在主星LVLH坐标系下的线性相对运动动力学
r_REL = y(7:9); % 副星在标称轨道的LVLH坐标系下的位置与速度
v_REL = y(10:12);

% LVLH的单位向量在M坐标系中的表示
h_C = cross(r_C,v_C); h_C_norm = norm(h_C);
k_LinC = h_C/h_C_norm;
i_LinC = r_C/r_M_C_norm;
j_LinC = cross(k_LinC,i_LinC);
M_C2L = [i_LinC,j_LinC,k_LinC]'; % 从会合坐标系至LVLH的转换矩阵

% 位置的模 与 角动量的模 的导数
r_norm_dot = dot(v_C, i_LinC);
h_C_dot = cross(r_C,v_C_dot);
h_norm_dot = dot(h_C_dot, k_LinC);

omega_lm_z = h_C_norm/r_M_C_norm^2;
omega_lm_x = r_M_C_norm/h_C_norm^2*dot(h_C,v_C_dot);
omega_lm = [omega_lm_x; 0; omega_lm_z];
omega_lm_z_dot = 1/r_M_C_norm*(h_norm_dot/r_M_C_norm - 2*r_norm_dot*omega_lm_z);
omega_lm_x_dot = (r_norm_dot/r_M_C_norm - 2*h_norm_dot/h_C_norm)*omega_lm_x + r_M_C_norm/h_C_norm^2*dot(h_C,v_C_ddot);
omega_lm_dot = [omega_lm_x_dot; 0; omega_lm_z_dot];

omega_li = omega_lm + M_C2L*omega_mi;
omega_li_dot = omega_lm_dot - cross(omega_lm,M_C2L*omega_mi);
Omega_li = [0, -omega_li(3), omega_li(2);
            omega_li(3), 0, -omega_li(1);
            -omega_li(2), omega_li(1), 0];
Omega_li_dot = [0, -omega_li_dot(3), omega_li_dot(2);
                omega_li_dot(3), 0, -omega_li_dot(1);
                -omega_li_dot(2), omega_li_dot(1), 0];

r_LVLH = M_C2L*r_C;
r_E_LVLH = M_C2L*r_E_C;

r_REL_dot = v_REL;
Arho = - (Omega_li_dot +  Omega_li^2 ...
     + mu/r_M_C3*(I-3*r_LVLH*r_LVLH'/r_M_C_norm^2) ...
     + (1-mu)/r_E_C3*(I-3*r_E_LVLH*r_E_LVLH'/r_E_C_norm^2));
v_REL_dot = -2*Omega_li*v_REL ...
     +Arho*r_REL;

A = [zeros(3),eye(3); Arho, -2*Omega_li];

%% phi
phi = reshape(y(13:end),6,6);
phidot = A*phi;

ydot = [r_C_dot; v_C_dot; r_REL_dot; v_REL_dot; phidot(:)];


