function ydot=sysSolveCRTBP_RelaDyn(t,y,mu)
%function ydot=sysSolveCRTBP_RelaDyn(t,y,options,flag, mu)
% function ydot=sysSolveCRTBP(t,y,options,flag, mu)


%This file contains the right hand side for the system which defines the
%state transition matrix for the spatial CRTBP.  
%状态转移矩阵通过积分获得
%-------------------------------------------------------------------------
%The structure of the system is 
%
%          A'=Df*A, x'=f(x)    
%
%Where A,Df, are 6X6 matrices and A' is the derivative of A.
%Df is the derivative of the N-body vector field f and has the form
%
%                    |  0   I  |  
%               Df = |         |
%                    | G    K  | 
%
%where all submatrices are 3X3 and G has depends only on the position
%vectors of the bodies (all be it in a complicated way), and x'=f(x)
% is the spatial CRTBP.  All of this is combined onto one 
%system of 1st order ODEs. The right hand side for the system is coded 
%in the remainder of the file
%------------------------------------------------------------------------
%% 主星/标称轨道在月球会合坐标系M下的动力学
% 原点位于月球旋转坐标系M，x轴从月球指向地球，y轴地球绕月球旋转的方向
r_C = y(37:39); % 主星/标称轨道在月球旋转坐标系M下的位置速度
v_C = y(40:42);
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
r_REL = y(43:45); % 副星在标称轨道的LVLH坐标系下的位置与速度
v_REL = y(46:48);

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
Arho1= -(Omega_li_dot +  Omega_li^2) ;
Arho2=((1-mu)*r_E_LVLH/(norm(r_E_LVLH))^3-(1-mu)*(r_E_LVLH+r_REL)/(norm(r_E_LVLH+r_REL))^3)+(mu*r_LVLH/(norm(r_LVLH))^3-mu*(r_LVLH+r_REL)/(norm(r_LVLH+r_REL))^3);

dArho2_dr=(1-mu)*3*r_E_LVLH*r_E_LVLH'/(norm(r_E_LVLH)^5)-(1-mu)/(norm(r_E_LVLH))^3+mu*(r_LVLH+r_REL)*(r_LVLH+r_REL)'/(norm((r_LVLH+r_REL)))^5-mu/(norm((r_LVLH+r_REL)))^3;

v_REL_dot = -2*Omega_li*v_REL +Arho1*r_REL+Arho2;

ydot = [r_C_dot; v_C_dot; r_REL_dot; v_REL_dot];
%--------------------------------------------------------------------------

%FIRST
%The matrix G is computed
%This matrix depends on the positions of the N-body problem.
%These positions are contained in the 6 entries y(145:150).
    
   
%Now compute the matrix G.  Since 'G' already denotes the gravatational 
%constant call the matrix G 'GMatrix'. 
%This is done by calling 'G_CRTBP.m'
%GMatrix为a对r求偏导
GMatrix=Arho1+dArho2_dr;

%SECOND
%The right hand side for the state transition system is computed
%To do this construct a matrices 'A' 'I' 'K' and 'O', where 'A' contains the 
%variables and 'O' is a 3X3 matrix of zeros and 'I' is the 3X3 identity. 
%K为a对v求偏导
O=zeros(3);
I=eye(3);
K= -2*Omega_li;


%Now the complete Jacobian of f is assembled
Df=0;
Df=[ O,       I;
     GMatrix, K];

%Make A
for i=1:6
    for j=1:6
        A(i,j)=y(6*(i-1)+j);
    end
end

%Then compute the 6X6 matrix Df*A, is named DfA.
DfA=0;
DfA=Df*A;

%This has to be put into vector format. We temporaly place the results in 
%the 36-vector 'a'.  Later this will be the first 36 components of ydot.
a=0;
for i=1:6
    for j=1:6
        a(6*(i-1)+j)=DfA(i,j);
    end
end    

%THIRD
%The last 6 entries are the vector field for the CRTBP,
%These are stored in 'c'.


%constructs a vector whose first 3 entries are the velocities and whose
%last 3 entries are the accelerations

c=ydot;

%Put it all toghether and pass back to integrator

ydot=[a';c];