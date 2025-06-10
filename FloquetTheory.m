function [Sol_linear]=FloquetTheory(x0_DRO_M_3d)
mu= 0.012150584269940; 
para.T0=0.5*2*pi;
opts_ode=odeset('RelTol',1e-13,'AbsTol',1e-20);
%% 线性化相对运动的周期状态转移矩阵（φ(0)^-1 φ(T)）
sol_temp = ode113(@(t,x)CRTBP_RelaDynPhi(t,x,mu),[0 para.T0], [x0_DRO_M_3d, zeros(1,6),reshape(eye(6),1,36)], opts_ode);
M_REL_lin = reshape(sol_temp.y(13:end,end),6,6);
[Meigve, Meigva] = eig(M_REL_lin);%特征值
% remove small numbers in Meigve and obtain the accurate eigen vector
Meigve_imga = imag(Meigve); Meigve_imga(abs(Meigve_imga)<1e-5) = 0;
Meigve_real = real(Meigve); Meigve_real(abs(Meigve_real)<1e-5) = 0;
Meigve = Meigve_real+Meigve_imga*1i;

% 判断单位特征值在第1-2行，还是第3-4行
[a,sort_index] = sort(real(diag(Meigva(1:4,1:4))));
index_temp1 = 3*(sort_index(1)>2) + 1*(sort_index(1)<=2);
index_temp2 = 1*(sort_index(1)>2) + 3*(sort_index(1)<=2);
H = blkdiag([real(Meigva(index_temp1,index_temp1)),-imag(Meigva(index_temp1,index_temp1)); ...
    imag(Meigva(index_temp1,index_temp1)),real(Meigva(index_temp1,index_temp1))],...
    [1,1; 0,1],...
    [real(Meigva(5,5)),-imag(Meigva(5,5)); imag(Meigva(5,5)),real(Meigva(5,5))]);
p2_real = real(Meigve(:,index_temp1)); p2_imag = imag(Meigve(:,index_temp1+1));
p3 = real(Meigve(:,index_temp2)); p3 = sign(p3(2)).*p3;
p4 = pinv(M_REL_lin-eye(6))*p3; % 由于M_REL_lin-eye(6)奇异，采用广义逆解
try
    p4 = p4-p4(2)/p3(2)*p3;
catch
end
p3 = p3/norm(p3); 
if norm(p4)>10 % 精度太高有时候算出来的p4超级大
    p4 = pinv(M_REL_lin-eye(6),1e-10)*p3;
end
p4 = p4/norm(p3);

p6_real = real(Meigve(:,6)); p6_imag = imag(Meigve(:,6));

S = [p2_real,p2_imag,p3,p4,p6_real,p6_imag];%特征向量
flag_Jordan = norm(S*H*S^-1-M_REL_lin);

J = 1/para.T0*logm(H);
J = J.*(abs(J)>1e-8);
omega1 = J(2,1);
omega2 = J(3,4); % 正好是1/T0
omega3 = J(6,5);

expJt = @(t)blkdiag([cos(omega1*t),-sin(omega1*t);sin(omega1*t),cos(omega1*t)],...
    [1,omega2*t;0,1],...
    [cos(omega3*t),-sin(omega3*t);sin(omega3*t),cos(omega3*t)]);
flag_Floquet = norm(M_REL_lin * S * expJt(-para.T0) - S);

% resort data
Meigva_diag = diag(Meigva)';
Meigva_diag = Meigva_diag([index_temp1,index_temp1+1,index_temp2,index_temp2+1,5,6]);%特征向量
Meigve = Meigve(:,[index_temp1,index_temp1+1,index_temp2,index_temp2+1,5,6]);

% alpha其实就是算拟周期的时候的旋转量，计算的时候取正虚部的特征值计算
alpha = [atan2(imag(Meigva_diag(2)),real(Meigva_diag(2))),...
        atan2(imag(Meigva_diag(6)),real(Meigva_diag(6)))];
T12_propotion = 2*pi./alpha;
% save data
Sol_linear.vec1 = p2_real;
Sol_linear.vec2 = p2_imag;
Sol_linear.vec3 = p3;
Sol_linear.vec4 = p4;
Sol_linear.vec5 = p6_real;
Sol_linear.vec6 = p6_imag;

% 这里采用第二个和第六个特征值，因为其虚部是正的，如此可以确保alpha大于0
EigenVector.p1_real = p2_real;
EigenVector.p1_imag = p2_imag;
EigenVector.p3 = p3;
EigenVector.p4 = p4;
EigenVector.p5_real = p6_real;
EigenVector.p5_imag = p6_imag;
save ../FloquetEig12 Meigve Meigva_diag EigenVector Sol_linear x0_DRO_M_3d