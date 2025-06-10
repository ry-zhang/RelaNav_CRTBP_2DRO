%--------------------------------------------------------------------------
%
%        CRTBP Autonomous Orbit Determination application 
%               
%        EKF estimation function
%
%        Use together with CRTBP_OD_EKF_main.m
%
% Inputs:
% P0            Initial State Covariance matrix m and m/s diag(12by12)
% Yerr          Initial State Error m and m/s (1by12)
% SC1           First S/C initial states from the database (1-16)
% SC2           Second S/C initial states from the database (1-16)
% simdur        Simulation duration in days
% simstep       
% rm            Satellite-to-satellite range measurement logical 1 or 0 
% sigma_range   SST range measurement error in Simluation time step in seconds meter (1 sigma)
% bias_range        SST range measurement bias in meter 
% 
% Output:
% rsspos        RSS position estimation state vector
% rssvel        RSS velocity estimation state vector
% rssposunc     RSS position uncertainty vector
% rssvelunc     RSS velocity uncertainty vector
% ErrYfull      Full Estimation error vector
% SigUncfull    Full Uncertainty vector
% All mean values are displayed in the end including plots
%
% Note: This script is used to run CRTBP_OD_EKF_main.m
%       

function [structRss, structState, distance] = ...
    CRTBP_CD_EKF_func(P0, Cerr, SC1, SC2, rm, simdur, simstep, sigma_range, bias_range)


%% 主要部分
global AuxParam

AuxParam.rangemeas      = 2;%R1,RA2
AuxParam.AdapFilter     = 0;
AuxParam.writeLog       =0;
AuxParam.CRLB           =0;
AuxParam.infowindow     =0;

AuxParam.estmode=1;%estRelaSta==1； estAbsSta_OneKnow==2；estAbsSta_twoEst=3;

if AuxParam.writeLog == 1
    %% 日志文件
    % 定义日志文件的名称，这里使用.log扩展名
    logFileName = 'EKF.log';
    
    % 打开文件以进行写入，如果文件不存在则创建它
    fileID = fopen(logFileName, 'wt');  % 'wt' 表示写入文本文件，如果文件存在则覆盖
    
    % 检查文件是否成功打开
    if fileID == -1
        error('无法打开日志文件进行写入。');
    end
end

% 加载常数
const


% 模拟测量数据的噪声
sigma_range = (1/LU)*sigma_range;  %range measurement in LU

% 测量噪声矩阵R
sigma_range_R = sigma_range*sqrt(1);
sigma_angle_R=0.0115*pi/180;
%step size in TU
step=simstep/TU;

%load initial states
IC_CRTBP
 
%simulation duration in days
daysnum=simdur;
n_obs = 4*(daysnum/(T0/86400))/step;
% n_obs = daysnum*86400/TU/step;

%主星（参考星）
Sc_1st=SC1;%dro

%相对轨道
load('Sol_linear.mat', 'Sol_linear');
load('Meigva_diag.mat','Meigva_diag');
%  Sol_linear=FloquetTheory(IC(Sc_1st,:));

x0_REL= [];
for jj_index = 1:6
        eval(['x0 = Sol_linear.vec',num2str(jj_index)])%,'*1e-5;'
        x0_REL(jj_index,:) =x0';
end
% x0_REL=[0.187399594815906	0.322457023672098	0	0.294610632214807	-0.559906894515994	0
% -0.0890406488847317	0.0552672741836855	0	-0.0695585398219571	0.197401932212126	0
% 0.147978004610028	0.480211879725382	0	0.0498699388446720	-0.436584298108024	0
% -0.0673423506744282	0.0306891039216372	0	-0.0729108447457729	0.116638296962979	0
% 0	0	0.237015932273651	0	0	0.730713560055757
% 0	0	-0.443366543113717	0	0	0.518191958845804];


x0_REL_true=x0_REL;


% Initial States 初始状态真值
Y0_true=[];
YL_true=[];
if AuxParam.estmode==1
%     Y0_rel_1=[-0.0562,0 ,0,0,0.0972,0]'*10^-05;
%     Y0_rel_2=[-0.079339635806951,0,0,0,0.186074431780780,0]'*10^-05;
%     Y0_rel_3=[0,0.576120458608156,0,0.817364800546933,0,0]'*10^-05;
%     Y0_rel_4=[-0.056189340713224,0,0,0,0.097258792961631,0]'*10^-05;
%     Y0_rel_5=[0,0,0,0,0,0.851143133756736]'*10^-05;
%     Y0_rel_6=[0,0,-0.524933677581048,0,0,0]'*10^-05;

    %相对状态初值
    Y0_rel_1=x0_REL(1,:)';%拟周期1
    Y0_rel_2=x0_REL(2,:)';%拟周期2
    Y0_rel_3=x0_REL(3,:)';%周期
    Y0_rel_4=x0_REL(4,:)';%发散
    Y0_rel_5=x0_REL(5,:)';%法向1
    Y0_rel_6=x0_REL(6,:)';%法向2

%     N1=50000/ Y0_rel_3(2)/LU;
%     Y0_true= N1*Y0_rel_3;
%     C0_true=[10^-5 10^-5 10^-5 0 10^-5 10^-5];%模态系数真值
     C0_true=[0 0 10^-5 0  0 0];%模态系数真值
    Y0_true=(C0_true*x0_REL)';%相对状态真值
    YL_true = [IC(Sc_1st,:)]';%dro主星初值 

%     YL_true=[-0.148303424555692	0.139644710773794	0	0.198220896282056	0.403698312959305	0]';
    C0_true=C0_true';
elseif AuxParam.estmode==2
    Y0_true = [IC(Sc_1st,:)]';%dro主星初值
    YL_true = [IC(Sc_2nd,:)]';%LEO真值
elseif AuxParam.estmode==3
    Y0_true = [IC(Sc_1st,:),IC(Sc_2nd,:)]';%dro+LEO初值
end


%Gravational constant
G=1;
mu=0.012150584269940;

% Initialization
t = 0;

%Initial state 状态初值=状态真值+初始误差
Cerr(:) =normrnd(0,Cerr(:));%C初始误差
Yerr =Cerr*x0_REL;%Y初始误差

C=C0_true+Cerr';
Y = Y0_true+Yerr';


%Initial covariance matrix 初始状态协方差矩阵
P = (2*10^-6)^2*eye(6,6);

% P = zeros(6);
% for i= 1:3
%     P(i,i) = Yerr(i);
% end
% for i=4:6
%     P(i,i) = Yerr(i);
% end

% 费雪信息矩阵的初值
J = inv(P);
Phi_CRLB=eye(6);%从开始时刻到结束的状态转移矩阵
Gramian=zeros(6,6);

% Initial Measurement Covariance 初始过程噪声协方差矩阵
sigm1 = 1e-9;
sigm2 = 1e-9;
sigm3 = 1e-9;
% delt = step;
% Qdt1=[((delt^4)*sigm1)/4 0 0 ((delt^3)*sigm1)/2 0 0;
%         0 ((delt^4)*sigm2)/4 0 0 ((delt^3)*sigm2)/2 0
%         0 0 ((delt^4)*sigm3)/4 0 0 ((delt^3)*sigm3)/2
%         ((delt^3)*sigm1)/2 0 0 ((delt^2)*sigm1) 0 0
%         0 ((delt^3)*sigm2)/2 0 0 ((delt^2)*sigm2) 0
%         0 0 ((delt^3)*sigm3)/2 0 0 ((delt^2)*sigm3)];%6*6
% Qdt1=[((delt^4)*sigm1)/4 0 0 ((delt^3)*sigm1)/2 0 0;
%         0 ((delt^4)*sigm2)/4 0 0 ((delt^3)*sigm2)/2 0
%         0 0 ((delt^4)*sigm3)/4 0 0 ((delt^3)*sigm3)/2
%         ((delt^3)*sigm1)/2 0 0 ((delt^2)*sigm1) 0 0
%         0 ((delt^3)*sigm2)/2 0 0 ((delt^2)*sigm2) 0
%         0 0 ((delt^3)*sigm3)/2 0 0 ((delt^2)*sigm3)];%6*6
% 
% tao=[delt^2/2 0 0 ;
%         0 delt^2/2 0; 
%         0 0 delt^2/2 ;
%         delt 0 0 ;
%         0 delt 0 ;
%         0 0 delt ];%6*3
% Q=Qdt1;

tao=[zeros(3,3);eye(3,3)];
Qc=(sigm1)^2*eye(6,6);
Q=Qc;


%Adaptive parameter (遗忘因子)
if AuxParam.AdapFilter == 1
    alpha=0.4;
else
    alpha=1;
end

% 测量误差的期望
bias_angle=bias_range;

k=0;
distance = [];

% ODE option
options=odeset('RelTol',1e-13,'AbsTol',1e-20);

pres=[];
vres=[];
% 进入导航
for i=1:n_obs

    k=k+1;
    % Previous step
    t_old = t;
    Y_old = Y;
    C_old = C;
    YL_true_old = YL_true;
    x0_REL_true_old=x0_REL_true;

    % Propagation to measurement epoch
    t = i*step;     % Time since epoch
    tspan = linspace(t_old, t, 2);

    %True trajectory真实轨道
    Y0_true_old = Y0_true; % 用于计算克拉美罗下限
    % 主星轨道和相对轨道
    [Y2,Y1,Sia_true]=LT(C0_true,YL_true,x0_REL_true_old,tspan,t);
    Y0_true=Y1;
    YL_true=Y2;
    x0_REL_true=Sia_true';




%     %% 时间更新
%     %%%时间更新Time Update
%     % 时间更新step1: 一步状态预测
%     %Estimated trajectory估计的轨迹
%     % 1st S/C
    [Y2,Y1,Sia]=LT(C,YL_true_old,x0_REL,tspan,t);
%     [~,xx]=ode113(@(t,x)CRTBP_RelaDyn_Noliner(t,x,mu),tspan, [YL_true;Y],options);
%     Y1=xx(end,7:12); 
%     Y2=xx(end,1:6);
     Y=Y1;
     x0_REL=Sia';

        
    % State-Transition Matrix 状态转移矩阵

%     STM=STM_CRTBP_RelaDyn_Noliner(t_old,t,[YL_true_old;Y_old],mu);
%      Phi=STM;
%      whitenoise=[normrnd(0,sigm1);normrnd(0,sigm2);normrnd(0,sigm3)];
%      Y=Y+tao*whitenoise;
%     
%      yres=Phi*Y_old-Y;
%      pres(k,:)=norm(yres(1:3))*LU;
%      vres(k,:)=norm(yres(4:6))*SU;
%     
%     % 时间更新step2: 先验估计的协方差矩阵的递推公式
%     P = Phi*P*Phi'+Q;
%       P=P+Q;
  
      
      whitenoise=[normrnd(0,sigm1);normrnd(0,sigm2);normrnd(0,sigm3)];
      C=C+tao*whitenoise;
%       P=Sia*P*Sia'+Q;
      P=P+Q;
    

    %% 测量更新
   bConsiderSignalBlock = 0;
   %  bConsiderSignalBlock = LinkIsBlocked(Y1', Y2');

    % 设计矩阵H以及测量噪声协方差矩阵R
    H = [];
    R = [];
    if AuxParam.rangemeas==1
        R = diag([sigma_range_R.^2]);
    elseif AuxParam.rangemeas==2
        R = diag([sigma_range_R.^2,sigma_angle_R.^2,sigma_angle_R.^2]);
      
%         R = diag([sigma_range_R.^2,sigma_angle_R.^2]);
    end
    % 若被遮挡，则只做时间更新；不被遮挡，进行时间和测量更新
    if bConsiderSignalBlock == true
        value_OC = NaN; % 赋空值
        distance = [distance;NaN];

    else
        % Measurements from models and their partials
%         dr        = Y(1:3)-Y(7:9); 
%         dr= Y(1:3)-YL_true(1:3); 
        dr= Y(1:3); %lvlh 相对状态b-a
%         dr= Y0_true(1:3); %lvlh 相对状态
        [udr rho] = unit(dr);
        %测距range
        obs_range = norm(Y0_true(1:3))+normrnd(bias_range,sigma_range_R);
        obs_azi=atan2(Y0_true(2),Y0_true(1))+normrnd(bias_angle,sigma_angle_R);
        if(abs(obs_azi)>pi)
            if(obs_azi>pi)
                obs_azi=obs_azi-2*pi;
            else(obs_azi<-pi)
                obs_azi=obs_azi+2*pi;
            end
        end
        obs_ele=asin(Y0_true(3)/norm(Y0_true(1:3)))+normrnd(bias_angle,sigma_angle_R);

        %obs from model
        rangeHat = norm(dr);
        AziHat=atan2(Y(2),Y(1));
        EleHat=asin(Y(3)/norm(Y(1:3)));
        if(abs(obs_azi-AziHat)>=359*pi/180)
           if ((obs_azi <= 0 && obs_azi >= -pi))
					obs_azi = 2 * pi + obs_azi;
           end
		   if (AziHat <= 0 && AziHat >= -pi)
					AziHat = 2 * pi + AziHat;
           end
        end
     

       
        %design mat
        Hrr = shiftdim(udr,-1); %将udr的维度向右移1位，比如1*3的矩阵，大小被reshape为3*1
        Hrv = zeros(size(Hrr));
        He_x=-2*Y(1)*Y(3)/(2 * sqrt(1 - Y(3) * Y(3)  / (Y(1)*Y(1) +Y(2) * Y(2) +Y(3) * Y(3)))*(Y(1)*Y(1) +Y(2) * Y(2) +Y(3) * Y(3))^1.5);
        He_y=-2*Y(2)*Y(3)/(2 * sqrt(1 - Y(3) * Y(3)  / (Y(1)*Y(1) +Y(2) * Y(2) +Y(3) * Y(3)))*(Y(1)*Y(1) +Y(2) * Y(2) +Y(3) * Y(3))^1.5 );
        He_z=((Y(1)*Y(1)) + (Y(2)*Y(2))) / (sqrt(1 - Y(3) * Y(3)/(Y(1)*Y(1) +Y(2) * Y(2) +Y(3) * Y(3)))*(Y(1)*Y(1) +Y(2) * Y(2) +Y(3) * Y(3))^ 1.5);
        Ha_x = -Y(2) / (Y(1)^ 2*((Y(2)/Y(1))^2 + 1));
		Ha_y = 1 / (Y(1)*((Y(2)/Y(1))^ 2 + 1));
		Ha_z = 0;
   
        Hr  = [Hrr, Hrv]; 
        He=[He_x He_y He_z 0 0 0];
        Ha=[Ha_x Ha_y Ha_z 0 0 0];

        Hat=[];
        obs=[];
        H=[];
        if AuxParam.rangemeas==1
            Hat=rangeHat;
            obs=obs_range;
            H = Hr;

        elseif AuxParam.rangemeas==2
            Hat=[rangeHat;EleHat;AziHat];
            obs=[obs_range;obs_ele;obs_azi];
            Hx = [Hr;He;Ha];
            H=Hx*Sia;

            
        end


        % 先验残差及设计矩阵H

        d = (obs - Hat);
       
        %%%测量更新Measurement update
        % step3: 更新卡尔曼增益
        K = P*H'*inv(R+H*P*H');
    
        % step4: 状态校正
%         Y = Y + K*d;
        C=C+K*d;
        corrC=K*d;
    
        % step5: 状态协方差矩阵更新
%         P = (eye(12)-K*H)*P;
%          P = (eye(6)-K*H)*P;
        P = (eye(6)-K*H)*P*(eye(6)-K*H)'+K*R*K';
    
        % Q矩阵更新（自适应卡尔曼滤波，alpha为1时，等同于EKF）
        % Update Q
        Q = alpha*Q + (1-alpha)*(K*d*d'*K');

%         distance = [distance;norm(Y0_true(1:3)-Y0_true(7:9));];
%         value_OC = obs_range - norm(Y(1:3)-Y(7:9));
%         distance = [distance;norm(Y0_true(1:3)-YL_true(1:3));];
%         value_OC = obs_range - norm(Y(1:3)-YL_true(1:3));
       % distance = [distance;norm(Y0_true(1:3));];
        value_OC = d- H*K*d;
        
    end

    %% 克拉美罗下限（CRLB），用的是真实轨道
    if AuxParam.CRLB == 1
        % State-Transition Matrix 状态转移矩阵
        % 1st S/C
        STM_CRLB=STM_CRTBP_RelaDyn_Noliner(t_old,t,[YL_true_old;Y0_true_old],mu);

      
        Phi_CRLB=STM_CRLB*Phi_CRLB;
        


        dr = Y0_true(1:3); 
        [udr rho] = unit(dr);

        He_x=2*Y0_true(1)*Y0_true(3)/(2 * sqrt(1 - Y0_true(3) * Y0_true(3)  / (Y0_true(1)*Y0_true(1) +Y0_true(2) * Y0_true(2) +Y0_true(3) * Y0_true(3)))*(Y0_true(1)*Y0_true(1) +Y0_true(2) * Y0_true(2) +Y0_true(3) * Y0_true(3))^1.5);
        He_y=2*Y0_true(2)*Y0_true(3)/(2 * sqrt(1 - Y0_true(3) * Y0_true(3)  / (Y0_true(1)*Y0_true(1) +Y0_true(2) * Y0_true(2) +Y0_true(3) * Y0_true(3)))*(Y0_true(1)*Y0_true(1) +Y0_true(2) * Y0_true(2) +Y0_true(3) * Y0_true(3))^1.5 );
        He_z=-((Y0_true(1)*Y0_true(1)) + (Y0_true(2)*Y0_true(2))) / (sqrt(1 - Y0_true(3) * Y0_true(3)/(Y0_true(1)*Y0_true(1) +Y0_true(2) * Y0_true(2) +Y0_true(3) * Y0_true(3)))*(Y0_true(1)*Y0_true(1) +Y0_true(2) * Y0_true(2) +Y0_true(3) * Y0_true(3))^ 1.5);
        Ha_x = -Y0_true(2) / (Y0_true(1)^ 2*((Y0_true(2)/Y0_true(1))^2 + 1));
		Ha_y = 1 / (Y0_true(1)*((Y0_true(2)/Y0_true(1))^ 2 + 1));
		Ha_z = 0;
   

        H_CRLB (1,:)= [udr', zeros(1,3)];
        H_CRLB (2,:)= [He_x,He_y,He_z, zeros(1,3)];
        H_CRLB (3,:)= [Ha_x,Ha_y,Ha_z, zeros(1,3)];

       


        % 费雪信息矩阵
        if isempty(R)
            J = inv(Phi_CRLB*inv(J)*Phi_CRLB');
        else
            J = inv(Phi_CRLB*inv(J)*Phi_CRLB') + H_CRLB'*inv(R)*H_CRLB;
        end


         Gramian=Gramian+Phi_CRLB'*H_CRLB'*H_CRLB*Phi_CRLB;
         condnum = cond(Gramian);
         eigenvalues = eig(Gramian);
         min_eigenvalue = min(eigenvalues);        % 获取最小特征值

    
        % P >= inv(J)
        P_CRLB = inv(J);
        CRLB = zeros(6,1);
        for ii=1:6
            CRLB(ii) = sqrt((P_CRLB(ii,ii)));
        end


        %RSS Pos and Vel errors 位置速度残差(CRLB)
        rssCRLBpos(k,:)=norm(CRLB(1:3));
        rssCRLBvel(k,:)=norm(CRLB(4:6));
        fullCRLB(k,:)=CRLB;
        rssCRLBCN(k,:)=condnum;
        rssCRLBOI(k,:)=min_eigenvalue;

    end

    %% 残差统计
    % OC后验残差
    rssposOC(k,:) = value_OC(1,:);
    rsseleOC(k,:) = value_OC(2,:)*180/pi;
    rssaziOC(k,:) = value_OC(3,:)*180/pi;
    % OC先验残差
%     rssposOC(k,:) = d(1,:);
%     rsseleOC(k,:) = d(2,:)*180/pi;
%     rssaziOC(k,:) = d(3,:)*180/pi;
%     rssposOC(k,:) =d(1,:);
%     rsseleOC(k,:) = AziHat*180/pi;
%     rssaziOC(k,:) = obs_azi*180/pi;
    
    % 状态残差
    ErrY = Y-Y0_true;
    ErrC = C-C0_true;
%     Sigma = zeros(12,1);
    Sigma = zeros(6,1);
%     for ii=1:12
    for ii=1:6
        Sigma(ii) = sqrt(P(ii,ii));
    end

    %RSS Pos and Vel errors 位置速度残差以及不确定度残差
%     rsspos(k,:)=[norm(ErrY(1:3)),norm(ErrY(7:9))];
%     rssvel(k,:)=[norm(ErrY(4:6)),norm(ErrY(10:12))];
%     rssposunc(k,:)=[norm(Sigma(1:3)),norm(Sigma(7:9))];
%     rssvelunc(k,:)=[norm(Sigma(4:6)),norm(Sigma(10:12))];
    rsspos(k,:)=[norm(ErrY(1:3))];
    rssvel(k,:)=[norm(ErrY(4:6))];
    rssposunc(k,:)=[norm(Sigma(1:3))];
    rssvelunc(k,:)=[norm(Sigma(4:6))];
    

    % Full history
    rssErrYfull(k,:) = ErrY;
    rssSigUncfull(k,:) = Sigma;
    estState(k,:) = Y;
    trueState(k,:) = Y0_true;
    CSia=zeros(6,6);
    if(i==n_obs)
        CSia=x0_REL_true;
    end

    rssErrCfull(k,:) = ErrC;
    rssSigUncCfull(k,:) = Sigma;
    estC(k,:) = C;
    trueC(k,:) = C0_true;
    YAState(k,:) = YL_true;

    if any(i == [0:1000:n_obs])
        formatSpec = 'Simulation is running: %4.2f ';
        disp([sprintf(formatSpec,(i/n_obs)*100) '%'])

        if AuxParam.infowindow == 1
            if getappdata(f,'canceling')
                break
            end
                
            % Update waitbar and message
            waitbar((i/n_obs),f,sprintf('Simulation is running: %4.2f %%',...
                (i/n_obs)*100))
        end  
    end


    %% 写入日志内容
    if AuxParam.writeLog == 1
        % 确定输出宽度（例如：10字符宽度，保留6位小数）
        formatSpec = '%10.6f\t';  % 每个元素占10个字符，6位小数
        fprintf(fileID, 'Epoch %d\n', i);
        fprintf(fileID, 'Hx设计矩阵：\n');
        for ii = 1:size(Hx, 1)
            fprintf(fileID, formatSpec, Hx(ii, :));
            fprintf(fileID, '\n');  % 换行符
        end
        fprintf(fileID, '设计矩阵：\n');
        for ii = 1:size(H, 1)
            fprintf(fileID, formatSpec, H(ii, :));
            fprintf(fileID, '\n');  % 换行符
        end
        fprintf(fileID, '后验状态协方差矩阵：\n');
        for ii = 1:size(P, 1)
            fprintf(fileID, formatSpec, P(ii, :)/10^-12);
            fprintf(fileID, '\n');  % 换行符
        end
        fprintf(fileID, '卡尔曼增益矩阵：\n');
        for ii = 1:size(K, 1)
            fprintf(fileID, formatSpec, K(ii,:));
            fprintf(fileID, '\n');  % 换行符
        end
        fprintf(fileID, '\n');  % 换行符

         fprintf(fileID, '模态矩阵：\n');
        for ii = 1:size(Sia, 1)
            fprintf(fileID, formatSpec, Sia(ii,:));
            fprintf(fileID, '\n');  % 换行符
        end
        fprintf(fileID, '\n');  % 换行符
         fprintf(fileID, '矫正量：\n');
        for ii = 1:3
            fprintf(fileID, formatSpec, corrC(ii, :)*10^6);          
        end
        for ii = 4:6
            fprintf(fileID, formatSpec, corrC(ii, :)*10^6);
        end
        fprintf(fileID, '\n');  % 换行符
    end

end

% % % %% 存储至结构体
% % % %状态序列
structState.estState = estState;
structState.trueState = trueState;
structState.estC= estC;
structState.trueC= trueC;
structState.yA=YAState;
structState.CSia=CSia;
% 残差序列
structRss.rsspos=rsspos;
structRss.rssvel=rssvel;
structRss.pres=pres;
structRss.vres=vres;
structRss.rssposunc=rssposunc;
structRss.rssvelunc=rssvelunc;
structRss.rssErrYfull=rssErrYfull;
structRss.rssSigUncfull=rssSigUncfull;
structRss.rssErrCfull=rssErrCfull;
structRss.rssSigUncCfull=rssSigUncCfull;
structRss.rssposOC=rssposOC;
structRss.rsseleOC=rsseleOC;
structRss.rssaziOC=rssaziOC;
structRss.CRLB=AuxParam.CRLB;
if AuxParam.CRLB == 1
    structRss.rssCRLBpos=rssCRLBpos;
    structRss.rssCRLBvel=rssCRLBvel;
    structRss.rssCRLBCN=rssCRLBCN;
    structRss.rssCRLBOI=rssCRLBOI;
end

figure()
plot3(estState(:,1)*LU,estState(:,2)*LU,estState(:,3)*LU,'b');
hold on;
plot3(trueState(:,1)*LU,trueState(:,2)*LU,trueState(:,3)*LU,'r');
hold on;
% plot3(trueState(:,1)*LU,trueState(:,2)*LU,trueState(:,3)*LU,'r');
% hold on;
%  plot3(posEarth(1)*1,posEarth(2)*1,posEarth(3)*1, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
plot3(posMoon(1)*LU, posMoon(2)*LU,posMoon(3)*LU,'o', 'MarkerSize', 2, 'MarkerFaceColor', 'y');
hold on;
legend('估计轨迹','真实轨迹','月球位置','Location', 'northeast');
xlabel('x[m]', 'fontsize',14,'interpreter','latex')
ylabel('y[m]', 'fontsize',14,'interpreter','latex')
zlabel('z[m]', 'fontsize',14,'interpreter','latex')
title('月心LVLH下相对轨迹 ', 'fontsize',16,'interpreter','latex')
save('structState.mat','structState');
save('structRss.mat','structRss');



%% 关闭
if AuxParam.infowindow == 1
    delete(f)
end

if AuxParam.writeLog == 1
    % 关闭文件
    fclose(fileID);
end






