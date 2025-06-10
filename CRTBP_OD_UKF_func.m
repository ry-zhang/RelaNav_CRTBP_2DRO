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
    CRTBP_OD_EKF_func(P0, Yerr, SC1, SC2, rm, simdur, simstep, sigma_range, bias_range)


%% 主要部分
global AuxParam

AuxParam.rangemeas      = 2;%R1,RA2
AuxParam.AdapFilter     = 0;
AuxParam.writeLog       = 0;
AuxParam.CRLB           = 0;
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


Sc_1st=SC1;%dro
%Sc_2nd=SC2;%leo

% Initial States 初始状态真值
Y0_true=[];
YL_true=[];
if AuxParam.estmode==1
    Y0_rel_3=[0,0.576120458608156,0,0.817364800546933,0,0]'*10^-05;%相对状态初值
    Y0_rel_5=[0,0,0,0,0,0.851143133756736]'*10^-05;%相对状态初值
    N=50000/ Y0_rel_3(2)/LU;
    Y0_true=N*(1*Y0_rel_3+1*Y0_rel_5);
    YL_true = [IC(Sc_1st,:)]';%dro主星初值 
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
Yerr(:) =normrnd(0,Yerr(:));
Yerr=Yerr.*[ones(1,3)*(1/LU) ones(1,3)*(1/SU)]; %初始误差
Y = Y0_true+Yerr';

%Initial covariance matrix 初始状态协方差矩阵
P = zeros(6);

for i= 1:3
    P(i,i) = P0(i)*((1/LU)^2);
end
for i=4:6
    P(i,i) = P0(i)*((1/SU)^2);
end

% 费雪信息矩阵的初值
J = inv(P);
Phi_CRLB=eye(6);%从开始时刻到结束的状态转移矩阵
Gramian=zeros(6,6);

% Initial Measurement Covariance 初始过程噪声协方差矩阵
sigm1 = 1e-8;
sigm2 = 1e-8;
sigm3 = 1e-8;
delt = step;

Qdt1=[((delt^4)*sigm1)/4 0 0 ((delt^3)*sigm1)/2 0 0;
        0 ((delt^4)*sigm2)/4 0 0 ((delt^3)*sigm2)/2 0;
        0 0 ((delt^4)*sigm3)/4 0 0 ((delt^3)*sigm3)/2;
        ((delt^3)*sigm1)/2 0 0 ((delt^2)*sigm1) 0 0;
        0 ((delt^3)*sigm2)/2 0 0 ((delt^2)*sigm2) 0;
        0 0 ((delt^3)*sigm3)/2 0 0 ((delt^2)*sigm3)];%6*6

tao=[delt^2/2 0 0 ;
        0 delt^2/2 0; 
        0 0 delt^2/2 ;
        delt 0 0 ;
        0 delt 0 ;
        0 0 delt ];%6*3
Q=Qdt1;

%Adaptive parameter (遗忘因子)
if AuxParam.AdapFilter == 1
    alpha=0.8;
else
    alpha=1;
end

% 测量误差的期望
bias_angle=bias_range;

k=0;
distance = [];

% ODE option
options=odeset('RelTol',1e-16,'AbsTol',1e-22);

pres=[];
vres=[];
% 进入导航
for i=1:n_obs

    k=k+1;
    % Previous step
    t_old = t;
    Y_old = Y;
    YL_true_old = YL_true;

    % Propagation to measurement epoch
    t = i*step;     % Time since epoch
    tspan = linspace(t_old, t, 45);

    %True trajectory真实轨道
    Y0_true_old = Y0_true; % 用于计算克拉美罗下限
    % 主星轨道和相对轨道
    [~,xx]=ode113(@(t,x)CRTBP_RelaDyn(t,x,mu),tspan, [YL_true;Y0_true],options);
    Y1=xx(end,7:12); %相对轨道(lvlh)

    Y2=xx(end,1:6); %主星轨道(M)

     Y0_true=Y1';
     YL_true=Y2';

    %% 时间更新
    %%%时间更新Time Update
    n = numel(Y);%状态维数
    m =3;%观测维数
    alpha=0.05;%默认系数
    ki=0;    %默认系数：3-n
    beta=2;%默认系数
    lambda=alpha^2*(n+ki)-n;%默认系数
    
    %%Step1：利用无迹变换公式获得一组采样点(称为Sigma点集)并计算相对于的权值
    %利用无迹变换公式获得一组采样点
    Xsigmaset=sigmas(Y,P,lambda,n);
    %计算相对应的权值
    Wm = lambda/(n+lambda);
    Wc = lambda/(n+lambda)+(1 - alpha^2 + beta); 
    for i = 1:2*n
        Wm = [Wm 0.5/(n+lambda)];
        Wc = [Wc 0.5/(n+lambda)];
    end
   
    
    %%Step2，Step3，Step4：对Sigma点集进行一步预测，得到均值X1means和方差P1和新的Sigma点集X1
    [X1means,X1,P1,X2]=ut(Xsigmaset,Wm,Wc,n,Q,tspan, [YL_true;Y],options,mu);%对状态进行UT变换%对状态进行UT变换
    Y=X1means;
    P=P1;

  

%     % 时间更新step1: 一步状态预测
%     %Estimated trajectory估计的轨迹
%     % 1st S/C
%     [~,xx]=ode113(@(t,x)CRTBP_RelaDyn(t,x,mu),tspan, [YL_true;Y],options);
%     Y1=xx(end,7:12); 
%     Y2=xx(end,1:6);
%      Y=Y1';
%         
%     % State-Transition Matrix 状态转移矩阵
%     % 1st S/C
%     STM=STM_CRTBP_RelaDyn(t_old,t,[YL_true_old;Y_old],mu);
%      Phi=STM;
%     
%      whitenoise=[normrnd(0,sigm1);normrnd(0,sigm2);normrnd(0,sigm3)];
%      Y=Y+tao*whitenoise;
% 
%     
% %      yres=Phi*Y_old-Y;
% %      pres(k,:)=norm(yres(1:3))*LU;
% %      vres(k,:)=norm(yres(4:6))*SU;
%     
%     % 时间更新step2: 先验估计的协方差矩阵的递推公式
%     P = Phi*P*Phi'+Q;

    

    %% 测量更新
   bConsiderSignalBlock = 0;
  %%Step5，Step6：得到观测预测，Z1为X1集合的预测，Zpre为Z1的均值
%   ob_temp=zeros(3);
%   ob_tempv=zeros(3);
%   ob_temp2=zeros(3);
%   ob_temp3=zeros(3);

  L=size(Xsigmaset,2);

  for k=1:L
%    for i=1:3
%       ob_temp(i)=X1(i,k);
%    end
  
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
        Y=X1(:,k);
        dr= Y(1:3); %lvlh 相对状态b-a

        [udr rho] = unit(dr);
        %测距range
        obs_range = norm(Y0_true(1:3))+normrnd(bias_range,sigma_range_R);
        obs_azi=atan2(Y0_true(2),Y0_true(1))+normrnd(bias_angle,sigma_angle_R);
        if(obs_azi<0)
            obs_azi=obs_azi+2*pi;
        end
        obs_ele=asin(Y0_true(3)/norm(Y0_true(1:3)))+normrnd(bias_angle,sigma_angle_R);


   
        %obs from model
        rangeHat = norm(dr);
        AziHat=atan2(Y(2),Y(1));
        EleHat=asin(Y(3)/norm(Y(1:3)));
        if(AziHat<0)
            AziHat=AziHat+2*pi;
        end

       
        %design mat
        Hrr = shiftdim(udr,-1); %将udr的维度向右移1位，比如1*3的矩阵，大小被reshape为3*1
        Hrv = zeros(size(Hrr));
        He_x=2*Y(1)*Y(3)/(2 * sqrt(1 - Y(3) * Y(3)  / (Y(1)*Y(1) +Y(2) * Y(2) +Y(3) * Y(3)))*(Y(1)*Y(1) +Y(2) * Y(2) +Y(3) * Y(3))^1.5);
        He_y=2*Y(2)*Y(3)/(2 * sqrt(1 - Y(3) * Y(3)  / (Y(1)*Y(1) +Y(2) * Y(2) +Y(3) * Y(3)))*(Y(1)*Y(1) +Y(2) * Y(2) +Y(3) * Y(3))^1.5 );
        He_z=-((Y(1)*Y(1)) + (Y(2)*Y(2))) / (sqrt(1 - Y(3) * Y(3)/(Y(1)*Y(1) +Y(2) * Y(2) +Y(3) * Y(3)))*(Y(1)*Y(1) +Y(2) * Y(2) +Y(3) * Y(3))^ 1.5);
        Ha_x = -Y(2) / (Y(1)^ 2*((Y(2)/Y(1))^2 + 1));
		Ha_y = 1 / (Y(1)*((Y(2)/Y(1))^ 2 + 1));
		Ha_z = 0;
   
        Hr  = [Hrr Hrv]; 
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
            H = [Hr;He;Ha];
         
        end
    end
    Z_dObservation=obs;
    Hx_dObservation=Hat;
    Z_SigmaPoint(:,k)=Hx_dObservation;
  end

    Z_eXpect=zeros(3,1);%均值
    Pzz=zeros(6,6);
    Zdiv=zeros(m,L);
    for k=1:L
        Z_eXpect=Z_eXpect+Wm(k)*Z_SigmaPoint(:,k); 
    end
    for k=1:L
        Zdiv(:,k)=Z_SigmaPoint(:,k)-Z_eXpect(:);%预测减去均值
    end
    Pzz=Zdiv*diag(Wc)*Zdiv'+R;%协方差

    Pxz=X1*diag(Wc)*Zdiv';%计算交叉协方差矩阵
    K=Pxz*inv(Pzz);%计算kalman增益
 

    % 先验残差
    d = (Z_dObservation - Z_eXpect);
   


    % step4: 状态校正
    Y = Y +  K*d;
    corrY= K*d;

    % step5: 状态协方差矩阵更新
   

    P = P-K*Pxz';

    value_OC = d- H*K*d;
        
%     end

    %% 克拉美罗下限（CRLB），用的是真实轨道
    if AuxParam.CRLB == 1
        % State-Transition Matrix 状态转移矩阵
        % 1st S/C
        STM_CRLB=STM_CRTBP_RelaDyn(t_old,t,[YL_true_old;Y0_true_old],mu);

      
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
%     rsseleOC(k,:) = EleHat*180/pi;
%     rssaziOC(k,:) = obs_ele*180/pi;
    
    % 状态残差
    ErrY = Y-Y0_true;
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
    YAState(k,:) = YL_true;

    if any(i == [0:100:n_obs])
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
        fprintf(fileID, '状态转移矩阵：\n');
        for ii = 1:size(Phi, 1)
            fprintf(fileID, formatSpec, Phi(ii, :));
            fprintf(fileID, '\n');  % 换行符
        end
        fprintf(fileID, '后验状态协方差矩阵：\n');
        
        for ii = 1:3
            for iii= 1:size(Phi, 2)
            if iii<=3&&iii>=1
                fprintf(fileID, formatSpec, P(ii, iii)*LU^2);
            else
                fprintf(fileID, formatSpec, P(ii, iii)*LU*SU);
            end
            end
             fprintf(fileID, '\n');  % 换行符
        end
        for ii = 4:6
            for iii= 1:size(Phi, 2)
            if iii<=3&&iii>=1
                fprintf(fileID, formatSpec, P(ii, iii)*LU*SU);
            else
                fprintf(fileID, formatSpec, P(ii, iii)*SU^2);
            end
%             fprintf(fileID, formatSpec, P(ii, :)*SU^2);
            end
            fprintf(fileID, '\n');  % 换行符
        end
        fprintf(fileID, '卡尔曼增益矩阵：\n');
        for ii = 1:size(K, 1)
            fprintf(fileID, formatSpec, K(ii,:));
            fprintf(fileID, '\n');  % 换行符
        end
        fprintf(fileID, '\n');  % 换行符

        fprintf(fileID, '矫正量：\n');
        for ii = 1:3
            fprintf(fileID, formatSpec, corrY(ii, :)*LU);          
        end
        for ii = 4:6
            fprintf(fileID, formatSpec, corrY(ii, :)*SU);
        end
        fprintf(fileID, '\n');  % 换行符
    end

end

%% 存储至结构体
%状态序列
structState.estState = estState;
structState.trueState = trueState;
structState.yA=YAState;
%残差序列
structRss.rsspos=rsspos;
structRss.rssvel=rssvel;
structRss.pres=pres;
structRss.vres=vres;
structRss.rssposunc=rssposunc;
structRss.rssvelunc=rssvelunc;
structRss.rssErrYfull=rssErrYfull;
structRss.rssSigUncfull=rssSigUncfull;
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
% plot3(YAState(:,1)*LU,YAState(:,2)*LU,YAState(:,3)*LU,'g');
% hold on;
%  plot3(posEarth(1)*1,posEarth(2)*1,posEarth(3)*1, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
plot3(posMoon(1)*LU, posMoon(2)*LU,posMoon(3)*LU,'o', 'MarkerSize', 2, 'MarkerFaceColor', 'g');
hold on;
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






