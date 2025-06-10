%%%%  第7段2 定义无迹变换子函数
function [Xmeans,Xsigma_pre,P,Xdiv]=ut(Xsigma,Wm,Wc,n,cov,tspan,YY,options,mu)

L=size(Xsigma,2);%得到Xsigma样本个数
Xmeans=zeros(n,1);%均值
Xsigma_pre=zeros(n,L);
for k=1:L
     [~,xx]=ode113(@(t,x)CRTBP_RelaDyn(t,x,mu),tspan, YY,options);
     Y1=xx(end,7:12); 
     Y2=xx(end,1:6);
     Xsigma_pre(:,k)=Y1;%计算一步预测
     Xmeans=Xmeans+Wm(k)*Xsigma_pre(:,k);
end
Xdiv=Xsigma_pre-Xmeans(:,ones(1,L));%预测减去均值
P=Xdiv*diag(Wc)*Xdiv'+cov;%协方差
end