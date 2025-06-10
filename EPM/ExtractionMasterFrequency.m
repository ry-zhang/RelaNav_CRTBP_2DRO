% 1. 读取星历数据
load('XX_A_MCR_Eph_70d');  % 读取数据，假设每行是一个时刻的状态
data = xx_MCR_target;
time = (0:size(data, 1)-1) * 60;  % 假设每个数据点的时间间隔为60秒
x = data(:, 1);  % 提取x位置
y = data(:, 2);  % 提取y位置
z = data(:, 3);  % 提取z位置
Fs = 1 / 60;  % 采样频率（60秒一个采样点，即每分钟一个点）

% 设计高通滤波器，去除低频部分，设定截止频率为1e-3Hz（可以根据需要调整）
d = designfilt('highpassfir', 'FilterOrder', 8, 'CutoffFrequency', 1e-6, 'SampleRate', Fs);

% 对x, y, z分量应用高通滤波器
x_filtered = filter(d, x);
y_filtered = filter(d, y);
z_filtered = filter(d, z);

% 2. 对XY平面位置进行傅里叶变换以提取主频
N = length(x);  % 数据点数
f = Fs * (0:(N/2)) / N;  % 频率范围（正频率部分）

% 去除均值（去趋势）
x_detrended = x_filtered - mean(x_filtered);  % 去除x分量的均值
y_detrended = y_filtered - mean(y_filtered);  % 去除y分量的均值
z_detrended = z_filtered - mean(z_filtered);  % 去除z分量的均值

% 对去趋势后的数据进行傅里叶变换
XY_signal_detrended = x_detrended + 1i * y_detrended;  % 合成复数信号，方便做傅里叶变换
XY_f_detrended = fft(XY_signal_detrended);

% 计算频谱（取正频率部分）
P_XY_detrended = abs(XY_f_detrended(1:N/2+1));  % XY 平面频谱

% 寻找主频（最大峰值位置对应的频率）
[~, idx_XY1] = max(P_XY_detrended);  % XY 平面主频位置
f_XY1 = f(idx_XY1);  % XY 平面方向的主频

% 寻找第二大频率
P_XY_detrended_filtered = P_XY_detrended;
P_XY_detrended_filtered(idx_XY1) = 0;  % 排除最大频率位置
[~, idx_XY2] = max(P_XY_detrended_filtered);  % 第二大频率位置在 XY 平面
f_XY2 = f(idx_XY2);  % XY 平面方向的第二大频率

% 输出结果
disp(['主频（XY平面）：', num2str(f_XY1), ' Hz']);
disp(['第二大频率（XY平面）：', num2str(f_XY2), ' Hz']);

% 3. 对z方向进行傅里叶变换
Z_f_filtered = fft(z_detrended);
P_z_filtered = abs(Z_f_filtered(1:N/2+1));  % 位置 z 的频谱

% 寻找z方向的主频和第二大频率
[~, idx_z1] = max(P_z_filtered);  % 主频位置在 z 分量
P_z_filtered_filtered = P_z_filtered;
P_z_filtered_filtered(idx_z1) = 0;  % 排除最大频率位置
[~, idx_z2] = max(P_z_filtered_filtered);  % 第二大频率位置在 z 分量

% 主频
f_z1 = f(idx_z1);  % z 方向的主频
f_z2 = f(idx_z2);  % z 方向的第二大频率

% 输出结果
disp(['主频（z方向）：', num2str(f_z1), ' Hz']);
disp(['第二大频率（z方向）：', num2str(f_z2), ' Hz']);

% 4. 绘制XY平面频谱图
figure;
subplot(2,1,1);
plot(f, P_XY_detrended);
title('XY平面频谱');
xlabel('频率（Hz）');
ylabel('幅度');

% 绘制z方向频谱图
subplot(2,1,2);
plot(f, P_z_filtered);
title('z方向频谱');
xlabel('频率（Hz）');
ylabel('幅度');
