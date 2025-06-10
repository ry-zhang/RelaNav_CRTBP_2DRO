%读取txt轨道文件

files_A = 'orbit-a.txt';
% files_AB = 'deltaX-ab.txt';
files_B = 'orbit-b.txt';
[tX_B_MCR_Eph]=ReadOrbit(files_B);

% XX_A_MCR_Eph=tX_A_MCR_Eph(:,2:7);
XX_B_MCR_Eph=tX_B_MCR_Eph(:,2:7);
% XX_AB_MCLVLH_Eph=tX_AB_MCLVLH_Eph(:,2:7);

% save('XX_A_MCR_Eph.mat','XX_A_MCR_Eph');
save('XX_B_MCR_Eph.mat','XX_B_MCR_Eph');
% save('XX_AB_MCLVLH_Eph.mat','XX_AB_MCLVLH_Eph');