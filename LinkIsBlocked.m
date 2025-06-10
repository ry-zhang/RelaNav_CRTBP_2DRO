function bConsiderSignalBlock = LinkIsBlocked(state1, state2)

    % 考虑信号遮挡（月球和地球）
    % input:
    % state1为航天器1的轨道状态；state2为航天器2的轨道状态
    % output:
    % bConsiderSignalBlock=false表示未被遮挡；=true表示被遮挡

    bConsiderSignalBlock = false;
    const
    % 1.计算地球对星间测距信号的遮挡，参考百度文库《第五章 星间链路及星间组网技术》
	dHp = 1e5/LU;    %距离地球表面100公里，余隙 
    vecEarthToSat1 = state1(1:3)-posEarth;
    vecEarthToSat2 = state2(1:3)-posEarth;
    dRangeSat1 = norm(vecEarthToSat1);
	dRangeSat2 = norm(vecEarthToSat2);
	dMaxAlpha = acos((R_Earth + dHp) / dRangeSat1) + acos((R_Earth + dHp) / dRangeSat2);
	dRealAlpha = acos(dot(vecEarthToSat1, vecEarthToSat2) / norm(vecEarthToSat1) / norm(vecEarthToSat2));
	if dRealAlpha > dMaxAlpha  %连线穿过地球了，大气高度1000km---->100km
        bConsiderSignalBlock = true; 
    end

    % 2.计算月球对星间测距信号的遮挡，参考百度文库《第五章 星间链路及星间组网技术》
    dHp = 0;    %月球没有大气，所以不考虑余隙 
    vecMoonToSat1 = state1(1:3)-posMoon;
    vecMoonToSat2 = state2(1:3)-posMoon;
    dRangeSat1 = norm(vecMoonToSat1);
	dRangeSat2 = norm(vecMoonToSat2);
	dMaxAlpha = acos((R_Moon + dHp) / dRangeSat1) + acos((R_Moon + dHp) / dRangeSat2);
	dRealAlpha = acos(dot(vecMoonToSat1, vecMoonToSat2) / norm(vecMoonToSat1) / norm(vecMoonToSat2));
	if dRealAlpha > dMaxAlpha  %连线穿过地球了，大气高度1000km---->100km
        bConsiderSignalBlock = true; 
    end

end


