function data = read_data_from_file(filename)
    % 从文本文件中读取数据
    % 假设文件格式为：儒略日时间  轨道位置 轨道速度
    fileID = fopen(filename, 'r');
    data = [];
        line = fgetl(fileID);
        while ischar(line)
            parsedData = sscanf(line,'%f %f %f %f %f %f %f', 7)';
            data = [data; parsedData];
            line = fgetl(fileID);  % 读取下一行
        end
    fclose(fileID);
end
