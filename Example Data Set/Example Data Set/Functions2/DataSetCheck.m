function DataSetTimings = DataSetCheck(fuel)


DataSetTimingsDiesel = [9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20]; % Dataset timings inside Diesel Folder
DataSetTimingsHVO50 = [14, 16, 18, 19, 20];                         % Dataset timings inside HVO50 folder
DataSetTimingsOther = [14, 15, 16, 17, 18, 19, 20];         % Dataset timings in other folders

if strcmp(fuel, 'Diesel100')
    DataSetTimings = DataSetTimingsDiesel;
elseif strcmp(fuel, 'HVO50')
    DataSetTimings = DataSetTimingsHVO50
else
   DataSetTimings = DataSetTimingsOther;
end
end

