clear

%% 载入数据
arraySize = 4;
% runID = 0;
NDLName = "ENDF-VIII.0/";
% NDLName = "./";
% NDLName = "USE_ONLY_PHOTO_EVAPORATION/";
array = "_" + num2str(arraySize) + "x" + num2str(arraySize);
runCondition = array + '_1e+06_NEUTRINO_Random';
% runCondition = array + '_1e+08_CRY_NEUTRON';
% runCondition = array + '_1e+06_CRY';
% runCondition = array + '_1e+06_Li9';
% runCondition = array + '_1e+06_He8';
capTimeGdOrg = [];
capTimeHOrg = [];
for runID = 0:0
    dirName = NDLName + num2str(runID) + runCondition + "/";
    fileName = dirName + 'moduleCapTimeGd' + runCondition + ".data";
    capTimeGdOrg = cat(3, capTimeGdOrg, ReadBinaryFile(fileName, arraySize));
    fileName = dirName + 'moduleCapTimeH' + runCondition + ".data";
    capTimeHOrg = cat(3, capTimeHOrg, ReadBinaryFile(fileName, arraySize));
end
capTime_org = cat(3, capTimeGdOrg, capTimeHOrg);

sliceNum = 100;
nTotGd = size(capTimeGdOrg, 3);
incGd = ceil(nTotGd ./ sliceNum);

nTotH = size(capTimeHOrg, 3);
incH = ceil(nTotH ./ sliceNum);

nn = 1;

eff_capT = zeros(sliceNum, 1);
eff_capT300 = zeros(sliceNum, 1);

for kk = 1:sliceNum
    if kk * incGd < nTotGd
        capTimeGd = capTimeGdOrg(:,:,1 + (kk - 1) * incGd:kk * incGd);
    else
        capTimeGd = capTimeGdOrg(:,:,1 + (kk - 1) * incGd:end);
    end
    if kk * incH < nTotH
        capTimeH = capTimeHOrg(:,:,1 + (kk - 1) * incH:kk * incH);
    else
        capTimeH = capTimeHOrg(:,:,1 + (kk - 1) * incH:end);
    end
    capTime = cat(3, capTimeGd, capTimeH);
    capTime(capTime == 0) = [];
    eff_capT(nn) = sum(capTime >= 8 & capTime <= 300) ./ sum(capTime > 0);
    eff_capT300(nn) = sum(capTime > 300) ./ sum(capTime > 0);
    nn = nn + 1;
end
meanCapT = mean(eff_capT);
stdCapT = std(eff_capT);

meanCapT300 = mean(eff_capT300);
stdCapT300 = std(eff_capT300);
