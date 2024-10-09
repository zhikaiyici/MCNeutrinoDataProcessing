clear

%% 载入数据
numFile = 11;
eff_delayed = zeros(numFile, 1);
eff_capN = zeros(numFile, 1);
arraySize = 4;
NDLName = "ENDF-VIII.0/";
% NDLName = "./";
% NDLName = "USE_ONLY_PHOTO_EVAPORATION/";
array = "_" + num2str(arraySize) + "x" + num2str(arraySize);
runCondition = array + '_1e+06_NEUTRINO_Random';
runCondition = array + '_1e+09_CRY_NEUTRON';
% runCondition = array + '_1e+06_CRY';
% runCondition = array + '_1e+06_Li9';
% runCondition = array + '_1e+06_He8';
for runID = 0:numFile - 1
    dirName = NDLName + num2str(runID) + runCondition + "/";

    fileName = dirName + 'moduleCapTimeGd' + runCondition + ".data";
    capTimeGd = ReadBinaryFile(fileName, arraySize);
    fileName = dirName + 'moduleCapTimeH' + runCondition + ".data";
    capTimeH = ReadBinaryFile(fileName, arraySize);

    fileName = dirName + 'moduleNeutronTrack' + runCondition + ".data";
    neutronTrack = ReadBinaryFile(fileName, arraySize);

    capTime = cat(3, capTimeGd, capTimeH);
    capTime(capTime == 0) = [];
    eff_capN(runID + 1) = length(capTime) ./ length(neutronTrack);

    fileName = dirName + 'moduleEdepDelayGd' + runCondition + ".data";
    delayedGd = ReadBinaryFile(fileName, arraySize);
    fileName = dirName + 'moduleEdepDelayH' + runCondition + ".data";
    delayedH = ReadBinaryFile(fileName, arraySize);

    delayed = cat(3, delayedGd, delayedH);
    eff_delayed(runID + 1) = length(delayed) ./ length(capTime);

end

meanCapN = mean(eff_capN);
stdCapN = std(eff_capN);

meanDelayed = mean(eff_delayed);
stdDelayed = std(eff_delayed);

