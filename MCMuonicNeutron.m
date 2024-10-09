clear;

arraySize = 4;
NDLName = "ENDF-VIII.0/";
NDLName = "./";
array = "_" + num2str(arraySize) + "x" + num2str(arraySize);
runCondition = array + '_2e+09_CRY_MUON';
runCondition = array + '_1e+08_CRY';
% runCondition = array + '_1e+09_MUON';
% runCondition = array + '_2e+09_CRY_WOROOF';

numFile = 11;
neutronCapEff = zeros(numFile,1);
neutronYield = zeros(numFile,1);
neutronYieldCorrected = zeros(numFile,1);
for runID = 0:numFile - 1
    dirName = NDLName + num2str(runID) + runCondition + "/";

    % fileName = dirName + 'neutronKineticEnergy' + runCondition + ".data";
    % neutronKineticEnergy = ReadBinaryFile(fileName, arraySize) ./ 1000;
    % fileName =  dirName + 'neutronGenicTime' + runCondition + ".data";
    % neutronGenicTime = ReadBinaryFile(fileName, arraySize) ./ 1000;

    fileName =  dirName + 'moduleMuTrackLength' + runCondition + ".data";
    moduleMuTrackLength = ReadBinaryFile(fileName, arraySize);

    fileName = dirName + 'neutronKE' + runCondition + ".data";
    neutronKE = ReadBinaryFile(fileName, 1) ./ 1000;
    fileName =  dirName + 'neutronGT' + runCondition + ".data";
    neutronGT = ReadBinaryFile(fileName, 1) .* 1000;

    fileName = dirName + 'capTimeH' + runCondition + ".data";
    capTimeH = ReadBinaryFile(fileName, 1);
    fileName =  dirName + 'capTimeGd' + runCondition + ".data";
    capTimeGd = ReadBinaryFile(fileName, 1);
    capTime = [capTimeH; capTimeGd];

    % nKE(nKE == 0) = [];
    fileName = dirName + 'moduleNumNeutron' + runCondition + ".txt";
    numNeutron = importdata(fileName);
    nNeutron = sum(numNeutron, "all");

    fileName = dirName + 'moduleNumNeutronCorrected' + runCondition + ".txt";
    numNeutronCorrected = importdata(fileName);
    nNeutronCorrected = sum(numNeutronCorrected, "all");

    neutronYield = nNeutron ./ (sum(moduleMuTrackLength, "all") ./ 10);
    neutronYieldCorrected = nNeutronCorrected ./ (sum(moduleMuTrackLength, "all") ./ 10);

    neutronCapEff(runID + 1) = length(capTime) ./ nNeutronCorrected;
end
meanNeutronCapEff = mean(neutronCapEff);
stdNeutronCapEff  =std(neutronCapEff);
