clear;

NDLName = "ENDF-VIII.0/";
% NDLName = "./";
runCondition = '_4x4_1e+07_NEUTRINO_Random';
runCondition = '_4x4_1e+07_Li9';
runCondition = '_4x4_2e+09_CRY';
runID = 0;
dirName = NDLName + num2str(runID) + runCondition + "/";
fileName = dirName + 'neutronKE' + runCondition + ".data";
nKE = ReadBinaryFile(fileName, 4) ./ 1000;
fileName =  dirName + 'neutronGenicTime' + runCondition + ".data";
nGT = ReadBinaryFile(fileName, 4) ./ 1000;
% nKE(nKE == 0) = [];
fileName = dirName + 'numNeutron' + runCondition + ".txt";
numNeutron = importdata(fileName);
nNeutron = sum(numNeutron, "all");

% figure;
% histogram(nKE);
% set(gca, 'YScale', 'log');