clear;

arraySize = 4;
NDLName = "ENDF-VIII.0/";
% NDLName = "./";
array = "_" + num2str(arraySize) + "x" + num2str(arraySize);
runCondition = array + '_1e+09_CRY_NEUTRON';
runCondition = array + "_1e+07_COSMICNEUTRON";
runCondition = array + '_2e+07_CRY';
runCondition = array + '_1e+09_CRY_NEUTRON2.5MeV';

numFile = 11;
eff_capN = zeros(numFile,1);
eff_delayed = zeros(numFile, 1);
eff_prompt = zeros(numFile, 1);
neutronTL = zeros(numFile,1);
neutronTLstd = zeros(numFile,1);
for runID = 0:numFile - 1
    dirName = NDLName + num2str(runID) + runCondition + "/";

    fileName = dirName + 'neutronKE' + runCondition + ".data";
    neutronKE = ReadBinaryFile(fileName, 1) ./ 1000;

    fileName = dirName + 'neutronKEPrimary' + runCondition + ".data";
    neutronKEPrimary = ReadBinaryFile(fileName, 1) .* 1000;
    eCenter = logspace(-11, 6, 180);
    % edges =
    nCounts = histcounts(neutronKEPrimary, eCenter);
    nCounts_ = nCounts .* eCenter(2:end);

    fileName = dirName + 'moduleEdepPrompt' + runCondition + ".data";
    prompt = ReadBinaryFile(fileName, arraySize);
    eff_prompt(runID + 1) = length(prompt) ./ length(neutronKEPrimary);

    fileName = dirName + 'capTimeH' + runCondition + ".data";
    capTimeH = ReadBinaryFile(fileName, 1);
    fileName =  dirName + 'capTimeGd' + runCondition + ".data";
    capTimeGd = ReadBinaryFile(fileName, 1);
    capTime = [capTimeH; capTimeGd];

    eff_capN(runID + 1) = length(capTime) ./ length(neutronKEPrimary);

    fileName = dirName + 'moduleEdepDelayGd' + runCondition + ".data";
    delayedGd = ReadBinaryFile(fileName, arraySize);
    fileName = dirName + 'moduleEdepDelayH' + runCondition + ".data";
    delayedH = ReadBinaryFile(fileName, arraySize);

    delayed = cat(3, delayedGd, delayedH);
    eff_delayed(runID + 1) = length(delayed) ./ length(capTime);

    fileName =  dirName + 'moduleNeutronTrack' + runCondition + ".data";
    moduleNeutronTrackLength = ReadBinaryFile(fileName, arraySize);
    neutronTL(runID + 1) = mean(sum(sum(moduleNeutronTrackLength ./ 10)));
    neutronTLstd(runID + 1) = std(sum(sum(moduleNeutronTrackLength ./ 10)));

    fileName = dirName + 'moduleNumNeutron' + runCondition + ".txt";
    numNeutron = importdata(fileName);
    nNeutron = sum(numNeutron, "all");

    fileName = dirName + 'moduleNumNeutronCorrected' + runCondition + ".txt";
    numNeutronCorrected = importdata(fileName);
    nNeutronCorrected = sum(numNeutronCorrected, "all");

end
meanNeutronCapEff = mean(eff_capN);
stdNeutronCapEff  =std(eff_capN);

meanPromptEff = mean(eff_prompt);
stdPromptEff  =std(eff_prompt);

meanDelayedEff = mean(eff_delayed);
stdDelayedEff  =std(eff_delayed);

meanNeutronTL = mean(neutronTL);
stdNeutronTL = std(neutronTL);

