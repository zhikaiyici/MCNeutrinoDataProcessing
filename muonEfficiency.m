clear
%% 载入数据
arraySize = 4;
NDLName = "ENDF-VIII.0/";
% NDLName = "./";
array = "_" + num2str(arraySize) + "x" + num2str(arraySize);
% NDLName = "USE_ONLY_PHOTO_EVAPORATION/";
runCondition = array + '_2e+09_CRY_MUON';
spectraData_org = [];
for runID = 0:0
    dirName = NDLName + num2str(runID) + runCondition + "/";
    fileName = dirName + 'moduleMuEdep' + runCondition + ".data";
    spectraData_org = cat(3, spectraData_org, ReadBinaryFile(fileName, arraySize, 0));
end

spectraData = spectraData_org;

%%
dscrTh = 0.2;
sliceNum = 100;
nTot = size(spectraData, 3);
inc = ceil(nTot ./ sliceNum);
nE = zeros(sliceNum, 1);
nMu = zeros(sliceNum, 1);
e = zeros(sliceNum, 1);
jj = 1;
for ii = 0:inc:nTot
    if ii + inc < nTot
        data = spectraData(:,:,ii + 1:ii + inc);
    else
        data = spectraData(:,:,ii + 1:end);
    end
    nE(jj, 1)  = size(data, 3);

    logic_dscr = data >= dscrTh;
    data(~logic_dscr) = 0;
    hit = sum(sum(logic_dscr));
    hit = permute(hit, [3,1,2]);
    droped = hit == 0;
    hit(droped,:) = [];
    data(:,:,droped) = [];
    edep = sum(sum(data));
    edep = permute(edep, [3,1,2]);

    logicMu = hit > 1 & edep > 10; 
    
    nMu(jj, 1) = sum(logicMu);

    e(jj, 1) = sum(logicMu) ./ size(data, 3);
    jj = jj + 1;
end
%%
effMu = sum(nMu(1:sliceNum)) ./ sum(nE(1:sliceNum));
stdnMu = std(nMu(1:sliceNum));
stdEffMu = stdnMu .* sliceNum ./ sum(nE(1:sliceNum));
%%
nLost = nE - nMu;
effLost = sum(nLost(1:sliceNum)) ./ sum(nE(1:sliceNum));
stdLost = std(nLost(1:sliceNum));
stdEffLost = stdLost .* sliceNum ./ sum(nE(1:sliceNum));

