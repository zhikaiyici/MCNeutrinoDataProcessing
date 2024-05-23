clear

%% 载入数据
arraySize = 4;
% runID = 0;
NDLName = "ENDF-VIII.0/";
% NDLName = "USE_ONLY_PHOTO_EVAPORATION/";
NDLName = "./";
array = "_" + num2str(arraySize) + "x" + num2str(arraySize);
runCondition = array + '_1e+06_NEUTRINO_Random';
% runCondition = array + '_1e+09_CRY_NEUTRON';
runCondition = array + '_2e+07_CRY';
% runCondition = array + '_1e+06_Li9';
% runCondition = array + '_1e+06_He8';
% runCondition = array + '_2e+09_CRY_MUON';
spectraData_org = [];
for runID = 0:0
    dirName = NDLName + num2str(runID) + runCondition + "/";
    fileName = dirName + 'moduleCalPhPromptLeft' + runCondition + ".data";
    fileName = dirName + 'moduleEdepPrompt' + runCondition + ".data";
    spectraData_org = cat(3, spectraData_org, ReadBinaryFile(fileName, arraySize, 0));
end
spectraDataEdep = spectraData_org .* 19.27 ./ 1000;
spectraDataEdep = spectraData_org;
%%
dscrTh = 0.2;
sliceNum = 10;
nTot = size(spectraDataEdep, 3);
inc = ceil(nTot ./ sliceNum);
nn = 1;
eff_dscr = zeros(sliceNum, 1);
eff_trig = zeros(sliceNum, 1);
eff_et = zeros(sliceNum, 1);
eff_1st = zeros(sliceNum, 1);
eff_2nd = zeros(sliceNum, 1);
eff_rest = zeros(sliceNum, 1);
eff_pattern = zeros(sliceNum, 1);
eff_total = zeros(sliceNum, 1);

for kk = 0:inc:nTot
    if kk + inc < nTot
        spectraData = spectraDataEdep(:,:,kk + 1:kk + inc);
    else
        spectraData = spectraDataEdep(:,:,kk + 1:end);
    end

    logic_dscr = spectraData >= dscrTh;
    spectraData(~logic_dscr) = 0;

    % 总谱
    if arraySize > 1
        totalEdep = sum(sum(spectraData));
    else
        totalEdep = spectraData;
    end

    logic_et0 = totalEdep == 0;
    totalEdep(logic_et0) = [];
    spectraData(:,:,logic_et0) = [];
    logic_dscr = spectraData >= dscrTh;

    % hit统计
    if arraySize > 1
        triggerEvents = sum(sum(logic_dscr));
    else
        triggerEvents = logic_dscr;
    end

    [triggerCounts, triggerNum] = histcounts(triggerEvents, 0.5:arraySize * arraySize + 0.5);
    triggerNum = triggerNum - 0.5;
    triggerNum(triggerNum == 0) = [];

    tempSpec = spectraData;
    tempSpec1 = reshape(permute(tempSpec, [2, 1, 3]), [arraySize, size(tempSpec, 3) * arraySize]);
    tempSpec1 = reshape(tempSpec1, [arraySize * arraySize, numel(tempSpec1) ./ (arraySize * arraySize)])';
    sortedSpec = sort(tempSpec1, 2, 'descend');
    % maxEdep = zeros(size(spectraData, 3), arraySize * arraySize);
    % iEmax = zeros(size(spectraData, 3), arraySize * arraySize);
    % % specData = cell(arraySize * arraySize, 1);
    % for ii = 1:arraySize * arraySize
    %     [maxEdep(:,ii), iEmax(:,ii)] = max(tempSpec, [], [1, 2], 'linear');
    %     tempSpec(tempSpec == permute(maxEdep(:,ii), [2, 3, 1])) = 0;
    %     % temp = maxEdep(:,ii);
    %     % temp(temp < dscrTh) = [];
    %     % h = histogram(temp, 'BinWidth', channelWidth, 'Visible', 'off');
    %     % counts = h.Values;
    %     % energy = (h.BinEdges(1:end - 1) + h.BinEdges(2:end)) ./ 2;
    %     % specData{ii} = [energy, counts];
    % end

    % 甄别效率
    trigminTh = 2;
    trigmaxTh = 16;
    etminTh = 2.5;
    etmaxTh = 7;
    e1minTh = 1;
    e1maxTh = 6;
    e2minTh = 0;
    e2maxTh = 0.52;
    ermaxTh = 0.7;
    eff_dscr(nn, 1) = sum((sum(sum(logic_dscr)) ~= 0)) ...
        ./ (ceil(length(spectraData_org) ./ sliceNum));
    eff_trig(nn, 1) = sum(triggerCounts(trigminTh:trigmaxTh)) ./ sum(triggerCounts);
    logic_trig = triggerEvents >= trigminTh & triggerEvents <= trigmaxTh;
    tempTE = totalEdep;
    tempTE(:,:,~logic_trig) = [];
    % tempiE = iEmax;
    % tempiE = tempiE(:,1:2) - arraySize .* arraySize ...
    %     .* ones(size(tempiE, 1), 1) .* (0:size(tempiE, 1) - 1)';
    % tempiE(~logic_trig,:) = [];
    logic_et = tempTE > etminTh & tempTE <= etmaxTh;
    eff_et(nn, 1) = sum(logic_et) ./ size(tempTE, 3);
    % tempiE(~logic_et,:) = [];
    tempME = sortedSpec;
    tempME(~logic_trig,:) = [];
    tempME(~logic_et,:) = [];
    logic_1st = tempME(:,1) > e1minTh & tempME(:,1) <= e1maxTh;
    eff_1st(nn, 1) = sum(logic_1st) ./ sum(tempME(:,1) > 0);
    tempME(~logic_1st,:) = [];
    % tempiE(~logic_1st,:) = [];
    logic_2nd = tempME(:,2) >= e2minTh & tempME(:,2) <= e2maxTh;
    eff_2nd(nn, 1) = sum(logic_2nd) ./ sum(tempME(:,2) > 0);
    tempME(~logic_2nd,:) = [];
    logic_rest = (sum(tempME, 2) - tempME(:,1) - tempME(:,2)) < ermaxTh;
    eff_rest(nn, 1) = sum(logic_rest) ./ size(logic_rest, 1);
    % tempiE(~logic_2nd,:) = [];
    % [r1st, c1st] = ind2sub(arraySize, tempiE(:,1));
    % [r2nd, c2nd] = ind2sub(arraySize, tempiE(:,2));
    % logic_pattern = ~(abs(r1st - r2nd) > 2 | abs(c1st - c2nd) > 2);
    % eff_pattern(nn, 1) = sum(logic_pattern) ./ size(logic_pattern, 1);
    eff_total(nn, 1) = eff_dscr(nn, 1) .* eff_trig(nn, 1) .* ...
        eff_et(nn, 1) .* eff_1st(nn, 1) .* eff_2nd(nn, 1); % .* eff_pattern(jj, 1);

    nn = nn + 1;
end
%%
syms effTotal(effDiscr, effTrigg, effET, eff1st, eff2nd, effRest);
effTotal(effDiscr, effTrigg, effET, eff1st, eff2nd, effRest) = ...
    effDiscr * effTrigg * effET * eff1st * eff2nd * effRest;
dEffDiscr = diff(effTotal, effDiscr);
dEffTrigg = diff(effTotal, effTrigg);
dEffET = diff(effTotal, effET);
dEff1st = diff(effTotal, eff1st);
dEff2nd = diff(effTotal, eff2nd);
dEffRest = diff(effTotal, effRest);
stdDiscr = std(eff_dscr);
stdTrigg = std(eff_trig);
stdET = std(eff_et);
std1st = std(eff_1st);
std2nd = std(eff_2nd);
stdRest = std(eff_rest);
stdTotal = std(eff_total);
meanDiscr = mean(eff_dscr);
meanTrigg = mean(eff_trig);
meanET = mean(eff_et);
mean1st = mean(eff_1st);
mean2nd = mean(eff_2nd);
meanRest = mean(eff_rest);
meanTotal = mean(eff_total);
stdTotal1 = double(sqrt( ...
    (dEffDiscr(meanDiscr, meanTrigg, meanET, mean1st, mean2nd, meanRest) .* stdDiscr) .^ 2 + ...
    (dEffTrigg(meanDiscr, meanTrigg, meanET, mean1st, mean2nd, meanRest) .* stdTrigg) .^ 2 + ...
    (dEffET(meanDiscr, meanTrigg, meanET, mean1st, mean2nd, meanRest) .* stdET) .^ 2 + ...
    (dEff1st(meanDiscr, meanTrigg, meanET, mean1st, mean2nd, meanRest) .* std1st) .^ 2 + ...
    (dEff2nd(meanDiscr, meanTrigg, meanET, mean1st, mean2nd, meanRest) .* std2nd) .^ 2 + ...
    (dEffRest(meanDiscr, meanTrigg, meanET, mean1st, mean2nd, meanRest) .* stdRest) .^ 2));
