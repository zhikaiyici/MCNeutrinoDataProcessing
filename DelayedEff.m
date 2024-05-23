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
% runCondition = array + '_2e+07_CRY';
spectraData_org_gd = [];
spectraData_org_h = [];
for runID = 0:0
    dirName = NDLName + num2str(runID) + runCondition + "/";
    fileName = dirName + 'moduleCalPhDelayGdLeft' + runCondition + ".data";
    fileName = dirName + 'moduleEdepDelayGd' + runCondition + ".data";
    spectraData_org_gd = cat(3, spectraData_org_gd, ReadBinaryFile(fileName, arraySize, 0));
    fileName = dirName + 'moduleCalPhDelayHLeft' + runCondition + ".data";
    fileName = dirName + 'moduleEdepDelayH' + runCondition + ".data";
    spectraData_org_h = cat(3, spectraData_org_h, ReadBinaryFile(fileName, arraySize, 0));
end
spectraDataEdepGd = spectraData_org_gd .* 19.27 ./ 1000;
spectraDataEdepH = spectraData_org_h .* 19.27 ./ 1000;
spectraDataEdepGd = spectraData_org_gd;
spectraDataEdepH = spectraData_org_h;

spectraData_org = cat(3, spectraData_org_gd, spectraData_org_h);
spectraDataEdep = spectraData_org .* 19.27 ./ 1000;
spectraDataEdep = spectraData_org;

%%
dscrTh = 0.2;
sliceNum = 10;
nTotGd = size(spectraDataEdepGd, 3);
incGd = ceil(nTotGd ./ sliceNum);

nTotH = size(spectraDataEdepH, 3);
incH = ceil(nTotH ./ sliceNum);

nn = 1;
eff_dscr = zeros(sliceNum, 1);
eff_trig = zeros(sliceNum, 1);
eff_et = zeros(sliceNum, 1);
eff_1st = zeros(sliceNum, 1);
eff_2nd = zeros(sliceNum, 1);
eff_pattern = zeros(sliceNum, 1);
eff_total = zeros(sliceNum, 1);

for kk = 1:sliceNum
    if kk * incGd < nTotGd
        spectraDataGd = spectraDataEdepGd(:,:,1 + (kk - 1) * incGd:kk * incGd);
    else
        spectraDataGd = spectraDataEdepGd(:,:,1 + (kk - 1) * incGd:end);
    end
    if kk * incH < nTotH
        spectraDataH = spectraDataEdepH(:,:,1 + (kk - 1) * incH:kk * incH);
    else
        spectraDataH = spectraDataEdepH(:,:,1 + (kk - 1) * incH:end);
    end
    spectraData = cat(3, spectraDataGd, spectraDataH);

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
    etminTh = 2.8;
    etmaxTh = 8;
    e1minTh = 0.5;
    e1maxTh = 6;
    e2minTh = 0;
    e2maxTh = 3;
    e3minTh = 0;
    e3maxTh = 2;
    e4minTh = 0;
    e4maxTh = 1;
    eff_dscr(nn, 1) = sum((sum(sum(logic_dscr)) ~= 0)) ...
        ./ (ceil(length(spectraData_org) ./ sliceNum));
    eff_trig(nn, 1) = sum(triggerCounts(trigminTh:trigmaxTh)) ./ sum(triggerCounts);
    logic_trig = triggerEvents >= trigminTh & triggerEvents <= trigmaxTh;
    tempTE = totalEdep;
    tempTE(:,:,~logic_trig) = [];
    logic_et = tempTE > etminTh & tempTE <= etmaxTh;
    eff_et(nn, 1) = sum(logic_et) ./ size(tempTE, 3);
    tempME = sortedSpec;
    tempME(~logic_trig, :) = [];
    tempME(~logic_et, :) = [];
    logic_1st = tempME(:,1) > e1minTh & tempME(:,1) <= e1maxTh;
    eff_1st(nn, 1) = sum(logic_1st) ./ sum(tempME(:,1) > 0);
    tempME(~logic_1st, :) = [];
    logic_2nd = tempME(:,2) >= e2minTh & tempME(:,2) <= e2maxTh;
    eff_2nd(nn, 1) = sum(logic_2nd) ./ sum(tempME(:,2) > 0);
    tempME(~logic_2nd,:) = [];
    logic_3rd = tempME(:,3) >= e3minTh & tempME(:,3) <= e3maxTh;
    eff_3rd = sum(logic_3rd) ./ sum(tempME(:,3) >= 0);
    tempME(~logic_3rd,:) = [];
    logic_4th = tempME(:,4) >= e4minTh & tempME(:,4) <= e4maxTh;
    eff_4th = sum(logic_4th) ./ sum(tempME(:,4) >= 0);
    eff_total(nn, 1) = eff_dscr(nn, 1) .* eff_trig(nn, 1) .* ...
        eff_et(nn, 1) .* eff_1st(nn, 1) .* eff_2nd(nn, 1); % .* eff_pattern(jj, 1);

    nn = nn + 1;
end
%%
syms effTotal(effDiscr, effTrigg, effET, eff1st, eff2nd, eff3rd, eff4th);
effTotal(effDiscr, effTrigg, effET, eff1st, eff2nd, eff3rd, eff4th) = ...
    effDiscr * effTrigg * effET * eff1st * eff2nd * eff3rd * eff4th;
dEffDiscr = diff(effTotal, effDiscr);
dEffTrigg = diff(effTotal, effTrigg);
dEffET = diff(effTotal, effET);
dEff1st = diff(effTotal, eff1st);
dEff2nd = diff(effTotal, eff2nd);
dEff3rd = diff(effTotal, eff3rd);
dEff4th = diff(effTotal, eff4th);
stdDiscr = std(eff_dscr);
stdTrigg = std(eff_trig);
stdET = std(eff_et);
std1st = std(eff_1st);
std2nd = std(eff_2nd);
std3rd = std(eff_3rd);
std4th = std(eff_4th);
stdTotal = std(eff_total);
meanDiscr = mean(eff_dscr);
meanTrigg = mean(eff_trig);
meanET = mean(eff_et);
mean1st = mean(eff_1st);
mean2nd = mean(eff_2nd);
mean3rd = mean(eff_3rd);
mean4th = mean(eff_4th);
meanTotal = mean(eff_total);
stdTotal1 = double(sqrt( ...
    (dEffDiscr(meanDiscr, meanTrigg, meanET, mean1st, mean2nd, mean3rd, mean4th) .* stdDiscr) .^2 + ...
    (dEffTrigg(meanDiscr, meanTrigg, meanET, mean1st, mean2nd, mean3rd, mean4th) .* stdTrigg) .^2 + ...
    (dEffET(meanDiscr, meanTrigg, meanET, mean1st, mean2nd, mean3rd, mean4th) .* stdET) .^2 + ...
    (dEff1st(meanDiscr, meanTrigg, meanET, mean1st, mean2nd, mean3rd, mean4th) .* std1st) .^2 + ...
    (dEff2nd(meanDiscr, meanTrigg, meanET, mean1st, mean2nd, mean3rd, mean4th) .* std2nd) .^2 + ...
    (dEff3rd(meanDiscr, meanTrigg, meanET, mean1st, mean2nd, mean3rd, mean4th) .* std3rd) .^2 + ...
    (dEff4th(meanDiscr, meanTrigg, meanET, mean1st, mean2nd, mean3rd, mean4th) .* std4th) .^2));
