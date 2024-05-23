clear
%%
arraySize = 4;
channelWidth = 0.02; % channalwidth MeV / pheNum
dscrTh = 0.2;
%% 载入数据
NDLName = "ENDF-VIII.0/";
% NDLName = "USE_ONLY_PHOTO_EVAPORATION/";
NDLName = "./";
array = "_" + num2str(arraySize) + "x" + num2str(arraySize);
runCondition = array + '_1e+06_NEUTRINO_Random';
runCondition = array + '_1e+09_CRY_NEUTRON';
% runCondition = array + '_1e+08_CRY_NEUTRON_WOROOF';
% runCondition = array + '_1e+06_He8';
% runCondition = array + '_1e+06_Li9';
runCondition = array + '_1e+07_CRY';
spectraData_org = [];
for runID = 0:0
    dirName = NDLName + num2str(runID) + runCondition + "/";
    fileName = dirName + 'moduleCalPhPromptLeft' + runCondition + ".data";
    fileName = dirName + 'moduleEdepPrompt' + runCondition + ".data";
    spectraData_org = cat(3, spectraData_org, ReadBinaryFile(fileName, arraySize, 0));
end
spectraData = spectraData_org .* 19.27 ./ 1000;
spectraData = spectraData_org;

%% 能量展宽
% FWHM = @(edep)0.0981 ./ sqrt(edep) + 0.052;
% sigma = spectraData .* FWHM(spectraData) ./ 2.355;
% sigma(isnan(sigma)) = 0;
% parfor ii = 1:numel(spectraData)
%     spectraData(ii) = normrnd(spectraData(ii),sigma(ii));
% end
%%
logic_dscr = spectraData >= dscrTh;
spectraData(~logic_dscr) = 0;
%% 总谱
if arraySize > 1
    totalEdep = sum(sum(spectraData));
else
    totalEdep = spectraData;
end
%%
logic_et0 = totalEdep == 0;
totalEdep(logic_et0) = [];
spectraData(:,:,logic_et0) = [];
logic_dscr = spectraData >= dscrTh;
%% hit统计
if arraySize > 1
    triggerEvents = sum(sum(logic_dscr));
else
    triggerEvents = logic_dscr;
end
figure('Name', 'Trigger_' + fileName);
h = histogram(triggerEvents, 'BinEdges', 0.5:arraySize * arraySize + 0.5, 'Visible', 'on');
triggerCounts = h.Values;
triggerNum = h.BinEdges - 0.5;
triggerNum(triggerNum == 0) = []; 
bar(triggerNum, triggerCounts, 1, 'w');
xlabel("Number of Triggered Modules");
ylabel("Counts");
set(gca, 'yscale', 'log');
set(gca, 'xtick', 0:16);
set(gca, 'fontname', 'times new roman', 'xgrid', 'off', 'ygrid', 'off');

%%
% spectraData(:,:,~logic_trig) = [];
%% 分立谱
subNum = 0;
moduleEdep = cell(arraySize);
moduleTotal = zeros(arraySize);
for ii = 1:arraySize
    for jj = 1:arraySize
        moduleEdep{ii, jj} = spectraData(ii, jj, :);
        moduleTotal(ii, jj) = sum(moduleEdep{ii,jj});
        temp = moduleEdep{ii, jj};
        temp(temp < dscrTh) = [];
    end
end
%% 非按模块统计能量排序
totalAxes = axes(figure('Name', 'Total_' + fileName));
sortAxes = axes(figure('name','1st to 16th by energy' + fileName));
PlotSpectrum(totalEdep, channelWidth, 'e', totalAxes);
hold(totalAxes,'on');
set(totalAxes, 'yscale', 'linear');
% PlotSpectrum(totalEdep, channelWidth, 'e', sortAxes);
hold(sortAxes,'on');
maxEdep = zeros(size(spectraData, 3), arraySize * arraySize);
iEmax = zeros(size(spectraData, 3), arraySize * arraySize);
tempSpec = spectraData;
specData = cell(arraySize * arraySize, 1);
for ii = 1:arraySize * arraySize
    [maxEdep(:,ii), iEmax(:,ii)] = max(tempSpec,[],[1, 2], 'linear');
    tempSpec(tempSpec == permute(maxEdep(:,ii), [2, 3, 1])) = 0;
    temp = maxEdep(:,ii);
    temp(temp < dscrTh) = [];
    if ii == 1
        [c, e] = PlotSpectrum(temp, channelWidth, 'e', totalAxes);
    else
        [c, e] = PlotSpectrum(temp, channelWidth, 'e', totalAxes);
    end
    specData{ii} = [e, c];
end
% set(sortAxes, 'yscale', 'log');
legend(sortAxes, '{\itE}_{2nd}', '{\itE}_{3rd}', '{\itE}_{4th}', ...
    'FontName', 'Times New Roman','Box','off');
hold(totalAxes,'off');
hold(sortAxes,'off');
legend(totalAxes, ['{\itE}_{total}'; '{\itE}_{1st  }'],'FontName', 'Times New Roman','box','off');
totalAxes.XLim = [0,8];
sortAxes.XLim = [0,4];
sortAxes.YScale = 'log';
%% 阵列内最大能量位置统计
maxCount = zeros(arraySize);
for ii = 1:length(spectraData)
    [r, c] = find(spectraData(:,:,ii) == maxEdep(ii, 1));
    maxCount(r, c) = maxCount(r, c) + 1;
end
%% 按模块统计能量排序
axesT = axes(figure('name','Total and 1st by module'  + fileName));
PlotSpectrum(totalEdep, channelWidth, 'e', axesT);
hold(axesT, 'on');
axesS = axes(figure('name','2nd to 16th by module' + fileName));
hold(axesS, 'on');
tempCount = maxCount;
for ii = 1:arraySize * arraySize
    [r, c] = find(tempCount == max(max(tempCount)));
    temp = moduleEdep{r, c};
    temp(temp < dscrTh) = [];
    if ii == 1
        PlotSpectrum(temp, channelWidth, 'e', axesT);
    else
        PlotSpectrum(temp, channelWidth, 'e', axesS);
    end
    tempCount(r, c) = 0;
end
% set(axesS, 'yscale', 'log');
legend(axesS, ['E_{1st}';'E_{2nd}'; 'E_{3rd}'; 'E_{4th}';'E_{5th}';'E_{6th}'], ...
    'FontName', 'Times New Roman','Box','off');
hold(axesT,'off');
hold(axesS,'off');
legend(axesT, ['E_{total}'; 'E_{1st  }'],'FontName', 'Times New Roman','box','off');
set(axesT,'xlim',[0,8]);
set(axesS,'xlim',[0,8]);

%% 甄别效率
trigminTh = 2;
trigmaxTh = 16;
etminTh = 2.5;
etmaxTh = 7;
e1minTh = 1;
e1maxTh = 6.0;
e2minTh = 0;
e2maxTh = 0.52;
eff_dscr = sum(~(sum(sum(logic_dscr))==0)) ./ (length(spectraData_org));
eff_trig = sum(triggerCounts(2:16)) ./ sum(triggerCounts);
logic_trig = triggerEvents >= trigminTh & triggerEvents <= trigmaxTh;
tempTE = totalEdep;
tempTE(:,:,~logic_trig) = [];
tempiE = iEmax;
tempiE = tempiE(:,1:2) - arraySize .* arraySize .* ...
    ones(size(tempiE, 1), 1) .* (0:size(tempiE, 1) - 1)';
tempiE(~logic_trig,:) = [];
logic_et = tempTE > etminTh & tempTE <= etmaxTh;
eff_et = sum(logic_et) ./ size(tempTE, 3);
tempiE(~logic_et,:) = [];
tempME = maxEdep;
tempME(~logic_trig,:) = [];
tempME(~logic_et,:) = [];
logic_1st = tempME(:,1) > e1minTh & tempME(:,1) <= e1maxTh;
eff_1st = sum(logic_1st) ./ sum(tempME(:,1) > 0);
tempME(~logic_1st,:) = [];
tempiE(~logic_1st,:) = [];
logic_2nd = tempME(:,2) >= e2minTh & tempME(:,2) <= e2maxTh;
eff_2nd = sum(logic_2nd) ./ sum(tempME(:,2) > 0);
tempiE(~logic_2nd,:) = [];
[r1st, c1st] = ind2sub(arraySize, tempiE(:,1));
[r2nd, c2nd] = ind2sub(arraySize, tempiE(:,2));
logic_pattern = abs(r1st - r2nd) <= 2 & abs(c1st - c2nd) <= 2;
% logic_pattern = (abs(r1st - r2nd) + abs(c1st - c2nd)) <= 1;
eff_pattern = sum(logic_pattern) ./ size(logic_pattern, 1);
eff_total = eff_dscr .* eff_trig .* eff_et .* eff_1st .* eff_2nd; % .* eff_pattern;
