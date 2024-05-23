clear
%%
arraySize = 4;
channelWidth = 0.02; % channalwidth MeV / pheNum
dscrTh = 0.2;
%% 载入数据
NDLName = "ENDF-VIII.0/";
% NDLName = "USE_ONLY_PHOTO_EVAPORATION/";
% NDLName = "./";

array = "_" + num2str(arraySize) + "x" + num2str(arraySize);
runCondition = array + '_1e+06_NEUTRINO_Random';
% runCondition = array + '1e+06_neutron_1 MeV_INSIDE';
% runCondition = array + '_1e+08_CRY_NEUTRON';
% runCondition = array + '_1e+06_He8';
% runCondition = array + '_1e+06_Li9';

spectraData_org_gd = [];
spectraData_org_h = [];
for runID = 0:0
    dirName = NDLName + num2str(runID) + runCondition + "/";
    fileName = dirName + 'moduleCalPhDelayGdLeft' + runCondition + ".data";
    % fileName = dirName + 'moduleEdepDelayGd' + runCondition + ".data";
    spectraData_org_gd = cat(3, spectraData_org_gd, ReadBinaryFile(fileName, arraySize, 1));
    fileName = dirName + 'moduleCalPhDelayHLeft' + runCondition + ".data";
    % fileName = dirName + 'moduleEdepDelayH' + runCondition + ".data";
    spectraData_org_h = cat(3, spectraData_org_h, ReadBinaryFile(fileName, arraySize, 1));
end
spectraData_org = cat(3, spectraData_org_gd, spectraData_org_h);

spectraData = spectraData_org .* 19.27 ./ 1000;
% spectraData = spectraData_org;

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
% tCap(~logic_dscr) = 0;
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
%         moduleEdep{ii, jj} = temp;
%         subNum = subNum + 1;
%         subplot(arraySize,arraySize, subNum);
%         figure('Name', [num2str(ii),', ', num2str(jj)]);
%         PlotSpectrum(temp, channelWidth, 'e');
%         set(gca, 'yscale', 'log');
%         hold on;
    end
end
%% 非按模块统计能量排序
axesE = axes(figure('Name', fileName));
PlotSpectrum(totalEdep, channelWidth, 'e', axesE);
hold(axesE,'on');
% set(axesE, 'yscale', 'log');
maxEdep = zeros(size(spectraData, 3), arraySize * arraySize);
tempSpec = spectraData;
specData = cell(arraySize * arraySize, 1);
for ii = 1:arraySize * arraySize
    maxEdep(:,ii) = max(tempSpec,[],[1, 2]);
    tempSpec(tempSpec == permute(maxEdep(:,ii), [2, 3, 1])) = 0;
    temp = maxEdep(:,ii);
    temp(temp < dscrTh) = [];
    [c, e] = PlotSpectrum(temp, channelWidth, 'e', axesE);
    specData{ii} = [e, c];
end
legend(axesE, '{\itE}_{total}', '{\itE}_{1st}', '{\itE}_{2nd}', '{\itE}_{3rd}', '{\itE}_{4th}', ...
    'FontName', 'Times New Roman','Box','off');
hold(axesE,'off');
axesE.XLim = [0,10];
set(axesE, 'yscale', 'log');

%% 阵列内最大能量位置统计
maxCount = zeros(arraySize);
for ii = 1:length(spectraData)
    [r, c] = find(spectraData(:,:,ii) == maxEdep(ii, 1));
    maxCount(r, c) = maxCount(r, c) + 1;
end
%% 按模块统计能量排序
axesM = axes(figure('name','Total and all module'));
% PlotSpectrum(totalEdep, channelWidth, 'e', axesM);
hold(axesM, 'on');
tempCount = maxCount;
for ii = 1:arraySize * arraySize
    [r, c] = find(tempCount == max(max(tempCount)));
    temp = moduleEdep{r, c};
    temp(temp < dscrTh) = [];
    PlotSpectrum(temp, channelWidth, 'e', axesM);
    tempCount(r, c) = 0;
end
% set(axesM, 'yscale', 'log');
legend(axesM, 'E_{1st}', 'E_{2nd}', 'E_{3rd}', 'E_{4th}',...
    'E_{5th}','E_{6th}', 'E_{7th}','E_{8th}','E_{9th}','E_{10th}',...
    'E_{11st}', 'E_{12nd}','E_{13rd}', 'E_{14th}','E_{15th}','E_{16th}',...
    'FontName', 'Times New Roman','Box','off');
hold(axesM,'off');
set(axesM,'xlim',[0,10]);
set(axesM, 'yscale', 'log');

%%
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
eff_dscr = sum(~(sum(sum(logic_dscr))==0)) ./ (length(spectraData_org));
eff_trig = sum(triggerCounts(2:16)) ./ sum(triggerCounts);
logic_trig = triggerEvents >= trigminTh & triggerEvents <= trigmaxTh;
tempTE = permute(totalEdep, [3, 1, 2]);
tempTE(~logic_trig,:) = [];
logic_et = tempTE > etminTh & tempTE <= etmaxTh;
eff_et = sum(logic_et) ./ size(tempTE, 1);
tempME = maxEdep;
tempME(~logic_trig,:) = [];
tempME(~logic_et,:) = [];
tempTE(~logic_et,:) = [];
e3 = tempME(:,3) ./ tempTE;
e2 = tempME(:,2) ./ tempTE;
e1 = tempME(:,1) ./ tempTE;
logic_ratio = e3 >= (e1 - 0.5) ./ 5;
eff_ratio = sum(logic_ratio) ./ size(tempTE, 1);
eff_PANDA = eff_dscr .* eff_trig .* eff_et .* eff_ratio;
logic_1st = tempME(:,1) > e1minTh & tempME(:,1) <= e1maxTh;
eff_1st = sum(logic_1st) ./ sum(tempME(:,1) > 0);
tempME(~logic_1st,:) = [];
logic_2nd = tempME(:,2) >= e2minTh & tempME(:,2) <= e2maxTh;
eff_2nd = sum(logic_2nd) ./ sum(tempME(:,2) > 0);
tempME(~logic_2nd,:) = [];
logic_3rd = tempME(:,3) >= e3minTh & tempME(:,3) <= e3maxTh;
eff_3rd = sum(logic_3rd) ./ sum(tempME(:,3) >= 0);
tempME(~logic_3rd,:) = [];
logic_4th = tempME(:,4) >= e4minTh & tempME(:,4) <= e4maxTh;
eff_4th = sum(logic_4th) ./ sum(tempME(:,4) >= 0);
eff_total = eff_dscr .* eff_trig .* eff_et .* eff_1st .* eff_2nd .* eff_3rd .* eff_4th;
