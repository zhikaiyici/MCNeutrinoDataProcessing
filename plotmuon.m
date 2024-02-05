clear
%% 载入数据
fileName = 'moduleEdepPrompt(4x4_1e+08_CRY).txt';
trackFileName = 'moduleMuTrackLength(4x4_1e+08_CRY).txt';
edepData_org = load(fileName);
trackData_org = load(trackFileName);
arraySize = 4;
channelWidth = 0.5; % channalwidth MeV / pheNum
dscrTh = 0.2;
edepData = ReshapeDataMatrix(arraySize, edepData_org);
trackData = ReshapeDataMatrix(arraySize, trackData_org);
spectraData = trackData;
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
%% hit统计
if arraySize > 1
    triggerEvents = sum(sum(logic_dscr));
else
    triggerEvents = logic_dscr;
end
figure('Name', ['Trigger_', fileName]);
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
trigTh = 16;
logic_trig = triggerEvents >= 2 & triggerEvents <= trigTh;
%%
spectraData(:,:,~logic_trig) = [];
%% 总谱
if arraySize > 1
    totalEdep = sum(sum(spectraData));
else
    totalEdep = spectraData;
end
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
set(axesE, 'yscale', 'log');
maxEdep = cell(arraySize);
tempSpec = spectraData;
for ii = 1:4
    maxEdep{ii} = max(tempSpec,[],[1, 2]);
    tempSpec(tempSpec == maxEdep{ii}) = 0;
    temp = maxEdep{ii};
    temp(temp < dscrTh) = [];
    PlotSpectrum(temp, channelWidth, 'e', axesE);
end
legend(axesE, 'E_{total}', 'E_{1st}', 'E_{2nd}', 'E_{3rd}', 'E_{4th}', ...
    'FontName', 'Times New Roman','Box','off');
hold(axesE,'off');
axesE.XLim = [0,600];

%% 阵列内最大能量位置统计
maxCount = zeros(arraySize);
for ii = 1:length(spectraData)
    [r, c] = find(spectraData(:,:,ii) == maxEdep{1}(ii));
    maxCount(r, c) = maxCount(r, c) + 1;
end
%% 按模块统计能量排序
axesM = axes(figure('name',fileName));
PlotSpectrum(totalEdep, channelWidth, 'e', axesM);
hold(axesM, 'on');
tempCount = maxCount;
for ii = 1:16
    [r, c] = find(tempCount == max(max(tempCount)));
    temp = moduleEdep{r, c};
    temp(temp < dscrTh) = [];
    PlotSpectrum(temp, channelWidth, 'e', axesM);
    tempCount(r, c) = 0;
end
set(axesM, 'yscale', 'log');
legend(axesM,'E_{total}',  'E_{1st}', 'E_{2nd}', 'E_{3rd}', 'E_{4th}',...
    'E_{5th}','E_{6th}', 'E_{7th}','E_{8th}','E_{9th}','E_{10th}',...
    'E_{11st}', 'E_{12nd}','E_{13rd}', 'E_{14th}','E_{15th}','E_{16th}',...
    'FontName', 'Times New Roman','Box','off');
hold(axesM,'off');
set(axesM,'xlim',[0,600]);

%%
eff_dscr = sum(~(sum(sum(logic_dscr))==0)) ./ (length(edepData_org) ./ arraySize);
eff_trig = sum(triggerCounts(2:trigTh)) ./ sum(triggerCounts);
eff_et = sum(totalEdep > 1 & totalEdep <= 10) ./ length(totalEdep);
eff_1st = sum(maxEdep{1} > 0.5 & maxEdep{1} <= 8) ./ sum(maxEdep{1} > 0);
eff_2nd = sum(maxEdep{2} > 0 &maxEdep{2} <= 3.5) ./ sum(maxEdep{2} > 0);
eff_3rd = sum(maxEdep{3} > 0 &maxEdep{3} <= 2) ./ sum(maxEdep{3} > 0);
eff_4th = sum(maxEdep{4} > 0 &maxEdep{4} <= 1) ./ sum(maxEdep{4} > 0);
eff_t = eff_dscr .* eff_trig .* eff_et .* eff_1st .* eff_2nd .* eff_3rd .* eff_4th;