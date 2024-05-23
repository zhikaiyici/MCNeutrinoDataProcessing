clear
%%
arraySize = 1;
channelWidth = 0.5; % channalwidth MeV / pheNum
dscrTh = 0.2;
%% 载入数据
NDLName = "./";
NDLName = "ENDF-VIII.0/";
% runCondition = '_4x4_1e+06_Cs137_CENTER';
runConditionCs = '_1x1_1e+06_Cs137_CENTER';
runConditionCo = '_1x1_1e+06_Co60_CENTER';
% fileName = NDLName + '0' + runCondition + '/' + 'moduleEdepPrompt' + runCondition + ".data";
fileNameCs = NDLName + '0' + runConditionCs + '/' + 'moduleCalPhPromptRight' + runConditionCs + ".data";
fileNameCo = NDLName + '1' + runConditionCo + '/' + 'moduleCalPhPromptRight' + runConditionCo + ".data";
spectraDataCs_org = ReadBinaryFile(fileNameCs, arraySize, 1);
spectraDataCo_org = ReadBinaryFile(fileNameCo, arraySize, 1);

spectraDataCs = spectraDataCs_org;
spectraDataCo = spectraDataCo_org;

%%
logic_dscr = spectraDataCs >= dscrTh;
spectraDataCs(~logic_dscr) = 0;

totalEdep = spectraDataCs;
logic_et0 = totalEdep == 0;
totalEdep(logic_et0) = [];

axesCs = axes(figure('Name', fileNameCs));
PlotSpectrum(totalEdep, channelWidth, 'p', axesCs);
legend(axesCs, '{\itE}_{total}', 'FontName', 'Times New Roman','Box','off');
axesCs.XLim = [0,70];
set(axesCs, 'yscale', 'linear');
%%
logic_dscr = spectraDataCo >= dscrTh;
spectraDataCo(~logic_dscr) = 0;

totalEdep = spectraDataCo;
logic_et0 = totalEdep == 0;
totalEdep(logic_et0) = [];

axesCo = axes(figure('Name', fileNameCo));
PlotSpectrum(totalEdep, channelWidth, 'p', axesCo);
legend(axesCo, '{\itE}_{total}', 'FontName', 'Times New Roman','Box','off');
axesCo.XLim = [0,150];
set(axesCo, 'yscale', 'linear');



