clear
%%
arraySize = 1;
arraySize = 4;
% runID = 0;
NDLName = "ENDF-VIII.0/";
% NDLName = "USE_ONLY_PHOTO_EVAPORATION/";
array = "_" + num2str(arraySize) + "x" + num2str(arraySize);
runCondition = array + '_1e+06_NEUTRINO_Random';
% runCondition = array + '_1e+08_CRY_NEUTRON';
% runCondition = array + '_1e+06_He8';
% runCondition = array + '_1e+06_Li9';
capTimeGdOrg = [];
capTimeHOrg = [];
for runID = 0:0
    dirName = NDLName + num2str(runID) + runCondition + "/";
    fileName = dirName + 'moduleCapTimeGd' + runCondition + ".data";
    capTimeGdOrg = cat(3, capTimeGdOrg, ReadBinaryFile(fileName, arraySize));
    fileName = dirName + 'moduleCapTimeH' + runCondition + ".data";
    capTimeHOrg = cat(3, capTimeHOrg, ReadBinaryFile(fileName, arraySize));
end
capTime_org = cat(3, capTimeGdOrg, capTimeHOrg);

t = capTime_org;
% %%
% t = [th;tgd];
% t = ReshapeDataMatrix(arraySize, t);
% %%
% t(~logic_dscr) = 0;
% t(:,:,~logic_trig) = [];
%%
capTimeGdOrg(capTimeGdOrg == 0) = [];
capTimeHOrg(capTimeHOrg == 0) = [];
figure('Name','CapTime');
hgd = histogram(capTimeGdOrg,'BinWidth', 1, 'DisplayStyle', 'stairs');
hold on;
hh = histogram(capTimeHOrg,'BinWidth', 1, 'DisplayStyle', 'stairs');

xlabel("{\it\DeltaT} (\mus)");    
ylabel("Events / \mus");
% ylabel('\fontname{宋体}计数');
legend('Captured by Gd','Captured by H','Box','off');
% legend('\fontname{宋体}被\fontname{Times new roman}Gd\fontname{宋体}俘获',...
%     '\fontname{宋体}被\fontname{Times new roman}H\fontname{宋体}俘获','Box','off')
set(gca, 'xlim',[0, 400]);
set(gca, 'fontname', 'times new roman', 'xgrid', 'off', 'ygrid', 'off');

%%
cth = hh.Values;
tth = hh.BinEdges - 0.5;
tth(tth < 0) = [];
capTimeH = [tth',cth'];

modelFun = @(beta0, x)(beta0(1) / beta0(2) * exp(-x / beta0(2)) + beta0(3));
start = [3e4; 10; 1e2];
nlm = fitnlm(tth, cth, modelFun, start);

ctgd = hgd.Values;
ttgd = hgd.BinEdges - 0.5;
ttgd(ttgd < 0) = [];

nlm2 = fitnlm(ttgd, ctgd, modelFun, start);

capTimeGd = [ttgd',ctgd'];
temp = zeros(1,abs(length(cth) - length(ctgd)));
counts = [cth,temp]' + ctgd';
counts = counts ./ sum(counts);
cdft = zeros(length(counts),1);
cdft(1) = counts(1);
for ii = 2:length(counts)
    cdft(ii) = counts(ii) + cdft(ii-1);
end
cdfT = [ttgd',cdft];
% save('cdfT.txt', 'cdfT', '-ascii');
tt = [capTimeGdOrg,capTimeHOrg];
eff_tcap = sum(tt > 8 & tt < 300) ./ sum(tt > 0);
%%
figure;

hh = histogram(tt,'BinWidth', 1, 'DisplayStyle', 'stairs');

xlabel("{\it\DeltaT} (\mus)");    
ylabel("Events / 1\mus");
% ylabel('\fontname{宋体}计数');
% legend('Captured by Gd','Captured by H','Box','off');
% legend('\fontname{宋体}被\fontname{Times new roman}Gd\fontname{宋体}俘获',...
%     '\fontname{宋体}被\fontname{Times new roman}H\fontname{宋体}俘获','Box','off')
set(gca, 'xlim',[0, 400]);
set(gca, 'fontname', 'times new roman', 'xgrid', 'off', 'ygrid', 'off');

%%
hold on
fcolor = '#6279c1';
falpha = 0.3;
% modelFun = @(b,x) b(1).*(1+exp(b(2).*x));
modelFun = @(beta0, x)(beta0(1) / beta0(2) * exp(-x / beta0(2)) + beta0(3));
start = [3e4; 10; 1e2];
wnlm = fitnlm(hh.BinEdges(41:301) - 0.5, hh.BinCounts(40:300), modelFun, start);%, 'Weight', wei(2:4)');
xx = 2:0.1:400;
[ypred, ypredci] = predict(wnlm, xx', 'Simultaneous', false);
plot(gca, xx, ypred, Color = fcolor);
fill(gca, [xx, xx(end:-1:1)], [ypredci(:,1)', fliplr(ypredci(:,2)')] ,...
    'r', FaceAlpha = falpha, EdgeColor = 'none', FaceColor = fcolor);
