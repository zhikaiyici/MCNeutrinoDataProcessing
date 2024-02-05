clear
nKE = ["25.3 meV", "10 keV", "100 keV", "500 keV", "1 MeV", "2 MeV", "5 MeV", "10 MeV"];
axGd = axes(figure('Name','Gd'));
% hold(axGd, "on");
axH = axes(figure('Name','H'));
% hold(axH, "on");
axA = axes(figure('Name','All'));
hold([axGd, axH, axA], "on");
for ii = 1:8
    runCondition = '_4x4_1e+06_neutron_'  + nKE(ii) + "_INSIDE";
    dirName = num2str(ii - 1) + runCondition;
    cd(dirName)
    % tH = importdata("moduleCapTimeH" + dirName + ".txt");
    % tGd = importdata("moduleCapTimeGd" + dirName + ".txt");
    tH = ReadBinaryFile("moduleCapTimeH" + runCondition + ".data", 4);
    tGd = ReadBinaryFile("moduleCapTimeGd" + runCondition + ".data", 4);
    cd ..
    tH(tH == 0) = [];
    tGd(tGd == 0) = [];
    t = [tH'; tGd'];

    hGd = histogram(axGd, tGd,'BinWidth',1, 'DisplayStyle', 'stairs', Normalization = 'pdf');
    hH = histogram(axH, tH,'BinWidth',1, 'DisplayStyle', 'stairs', Normalization = 'pdf');
    h = histogram(axA, t,'BinWidth',1, 'DisplayStyle', 'stairs', Normalization = 'pdf');
end
legend(axGd, nKE);
legend(axH, nKE);
legend(axA, nKE);
set([axGd, axH, axA], 'xlim',[0, 400]);
hold([axGd, axH, axA], 'off');
%%
e = hH.BinEdges;
e = e(1:end-1) + e(2:end);
e = e ./ 2;
c = hH.BinCounts;
%%
e = hGd.BinEdges;
e = e(1:end-1) + e(2:end);
e = e ./ 2;
e = e(10:end);
c = hGd.BinCounts;
c = c(10:end);
%%
e = h.BinEdges;
e = e(1:end-1) + e(2:end);
e = e ./ 2;
c = h.BinCounts;
e = e(60:end);
c = c(60:end);
%%
modelFun = @(beta0, x)(beta0(1) / beta0(2) * exp(-x / beta0(2)) + beta0(3));
beta0 = [50000; 10; 100];
nlm = fitnlm(e, c, modelFun, beta0);

xx = 0:0.01:400;
[ypred, ypredci] = predict(nlm, xx', 'Simultaneous', false);
plot(xx, ypred);
fill([xx, xx(end:-1:1)], [ypredci(:,1)', fliplr(ypredci(:,2)')] ,...
    'r', FaceAlpha = 0.3, EdgeColor = 'none');

nlm.Coefficients.Estimate(2)
nlm.Rsquared.Adjusted


