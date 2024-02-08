clear
% nKE = ["25.3 meV", "10 keV", "100 keV", "500 keV", "1 MeV", "2 MeV", "5 MeV", "10 MeV"];
nKE = "100 keV";
axGd = axes(figure('Name','Gd'));
% hold(axGd, "on");
axH = axes(figure('Name','H'));
% hold(axH, "on");
axA = axes(figure('Name','All'));
hold([axGd, axH, axA], "on");
for ii = 1:length(nKE)
    runCondition = '_4x4_1e+07_neutron_'  + nKE(ii) + "_INSIDE";
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

    hGd = histogram(axA, tGd,'BinWidth',1, 'DisplayStyle', 'stairs', Normalization = 'count');
    hH = histogram(axA, tH,'BinWidth',1, 'DisplayStyle', 'stairs', Normalization = 'count');
    hA = histogram(axA, t,'BinWidth',1, 'DisplayStyle', 'stairs', Normalization = 'count');
    
    %%
    eH = hH.BinEdges;
    eH = eH(1:end-1) + eH(2:end);
    eH = eH ./ 2;
    cH = hH.BinCounts;
    %%
    eGd = hGd.BinEdges;
    eGd = eGd(1:end-1) + eGd(2:end);
    eGd = eGd ./ 2;
    cGd = hGd.BinCounts;
    eGd = eGd(10:end);
    cGd = cGd(10:end);
    %%
    eA = hA.BinEdges;
    eA = eA(1:end-1) + eA(2:end);
    eA = eA ./ 2;
    cA = hA.BinCounts;
    eA = eA(80:end);
    cA = cA(80:end);
    %%
    modelFun = @(beta0, x)(beta0(1) / beta0(2) * exp(-x / beta0(2)) + beta0(3));
    beta0 = [5000000; 100; 10];
    w = cA;
    w(w == 0) = 1;
    nlmA = fitnlm(eA, cA, modelFun, beta0, 'Weight', w);

    xx = 0:0.1:400;
    [ypredA, ypredciA] = predict(nlmA, xx', 'Simultaneous', false);
    % yyaxis(axA, 'right');
    plot(axA, xx, ypredA);
    fill(axA, [xx, xx(end:-1:1)], [ypredciA(:,1)', fliplr(ypredciA(:,2)')] ,...
        'r', FaceAlpha = 0.3, EdgeColor = 'none');

    coef = nlmA.Coefficients.Estimate(2);
    r2 = nlmA.Rsquared.Adjusted;
    ci = coefCI(nlmA);
end
legend(axGd, nKE);
legend(axH, nKE);
legend(axA, nKE);
set([axGd, axH, axA], 'xlim',[0, 400]);
hold([axGd, axH, axA], 'off');

