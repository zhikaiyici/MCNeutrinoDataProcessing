clear
nKE = ["25.3 meV", "10 keV", "100 keV", "500 keV", "1 MeV", "2 MeV", "5 MeV", "10 MeV"];
axGd = axes(figure('Name','Gd'));
% hold(axGd, "on");
axH = axes(figure('Name','H'));
% hold(axH, "on");
axA = axes(figure('Name','All'));
hold([axGd, axH, axA], "on");
for ii = 1:length(nKE)
    runCondition = '_4x4_1e+06_neutron_'  + nKE(ii) + "_INSIDE";
    dirName = num2str(ii - 1) + runCondition;
    cd(dirName)
    % tH = importdata("moduleCapTimeH" + dirName + ".txt");
    % tGd = importdata("moduleCapTimeGd" + dirName + ".txt");
    tH = ReadBinaryFile("moduleEdepDelayH" + runCondition + ".data", 4);
    tGd = ReadBinaryFile("moduleEdepDelayGd" + runCondition + ".data", 4);
    cd ..
    tH = permute(sum(sum(tH)), [3, 1, 2]);
    tGd = permute(sum(sum(tGd)), [3, 1, 2]);
    tH(tH == 0) = [];
    tGd(tGd == 0) = [];
    t = [tH; tGd];

    hGd = histogram(axGd, tGd,'BinWidth', 0.01, 'DisplayStyle', 'stairs', Normalization = 'pdf');
    hH = histogram(axH, tH,'BinWidth', 0.01, 'DisplayStyle', 'stairs', Normalization = 'pdf');
    h = histogram(axA, t,'BinWidth', 0.01, 'DisplayStyle', 'stairs', Normalization = 'pdf');
end
legend(axGd, nKE);
legend(axH, nKE);
legend(axA, nKE);
hold([axGd, axH, axA], 'off');
set([axGd, axH, axA], 'xlim', [0, 10]);
set([axGd, axH, axA], 'yscale', 'log');
