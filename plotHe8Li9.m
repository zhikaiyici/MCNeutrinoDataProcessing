clear

NDLName = "ENDF-VIII.0/";
arraySize = 4;
array = "_" + num2str(arraySize) + "x" + num2str(arraySize);
runCondition = array + '_1e+06_He8';
runID = 0;
dirName = NDLName + num2str(runID) + runCondition + "/";

% name = [dirName, '\betaKEHe8', dirName(2:end), '.data'];
name = dirName + 'betaKEHe8' + runCondition + '.data';
bkeHe8 = ReadBinaryFile(name, 1);
figure(Name = name)
h = histogram(bkeHe8, 'DisplayStyle', 'stairs', BinEdges = 0:50:14000, Normalization = 'pdf');
set(gca,'xlim',[0,14000])

he8 = reshape(bkeHe8, 100000, []);
for ii = 1:size(he8, 2)
    [cHe8(ii,:), eHe8(ii,:)] = histcounts(he8(:,ii), 0:50:14000);
end
scHe8 = std(cHe8);

name = [dirName, '\decayTimeHe8', dirName(2:end), '.data'];
name = dirName + 'decayTimeHe8' + runCondition + '.data';
dtHe8 = ReadBinaryFile(name, 1);
figure(Name = name)
histogram(dtHe8, 'DisplayStyle', 'stairs');

name = [dirName, '\neutronGenicTime', dirName(2:end), '.data'];
name = dirName + 'ModuleNeutronGenicTime' + runCondition + '.data';
gtHe8 = ReadBinaryFile(name, 4);
gtHe8(gtHe8 == 0) = [];
figure(Name = name)
histogram(gtHe8,'DisplayStyle','stairs');

runCondition = array + '_1e+06_Li9';
dirName = NDLName + num2str(runID) + runCondition + "/";

name = [dirName, '\betaKELi9', dirName(2:end), '.data'];
name = dirName + 'betaKELi9' + runCondition + '.data';
bkeLi9 = ReadBinaryFile(name, 1);
figure(Name = name)
histogram(bkeLi9,'DisplayStyle','stairs');
set(gca,'xlim',[0,14000])

li9 = reshape(bkeLi9, 100000, []);
for ii = 1:size(li9, 2)
    [cLi9(ii,:), eLi9(ii,:)] = histcounts(li9(:,ii), 0:50:14000);
end
scLi9 = std(cLi9);

% name = [dirName, '\decayTimeLi9', dirName(2:end), '.data'];
name = dirName + 'decayTimeLi9' + runCondition + '.data';
dtLi9 = ReadBinaryFile(name, 1);
figure(Name = name)
histogram(dtLi9,'DisplayStyle','stairs');

name = [dirName, '\neutronGenicTime', dirName(2:end), '.data'];
name = dirName + 'ModuleNeutronGenicTime' + runCondition + '.data';
gtLi9 = ReadBinaryFile(name, 4);
gtLi9(gtLi9 == 0) = [];
figure(Name = name)
histogram(gtLi9, 'DisplayStyle', 'stairs');

figure;
errorbar(25:50:14000,mean(cLi9),scLi9);
hold on;
errorbar(25:50:14000,mean(cHe8),scHe8);
hold off;

figure
histogram(bkeHe8, 'DisplayStyle', 'stairs', BinEdges = 0:50:14000, Normalization = 'pdf');
hold on
histogram(bkeLi9, 'DisplayStyle', 'stairs', BinEdges = 0:50:14000, Normalization = 'pdf');
hold off
