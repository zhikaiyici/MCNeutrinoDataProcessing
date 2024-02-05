
dirName = '0_4x4_1e+07_He8';
name = [dirName, '\betaKEHe8', dirName(2:end), '.data'];
bkeHe8 = ReadBinaryFile(name, 1);
figure(Name = name)
h = histogram(bkeHe8, 'DisplayStyle', 'stairs');
set(gca,'xlim',[0,14000])

he8 = reshape(bkeHe8, 100000, []);
for ii = 1:size(he8, 2)
    [c(ii,:), e(ii,:)] = histcounts(he8(:,ii), 0:50:14000);
end

name = [dirName, '\decayTimeHe8', dirName(2:end), '.data'];
dtHe8 = ReadBinaryFile(name, 1);
figure(Name = name)
histogram(dtHe8, 'DisplayStyle', 'stairs');

name = [dirName, '\neutronGenicTime', dirName(2:end), '.data'];
gtHe8 = ReadBinaryFile(name, 4);
gtHe8(gtHe8 == 0) = [];
figure(Name = name)
histogram(gtHe8,'DisplayStyle','stairs');

dirName = '1_4x4_1e+07_Li9';
name = [dirName, '\betaKELi9', dirName(2:end), '.data'];
bkeLi9 = ReadBinaryFile(name, 1);
figure(Name = name)
histogram(bkeLi9,'DisplayStyle','stairs');
set(gca,'xlim',[0,14000])

name = [dirName, '\decayTimeLi9', dirName(2:end), '.data'];
dtLi9 = ReadBinaryFile(name, 1);
figure(Name = name)
histogram(dtLi9,'DisplayStyle','stairs');

name = [dirName, '\neutronGenicTime', dirName(2:end), '.data'];
gtLi9 = ReadBinaryFile(name, 4);
gtLi9(gtLi9 == 0) = [];
figure(Name = name)
histogram(gtLi9, 'DisplayStyle', 'stairs');


