clear
%% 载入数据
fileName = 'moduleMuEdep(4x4_2e+09_CRY).txt';
spectraData_org = importdata(fileName);
%%
arraySize = 4;
spectraData = ReshapeDataMatrix(arraySize, spectraData_org);

%%
dscrTh = 0.2;
nTot = size(spectraData, 3);
inc = floor(nTot ./ 100);
jj = 1;
for ii = 1:inc:nTot
    if ii + inc < nTot
        data = spectraData(:,:,ii:ii + inc);
    else
        data = spectraData(:,:,ii:end);
    end
    nE(jj, 1)  = size(data, 3);

    logic_dscr = data >= dscrTh;
    data(~logic_dscr) = 0;
    hit = sum(sum(logic_dscr));
    hit = permute(hit, [3,1,2]);
    droped = hit == 0;
    hit(droped,:) = [];
    data(:,:,droped) = [];
    edep = sum(sum(data));
    edep = permute(edep, [3,1,2]);

    logicMu = hit > 1 & edep > 10; 
    
    nMu(jj, 1) = sum(logicMu);

    e(jj, 1) = sum(logicMu) ./ size(data, 3);
    jj = jj + 1;
end
%%
effMu = sum(nMu(1:100)) ./ sum(nE(1:100));
stdnMu = std(nMu(1:100));
stdEffMu = stdnMu .* 100 ./ sum(nE(1:100));
%%
nLost = nE - nMu;
effLost = sum(nLost(1:100)) ./ sum(nE(1:100));
stdLost = std(nLost(1:100));
stdEffLost = stdLost .* 100 ./ sum(nE(1:100));

