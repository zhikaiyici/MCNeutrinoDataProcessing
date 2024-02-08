clear
%% 载入数据
runCondition = "_4x4_2e+09_CRY";
dirName = "0" + runCondition + "/";
trackFileName = 'moduleMuTrackLength' + runCondition + '.data';

arraySize = 4;
trackData = ReadBinaryFile(dirName + trackFileName, arraySize);

%%
dscrTh = 0.1;
data = trackData ./ 10; 
logic_dscr = data >= dscrTh;
data(~logic_dscr) = 0;
hit = sum(sum(logic_dscr));
hit = permute(hit, [3,1,2]);
droped = hit == 0;
hit(droped,:) = [];
data(:,:,droped) = [];
trajectory = sum(sum(data));
trajectory = permute(trajectory, [3,1,2]);
trajectory1 = trajectory(hit == 1);

nTot = length(trajectory1);
inc = floor(nTot ./ 100);
jj = 1;
for ii = 1:inc:nTot
    if ii + inc < nTot
        data = trajectory1(ii:ii + inc);
    else
        data = trajectory1(ii:end);
    end
    meanTL1(jj,:) = mean(data);
    stdMTL1(jj,:) = std(data);
    jj = jj + 1;
end
%%
meanTL1 = mean(trajectory1);
stdMTL1 = std(stdMTL1(1:100)) .* 100;




