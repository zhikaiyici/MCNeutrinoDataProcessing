clear
%% 载入数据
arraySize = 4;
NDLName = "ENDF-VIII.0/";
array = "_" + num2str(arraySize) + "x" + num2str(arraySize);
runCondition = array + '_2e+09_CRY_MUON';

trackData = [];
for runID = 0:0
    dirName = NDLName + num2str(runID) + runCondition + "/";
    fileName = dirName + 'moduleMuTrackLength' + runCondition + ".data";
    trackData = cat(3, trackData, ReadBinaryFile(fileName, arraySize, 0));
end

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
sliceNum = 100;
meanTL1 = zeros(sliceNum, 1);
stdMTL1 = zeros(sliceNum, 1);
inc = ceil(nTot ./ sliceNum);
jj = 1;
for ii = 1:inc:nTot
    if ii + inc < nTot
        data = trajectory1(ii + 1:ii + inc);
    else
        data = trajectory1(ii + 1:end);
    end
    meanTL1(jj,:) = mean(data);
    stdMTL1(jj,:) = std(data);
    jj = jj + 1;
end
%%
meanTL1_ = mean(trajectory1);
stdMTL1_ = std(trajectory1);




