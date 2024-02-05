clear
%% 载入数据
fileName = 'moduleMuEdep(4x4_2e+09_CRY).txt';
% fileName = 'moduleEdepPrompt(4x4_2e+09_CRY).txt';
trackFileName = 'moduleMuTrackLength(4x4_2e+09_CRY).txt';
spectraData_org = importdata(fileName);
trackData_org = importdata(trackFileName);
%%
arraySize = 4;
spectraData = ReshapeDataMatrix(arraySize, spectraData_org);
trackData = ReshapeDataMatrix(arraySize, trackData_org);

%%
data = spectraData; yBinEdges = 0:0.5:220; dscrTh = 0.2; width = 8;
data = trackData ./ 10; yBinEdges = 0:0.5:140; dscrTh = 0.1; width = 4.5;

logic_dscr = data >= dscrTh;
data(~logic_dscr) = 0;
hit = sum(sum(logic_dscr));
hit = permute(hit, [3,1,2]);
droped = hit == 0;
hit(droped,:) = [];
data(:,:,droped) = [];
trajectory = sum(sum(data));
trajectory = permute(trajectory, [3,1,2]);
%%
trajectory1 = trajectory(hit == 1);
totalTL1 = sum(trajectory1);
sigmaTTL1 = sqrt(sum(trajectory1 .^ 2) - totalTL1 .^ 2 ./ length(trajectory1));
meanTL1 = mean(trajectory1);
sigmaTL1 = sigmaTTL1 ./ length(trajectory1);

%%
fcolor = '#6279c1';
falpha = 0.5;

pfig = figure('Units','normalized');
bo = 0.15;
he = 0.7;
pfig.Position = [bo .* 16 ./ 9, bo, he .* 9 ./ 16 ./ 0.75, he];
axM = axes(pfig);
gap = 0.0; ruf = 0.3;
mle = 0.12; mbo = mle;
mwi = (1 - gap - mle - 0.05) ./ (1 + ruf); mhe = mwi; 
axM.Position = [mle, mbo, mwi, mhe];

axR = axes(pfig);
axR.Position = [mwi + gap + mle, mbo, mwi .* ruf, mhe];

axU = axes(pfig);
axU.Position = [mle, mbo + mhe + gap, mwi, mhe .* ruf];

%%
ht = histogram2(axM, hit, trajectory, 0.5:1:16.5, yBinEdges, 'DisplayStyle', 'tile',...
    'ShowEmptyBins', 'off', 'LineStyle', 'none', 'FaceColor', 'flat', Normalization = 'probability');
c = colorbar(axM);
c.Label.String = 'Probability';
c.AxisLocation = "in";
% c.Location = 'north';
cwi = 0.015;
chef = 0.8;
c.Position = [mwi + mle - 2 .* cwi, mbo + mhe .* (1 - chef) ./ 2, cwi, mhe .* chef];
c.FontSize = 18;
colormap(axM, sky);

hh = histogram(axU, hit, 'Normalization', 'probability',...
    EdgeColor = fcolor, FaceColor = fcolor, FaceAlpha = falpha);

%%
% xx = linspace(min(hit), max(hit), 1000);
% [f, yi] = ksdensity(hit, xx);
% fill(axU, [xx, fliplr(xx)], [zeros(size(f)), fliplr(f)],...
%     [0, 0.447, 0.741], 'FaceAlpha', 0.3, 'EdgeColor','none');

vfig = figure('Units','normalized');
vfig.Position = pfig.Position;
axV= axes(vfig);
X = (ht.XBinEdges(1:end-1) + ht.XBinEdges(2:end)) ./ 2;
meanY = zeros(size(X));
stdY = zeros(size(X));
wei = zeros(size(X));

coloraxR = [
     49,  10,  11; 101,  22,  23; 153,  34,  36; 181,  71, 100;
    227,  98,  93; 239, 139, 103; 240, 194, 132; 128, 116, 200;
    120, 149, 193; 168, 203, 223; 214, 239, 244; 242, 250, 252; ];
coloraxR = coloraxR ./ 255;

boxwidth = 0.2;

for ii = 1:length(X)
    x = X(ii);
    y = trajectory(hit == ii);
    meanY(ii) = mean(y);
    stdY(ii) = std(y);
    wei(ii) = length(y) ./ length(trajectory);
    hold(axV, "on");
    hold(axR, "on");
    if length(y) > 10
        violinchart(axV, x, y, fcolor, 0.3, width, boxwidth);
        % boxchart(axV, x .* ones(size(y)), y, BoxWidth = 0.2,...
        %     BoxFaceColor = fcolor, MarkerColor = fcolor, MarkerStyle = '.');
        
        xx = linspace(min(y), max(y), 1000);
        [f, yi] = ksdensity(y, xx);
        fill(axR, [zeros(size(f)), fliplr(f)], [xx, fliplr(xx)],...
            coloraxR(ii,:), 'FaceAlpha', falpha, 'EdgeColor', coloraxR(ii,:));
        
        k = ii;
    else
    scatter(axV, ii, y, 6, 'filled', 'MarkerFaceColor', fcolor, 'MarkerEdgeColor', 'none');
    end
    hold(axV, "off");
    hold(axR, "off");
end
leg = legend(axR, {num2str((1:k)')});
leg.Box = "off";
leg.FontSize = 18;

%%
axfit = axes(pfig);
axfit.Position = axM.Position;
axfit.Box = "off";
axfit.Color = 'none';

falpha = 0.3;

hold(axfit, "on");

errorbar(axfit, X(1:k), meanY(1:k), stdY(1:k), 'o',...
    Color = fcolor, MarkerFaceColor = 'auto');

% [p, S] = polyfit(X(2:4), meanY(2:4), 1);
% xx = 2:0.1:4;
% [yy,ydelta] = polyconf(p, xx, S, 'predopt', 'curve');
% plot(ax2, xx, yy, Color = fitcolor);
% fill([xx, xx(end:-1:1)], [yy - ydelta, fliplr(yy + ydelta)] ,...
%     'r', FaceAlpha = 0.5, EdgeColor = 'none', FaceColor = fitcolor);

wei = hh.BinCounts ./ length(trajectory);
% wei = ones(1, 16);

modelFun = @(b,x) b(1).*(1+exp(b(2).*x));
start = [240; 0.5];
wnlm = fitnlm(X(2:4)', meanY(2:4)', modelFun, start, 'Weight', wei(2:4)');
xx = 2:0.001:4;
[ypred, ypredci] = predict(wnlm, xx', 'Simultaneous', false);
plot(axfit, xx, ypred, Color = fcolor);
fill(axfit, [xx, xx(end:-1:1)], [ypredci(:,1)', fliplr(ypredci(:,2)')] ,...
    'r', FaceAlpha = falpha, EdgeColor = 'none', FaceColor = fcolor);

modelFun = @(b,x) b(1).*(1-exp(-b(2).*x));
start = [240; 0.5];
wnlm1 = fitnlm(X(5:7)', meanY(5:7)', modelFun, start, 'Weight', wei(5:7)');
xx = 5:0.001:7;
[ypred, ypredci] = predict(wnlm1, xx', 'Simultaneous', false);
plot(axfit, xx, ypred, Color = fcolor);
fill(axfit, [xx, xx(end:-1:1)], [ypredci(:,1)', fliplr(ypredci(:,2)')] ,...
    'r', FaceAlpha = falpha, EdgeColor = 'none', FaceColor = fcolor);

modelFun = @(b,x) b(1).*(1+exp(b(2).*x));
start = [240; 0.5];
wnlm2 = fitnlm(X(8:11)', meanY(8:11)', modelFun, start, 'Weight', wei(8:11)');
xx = 8:0.001:11;
[ypred,ypredci] = predict(wnlm2 ,xx', 'Simultaneous', false);
plot(axfit, xx, ypred, Color = fcolor);
fill(axfit, [xx, xx(end:-1:1)], [ypredci(:,1)', fliplr(ypredci(:,2)')] ,...
    'r', FaceAlpha = falpha, EdgeColor = 'none', FaceColor = fcolor);

hold(axfit, "off");

save('hitfit-2-4.mat','wnlm');
save('hitfit-5-7.mat','wnlm1');
save('hitfit-8-11.mat','wnlm2');

%%
linkaxes([axM, axR], 'y');
linkaxes([axM, axU], 'x');
linkaxes([axM, axfit], 'xy');

% axR.YTickLabel = '';
axR.XAxis.Visible = "off";
axR.YAxis.Visible = "off";
% axU.XTickLabel = '';
axU.XAxis.Visible = "off";
axU.YAxis.Visible = "off";
axM.XAxis.Visible = "off";
axM.YAxis.Visible = "off";
xlabel(axfit, 'Hit number');
ylabel(axfit, 'Muon trajectory length / cm');
xlabel(axV, 'Hit number');
ylabel(axV, 'Muon trajectory length / cm');
% set([axM, axR, axU, axV, axfit],...
%     "XGrid","off", "YGrid", "off", "ZGrid", "off", "Box", "off", ...
%     'fontname', 'times new roman', 'fontsize', 20);
myfigstyle([axM, axR, axU, axV, axfit]);
% set(axV, 'fontsize', 12);
set([axM, axV], 'ylim', [0, yBinEdges(end)]);
set(axfit, 'xlim', [0, 16], 'XTick', 2:2:16);
set(axV, 'xlim', [0, 14], 'XTick', 1:14);

%%
orient(vfig, 'landscape');
print(vfig, 'vio','-dpdf');
orient(pfig, 'landscape');
print(pfig, 'hit','-dpdf');


