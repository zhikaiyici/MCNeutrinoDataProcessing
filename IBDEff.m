clear;

syms effTotal(effPrompt, effDelayed, effCapTime, effCapNeutron);
effTotal(effPrompt, effDelayed, effCapTime, effCapNeutron) = ...
    effPrompt * effDelayed * effCapTime * effCapNeutron;
dEffPrompt = diff(effTotal, effPrompt);
dEffDelayed = diff(effTotal, effDelayed);
dEffCapTime = diff(effTotal, effCapTime);
dEffCapNeutron = diff(effTotal, effCapNeutron);

stdCapNeutron = 0.0012;
stdPrompt = sqrt((0.0106 .* 0.0028) .^ 2 + (0.7261 .* 0.0014) .^ 2);
stdDelayed = sqrt((0.2426 .* 0.0010) .^ 2 + (0.8679 .* 0.0050) .^ 2);
stdCapTime = 0.0037;

meanCapNeutron = 0.1031;
meanPrompt = 0.0106 .* 0.7261;
meanDelayed = 0.2426 .* 0.8679;
meanCapTime = 0.8846;

meanTotal = double(effTotal(meanPrompt, meanDelayed, meanCapTime, meanCapNeutron));
stdTotal = double(sqrt( ...
    (dEffPrompt(meanPrompt, meanDelayed, meanCapTime, meanCapNeutron) .* stdPrompt) .^ 2 + ...
    (dEffDelayed(meanPrompt, meanDelayed, meanCapTime, meanCapNeutron) .* stdDelayed) .^ 2 + ...
    (dEffCapTime(meanPrompt, meanDelayed, meanCapTime, meanCapNeutron) .* stdCapTime) .^ 2 + ...
    (dEffCapNeutron(meanPrompt, meanDelayed, meanCapTime, meanCapNeutron) .* stdCapNeutron) .^ 2));
