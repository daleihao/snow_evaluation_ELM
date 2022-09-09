function [R,RMSE,Bias] = calculateR2(obs, sim)
filters = ~isnan(obs) & ~isnan(sim) & obs>=0 & sim>=0;
obs = obs(filters);
sim = sim(filters);

if(sum(filters)<10)
    R = nan;
    RMSE = nan;
    Bias = nan;
else
R = corrcoef(obs, sim);
R = R(2,1);
RMSE = sqrt(nanmean((obs-sim).^2));
Bias = nanmean(sim-obs);
end