function R = calculateR(obs, sim)
datanum = size(obs(~isnan(obs) & ~isnan(sim)),1);

filters = ~isnan(obs) & ~isnan(sim) & obs>0 & sim>0;
if((datanum<=0) || (sum(filters)/datanum <= 10/91))
    R = nan;
else
    filters = ~isnan(obs) & ~isnan(sim) & obs>=0 & sim>=0;
    
    obs = obs(filters);
    sim = sim(filters);
    
    R_tmp = corrcoef(obs, sim);
    R = R_tmp(2,1);
end