function R = calculateR(obs, sim,obs2,sim2)
datanum = size(obs(~isnan(obs) & ~isnan(sim)),1);

filters = ~isnan(obs) & ~isnan(sim) & obs>0 & sim>0 & ~isnan(obs2) & ~isnan(sim2) & obs2>0 & sim2>0;
if((datanum<=0) || (sum(filters)/datanum <= 10/91) || sum(filters)<=10)
    R = nan;
else
    filters = ~isnan(obs) & ~isnan(sim) & obs>0 & sim>0 & ~isnan(obs2) & ~isnan(sim2) & obs2>0 & sim2>0;
    
    obs = obs(filters);
    sim = sim(filters);
    
    R_tmp = corrcoef(obs, sim);
    R = R_tmp(2,1);
end