function alloutput = calStat(mod,obs,areas)
filters = mod > 0 & obs >0;

mod = mod(filters);
obs = obs(filters);
areas = areas(filters);

areas = areas./nansum(areas(:));
R = corrcoef(mod, obs);
R = R(1,2);

Bias = nansum((mod - obs).*areas(:));
rBias = Bias./nansum(obs.*areas(:))*100;

RMSE = sqrt(nansum(((mod - obs).^2).*areas(:)));
rRMSE = RMSE./nansum(obs.*areas(:))*100;

alloutput = [R, Bias, rBias, RMSE, rRMSE];
end

