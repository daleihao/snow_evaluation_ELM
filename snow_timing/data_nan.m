for year_i = 2001:2019
    
    daynum = days(datetime(year_i,12,31) - datetime(year_i,1,1) + 1);
    doys = 1:daynum;
    load(['../../../MODIS/daily_data/ELM_daily_' num2str(year_i) '.mat']);
    figure;
    subplot(131)
    imagesc(sum(isnan(snowdepths_daily),3)>1);
    colorbar
    title(year_i)
    
    subplot(132)
    imagesc(sum(isnan(fsnos_daily),3)>1);
    colorbar
    title(year_i)
    
    subplot(133)
    imagesc(sum(isnan(swes_daily),3)>1);
    colorbar
    title(year_i)
    
end

%%
for year_i = 2001:2019
    
    daynum = days(datetime(year_i,12,31) - datetime(year_i,1,1) + 1);
    doys = 1:daynum;
    load(['../../../MODIS/daily_data/MODSCAG_daily_' num2str(year_i) '.mat']);
    figure;
    subplot(221)
    imagesc(sum(isnan(fsnos_daily),3)>0);
    colorbar
    title(year_i)
    
    subplot(222)
    imagesc(sum(isnan(fsnos_daily),3)>1);
    colorbar
    title(year_i)
    
    load(['../../../MODIS/daily_data/SPIRES_daily_' num2str(year_i) '.mat']);
    
    subplot(223)
    imagesc(sum(isnan(fsnos_daily),3)>0);
    colorbar
    title(year_i)
    subplot(224)
    imagesc(sum(isnan(fsnos_daily),3)>1);
    colorbar
    title(year_i)
    
end



%%
for year_i = 2001:2019
    
    daynum = days(datetime(year_i,12,31) - datetime(year_i,1,1) + 1);
    doys = 1:daynum;
    load(['../../../MODIS/daily_data/UA_daily_' num2str(year_i) '.mat']);
    figure;
    subplot(221)
    imagesc(sum(isnan(swes_daily),3)>0);
    colorbar
    title(year_i)
    
    subplot(222)
    imagesc(sum(isnan(swes_daily),3)>1);
    colorbar
    title(year_i)
    
    
    subplot(223)
    imagesc(sum(isnan(snowdepths_daily),3)>0);
    colorbar
    title(year_i)
    subplot(224)
    imagesc(sum(isnan(snowdepths_daily),3)>1);
    colorbar
    title(year_i)
    
end

%%
for year_i = 2001:2019
    
    daynum = days(datetime(year_i,12,31) - datetime(year_i,1,1) + 1);
    doys = 1:daynum;
    load(['../../../MODIS/daily_data/SNODAS_daily_' num2str(year_i) '.mat']);
    figure;
    subplot(221)
    imagesc(sum(isnan(swes_daily),3));
    colorbar
    title(year_i)
    
    subplot(222)
    imagesc(sum(isnan(swes_daily),3)>1);
    colorbar
    title(year_i)
    
    
    subplot(223)
    imagesc(sum(isnan(snowdepths_daily),3));
    colorbar
    title(year_i)
    subplot(224)
    imagesc(sum(isnan(snowdepths_daily),3)>1);
    colorbar
    title(year_i)
    
end