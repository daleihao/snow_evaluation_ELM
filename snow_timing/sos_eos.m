
for year_i = 2014:2014
    
    daynum = days(datetime(year_i,12,31) - datetime(year_i,1,1) + 1);
    doys = 1:daynum;
    load(['../daily_data/ELM_daily_' num2str(year_i) '.mat']);
    ELM_end_of_snow = nan(144, 168);
    ELM_start_of_snow = nan(144, 168);
    num_of_nonon = zeros(144, 168);
    for row_i = 1:144
        for col_j = 1:168
            fsnos_all = squeeze(fsnos_daily(row_i, col_j, :));
            
            num_of_nonon(row_i, col_j) = sum(fsnos_all>0.15);
            fsnos_all(250:end) = 0;
            is_snowcover = fsnos_all>0.15;
            tmp = conv([1 1 1 1 1],is_snowcover);
            tmp = tmp(5:(end-4));
            
            if(~isempty(find(tmp == 5, 1, 'last')))
                ELM_end_of_snow(row_i, col_j) = find(tmp == 5, 1, 'last')+4;
            end
            
            fsnos_all = squeeze(fsnos_daily(row_i, col_j, :));
            
            fsnos_all(1:250) = 0;
            is_snowcover = fsnos_all>0.15;
            tmp = conv([1 1 1 1 1],is_snowcover);
            tmp = tmp(5:(end-4));
            
            if(~isempty(find(tmp == 5, 1, 'first')))
                ELM_start_of_snow(row_i, col_j) = find(tmp == 5, 1, 'first') + 1;
            end
        end
    end
    nanmean(ELM_start_of_snow(:))
    nanmean(ELM_end_of_snow(:))
    
    figure;
    subplot(121)
    imagesc(ELM_start_of_snow, [250 365])
    colorbar
    subplot(122)
    imagesc(ELM_end_of_snow, [1 200])
    colorbar
end

%% plot

for year_i = 2019:2019
    
    daynum = days(datetime(year_i,12,31) - datetime(year_i,1,1) + 1);
    doys = 1:daynum;
    load(['../daily_data/MODSCAG_daily_' num2str(year_i) '.mat']);
    MODSCAG_end_of_snow = nan(144, 168);
    MODSCAG_start_of_snow = nan(144, 168);
    num_of_nonon = zeros(144, 168);
    for row_i = 1:144
        for col_j = 1:168
            fsnos_all = squeeze(fsnos_daily(row_i, col_j, :));
            
            num_of_nonon(row_i, col_j) = sum(fsnos_all>0.15);
            fsnos_all(250:end) = 0;
            is_snowcover = fsnos_all>0.15;
            tmp = conv([1 1 1 1 1],is_snowcover);
            tmp = tmp(5:(end-4));
            
            if(~isempty(find(tmp == 5, 1, 'last')))
                MODSCAG_end_of_snow(row_i, col_j) = find(tmp == 5, 1, 'last')+4;
            end
            
            fsnos_all = squeeze(fsnos_daily(row_i, col_j, :));
            
            fsnos_all(1:250) = 0;
            is_snowcover = fsnos_all>0.15;
            tmp = conv([1 1 1 1 1],is_snowcover);
            tmp = tmp(5:(end-4));
            
            if(~isempty(find(tmp == 5, 1, 'first')))
                MODSCAG_start_of_snow(row_i, col_j) = find(tmp == 5, 1, 'first') + 1;
            end
        end
    end
    nanmean(MODSCAG_start_of_snow(:))
    nanmean(MODSCAG_end_of_snow(:))
    
    figure;
    subplot(121)
    imagesc(MODSCAG_start_of_snow, [250 365])
    colorbar
    subplot(122)
    imagesc(MODSCAG_end_of_snow, [1 200])
    colorbar
end

%% plot

for year_i = 2014:2014
    
    daynum = days(datetime(year_i,12,31) - datetime(year_i,1,1) + 1);
    doys = 1:daynum;
    load(['../daily_data/SPIRES_daily_' num2str(year_i) '.mat']);
    SPIRES_end_of_snow = nan(144, 168);
    SPIRES_start_of_snow = nan(144, 168);
    num_of_nonon = zeros(144, 168);
    for row_i = 1:144
        for col_j = 1:168
            fsnos_all = squeeze(fsnos_daily(row_i, col_j, :));
            
            num_of_nonon(row_i, col_j) = sum(fsnos_all>0.15);
            fsnos_all(250:end) = 0;
            is_snowcover = fsnos_all>0.15;
            tmp = conv([1 1 1 1 1],is_snowcover);
            tmp = tmp(5:(end-4));
            
            if(~isempty(find(tmp == 5, 1, 'last')))
                SPIRES_end_of_snow(row_i, col_j) = find(tmp == 5, 1, 'last')+4;
            end
            
            fsnos_all = squeeze(fsnos_daily(row_i, col_j, :));
            
            fsnos_all(1:250) = 0;
            is_snowcover = fsnos_all>0.15;
            tmp = conv([1 1 1 1 1],is_snowcover);
            tmp = tmp(5:(end-4));
            
            if(~isempty(find(tmp == 5, 1, 'first')))
                SPIRES_start_of_snow(row_i, col_j) = find(tmp == 5, 1, 'first') + 1;
            end
        end
    end
    nanmean(SPIRES_start_of_snow(:))
    nanmean(SPIRES_end_of_snow(:))
    
    figure;
    subplot(121)
    imagesc(SPIRES_start_of_snow, [250 365])
    colorbar
    subplot(122)
    imagesc(SPIRES_end_of_snow, [1 200])
    colorbar
end

save('start_end_snow.mat','SPIRES_start_of_snow','SPIRES_end_of_snow',...
    'MODSCAG_start_of_snow','MODSCAG_end_of_snow',...
    'ELM_start_of_snow','ELM_end_of_snow'...
    );