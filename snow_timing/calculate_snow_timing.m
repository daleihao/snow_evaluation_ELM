
for year_i = 2015:2015
    
    daynum = days(datetime(year_i,12,31) - datetime(year_i,1,1) + 1);
    doys = 1:daynum;
    load(['../daily_data/ELM_daily_' num2str(year_i) '.mat']);
    end_of_snow = nan(144, 168);
    start_of_snow = nan(144, 168);
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
                end_of_snow(row_i, col_j) = find(tmp == 5, 1, 'last')+4;
            end
            
            fsnos_all = squeeze(fsnos_daily(row_i, col_j, :));
            
            fsnos_all(1:250) = 0;
            is_snowcover = fsnos_all>0.15;
            tmp = conv([1 1 1 1 1],is_snowcover);
            tmp = tmp(5:(end-4));
            
            if(~isempty(find(tmp == 5, 1, 'first')))
                start_of_snow(row_i, col_j) = find(tmp == 5, 1, 'first') + 1;
            end
        end
    end
    nanmean(start_of_snow(:))
    nanmean(end_of_snow(:))
    
    figure;
    subplot(121)
    imagesc(start_of_snow, [250 365])
    colorbar
    subplot(122)
    imagesc(end_of_snow, [1 200])
    colorbar
end

%% ELM
means_sos=nan(19,1);
means_eos=nan(19,1);

for year_i = 2002:2002
    
    daynum = days(datetime(year_i,12,31) - datetime(year_i,1,1) + 1);
    doys = 1:daynum;
    load(['../../../MODIS/daily_data/ELM_daily_' num2str(year_i) '.mat']);
    end_of_snow = nan(144, 168);
    start_of_snow = nan(144, 168);
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
                end_of_snow(row_i, col_j) = find(tmp == 5, 1, 'last')+4;
            end
            
            fsnos_all = squeeze(fsnos_daily(row_i, col_j, :));
            
            fsnos_all(1:250) = 0;
            is_snowcover = fsnos_all>0.15;
            tmp = conv([1 1 1 1 1],is_snowcover);
            tmp = tmp(5:(end-4));
            
            if(~isempty(find(tmp == 5, 1, 'first')))
                start_of_snow(row_i, col_j) = find(tmp == 5, 1, 'first') + 1;
            end
        end
    end
    means_sos(year_i - 2000) = nanmean(start_of_snow(:));
    means_eos(year_i - 2000) = nanmean(end_of_snow(:));
    
    figure;
    subplot(121)
    imagesc(start_of_snow, [250 365])
    colorbar
    subplot(122)
    imagesc(end_of_snow, [1 200])
    colorbar
end

figure;
plot(means_sos)
figure;
plot(means_eos)
%% end

for year_i = 2004:2004
    
    daynum = days(datetime(year_i,12,31) - datetime(year_i,1,1) + 1);
    doys = 1:daynum;
    load(['../../../MODIS/daily_data/ELM_daily_' num2str(year_i) '.mat']);
    end_of_snow2 = nan(144, 168);
    start_of_snow2 = nan(144, 168);
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
                end_of_snow2(row_i, col_j) = find(tmp == 5, 1, 'last')+4;
            end
            
            fsnos_all = squeeze(fsnos_daily(row_i, col_j, :));
            
            fsnos_all(1:250) = 0;
            is_snowcover = fsnos_all>0.15;
            tmp = conv([1 1 1 1 1],is_snowcover);
            tmp = tmp(5:(end-4));
            
            if(~isempty(find(tmp == 5, 1, 'first')))
                start_of_snow2(row_i, col_j) = find(tmp == 5, 1, 'first') + 1;
            end
        end
    end
    figure;
    subplot(121)
    imagesc(start_of_snow2, [250 365])
    colorbar
    subplot(122)
    imagesc(end_of_snow2, [1 150])
    colorbar
end



for year_i = 2004:2004
    
    daynum = days(datetime(year_i,12,31) - datetime(year_i,1,1) + 1);
    doys = 1:daynum;
    load(['../../../MODIS/daily_data/SPIRES_daily_' num2str(year_i) '.mat']);
    end_of_snow2 = nan(144, 168);
    start_of_snow2 = nan(144, 168);
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
                end_of_snow2(row_i, col_j) = find(tmp == 5, 1, 'last')+4;
            end
            
            fsnos_all = squeeze(fsnos_daily(row_i, col_j, :));
            
            fsnos_all(1:250) = 0;
            is_snowcover = fsnos_all>0.15;
            tmp = conv([1 1 1 1 1],is_snowcover);
            tmp = tmp(5:(end-4));
            
            if(~isempty(find(tmp == 5, 1, 'first')))
                start_of_snow2(row_i, col_j) = find(tmp == 5, 1, 'first') + 1;
            end
        end
    end
    figure;
    subplot(121)
    imagesc(start_of_snow2, [250 365])
    colorbar
    subplot(122)
    imagesc(end_of_snow2, [1 150])
    colorbar
end
%% snow melt


for year_i = 2004:2004
    
    
    daynum = days(datetime(year_i,12,31) - datetime(year_i,1,1) + 1);
    doys = 1:daynum;
    load(['daily_data/MCD43_daily_' num2str(year_i) '.mat']);
    melt_of_snow = nan(144, 168);
    num_of_nonon = zeros(144, 168);
    for row_i = 2:144
        for col_j = 24:168
            BSAs_daily2 = squeeze(fsnos_daily(row_i, col_j, :));
            
            num_of_nonon(row_i, col_j) = sum(~isnan(fsnos_all));
            %             fsnos_all(250:end) = 0;
            %             is_snowcover = fsnos_all>0.01;
            %             tmp = conv([1 1 1 1 1],is_snowcover);
            %             tmp = tmp(5:(end-4));
            %             melt_of_snow(row_i, col_j) = tmp;
            %
            opts = fitoptions( 'Method', 'NonlinearLeastSquares');
            ft =  fittype('a + b/(c + exp((x-e)/4))'); %fit function
            mdl = fit(doys(20:250)',BSAs_daily2(20:250)/1000,ft);
            figure;
            plot(mdl)
            hold on
            plot(20:250,BSAs_daily2(20:250)/1000,'*')
        end
    end
    figure;
    subplot(121)
    imagesc(melt_of_snow, [1 365])
    subplot(122)
    imagesc(num_of_nonon)
end

%% code

for year_i = [2001:2019]
    
    
    daynum = days(datetime(year_i,12,31) - datetime(year_i,1,1) + 1);
    doys = 1:daynum;
    load(['../../../MODIS/daily_data/ELM_daily_' num2str(year_i) '.mat']);
    melt_of_snow = nan(144, 168);
    midpoint_of_snow = nan(144, 168);

    for row_i = 1:144
        for col_j = 1:168
            BSAs_daily2 = squeeze(fsnos_daily(row_i, col_j, :));
            
            
            opts = fitoptions( 'Method', 'NonlinearLeastSquares');
            ft =  fittype('a+b/(1 + exp(c*(x-d)))'); %fit function
            x = doys(32:243)';
            y = BSAs_daily2(32:243);
            filter = y>=0;
            if(sum(filter)>30)
                [mdl,gof,out] = fit(x(filter),y(filter),ft, 'StartPoint', [0, 1, 0.2, 100],'Lower',[0,0,0,32],'Upper',[1,1,1,243]);
                estimates = mdl(32:243);
                
                if((gof.rsquare > 0.95) && (gof.rmse < 0.2) && ((estimates(1) - estimates(end))>0.05))
                    
                    dif_tmp = abs((estimates(1) - estimates(212))*0.99 + estimates(212) - estimates);
                    snowmelt_doy = find(dif_tmp == min(dif_tmp), 1)+31;
                    
                    melt_of_snow(row_i, col_j) = snowmelt_doy;
                    midpoint_of_snow(row_i, col_j) = mdl.d;
                end
                
            end
        end
    end    
    save(['ELM_' num2str(year_i) '.mat'],'melt_of_snow','midpoint_of_snow');
end



for year_i = [2001:2019]
    
    
    daynum = days(datetime(year_i,12,31) - datetime(year_i,1,1) + 1);
    doys = 1:daynum;
    load(['../../../MODIS/daily_data/MODSCAG_daily_' num2str(year_i) '.mat']);
    melt_of_snow = nan(144, 168);
    midpoint_of_snow = nan(144, 168);

    for row_i = 1:144
        for col_j = 1:168
            BSAs_daily2 = squeeze(fsnos_daily(row_i, col_j, :));
            
            
            opts = fitoptions( 'Method', 'NonlinearLeastSquares');
            ft =  fittype('a+b/(1 + exp(c*(x-d)))'); %fit function
            x = doys(32:243)';
            y = BSAs_daily2(32:243);
            filter = y>=0;
            if(sum(filter)>30)
                [mdl,gof,out] = fit(x(filter),y(filter),ft, 'StartPoint', [0, 1, 0.2, 100],'Lower',[0,0,0,32],'Upper',[1,1,1,243]);
                estimates = mdl(32:243);
                
                if((gof.rsquare > 0.95) && (gof.rmse < 0.2) && ((estimates(1) - estimates(end))>0.05))
                    
                    dif_tmp = abs((estimates(1) - estimates(212))*0.99 + estimates(212) - estimates);
                    snowmelt_doy = find(dif_tmp == min(dif_tmp), 1)+31;
                    
                    melt_of_snow(row_i, col_j) = snowmelt_doy;
                    midpoint_of_snow(row_i, col_j) = mdl.d;
                end
                
            end
        end
    end    
    save(['MODSCAG_' num2str(year_i) '.mat'],'melt_of_snow','midpoint_of_snow');
end

    
 %% spires
for year_i = [2001:2019]
    
    
    daynum = days(datetime(year_i,12,31) - datetime(year_i,1,1) + 1);
    doys = 1:daynum;
    load(['../../../MODIS/daily_data/SPIRES_daily_' num2str(year_i) '.mat']);
    melt_of_snow = nan(144, 168);
    midpoint_of_snow = nan(144, 168);

    for row_i = 1:144
        for col_j = 1:168
            BSAs_daily2 = squeeze(fsnos_daily(row_i, col_j, :));
            
            
            opts = fitoptions( 'Method', 'NonlinearLeastSquares');
            ft =  fittype('a+b/(1 + exp(c*(x-d)))'); %fit function
            x = doys(32:243)';
            y = BSAs_daily2(32:243);
            filter = y>=0;
            if(sum(filter)>30)
                [mdl,gof,out] = fit(x(filter),y(filter),ft, 'StartPoint', [0, 1, 0.2, 100],'Lower',[0,0,0,32],'Upper',[1,1,1,243]);
                estimates = mdl(32:243);
                
                if((gof.rsquare > 0.95) && (gof.rmse < 0.2) && ((estimates(1) - estimates(end))>0.05))
                    
                    dif_tmp = abs((estimates(1) - estimates(212))*0.99 + estimates(212) - estimates);
                    snowmelt_doy = find(dif_tmp == min(dif_tmp), 1)+31;
                    
                    melt_of_snow(row_i, col_j) = snowmelt_doy;
                    midpoint_of_snow(row_i, col_j) = mdl.d;
                end
                
            end
        end
    end    
    save(['SPIRES_' num2str(year_i) '.mat'],'melt_of_snow','midpoint_of_snow');
end
   