%% code
parpool(18)
parfor year_i = [2001:2019]
    
    
    daynum = days(datetime(year_i,12,31) - datetime(year_i,1,1) + 1);
    doys = 1:daynum;
    load(['../daily_data_N_3000/ELM_daily_' num2str(year_i) '.mat']);
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



parfor year_i = [2001:2019]
    
    
    daynum = days(datetime(year_i,12,31) - datetime(year_i,1,1) + 1);
    doys = 1:daynum;
    load(['../daily_data_N_3000/MODSCAG_daily_' num2str(year_i) '.mat']);
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
parfor year_i = [2001:2019]
    
    
    daynum = days(datetime(year_i,12,31) - datetime(year_i,1,1) + 1);
    doys = 1:daynum;
    load(['../daily_data_N_3000/SPIRES_daily_' num2str(year_i) '.mat']);
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
   