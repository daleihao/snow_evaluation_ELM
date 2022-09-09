clc
clear all;

%% code

for year_i = [2015:2015]
    
    
    daynum = days(datetime(year_i,12,31) - datetime(year_i,1,1) + 1);
    
    % 32 85 95 45 94 35 90
    for row_i = 32:120
        for col_j = 85:103
            doys = 1:daynum;

            load(['../../../../MODIS/daily_data/ELM_daily_' num2str(year_i-1) '.mat']);
            fsnos = squeeze(fsnos_daily(row_i, col_j, :));
            fsnos_1 = fsnos(1:365);
            
            load(['../../../../MODIS/daily_data/ELM_daily_' num2str(year_i) '.mat']);
            fsnos = squeeze(fsnos_daily(row_i, col_j, :));
            fsnos_2 = fsnos(1:365);
            
            fsnos = [fsnos_1; fsnos_2];
            doys = [doys doys+365];
            
            
            fsnos_all = fsnos_2;
            
            num_of_nonon(row_i, col_j) = sum(fsnos_all>0.05);%0.15-0
            fsnos_all(250:end) = 0;
            is_snowcover = fsnos_all>0.05;%0.15-0
            tmp = conv([1 1 1 1 1],is_snowcover);
            tmp = tmp(5:(end-4));
            
            if(~isempty(find(tmp == 5, 1, 'last')))
                MODSCAG_end_of_snow = find(tmp == 5, 1, 'last')+4;
            end
            
            fsnos_all = fsnos_1;
            
            fsnos_all(1:250) = 0;
            is_snowcover = fsnos_all>0.05;%0.15-0
            tmp = conv([1 1 1 1 1],is_snowcover);
            tmp = tmp(5:(end-4));
            
            if(~isempty(find(tmp == 5, 1, 'first')))
                MODSCAG_start_of_snow = find(tmp == 5, 1, 'first') ;
            end
            
            opts = fitoptions( 'Method', 'NonlinearLeastSquares');
            ft =  fittype('a+b/(1 + exp(c*(x-d)))'); %fit function
            x = doys(32:243)';
            y = fsnos( 365 + (32:243));
            filter = y>=0;
            if(sum(filter)>30)
                [mdl,gof,out] = fit(x(filter),y(filter),ft, 'StartPoint', [0, 1, 0.2, 100],'Lower',[0,0,0,32],'Upper',[1,1,1,243]);
                estimates = mdl(32:243);
                
                %if((gof.rsquare > 0.95) && (gof.rmse < 0.2) && ((estimates(1) - estimates(end))>0.05))
                    
                    dif_tmp = abs((estimates(1) - estimates(212))*0.99 + estimates(212) - estimates);
                    snowmelt_doy = find(dif_tmp == min(dif_tmp), 1)+31;
                    
                    
                    figure;
                    set(gcf,'unit','normalized','position',[0.1,0.1,0.8,0.5]);
                    set(gca, 'Position', [0.08 0.15 0.9 0.8])
                    hold on
                    %plot(doys, fsnos,'o','MarkerFaceColor',[0.5 0.5
                    %0.5],'MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5,'Linewidth',1);
                    
                    pc_x = [MODSCAG_start_of_snow MODSCAG_end_of_snow+365 MODSCAG_end_of_snow+365 MODSCAG_start_of_snow];
                    pc_y = [0.04 0.04 0.06 0.06];
                    patch(pc_x,pc_y,[0.5 0.5 0.5],'FaceColor',[0 0 1],'FaceAlpha',.5,'EdgeColor',[0.5 0.5 0.5]);
                    
                    pc_x = [244 244+153 244+153 244];
                    pc_y = [0 0 1 1];
                    patch(pc_x,pc_y,[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor',[0.5 0.5 0.5]);
                    alpha(0.2)
                    plot(doys, fsnos,'-k','Linewidth',1);
                    plot(x+365, estimates,'-r','linewidth',6,'color',[0.8500 0.3250 0.0980 1])
                    plot([1 730],[0.05 0.05],'--k')
                    plot([1 730],[0.9 0.9],'-k')

                    plot([snowmelt_doy+365 snowmelt_doy+365], [0 0.9],'-','linewidth',6,'color',[0 0 1 0.4]);
                    plot([mdl.d+365 mdl.d+365], [0 0.9],'-','linewidth',6,'color',[0 0 1 0.4]);
                    plot([MODSCAG_start_of_snow MODSCAG_start_of_snow], [0 0.9],'-','linewidth',6,'color',[0 0 1 0.4])
                    plot([MODSCAG_end_of_snow+365 MODSCAG_end_of_snow+365], [0 0.9],'-','linewidth',6,'color',[0 0 1 0.4])

                                        
                    text(snowmelt_doy+4+365, 0.82, 'Depletion onset date','fontsize',12,'fontweight','bold');
                    text(mdl.d+4+365, 0.4, 'Midpoint date','fontsize',12,'fontweight','bold');
                    text(MODSCAG_start_of_snow+4, 0.12, 'Accumulation onset date','fontsize',12,'fontweight','bold');
                    text(MODSCAG_end_of_snow+4+365, 0.12, 'End date','fontsize',12,'fontweight','bold');
                    text(383, 0.1, 'Duration','fontsize',12,'fontweight','bold');
                    text(265, 0.95, 'Snow accumulation season','fontsize',16,'fontweight','bold');
                    text(470, 0.95, 'Snowmelt season','fontsize',16,'fontweight','bold');
 
                    
                    ylabel('{\it f}_{sno}')

                    xlim([244 243+365])
                    ylim([0 1])
                    
                    box on
                    set(gca,'linewidth',1,'fontsize',16,...
                        'xtick',[258 288 319 349 365+15 365+46 365+75 365+106 365+136 365+167 365+197 365+228],...
                        'xticklabel',{'Sept','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug'})
                    
                    %print(gcf, '-dtiff', '-r300', '../../figure_all_tif/sos_extraction_method.tif')
                    
                    close all
               % end
                
            end
        end
    end
end

