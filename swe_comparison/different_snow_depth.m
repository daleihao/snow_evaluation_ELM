
clc;
clear all;
close all;
%% lat lon
load('../monthly_data.mat');

load('../isUS.mat');


res_v = 0.125;
res_h = 0.125;
lon = (-125+res_h/2):res_h: (-104-res_h/2);
lat = (50-res_v/2):-res_v: (32 + res_v/2);
[lons,lats]=meshgrid(lon,lat);


indexs  = repmat(0:18,3,1);

GridNames = {'Winter', 'Spring', 'Summer', 'Autumn'};


%% plot
colors = flipud(brewermap(101, 'Spectral'));
colors_2 =  flipud(brewermap(101, 'RdBu'));

figure;
set(gcf,'unit','normalized','position',[0.1,0.1,0.45,0.55]);
set(gca, 'Position', [0 0 1 1])

for season_i = 1:4
    switch season_i
        case 1
            seasons_all = [12 1 2];
        case 2
            seasons_all = [3 4 5];
        case 3
            seasons_all = [6 7 8];
        case 4
            seasons_all = [9 10 11];
    end
    
    for month_i = seasons_all
        
        filters = indexs * 12 + seasons_all';
        filters = filters(:);
        filters = filters(filters>36); %2004
        ELM_snowdepths_tmp = nanmean(ELM_snowdepths(:,:,filters),3)*1000;
        UA_snowdepths_tmp = nanmean(UA_snowdepths(:,:,filters),3);
        SNODAS_snowdepths_tmp = nanmean(SNODAS_snowdepths(:,:,filters),3);
        
        ELM_snowdepths_tmp(isUS<1) = nan;
        UA_snowdepths_tmp(isUS<1) = nan;
        SNODAS_snowdepths_tmp(isUS<1) = nan;
        
        
        
        switch season_i
            case 1
                index = "(a) ";
            case 2
                index = "(b) ";
            case 3
                index = "(c) ";
            case 4
                index = "(d) ";
        end
        
        
        
        min_values = [0 0 0];
        max_values = [1000 1000 1000];
        
        %% plot 1
        ax1 = subplot('position', [0.06 + 0.215*(season_i-1) 0.65 0.21 0.28]);
        colormap(ax1, colors);
        hold on
        if(season_i==1)
            plot_global_map_other(lats, lons, ELM_snowdepths_tmp, min_values(1), max_values(1),strcat(index, GridNames{season_i}),1, 0);
        else
            plot_global_map_other(lats, lons, ELM_snowdepths_tmp, min_values(1), max_values(1),strcat( index, GridNames{season_i}),0, 0);
        end
        
        if(season_i==1)
            ylabel({'ELM'})
        end
        
        if(season_i==4)
            hcb = colorbar;
            hcb.Title.String = 'mm';
            hcb.Title.FontSize = 7;
            hcb.Title.FontWeight = 'Bold';
            x1=get(gca,'position');
            x=get(hcb,'Position');
            x(3)=0.012;
            x(4)=0.23;
            x(1)=0.94;
            x(2)=0.68;
            set(hcb,'Position',x)
        end
        set(gca,'fontsize',8,'fontname','time new roman')
        
        ax2 = subplot('position', [0.06 + 0.215*(season_i-1) 0.36 0.21 0.28]);
        colormap(ax2, colors);
        hold on
        if(season_i==1)
            plot_global_map_other(lats, lons, UA_snowdepths_tmp, min_values(2), max_values(2),"",1, 0);
        else
            plot_global_map_other(lats, lons, UA_snowdepths_tmp, min_values(2), max_values(2),"",0, 0);
        end
        
        if(season_i==1)
            ylabel({'UA'})
        end
        if(season_i==4)
            hcb = colorbar;
            hcb.Title.String = 'mm';
            hcb.Title.FontSize = 7;
            hcb.Title.FontWeight = 'Bold';
            x1=get(gca,'position');
            x=get(hcb,'Position');
            x(3)=0.012;
            x(4)=0.23;
            x(1)=0.94;
            x(2)=0.38;
            set(hcb,'Position',x)
        end
        
        ax3 = subplot('position', [0.06 + 0.215*(season_i-1) 0.07 0.21 0.28]);
        colormap(ax3, colors);
        hold on
        if(season_i==1)
            plot_global_map_other(lats, lons, SNODAS_snowdepths_tmp, min_values(3), max_values(3),"",1, 1);
        else
            plot_global_map_other(lats, lons, SNODAS_snowdepths_tmp, min_values(3), max_values(3),"",0, 1);
        end
        
        if(season_i==1)
            ylabel({'SNODAS'})
        end
        if(season_i==4)
            hcb = colorbar;
            hcb.Title.String = 'mm';
            hcb.Title.FontSize = 7;
            hcb.Title.FontWeight = 'Bold';
            x1=get(gca,'position');
            x=get(hcb,'Position');
            x(3)=0.012;
            x(4)=0.23;
            x(1)=0.94;
            x(2)=0.08;
            set(hcb,'Position',x)
        end
        set(gca,'fontsize',8,'fontname','time new roman')
        
    end
end

print(gcf, '-dtiff', '-r300', 'different_snowdepth.tif')

close all

