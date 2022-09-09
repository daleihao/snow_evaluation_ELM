
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
set(gcf,'unit','normalized','position',[0.1,0.1,0.45,0.4]);
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
    
 
        
        filters = indexs * 12 + seasons_all';
        filters = filters(:);
        %filters = filters(filters>36); %2004
        
        ELM_surfacealbedo_tmp = nanmean(ELM_Albedos(:,:,filters),3);
        
        MCD43_surfacealbedo = MCD43_BSAs.* (1 - ELM_Skyls) + MCD43_BSAs .* ELM_Skyls; 
        MCD43_surfacealbedo_tmp = nanmean(MCD43_surfacealbedo(:,:,filters),3)/1000;
        
        ELM_fsnos_tmp = nanmean(ELM_fsnos(:,:,filters),3);

                
        ELM_surfacealbedo_tmp(isUS<1 | ELM_fsnos_tmp<=0) = nan;
        MCD43_surfacealbedo_tmp(isUS<1| ELM_fsnos_tmp<=0) = nan;
        
        
        
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
        max_values = [0.8 0.8 0.8];
        
        %% plot 1
        ax1 = subplot('position', [0.06 + 0.215*(season_i-1) 0.48 0.21 0.4]);
        colormap(ax1, colors);
        hold on
        if(season_i==1)
            plot_global_map_other(lats, lons, ELM_surfacealbedo_tmp, min_values(1), max_values(1),strcat(index, GridNames{season_i}),1, 0);
        else
            plot_global_map_other(lats, lons, ELM_surfacealbedo_tmp, min_values(1), max_values(1),strcat( index, GridNames{season_i}),0, 0);
        end
        
        if(season_i==1)
            ylabel({'ELM'})
        end
        
        if(season_i==4)
            hcb = colorbar;
            hcb.Title.String = 'Unitless';
            hcb.Title.FontSize = 7;
            hcb.Title.FontWeight = 'Bold';
            x1=get(gca,'position');
            x=get(hcb,'Position');
%             x(3)=0.012;
             x(4)=0.33;
            x(1)=0.94;
            x(2)=0.5;
            set(hcb,'Position',x)
        end
        set(gca,'fontsize',8,'fontname','time new roman')
        
        ax2 = subplot('position', [0.06 + 0.215*(season_i-1) 0.07 0.21 0.4]);
        colormap(ax2, colors);
        hold on
        if(season_i==1)
            plot_global_map_other(lats, lons, MCD43_surfacealbedo_tmp, min_values(2), max_values(2),"",1, 1);
        else
            plot_global_map_other(lats, lons, MCD43_surfacealbedo_tmp, min_values(2), max_values(2),"",0, 1);
        end
        
        if(season_i==1)
            ylabel({'MCD43'})
        end
        if(season_i==4)
            hcb = colorbar;
            hcb.Title.String = 'Unitless';
            hcb.Title.FontSize = 7;
            hcb.Title.FontWeight = 'Bold';
            x1=get(gca,'position');
            x=get(hcb,'Position');
%             x(3)=0.012;
             x(4)=0.33;
            x(1)=0.94;
            x(2)=0.09;
            set(hcb,'Position',x)
        end
        
      
end

print(gcf, '-dtiff', '-r300', 'different_surfacealbedo.tif')

close all

