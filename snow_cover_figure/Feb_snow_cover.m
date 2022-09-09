
clc;
clear all;
close all;
%% lat lon
load('../../../../data_processing/monthly_data_3000_v2.mat');

load('../isUS.mat');


res_v = 0.125;
res_h = 0.125;
lon = (-125+res_h/2):res_h: (-104-res_h/2);
lat = (50-res_v/2):-res_v: (32 + res_v/2);
[lons,lats]=meshgrid(lon,lat);


indexs  = repmat(0:18,3,1);

GridNames = {'Winter', 'Spring', 'Summer', 'Autumn'};


%% plot
colors = flipud(brewermap(10, 'Spectral'));
colors_2 =  flipud(brewermap(11, 'RdBu'));

figure;
set(gcf,'unit','normalized','position',[0.1,0.1,0.8,0.42]);
set(gca, 'Position', [0 0 1 1])

for season_i = 1
    switch season_i
        case 1
            seasons_all = [2];
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
        ELM_fsnos_tmp = nanmean(ELM_fsnos(:,:,filters),3);
        MODSCAG_fsnos_tmp = nanmean(MODSCAG_fsnos(:,:,filters),3);
        SPIRES_fsnos_tmp = nanmean(SPIRES_fsnos(:,:,filters),3);
        
        ELM_fsnos_tmp(isUS<1 | ELM_fsnos_tmp<=0) = nan;
        MODSCAG_fsnos_tmp(isUS<1 | MODSCAG_fsnos_tmp<=0) = nan;
        SPIRES_fsnos_tmp(isUS<1 | SPIRES_fsnos_tmp<=0) = nan;
        
          
        
        
        min_values = [0 0 0];
        max_values = [1 1 1];
        
        %% plot 1
        ax1 = subplot('position', [0.03 + 0.215*(season_i-1) 0.08 0.28 0.83]);
        colormap(ax1, colors);
        hold on
        if(season_i==1)
            plot_global_map_other3(lats, lons, ELM_fsnos_tmp, min_values(1), max_values(1),'(a) ELM',1, 1);
        else
            plot_global_map_other3(lats, lons, ELM_fsnos_tmp, min_values(1), max_values(1),'(a) ELM',0, 0);
        end
        
       
        
        if(season_i==4)
            hcb = colorbar;
            hcb.Title.String = 'Unitless';
            hcb.Title.FontSize = 7;
            hcb.Title.FontWeight = 'Bold';
            x1=get(gca,'position');
            x=get(hcb,'Position');
            x(3)=0.012;
            x(4)=0.23;
            x(1)=0.945;
            x(2)=0.68;
            set(hcb,'Position',x)
        end
        set(gca,'fontsize',14,'fontname','time new roman')
        
        ax2 = subplot('position', [0.33 + 0.215*(season_i-1) 0.08 0.28 0.83]);
        colormap(ax2, colors);
        hold on
        if(season_i==1)
            plot_global_map_other3(lats, lons, MODSCAG_fsnos_tmp, min_values(2), max_values(2),'(b) STC-MODSCAG',0, 1);
        else
            plot_global_map_other3(lats, lons, MODSCAG_fsnos_tmp, min_values(2), max_values(2),'(b) STC-MODSCAG',0, 0);
        end
        
        
        if(season_i==4)
            hcb = colorbar;
            hcb.Title.String = 'Unitless';
            hcb.Title.FontSize = 7;
            hcb.Title.FontWeight = 'Bold';
            x1=get(gca,'position');
            x=get(hcb,'Position');
            x(3)=0.012;
            x(4)=0.23;
            x(1)=0.945;
            x(2)=0.38;
            set(hcb,'Position',x)
        end
         set(gca,'fontsize',14,'fontname','time new roman')

        ax3 = subplot('position', [0.63 + 0.215*(season_i-1) 0.08 0.28 0.83]);
        colormap(ax3, colors);
        hold on
        if(season_i==1)
            plot_global_map_other3(lats, lons, SPIRES_fsnos_tmp, min_values(3), max_values(3),"(c) SPIReS",0, 1);
        else
            plot_global_map_other3(lats, lons, SPIRES_fsnos_tmp, min_values(3), max_values(3),"(c) SPIReS",0, 1);
        end
        
       
        if(season_i==1)
            hcb = colorbar;
            hcb.Title.String = 'Unitless';
            hcb.Title.FontSize = 12;
            hcb.Title.FontWeight = 'Bold';
            x1=get(gca,'position');
            x=get(hcb,'Position');
            x(3)=0.02;
            x(4)=0.8;
            x(1)=0.93;
            x(2)=0.08;
            set(hcb,'Position',x)
        end
        set(gca,'fontsize',14,'fontname','time new roman')
        
end

print(gcf, '-dtiff', '-r300', '../../figure_all_tif/Feb_spaital_pattern_fsno.tif')

close all

