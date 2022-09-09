
clc;
clear all;
close all;
%% lat lon

load('../isUS.mat');


res_v = 0.125;
res_h = 0.125;
lon = (-125+res_h/2):res_h: (-104-res_h/2);
lat = (50-res_v/2):-res_v: (32 + res_v/2);
[lons,lats]=meshgrid(lon,lat);

indexs  = repmat(0:18,3,1);

radius = 6378137;
height = radius * 0.125*pi/180;
width1 = radius*cos(lats*pi/180) * 0.125*pi/180;
width2 = radius*(cos((lats-0.125/2)*pi/180)+cos((lats+0.125/2)*pi/180))/2 * 0.125*pi/180;
Areas = width2.*height/1e6;

GridNames = {'Winter', 'Spring', 'Summer', 'Autumn'};


%% plot
colors_1 = flipud(brewermap(10, 'Spectral'));
colors_2 =  flipud(brewermap(11, 'RdBu'));
colors_3 =  flipud((brewermap(10, 'Spectral')));
colors_2(6,:) = [0.85 0.85 0.85];

figure;
set(gcf,'unit','normalized','position',[0.1,0.1,0.52,0.4]);
set(gca, 'Position', [0 0 1 1])

for season_i = 1:2
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
    
    
    load('../../../../data_processing/monthly_data_3000.mat');
    
    filters = indexs * 12 + seasons_all';
    filters = filters(:);
    %filters = filters(filters>36); %2004
    ELM_surfacealbedos_tmp = nanmean(ELM_Albedos(:,:,filters),3);
    
    MCD43_surfacealbedo = MCD43_BSAs.* (1 - ELM_Skyls) + MCD43_WSAs .* ELM_Skyls;
    MCD43_surfacealbedos_tmp = nanmean(MCD43_surfacealbedo(:,:,filters),3)/1000;
        
    
    ELM_surfacealbedos_tmp(isUS<1 | ELM_surfacealbedos_tmp<=0) = nan;
    MCD43_surfacealbedos_tmp(isUS<1| ELM_surfacealbedos_tmp<=0) = nan;
    
    filters = indexs * 12 + seasons_all';
    filters = filters(:);
    %filters = filters(filters>36); %2004
    ELM_fsnos_tmp = nanmean(ELM_fsnos(:,:,filters),3);
    MODSCAG_fsnos_tmp = nanmean(MODSCAG_fsnos(:,:,filters),3);
    SPIRES_fsnos_tmp = nanmean(SPIRES_fsnos(:,:,filters),3);
    
    ELM_fsnos_tmp(isUS<1) = nan;
    MODSCAG_fsnos_tmp(isUS<1) = nan;
    SPIRES_fsnos_tmp(isUS<1) = nan;
    

    min_values = [0 -0.2 -0.2 0 0];
    max_values = [0.8 0.2 0.2 1 1];
    %%spatial_pattern
    %% plot 1
    
     filters = ELM_fsnos_tmp<=0 | MODSCAG_fsnos_tmp<=0 | SPIRES_fsnos_tmp<=0 |...
        isnan(ELM_fsnos_tmp) | isnan(MODSCAG_fsnos_tmp) | isnan(SPIRES_fsnos_tmp) | ...
         ELM_surfacealbedos_tmp<=0 | MCD43_surfacealbedos_tmp<=0 | ...
        isnan(ELM_surfacealbedos_tmp) | isnan(MCD43_surfacealbedos_tmp);
    ELM_surfacealbedos_tmp(filters) = nan;
    MCD43_surfacealbedos_tmp(filters) = nan;
         Areas_tmp = Areas;

            alloutputs = nan(1,5);
    alloutputs(1,:) = calStat(ELM_surfacealbedos_tmp,MCD43_surfacealbedos_tmp,Areas_tmp);

    
     Areas_tmp = Areas;
     Areas_tmp(isnan(ELM_surfacealbedos_tmp)) = nan;
     Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));

     
    ax1 = subplot('position', [0.055 + 0.195*(season_i-1) 0.5 0.19 0.4]);
    colormap(ax1, colors_1);
    hold on
    if(season_i==1)
        plot_global_map_other(lats, lons, ELM_surfacealbedos_tmp, min_values(1), max_values(1),"(a)",1, 0);
    else
        plot_global_map_other(lats, lons, ELM_surfacealbedos_tmp, min_values(1), max_values(1),"(b)",0, 0);
    end
    
    m_text(-110,33.2,num2str(nansum(abs(ELM_surfacealbedos_tmp(:)).*Areas_tmp(:)),'%4.2f'),'fontsize',10, 'fontweight', 'bold');
    
    set(gca,'fontsize',8,'fontname','time new roman')
    
    if(season_i==1)
        ylabel({'ELM'})
        title("Winter",'fontsize',10, 'fontweight', 'bold')
    end
    
    if(season_i==2)
        title("Spring",'fontsize',10, 'fontweight', 'bold')
        
        hcb = colorbar;
        hcb.Title.String = 'Unitless';
        hcb.Title.FontSize = 6;
        hcb.Title.FontWeight = 'Bold';
        x1=get(gca,'position');
        x=get(hcb,'Position');
        x(3)=0.012;
        x(4)=0.36;
        x(1)=0.46;
        x(2)=0.51;
        set(hcb,'Position',x)
    end
    
    ax2 = subplot('position', [0.055 + 0.195*(season_i-1) 0.09 0.19 0.4]);
    colormap(ax2, colors_2);
    hold on
    if(season_i==1)
        plot_global_map_other(lats, lons, ELM_surfacealbedos_tmp- MCD43_surfacealbedos_tmp, min_values(2), max_values(2),"(c)",1, 1);
    else
        plot_global_map_other(lats, lons, ELM_surfacealbedos_tmp - MCD43_surfacealbedos_tmp, min_values(2), max_values(2),"(d)",0, 1);
    end
    m_text(-110,33.2,num2str(nansum((ELM_surfacealbedos_tmp(:)- MCD43_surfacealbedos_tmp(:)).*Areas_tmp(:)),'%4.2f'),'fontsize',10, 'fontweight', 'bold');

    set(gca,'fontsize',8,'fontname','time new roman')
    
    if(season_i==1)
        ylabel({'ELM - MCD43'})
    end
    if(season_i==2)
        hcb = colorbar;
        hcb.Title.String = 'Unitless';
        hcb.Title.FontSize = 6;
        hcb.Title.FontWeight = 'Bold';
        x1=get(gca,'position');
        x=get(hcb,'Position');
        x(3)=0.012;
        x(4)=0.36;
        x(1)=0.46;
        x(2)=0.10;
        set(hcb,'Position',x)
    end
    
%     ax3 = subplot('position', [0.055 + 0.195*(season_i-1) 0.07 0.19 0.28]);
%     colormap(ax3, colors_2);
%     hold on
%     if(season_i==1)
%         plot_global_map_other(lats, lons, ELM_surfacealbedos_tmp - SPIRES_surfacealbedos_tmp, min_values(3), max_values(3),"(g)",1, 1);
%     else
%         plot_global_map_other(lats, lons, ELM_surfacealbedos_tmp - SPIRES_surfacealbedos_tmp, min_values(3), max_values(3),"(h)",0, 1);
%     end
%     set(gca,'fontsize',8,'fontname','time new roman')
%     
%     if(season_i==1)
%         ylabel({'ELM - SPIReS'})
%     end
%     if(season_i==2)
%         hcb = colorbar;
%         hcb.Title.String = 'Unitless';
%         hcb.Title.FontSize = 6;
%         hcb.Title.FontWeight = 'Bold';
%         x1=get(gca,'position');
%         x=get(hcb,'Position');
%         x(3)=0.012;
%         x(4)=0.23;
%         x(1)=0.46;
%         x(2)=0.08;
%         set(hcb,'Position',x)
%     end
    
    %% plot 3
    load('Albedo_temporal_R.mat');
     filters = ELM_fsnos_tmp<=0 | MODSCAG_fsnos_tmp<=0 | SPIRES_fsnos_tmp<=0 |...
        isnan(ELM_fsnos_tmp) | isnan(MODSCAG_fsnos_tmp) | isnan(SPIRES_fsnos_tmp);
    R_tmp = squeeze(Rs_1(:,:,season_i));
    R_tmp(isUS<1 | filters) = nan;
    
    Areas_tmp = Areas;
     Areas_tmp(isnan(R_tmp)) = nan;
     Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));
     
    ax4 = subplot('position', [0.16 + 0.195*(season_i-1+2) 0.09 0.19 0.4]);
    colormap(ax4, colors_3);
    if(season_i==1)
        plot_global_map_other(lats, lons, R_tmp, min_values(4), max_values(4),"(e)",0, 1);
    else
        plot_global_map_other(lats, lons, R_tmp, min_values(4), max_values(4),"(f)",0, 1);
    end
    
        m_text(-110,33.2,num2str(nansum((R_tmp(:).*Areas_tmp(:))),'%4.2f'),'fontsize',10, 'fontweight', 'bold');

    if(season_i==1)
        title("Winter",'fontsize',10, 'fontweight', 'bold')
        
        label_y = ylabel({'R_{ELM & MCD43}'});
        label_y.Position(1) = -0.20; % change horizontal position of ylabel
        
    end
    if(season_i==2)
        title("Spring",'fontsize',10, 'fontweight', 'bold')
        
        hcb = colorbar;
        hcb.Title.String = 'Unitless';
        hcb.Title.FontSize = 6;
        hcb.Title.FontWeight = 'Bold';
        x1=get(gca,'position');
        x=get(hcb,'Position');
        x(3)=0.012;
        x(4)=0.36;
        x(1)=0.955;
        x(2)=0.10;
        set(hcb,'Position',x)
    end
    
%     load('surfacealbedo_temporal_R.mat');
%     ax5 = subplot('position', [0.16 + 0.195*(season_i-1+2) 0.07 0.19 0.28]);
%     colormap(ax5, colors_3);
%     if(season_i==1)
%         plot_global_map_other(lats, lons, squeeze(Rs_2(:,:,season_i)), min_values(5), max_values(5),"(i)",0, 1);
%     else
%         plot_global_map_other(lats, lons, squeeze(Rs_2(:,:,season_i)), min_values(5), max_values(5),"(j)",0, 1);
%     end
%     if(season_i==1)
%         label_y = ylabel({'R_{ELM & SPIReS}'});
%         label_y.Position(1) = -0.20; % change horizontal position of ylabel
%     end
%     
%     if(season_i==2)
%         hcb = colorbar;
%         hcb.Title.String = 'Unitless';
%         hcb.Title.FontSize = 6;
%         hcb.Title.FontWeight = 'Bold';
%         x1=get(gca,'position');
%         x=get(hcb,'Position');
%         x(3)=0.012;
%         x(4)=0.23;
%         x(1)=0.955;
%         x(2)=0.08;
%         set(hcb,'Position',x)
%     end
    
end

print(gcf, '-dtiff', '-r300', '../../figure_all_tif/surfacealbedo_spaital_pattern_temporal_correlation.tif')

close all

