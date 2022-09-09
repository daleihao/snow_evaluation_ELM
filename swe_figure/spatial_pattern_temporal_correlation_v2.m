
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

radius = 6378137;
height = radius * 0.125*pi/180;
width1 = radius*cos(lats*pi/180) * 0.125*pi/180;
width2 = radius*(cos((lats-0.125/2)*pi/180)+cos((lats+0.125/2)*pi/180))/2 * 0.125*pi/180;
Areas = width2.*height/1e6;

indexs  = repmat(0:18,3,1);

GridNames = {'Winter', 'Spring', 'Summer', 'Autumn'};


%% plot
colors_1 = flipud(brewermap(10, 'Spectral'));
colors_2 =  flipud(brewermap(11, 'RdBu'));
colors_3 =  flipud((brewermap(10, 'Spectral')));
colors_2(6,:) = [0.85 0.85 0.85];

figure;
set(gcf,'unit','normalized','position',[0.1,0.1,0.52,0.55]);
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
    ELM_swes_tmp = nanmean(ELM_swes(:,:,filters),3);
    UA_swes_tmp = nanmean(UA_swes(:,:,filters),3);
    SNODAS_swes_tmp = nanmean(SNODAS_swes(:,:,filters),3);
    
    ELM_swes_tmp(isUS<1) = nan;
    UA_swes_tmp(isUS<1) = nan;
    SNODAS_swes_tmp(isUS<1) = nan;
    
    ELM_fsnos_tmp = nanmean(ELM_fsnos(:,:,filters),3);
    ELM_fsnos_tmp = nanmean(ELM_fsnos(:,:,filters),3);
    MODSCAG_fsnos_tmp = nanmean(MODSCAG_fsnos(:,:,filters),3);
    SPIRES_fsnos_tmp = nanmean(SPIRES_fsnos(:,:,filters),3);
    
   ELM_fsnos_tmp(isUS<1) = nan;
    MODSCAG_fsnos_tmp(isUS<1) = nan;
    SPIRES_fsnos_tmp(isUS<1) = nan;

    
    min_values = [0 -200 -200 0 0];
    max_values = [600 200 200 1 1];
    %%spatial_pattern
    %% plot 1
       filters = isnan(ELM_swes_tmp) | isnan(UA_swes_tmp) | isnan(SNODAS_swes_tmp) |...
       ELM_swes_tmp<=0 | UA_swes_tmp<=0 | SNODAS_swes_tmp<=0 |...
          ELM_fsnos_tmp<=0 | MODSCAG_fsnos_tmp<=0 | SPIRES_fsnos_tmp<=0;

    ELM_swes_tmp(filters) = nan;
    UA_swes_tmp(filters) = nan;
    SNODAS_swes_tmp(filters) = nan;
         Areas_tmp = Areas;

                alloutputs = nan(2,5);
    alloutputs(1,:) = calStat(ELM_swes_tmp,UA_swes_tmp,Areas);
    alloutputs(2,:) = calStat(ELM_swes_tmp,SNODAS_swes_tmp,Areas);

    ax1 = subplot('position', [0.055 + 0.195*(season_i-1) 0.65 0.19 0.28]);
    
        Areas_tmp = Areas;
     Areas_tmp(isnan(ELM_swes_tmp)) = nan;
     Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));

    colormap(ax1, colors_1);
    hold on
    if(season_i==1)
        plot_global_map_other(lats, lons, ELM_swes_tmp, min_values(1), max_values(1),"(a)",1, 0);
        m_text(-110,33.2,num2str(nansum(ELM_swes_tmp(:).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
    else
        plot_global_map_other(lats, lons, ELM_swes_tmp, min_values(1), max_values(1),"(b)",0, 0);
        m_text(-110,33.2,num2str(nansum(ELM_swes_tmp(:).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
        
    end
    set(gca,'fontsize',8,'fontname','time new roman')
    
    if(season_i==1)
        ylabel({'ELM'})
        title("Winter",'fontsize',10, 'fontweight', 'bold')
    end
    
    if(season_i==2)
        title("Spring",'fontsize',10, 'fontweight', 'bold')
        
        hcb = colorbar;
        hcb.Title.String = 'mm';
        hcb.Title.FontSize = 6;
        hcb.Title.FontWeight = 'Bold';
        x1=get(gca,'position');
        x=get(hcb,'Position');
        x(3)=0.012;
        x(4)=0.23;
        x(1)=0.46;
        x(2)=0.66;
        set(hcb,'Position',x)
    end
    
    
    
    ax2 = subplot('position', [0.055 + 0.195*(season_i-1) 0.36 0.19 0.28]);
    colormap(ax2, colors_2);
    
            Areas_tmp = Areas;
     Areas_tmp(isnan(ELM_swes_tmp)) = nan;
     Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));

    hold on
    if(season_i==1)
        plot_global_map_other(lats, lons, ELM_swes_tmp- UA_swes_tmp, min_values(2), max_values(2),"(c)",1, 0);
        m_text(-110,33.2,num2str(nansum((ELM_swes_tmp(:)- UA_swes_tmp(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
        
    else
        plot_global_map_other(lats, lons, ELM_swes_tmp - UA_swes_tmp, min_values(2), max_values(2),"(d)",0, 0);
        m_text(-110,33.2,num2str(nansum((ELM_swes_tmp(:)- UA_swes_tmp(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
        
    end
    set(gca,'fontsize',8,'fontname','time new roman')
    
    if(season_i==1)
        ylabel({'ELM - UA'})
    end
    if(season_i==2)
        hcb = colorbar;
        hcb.Title.String = 'mm';
        hcb.Title.FontSize = 6;
        hcb.Title.FontWeight = 'Bold';
        x1=get(gca,'position');
        x=get(hcb,'Position');
        x(3)=0.012;
        x(4)=0.23;
        x(1)=0.46;
        x(2)=0.37;
        set(hcb,'Position',x)
    end
    
    ax3 = subplot('position', [0.055 + 0.195*(season_i-1) 0.07 0.19 0.28]);
    colormap(ax3, colors_2);
    hold on
    
            Areas_tmp = Areas;
     Areas_tmp(isnan(ELM_swes_tmp)) = nan;
     Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));

    if(season_i==1)
        plot_global_map_other(lats, lons, ELM_swes_tmp - SNODAS_swes_tmp, min_values(3), max_values(3),"(g)",1, 1);
        m_text(-110,33.2,num2str(nansum((ELM_swes_tmp(:)- SNODAS_swes_tmp(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
        
    else
        plot_global_map_other(lats, lons, ELM_swes_tmp - SNODAS_swes_tmp, min_values(3), max_values(3),"(h)",0, 1);
        m_text(-110,33.2,num2str(nansum((ELM_swes_tmp(:)- SNODAS_swes_tmp(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
        
    end
    set(gca,'fontsize',8,'fontname','time new roman')
    
    if(season_i==1)
        ylabel({'ELM - SNODAS'})
    end
    if(season_i==2)
        hcb = colorbar;
        hcb.Title.String = 'mm';
        hcb.Title.FontSize = 6;
        hcb.Title.FontWeight = 'Bold';
        x1=get(gca,'position');
        x=get(hcb,'Position');
        x(3)=0.012;
        x(4)=0.23;
        x(1)=0.46;
        x(2)=0.08;
        set(hcb,'Position',x)
    end
    
    %% plot 3
    load('swe_temporal_R.mat');
    
    Rtmp = squeeze(Rs_1(:,:,season_i));
    
    Rtmp(isUS<1) = nan;
    
    Areas_tmp = Areas;
     Areas_tmp(isnan(Rtmp)) = nan;
     Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));
     
    ax4 = subplot('position', [0.16 + 0.195*(season_i-1+2) 0.36 0.19 0.28]);
    colormap(ax4, colors_3);
    if(season_i==1)
        plot_global_map_other(lats, lons, Rtmp, min_values(4), max_values(4),"(e)",0, 0);
    else
        plot_global_map_other(lats, lons, Rtmp, min_values(4), max_values(4),"(f)",0, 0);
    end
    
    m_text(-109,33.2,num2str(nansum((Rtmp(:)).*Areas_tmp(:)),'%4.2f'),'fontsize',10, 'fontweight', 'bold');
    
    if(season_i==1)
        title("Winter",'fontsize',10, 'fontweight', 'bold')
        
        label_y = ylabel({'R_{ELM & UA}'});
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
        x(4)=0.23;
        x(1)=0.955;
        x(2)=0.37;
        set(hcb,'Position',x)
    end
    
    load('swe_temporal_R.mat');
    ax5 = subplot('position', [0.16 + 0.195*(season_i-1+2) 0.07 0.19 0.28]);
    
    Rtmp = squeeze(Rs_2(:,:,season_i));
    
    Rtmp(isUS<1) = nan;
    
        Areas_tmp = Areas;
     Areas_tmp(isnan(Rtmp)) = nan;
     Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));

    colormap(ax5, colors_3);
    if(season_i==1)
        plot_global_map_other(lats, lons, Rtmp, min_values(5), max_values(5),"(i)",0, 1);
    else
        plot_global_map_other(lats, lons, Rtmp, min_values(5), max_values(5),"(j)",0, 1);
    end
    
    m_text(-109,33.2,num2str(nansum((Rtmp(:)).*Areas_tmp(:)),'%4.2f'),'fontsize',10, 'fontweight', 'bold');
    
    if(season_i==1)
        label_y = ylabel({'R_{ELM & SNODAS}'});
        label_y.Position(1) = -0.20; % change horizontal position of ylabel
    end
    
    if(season_i==2)
        hcb = colorbar;
        hcb.Title.String = 'Unitless';
        hcb.Title.FontSize = 6;
        hcb.Title.FontWeight = 'Bold';
        x1=get(gca,'position');
        x=get(hcb,'Position');
        x(3)=0.012;
        x(4)=0.23;
        x(1)=0.955;
        x(2)=0.08;
        set(hcb,'Position',x)
    end
    
end

print(gcf, '-dtiff', '-r300', '../../figure_all_tif/swe_spaital_pattern_temporal_correlation.tif')

close all

