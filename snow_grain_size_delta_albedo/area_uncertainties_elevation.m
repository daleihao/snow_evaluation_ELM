
clc;
clear all;
close all;
%% lat lon
load('../../../../data_processing/monthly_data_3000_v2.mat');
load('../isUS.mat');

colors =  [0.45, 0.80, 0.69;...
    0.98, 0.40, 0.35;...
    0.55, 0.60, 0.79];

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

ELM_deltaAlbedos_yearly = nan(19,2);
MODSCAG_deltaAlbedos_yearly = nan(19,2);
SPIRES_deltaAlbedos_yearly = nan(19,2);
ELM_deltaAlbedos_yearly_std = nan(19,2);
MODSCAG_deltaAlbedos_yearly_std = nan(19,2);
SPIRES_deltaAlbedos_yearly_std = nan(19,2);


for year_i = 2001:2019
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
        
        filters = (year_i-2001) * 12 + seasons_all';
        %filters = filters(filters>36); %2004
        ELM_deltaAlbedos_tmp = nanmean(ELM_deltaAlbedos(:,:,filters),3);
        MODSCAG_deltaAlbedos_tmp = nanmean(MODSCAG_deltaAlbedos(:,:,filters),3);
        SPIRES_deltaAlbedos_tmp = nanmean(SPIRES_deltaAlbedos(:,:,filters),3);
        
        ELM_deltaAlbedos_tmp(isUS<1) = nan;
        MODSCAG_deltaAlbedos_tmp(isUS<1) = nan;
        SPIRES_deltaAlbedos_tmp(isUS<1) = nan;
        
        
        filters = ELM_deltaAlbedos_tmp>=0 & MODSCAG_deltaAlbedos_tmp>=0 & SPIRES_deltaAlbedos_tmp>=0;
        area_weight = Areas(filters)./sum(Areas(filters));
        
        ELM_deltaAlbedos_yearly(year_i-2000,season_i)= nansum(ELM_deltaAlbedos_tmp(filters).*area_weight);
        MODSCAG_deltaAlbedos_yearly(year_i-2000,season_i) = nansum(MODSCAG_deltaAlbedos_tmp(filters).*area_weight);
        SPIRES_deltaAlbedos_yearly(year_i-2000,season_i) = nansum(SPIRES_deltaAlbedos_tmp(filters).*area_weight);
        
        ELM_deltaAlbedos_yearly_std(year_i-2000,season_i)= std(ELM_deltaAlbedos_tmp(filters),area_weight,1,'omitnan');
        MODSCAG_deltaAlbedos_yearly_std(year_i-2000,season_i) = std(MODSCAG_deltaAlbedos_tmp(filters),area_weight,1,'omitnan');
        SPIRES_deltaAlbedos_yearly_std(year_i-2000,season_i) = std(SPIRES_deltaAlbedos_tmp(filters),area_weight,1,'omitnan');
        
    end
end

%% elevation distribution
elevations = imread('../WUS_mean_of_Elevation.tif');
elevations = squeeze(elevations(:,:,1));
indexs  = repmat(0:18,3,1);

ELM_fsnos_elevation = nan(7,2);
MODSCAG_fsnos_elevation = nan(7,2);
SPIRES_fsnos_elevation = nan(7,2);

ELM_fsnos_elevation_std = nan(7,2);
MODSCAG_fsnos_elevation_std = nan(7,2);
SPIRES_fsnos_elevation_std = nan(7,2);

ELM_deltaAlbedos_elevation = nan(7,2);
MODSCAG_deltaAlbedos_elevation = nan(7,2);
SPIRES_deltaAlbedos_elevation = nan(7,2);

ELM_deltaAlbedos_elevation_std = nan(7,2);
MODSCAG_deltaAlbedos_elevation_std = nan(7,2);
SPIRES_deltaAlbedos_elevation_std = nan(7,2);

ELM_grainsizes_elevation = nan(7,2);
MODSCAG_grainsizes_elevation = nan(7,2);
SPIRES_grainsizes_elevation = nan(7,2);

ELM_grainsizes_elevation_std = nan(7,2);
MODSCAG_grainsizes_elevation_std = nan(7,2);
SPIRES_grainsizes_elevation_std = nan(7,2);

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
    
    filters = indexs * 12 + seasons_all';
    filters = filters(:);
    %filters = filters(filters>36); %2004
    ELM_fsnos_tmp = nanmean(ELM_fsnos(:,:,filters),3);
    MODSCAG_fsnos_tmp = nanmean(MODSCAG_fsnos(:,:,filters),3);
    SPIRES_fsnos_tmp = nanmean(SPIRES_fsnos(:,:,filters),3);
    
    filters = ELM_fsnos_tmp<=0 | MODSCAG_fsnos_tmp<=0 | SPIRES_fsnos_tmp<=0 | ...
        isnan(ELM_fsnos_tmp) | isnan(MODSCAG_fsnos_tmp)  | isnan(SPIRES_fsnos_tmp);
    
    ELM_fsnos_tmp(isUS<1 | filters) = nan;
    MODSCAG_fsnos_tmp(isUS<1 | filters) = nan;
    SPIRES_fsnos_tmp(isUS<1 | filters) = nan;
    
    for elevation_band = 1:7
        filters = elevations >=((elevation_band-1)*500)  & elevations<(elevation_band*500);
        if elevation_band == 7
            filters = elevations >=((elevation_band-1)*500) ;
        end
        area_weight = Areas(filters)./sum(Areas(filters));
        
        ELM_fsnos_elevation(elevation_band,season_i) = nansum(ELM_fsnos_tmp(filters).*area_weight);
        MODSCAG_fsnos_elevation(elevation_band,season_i) = nansum(MODSCAG_fsnos_tmp(filters).*area_weight);
        SPIRES_fsnos_elevation(elevation_band,season_i) = nansum(SPIRES_fsnos_tmp(filters).*area_weight);
        
        ELM_fsnos_elevation_std(elevation_band,season_i) = std(ELM_fsnos_tmp(filters),area_weight,1,'omitnan');
        MODSCAG_fsnos_elevation_std(elevation_band,season_i) = std(MODSCAG_fsnos_tmp(filters),area_weight,1,'omitnan');
        SPIRES_fsnos_elevation_std(elevation_band,season_i) = std(SPIRES_fsnos_tmp(filters),area_weight,1,'omitnan');
    end
    
    filters = indexs * 12 + seasons_all';
    filters = filters(:);
    ELM_grainsizes_tmp = nanmean(ELM_grainsizes(:,:,filters),3);
    MODSCAG_grainsizes_tmp = nanmean(MODSCAG_grainsizes(:,:,filters),3);
    SPIRES_grainsizes_tmp = nanmean(SPIRES_grainsizes(:,:,filters),3);
    
    filters = ELM_grainsizes_tmp<0 | MODSCAG_grainsizes_tmp<0 | SPIRES_grainsizes_tmp<0 | ...
        isnan(ELM_grainsizes_tmp) | isnan(MODSCAG_grainsizes_tmp)  | isnan(SPIRES_grainsizes_tmp) | ...
        ELM_fsnos_tmp<=0 | MODSCAG_fsnos_tmp<=0 | SPIRES_fsnos_tmp<=0 | ...
        isnan(ELM_fsnos_tmp) | isnan(MODSCAG_fsnos_tmp)  | isnan(SPIRES_fsnos_tmp);
    
    ELM_grainsizes_tmp(isUS<1 | filters) = nan;
    MODSCAG_grainsizes_tmp(isUS<1 | filters) = nan;
    SPIRES_grainsizes_tmp(isUS<1 | filters) = nan;
    
    for elevation_band = 1:7
        filters = elevations >=((elevation_band-1)*500)  & elevations<(elevation_band*500);
        if elevation_band == 7
            filters = elevations >=((elevation_band-1)*500) ;
        end
        area_weight = Areas(filters)./sum(Areas(filters));
        
        ELM_grainsizes_elevation(elevation_band,season_i) = nansum(ELM_grainsizes_tmp(filters).*area_weight);
        MODSCAG_grainsizes_elevation(elevation_band,season_i) = nansum(MODSCAG_grainsizes_tmp(filters).*area_weight);
        SPIRES_grainsizes_elevation(elevation_band,season_i) = nansum(SPIRES_grainsizes_tmp(filters).*area_weight);
        
        ELM_grainsizes_elevation_std(elevation_band,season_i) = std(ELM_grainsizes_tmp(filters),area_weight,1,'omitnan');
        MODSCAG_grainsizes_elevation_std(elevation_band,season_i) = std(MODSCAG_grainsizes_tmp(filters),area_weight,1,'omitnan');
        SPIRES_grainsizes_elevation_std(elevation_band,season_i) = std(SPIRES_grainsizes_tmp(filters),area_weight,1,'omitnan');
    end
    
    filters = indexs * 12 + seasons_all';
    filters = filters(:);
    ELM_deltaAlbedos_tmp = nanmean(ELM_deltaAlbedos(:,:,filters),3);
    MODSCAG_deltaAlbedos_tmp = nanmean(MODSCAG_deltaAlbedos(:,:,filters),3);
    SPIRES_deltaAlbedos_tmp = nanmean(SPIRES_deltaAlbedos(:,:,filters),3);
    
    filters = ELM_deltaAlbedos_tmp<0 | MODSCAG_deltaAlbedos_tmp<0 | SPIRES_deltaAlbedos_tmp<0 | ...
        isnan(ELM_deltaAlbedos_tmp) | isnan(MODSCAG_deltaAlbedos_tmp)  | isnan(SPIRES_deltaAlbedos_tmp) | ...
        ELM_fsnos_tmp<=0 | MODSCAG_fsnos_tmp<=0 | SPIRES_fsnos_tmp<=0 | ...
        isnan(ELM_fsnos_tmp) | isnan(MODSCAG_fsnos_tmp)  | isnan(SPIRES_fsnos_tmp);
    
    ELM_deltaAlbedos_tmp(isUS<1 | filters) = nan;
    MODSCAG_deltaAlbedos_tmp(isUS<1 | filters) = nan;
    SPIRES_deltaAlbedos_tmp(isUS<1 | filters) = nan;
    
    for elevation_band = 1:7
        filters = elevations >=((elevation_band-1)*500)  & elevations<(elevation_band*500);
        if elevation_band == 7
            filters = elevations >=((elevation_band-1)*500) ;
        end
        area_weight = Areas(filters)./sum(Areas(filters));
        
        ELM_deltaAlbedos_elevation(elevation_band,season_i) = nansum(ELM_deltaAlbedos_tmp(filters).*area_weight);
        MODSCAG_deltaAlbedos_elevation(elevation_band,season_i) = nansum(MODSCAG_deltaAlbedos_tmp(filters).*area_weight);
        SPIRES_deltaAlbedos_elevation(elevation_band,season_i) = nansum(SPIRES_deltaAlbedos_tmp(filters).*area_weight);
        
        %         ELM_deltaAlbedos_elevation(elevation_band,season_i) = nanmean(ELM_deltaAlbedos_tmp(filters));
        %         MODSCAG_deltaAlbedos_elevation(elevation_band,season_i) = nanmean(MODSCAG_deltaAlbedos_tmp(filters));
        %         SPIRES_deltaAlbedos_elevation(elevation_band,season_i) = nanmean(SPIRES_deltaAlbedos_tmp(filters));
        
        
        ELM_deltaAlbedos_elevation_std(elevation_band,season_i) = std(ELM_deltaAlbedos_tmp(filters),area_weight,1,'omitnan');
        MODSCAG_deltaAlbedos_elevation_std(elevation_band,season_i) = std(MODSCAG_deltaAlbedos_tmp(filters),area_weight,1,'omitnan');
        SPIRES_deltaAlbedos_elevation_std(elevation_band,season_i) = std(SPIRES_deltaAlbedos_tmp(filters),area_weight,1,'omitnan');
    end
end

figure
set(gcf,'unit','normalized','position',[0.1,0.1,0.6,0.8]);
set(gca, 'Position', [0 0 1 1])


for season_i = 1:2
    %% plot 1
    subplot('position', [0.09+(season_i-1)*0.455 0.69 0.44 0.275]);
    hold on
    plot(1*[1:6],ELM_fsnos_elevation(2:end,season_i),'ko','MarkerSize',6,...
        'MarkerEdgeColor',colors(1,:),...
        'MarkerFaceColor',colors(1,:))
    if(season_i==1)
        errorbar([1:6]+0.2,SPIRES_fsnos_elevation(2:end,season_i),5.1/100/2*ones(6,1),'s','MarkerSize',6,...
            'MarkerEdgeColor',colors(3,:),'MarkerFaceColor',colors(3,:),'color',colors(3,:),'CapSize',10,'linewidth',1)
    else
        errorbar([1:6]-0.2,MODSCAG_fsnos_elevation(2:end,season_i),5.8/100/2*ones(6,1),'s','MarkerSize',6,...
            'MarkerEdgeColor',colors(2,:),'MarkerFaceColor',colors(2,:),'color',colors(2,:),'CapSize',10,'linewidth',1)
        errorbar([1:6]+0.2,SPIRES_fsnos_elevation(2:end,season_i),5.1/100/2*ones(6,1),'s','MarkerSize',6,...
            'MarkerEdgeColor',colors(3,:),'MarkerFaceColor',colors(3,:),'color',colors(3,:),'CapSize',10,'linewidth',1)
    end
    
    box on
    if(season_i==1)
        ylabel('\itf\rm_{sno}','fontsize',10)
    else
        set(gca,'yticklabel',[])
        
    end
    if(season_i==1)
        text(6, -0.05+1.05*0.1,'(a)','fontsize',12,'fontweight','bold')
        
        title('Winter')
    else
        text(6, -0.05+1.05*0.1,'(b)','fontsize',12,'fontweight','bold')
        legend('ELM','STC-MODSCAG/STC-MODDRFS','SPIReS','fontsize',8,'location','northwest')
        title('Spring')
    end
    
    xlim([0.5 6.5])
    ylim([-0.05 1])
    set(gca, 'linewidth',1,'fontsize',10,'xtick', [1:7],'xticklabel',[])
    
    %% plot 2
    subplot('position', [0.09+(season_i-1)*0.455 0.39 0.44 0.275]);
    hold on
    
    plot(1*[1:6],ELM_grainsizes_elevation(2:end,season_i),'ko','MarkerSize',6,...
        'MarkerEdgeColor',colors(1,:),...
        'MarkerFaceColor',colors(1,:))
    if(season_i==1)
        errorbar([1:6]+0.2,SPIRES_grainsizes_elevation(2:end,season_i),118/2*ones(6,1),'s','MarkerSize',6,...
            'MarkerEdgeColor',colors(3,:),'MarkerFaceColor',colors(3,:),'color',colors(3,:),'CapSize',10,'linewidth',1)
    else
        
        errorbar([1:6]-0.2,MODSCAG_grainsizes_elevation(2:end,season_i),118/2*ones(6,1),'s','MarkerSize',6,...
            'MarkerEdgeColor',colors(2,:),'MarkerFaceColor',colors(2,:),'color',colors(2,:),'CapSize',10,'linewidth',1)
        errorbar([1:6]+0.2,SPIRES_grainsizes_elevation(2:end,season_i),118/2*ones(6,1),'s','MarkerSize',6,...
            'MarkerEdgeColor',colors(3,:),'MarkerFaceColor',colors(3,:),'color',colors(3,:),'CapSize',10,'linewidth',1)
    end
    
    
    if(season_i==1)
        text(6, -100+800*0.1,'(c)','fontsize',12,'fontweight','bold')
    else
        text(6, -100+800*0.1,'(d)','fontsize',12,'fontweight','bold')
    end
    
    if(season_i==1)
        ylabel('\itS\rm_{sno} (\mum)','fontsize',10)
    else
        set(gca,'yticklabel',[])
    end
    box on
    xlim([0.5 6.5])
    ylim([-100 600])
    set(gca, 'linewidth',1,'fontsize',10,'xtick', [1:7],'xticklabel',{})
    
    
    %% plot 3
    subplot('position', [0.09+(season_i-1)*0.455 0.09 0.44 0.275]);
    hold on
    
    plot(1*[1:6],ELM_deltaAlbedos_elevation(2:end,season_i),'ko','MarkerSize',6,...
        'MarkerEdgeColor',colors(1,:),...
        'MarkerFaceColor',colors(1,:))
    if(season_i==1)
        errorbar([1:6]+0.2,SPIRES_deltaAlbedos_elevation(2:end,season_i),3.6/100/2*0.63*ones(6,1),'s','MarkerSize',6,...
            'MarkerEdgeColor',colors(3,:),'MarkerFaceColor',colors(3,:),'color',colors(3,:),'CapSize',10,'linewidth',1)
    else
        
        errorbar([1:6]-0.2,MODSCAG_deltaAlbedos_elevation(2:end,season_i),3.6/100/2*0.63*ones(6,1),'s','MarkerSize',6,...
            'MarkerEdgeColor',colors(2,:),'MarkerFaceColor',colors(2,:),'color',colors(2,:),'CapSize',10,'linewidth',1)
        errorbar([1:6]+0.2,SPIRES_deltaAlbedos_elevation(2:end,season_i),3.6/100/2*0.63*ones(6,1),'s','MarkerSize',6,...
            'MarkerEdgeColor',colors(3,:),'MarkerFaceColor',colors(3,:),'color',colors(3,:),'CapSize',10,'linewidth',1)
        
    end
    if(season_i==1)
        text(6, -0.02+0.07*0.1,'(e)','fontsize',12,'fontweight','bold')
    else
        text(6, -0.02+0.07*0.1,'(f)','fontsize',12,'fontweight','bold')
    end
    
    if(season_i==1)
        ylabel('\itR\rm_{sno}','fontsize',10)
    else
        set(gca,'yticklabel',[])
        
    end
    box on
    xlim([0.5 6.5])
    ylim([-0.02,0.05])
    
    set(gca, 'linewidth',1,'fontsize',10,'xtick', [1:7],'xticklabel',{'0.5-1.0', '1.0-1.5', '1.5-2.0', '2.0-2.5' , '2.5-3.0', '\geq3.0'})
    xlabel('Elevation (km)','fontsize',12)
    
    
end

print(gcf, '-dtiff', '-r300', '../../figure_all_tif/uncertainties_deltaalbedo.tif')

close all

