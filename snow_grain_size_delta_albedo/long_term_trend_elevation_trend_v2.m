
clc;
clear all;
close all;
%% lat lon
load('../../../../data_processing/monthly_data_3000_v2.mat');
load('../isUS.mat');

colors =  [0.45, 0.80, 0.69;...
    0.98, 0.40, 0.35;...
    0.55, 0.60, 0.79];
colors2 =  [0.45, 0.80, 0.69;...
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

ELM_grainsizes_yearly = nan(19,2);
MODSCAG_grainsizes_yearly = nan(19,2);
SPIRES_grainsizes_yearly = nan(19,2);
ELM_grainsizes_yearly_std = nan(19,2);
MODSCAG_grainsizes_yearly_std = nan(19,2);
SPIRES_grainsizes_yearly_std = nan(19,2);


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
        ELM_grainsizes_tmp = nanmean(ELM_grainsizes(:,:,filters),3);
        MODSCAG_grainsizes_tmp = nanmean(MODSCAG_grainsizes(:,:,filters),3);
        SPIRES_grainsizes_tmp = nanmean(SPIRES_grainsizes(:,:,filters),3);
        
        ELM_grainsizes_tmp(isUS<1) = nan;
        MODSCAG_grainsizes_tmp(isUS<1) = nan;
        SPIRES_grainsizes_tmp(isUS<1) = nan;
        
        ELM_fsnos_tmp = nanmean(ELM_fsnos(:,:,filters),3);
        MODSCAG_fsnos_tmp = nanmean(MODSCAG_fsnos(:,:,filters),3);
        SPIRES_fsnos_tmp = nanmean(SPIRES_fsnos(:,:,filters),3);
        
        ELM_fsnos_tmp(isUS<1) = nan;
        MODSCAG_fsnos_tmp(isUS<1) = nan;
        SPIRES_fsnos_tmp(isUS<1) = nan;
        
        
        filters = ELM_grainsizes_tmp>0 & MODSCAG_grainsizes_tmp>0 & SPIRES_grainsizes_tmp>0 & ...
            ELM_fsnos_tmp>0;
        area_weight = Areas(filters)./sum(Areas(filters));
        
        ELM_grainsizes_yearly(year_i-2000,season_i)= nansum(ELM_grainsizes_tmp(filters).*area_weight);
        MODSCAG_grainsizes_yearly(year_i-2000,season_i) = nansum(MODSCAG_grainsizes_tmp(filters).*area_weight);
        SPIRES_grainsizes_yearly(year_i-2000,season_i) = nansum(SPIRES_grainsizes_tmp(filters).*area_weight);
        
        ELM_grainsizes_yearly_std(year_i-2000,season_i)= std(ELM_grainsizes_tmp(filters),area_weight,1,'omitnan');
        MODSCAG_grainsizes_yearly_std(year_i-2000,season_i) = std(MODSCAG_grainsizes_tmp(filters),area_weight,1,'omitnan');
        SPIRES_grainsizes_yearly_std(year_i-2000,season_i) = std(SPIRES_grainsizes_tmp(filters),area_weight,1,'omitnan');
        
    end
end

%% elevation distribution
elevations = imread('../WUS_mean_of_Elevation.tif');
elevations = squeeze(elevations(:,:,1));
indexs  = repmat(0:18,3,1);
elevations_bands = nan(size(elevations));
for elevation_band = 1:7
    filters = elevations >=((elevation_band-1)*500)  & elevations<(elevation_band*500);
    if elevation_band == 7
        filters = elevations >=((elevation_band-1)*500) ;
    end
    elevations_bands(filters) = elevation_band;
end

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
    ELM_grainsizes_tmp = nanmean(ELM_grainsizes(:,:,filters),3);
    MODSCAG_grainsizes_tmp = nanmean(MODSCAG_grainsizes(:,:,filters),3);
    SPIRES_grainsizes_tmp = nanmean(SPIRES_grainsizes(:,:,filters),3);
    
        ELM_fsnos_tmp = nanmean(ELM_fsnos(:,:,filters),3);
        MODSCAG_fsnos_tmp = nanmean(MODSCAG_fsnos(:,:,filters),3);
        SPIRES_fsnos_tmp = nanmean(SPIRES_fsnos(:,:,filters),3);
        
        ELM_fsnos_tmp(isUS<1) = nan;
        MODSCAG_fsnos_tmp(isUS<1) = nan;
        SPIRES_fsnos_tmp(isUS<1) = nan;

        
    filters = ELM_grainsizes_tmp<0 | MODSCAG_grainsizes_tmp<0 | SPIRES_grainsizes_tmp<0 | ...
        isnan(ELM_grainsizes_tmp) | isnan(MODSCAG_grainsizes_tmp)  | isnan(SPIRES_grainsizes_tmp);
    
    ELM_grainsizes_tmp(isUS<1 | filters) = nan;
    MODSCAG_grainsizes_tmp(isUS<1 | filters) = nan;
    SPIRES_grainsizes_tmp(isUS<1 | filters) = nan;
    
    
    
    
    if season_i==1
        snow_data_all_winter = [];
        
        for elevation_band = 1:7
            filters = elevations_bands>=0 & ELM_grainsizes_tmp>0 & MODSCAG_grainsizes_tmp>0 & SPIRES_grainsizes_tmp>0 & ...
            ELM_fsnos_tmp>0;
            
            ELM_grainsizes_tmp2 = ELM_grainsizes_tmp(filters);
            MODSCAG_grainsizes_tmp2 = MODSCAG_grainsizes_tmp(filters);
            SPIRES_grainsizes_tmp2 = SPIRES_grainsizes_tmp(filters);
            ELM_grainsizes_tmp2(elevations_bands(filters)~=elevation_band) = nan;
            MODSCAG_grainsizes_tmp2(elevations_bands(filters)~=elevation_band) = nan;
            SPIRES_grainsizes_tmp2(elevations_bands(filters)~=elevation_band) = nan;
            
            tmp = [ELM_grainsizes_tmp2; MODSCAG_grainsizes_tmp2; SPIRES_grainsizes_tmp2];
            group_all_winter = [ones(size(tmp,1)/3, 1); 2*ones(size(tmp,1)/3, 1); 3*ones(size(tmp,1)/3, 1)];
            
            snow_data_all_winter =[snow_data_all_winter tmp];
        end
    else
        
        snow_data_all_spring = [];
        for elevation_band = 1:7
            filters = elevations_bands>=0 & ELM_grainsizes_tmp>0 & MODSCAG_grainsizes_tmp>0 & SPIRES_grainsizes_tmp>0 & ...
            ELM_fsnos_tmp>0;
        
            ELM_grainsizes_tmp2 = ELM_grainsizes_tmp(filters);
            MODSCAG_grainsizes_tmp2 = MODSCAG_grainsizes_tmp(filters);
            SPIRES_grainsizes_tmp2 = SPIRES_grainsizes_tmp(filters);
            ELM_grainsizes_tmp2(elevations_bands(filters)~=elevation_band) = nan;
            MODSCAG_grainsizes_tmp2(elevations_bands(filters)~=elevation_band) = nan;
            SPIRES_grainsizes_tmp2(elevations_bands(filters)~=elevation_band) = nan;
            
            tmp = [ELM_grainsizes_tmp2; MODSCAG_grainsizes_tmp2; SPIRES_grainsizes_tmp2];
            
            
            %tmp(elevations_bands(filters)~=elevation_band) = nan;
            group_all_spring = [ones(size(tmp,1)/3, 1); 2*ones(size(tmp,1)/3, 1); 3*ones(size(tmp,1)/3, 1)];
            
            snow_data_all_spring =[snow_data_all_spring tmp];
        end
    end
end


%% forest distribution
forest_cover = imread('../forest_cover.tif');

indexs  = repmat(0:18,3,1);
forest_cover(isUS<1) = nan;

fc_bands = nan(size(forest_cover));
for forest_band = 1:5
    
    if forest_band == 1
        filters = forest_cover >=0  & forest_cover<10;
    elseif forest_band<5
        filters = forest_cover >=((forest_band-2)*20+10)  & forest_cover<((forest_band-1)*20+10);
    else
        filters = forest_cover >=((forest_band-2)*20+10) ;
    end
    fc_bands(filters) = forest_band;
end


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
    ELM_grainsizes_tmp = nanmean(ELM_grainsizes(:,:,filters),3);
    MODSCAG_grainsizes_tmp = nanmean(MODSCAG_grainsizes(:,:,filters),3);
    SPIRES_grainsizes_tmp = nanmean(SPIRES_grainsizes(:,:,filters),3);
    
    ELM_fsnos_tmp = nanmean(ELM_fsnos(:,:,filters),3);
    MODSCAG_fsnos_tmp = nanmean(MODSCAG_fsnos(:,:,filters),3);
    SPIRES_fsnos_tmp = nanmean(SPIRES_fsnos(:,:,filters),3);
    
    ELM_fsnos_tmp(isUS<1) = nan;
    MODSCAG_fsnos_tmp(isUS<1) = nan;
    SPIRES_fsnos_tmp(isUS<1) = nan;

        
    filters = ELM_grainsizes_tmp<0 | MODSCAG_grainsizes_tmp<0 | SPIRES_grainsizes_tmp<0 | ...
        isnan(ELM_grainsizes_tmp) | isnan(MODSCAG_grainsizes_tmp)  | isnan(SPIRES_grainsizes_tmp);
    
    ELM_grainsizes_tmp(isUS<1 | filters) = nan;
    MODSCAG_grainsizes_tmp(isUS<1 | filters) = nan;
    SPIRES_grainsizes_tmp(isUS<1 | filters) = nan;
    
    
    
    
    if season_i==1
        snow_data_all2_winter = [];
        
        for forest_band = 1:5
            filters = fc_bands>=0 & ELM_grainsizes_tmp>0 & MODSCAG_grainsizes_tmp>0 & SPIRES_grainsizes_tmp>0 & ...
            ELM_fsnos_tmp>0;
            
            ELM_grainsizes_tmp2 = ELM_grainsizes_tmp(filters);
            MODSCAG_grainsizes_tmp2 = MODSCAG_grainsizes_tmp(filters);
            SPIRES_grainsizes_tmp2 = SPIRES_grainsizes_tmp(filters);
            ELM_grainsizes_tmp2(fc_bands(filters)~=forest_band) = nan;
            MODSCAG_grainsizes_tmp2(fc_bands(filters)~=forest_band) = nan;
            SPIRES_grainsizes_tmp2(fc_bands(filters)~=forest_band) = nan;
            
            tmp = [ELM_grainsizes_tmp2; MODSCAG_grainsizes_tmp2; SPIRES_grainsizes_tmp2];
            group_all2_winter = [ones(size(tmp,1)/3, 1); 2*ones(size(tmp,1)/3, 1); 3*ones(size(tmp,1)/3, 1)];
            
            snow_data_all2_winter =[snow_data_all2_winter tmp];
        end
    else
        
        snow_data_all2_spring = [];
        for forest_band = 1:7
            filters = fc_bands>=0 & ELM_grainsizes_tmp>0 & MODSCAG_grainsizes_tmp>0 & SPIRES_grainsizes_tmp>0 & ...
            ELM_fsnos_tmp>0;
        
            ELM_grainsizes_tmp2 = ELM_grainsizes_tmp(filters);
            MODSCAG_grainsizes_tmp2 = MODSCAG_grainsizes_tmp(filters);
            SPIRES_grainsizes_tmp2 = SPIRES_grainsizes_tmp(filters);
            ELM_grainsizes_tmp2(fc_bands(filters)~=forest_band) = nan;
            MODSCAG_grainsizes_tmp2(fc_bands(filters)~=forest_band) = nan;
            SPIRES_grainsizes_tmp2(fc_bands(filters)~=forest_band) = nan;
            
            tmp = [ELM_grainsizes_tmp2; MODSCAG_grainsizes_tmp2; SPIRES_grainsizes_tmp2];
            
            
            tmp(fc_bands(filters)~=forest_band) = nan;
            group_all2_spring = [ones(size(tmp,1)/3, 1); 2*ones(size(tmp,1)/3, 1); 3*ones(size(tmp,1)/3, 1)];
            
            snow_data_all2_spring =[snow_data_all2_spring tmp];
        end
    end
end


figure
set(gcf,'unit','normalized','position',[0.1,0.1,0.65,0.8]);
set(gca, 'Position', [0 0 1 1])

%% plot 1
subplot('position', [0.08 0.72 0.43 0.25]);
hold on

plot(2001:2019,ELM_grainsizes_yearly(:,1),'color',colors(1,:),'linewidth',1.5)
%plot(2001:2019,MODSCAG_grainsizes_yearly(:,1),'color',colors(2,:),'linewidth',1.5)
plot(2001:2019,SPIRES_grainsizes_yearly(:,1),'color',colors(3,:),'linewidth',1.5)


title('Winter')
box on
set(gca,'linewidth',1,'fontsize',10)
xlim([2001 2019])
ylim([150 500])
text(2001.5,500-450*0.06,'(a)','fontsize',12,'fontweight','bold')
ylabel('{\it S}_{sno} (\mum)','fontsize',12)
IAV1 = std(ELM_grainsizes_yearly(:,1)-nanmean(ELM_grainsizes_yearly(:,2)),1);
%IAV2 = std(MODSCAG_grainsizes_yearly(:,1)-nanmean(MODSCAG_grainsizes_yearly(:,2)),1);
IAV3 = std(SPIRES_grainsizes_yearly(:,1)-nanmean(SPIRES_grainsizes_yearly(:,2)),1);
text(2001.7,500-350*0.18,['IAV\_ELM: ' num2str(IAV1,'%4.1f')],'fontsize',11,'fontweight','bold','color',colors(1,:))
%text(2001.7,500-350*0.28,['IAV\_STC-MODSCAG: ' num2str(IAV2,'%4.1f')],'fontsize',11,'fontweight','bold','color',colors(2,:))
text(2001.7,500-350*0.28,['IAV\_SPIReS: ' num2str(IAV3,'%4.1f')],'fontsize',11,'fontweight','bold','color',colors(3,:))

xlabel('Year')

sen_ps1 = nan(3,2);
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' ELM_grainsizes_yearly(:,1)], 0.05, 0);
sen_ps1(1,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' MODSCAG_grainsizes_yearly(:,1)], 0.05, 0);
sen_ps1(2,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' SPIRES_grainsizes_yearly(:,1)], 0.05, 0);
sen_ps1(3,:) = [sen sig];

%% plot 2

subplot('position', [0.56 0.72 0.43 0.25]);
hold on


plot(2001:2019,ELM_grainsizes_yearly(:,2),'color',colors(1,:),'linewidth',1.5)
plot(2001:2019,MODSCAG_grainsizes_yearly(:,2),'color',colors(2,:),'linewidth',1.5)
plot(2001:2019,SPIRES_grainsizes_yearly(:,2),'color',colors(3,:),'linewidth',1.5)
title('Spring')
xlim([2001 2019])
ylim([200 800])
box on
set(gca,'linewidth',1,'fontsize',10)
text(2001.5,800-600*0.06,'(b)','fontsize',12,'fontweight','bold')

IAV1 = std(ELM_grainsizes_yearly(:,2)-nanmean(ELM_grainsizes_yearly(:,2)),1);
IAV2 = std(MODSCAG_grainsizes_yearly(:,2)-nanmean(MODSCAG_grainsizes_yearly(:,2)),1);
IAV3 = std(SPIRES_grainsizes_yearly(:,2)-nanmean(SPIRES_grainsizes_yearly(:,2)),1);
text(2001.7,800-600*0.18,['IAV\_ELM: ' num2str(IAV1,'%4.1f')],'fontsize',11,'fontweight','bold','color',colors(1,:))
text(2001.7,800-600*0.28,['IAV\_STC-MODSCAG: ' num2str(IAV2,'%4.1f')],'fontsize',11,'fontweight','bold','color',colors(2,:))
text(2001.7,800-600*0.38,['IAV\_SPIReS: ' num2str(IAV3,'%4.1f')],'fontsize',11,'fontweight','bold','color',colors(3,:))


xlabel('Year','fontsize',12)
sen_ps2 = nan(3,2);
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' ELM_grainsizes_yearly(:,1)], 0.05, 0);
sen_ps2(1,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' MODSCAG_grainsizes_yearly(:,1)], 0.05, 0);
sen_ps2(2,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' SPIRES_grainsizes_yearly(:,1)], 0.05, 0);
sen_ps2(3,:) = [sen sig];


%% plot 3
snow_data_all_winter_2 = snow_data_all_winter(group_all_winter~=2,:);
group_all_winter_2 = group_all_winter(group_all_winter~=2);

subplot('position', [0.08 0.4 0.43 0.25]);
hold on

group_names = {'Humans', 'Dogs' , 'God'};
condition_names = {'0.5-1.0', '1.0-1.5', '1.5-2.0', '2.0-2.5' , '2.5-3.0', '\geq3.0'};


h = daboxplot(snow_data_all_winter_2(:,2:end),'groups',group_all_winter_2,'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors2,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',13,'mean',1, 'outliers',0);


box on
set(gca,'linewidth',1,'fontsize',10)
xlim([0.5 6.5])
ylim([0 1050])
text(0.6,1050-1050*0.06,'(c)','fontsize',12,'fontweight','bold')
xlabel('Elevation (km)','fontsize',12)
ylabel('{\it S}_{sno} (\mum)','fontsize',12)

%% plot 4
subplot('position', [0.56 0.4 0.43 0.25]);
hold on

h = daboxplot(snow_data_all_spring(:,2:end),'groups',group_all_spring,'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',13,'mean',1, 'outliers',0,...
       'legend',{'ELM','STC-MODSCAG','SPIReS'});

legend('boxoff')


box on
set(gca,'linewidth',1,'fontsize',10)
xlim([0.5 6.5])
ylim([0 1200])
text(0.6,1200-1200*0.06,'(d)','fontsize',12,'fontweight','bold')


set(gca,'linewidth',1,'fontsize',10)
xlabel('Elevation (km)','fontsize',12)

%% forest
snow_data_all2_winter_2 = snow_data_all2_winter(group_all2_winter~=2,:);
group_all2_winter_2 = group_all2_winter(group_all2_winter~=2);

subplot('position', [0.08 0.08 0.43 0.25]);
hold on

group_names = {'Humans', 'Dogs' , 'God'};
condition_names = {'0-10', '10-30', '30-50', '50-70' , '\geq70'};

h = daboxplot(snow_data_all2_winter_2(:,1:end),'groups',group_all2_winter_2,'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors2,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',13,'mean',1, 'outliers',0);



box on
set(gca,'linewidth',1,'fontsize',10)
xlim([0.5 5.5])
ylim([0 1300])
text(0.6,1300-1300*0.06,'(e)','fontsize',12,'fontweight','bold')
xlabel('Forest cover (%)','fontsize',12)
ylabel('{\it S}_{sno} (\mum)','fontsize',12)



subplot('position', [0.56 0.08 0.43 0.25]);
hold on

group_names = {'Humans', 'Dogs' , 'God'};
condition_names = {'0-10', '10-30', '30-50', '50-70' , '\geq70'};

h = daboxplot(snow_data_all2_spring(:,1:end),'groups',group_all2_spring,'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',13,'mean',1, 'outliers',0);


box on
set(gca,'linewidth',1,'fontsize',10)
xlim([0.5 5.5])
ylim([0 1200])
text(0.6,1200 - 1200*0.06,'(f)','fontsize',12,'fontweight','bold');
xlabel('Forest cover (%)','fontsize',12)


print(gcf, '-dtiff', '-r300', '../../figure_all_tif/grainsize_trend_elevation_time.tif')

close all

