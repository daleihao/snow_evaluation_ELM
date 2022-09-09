
clc;
clear all;
close all;
%% lat lon
load('../../../../data_processing/monthly_data_3000.mat');
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

ELM_swes_yearly = nan(19,2);
UA_swes_yearly = nan(19,2);
SNODAS_swes_yearly = nan(19,2);
ELM_swes_yearly_std = nan(19,2);
UA_swes_yearly_std = nan(19,2);
SNODAS_swes_yearly_std = nan(19,2);


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
        ELM_swes_tmp = nanmean(ELM_swes(:,:,filters),3);
        UA_swes_tmp = nanmean(UA_swes(:,:,filters),3);
        SNODAS_swes_tmp = nanmean(SNODAS_swes(:,:,filters),3);
        
        ELM_swes_tmp(isUS<1) = nan;
        UA_swes_tmp(isUS<1) = nan;
        SNODAS_swes_tmp(isUS<1) = nan;
        
    ELM_fsnos_tmp = nanmean(ELM_fsnos(:,:,filters),3);
    MODSCAG_fsnos_tmp = nanmean(MODSCAG_fsnos(:,:,filters),3);
    SPIRES_fsnos_tmp = nanmean(SPIRES_fsnos(:,:,filters),3);
    
    ELM_fsnos_tmp(isUS<1) = nan;
    MODSCAG_fsnos_tmp(isUS<1) = nan;
    SPIRES_fsnos_tmp(isUS<1) = nan;
        
        
        filters = ELM_swes_tmp>0 & UA_swes_tmp>0 & SNODAS_swes_tmp>0 & ...
               ELM_fsnos_tmp>0 | MODSCAG_fsnos_tmp>0 | SPIRES_fsnos_tmp>0;
        area_weight = Areas(filters)./sum(Areas(filters));
        
        ELM_swes_yearly(year_i-2000,season_i)= nansum(ELM_swes_tmp(filters).*area_weight);
        UA_swes_yearly(year_i-2000,season_i) = nansum(UA_swes_tmp(filters).*area_weight);
        SNODAS_swes_yearly(year_i-2000,season_i) = nansum(SNODAS_swes_tmp(filters).*area_weight);
        
        ELM_swes_yearly_std(year_i-2000,season_i)= std(ELM_swes_tmp(filters),area_weight,1,'omitnan');
        UA_swes_yearly_std(year_i-2000,season_i) = std(UA_swes_tmp(filters),area_weight,1,'omitnan');
        SNODAS_swes_yearly_std(year_i-2000,season_i) = std(SNODAS_swes_tmp(filters),area_weight,1,'omitnan');
        
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
    ELM_swes_tmp = nanmean(ELM_swes(:,:,filters),3);
    UA_swes_tmp = nanmean(UA_swes(:,:,filters),3);
    SNODAS_swes_tmp = nanmean(SNODAS_swes(:,:,filters),3);
    
         ELM_fsnos_tmp = nanmean(ELM_fsnos(:,:,filters),3);
    MODSCAG_fsnos_tmp = nanmean(MODSCAG_fsnos(:,:,filters),3);
    SPIRES_fsnos_tmp = nanmean(SPIRES_fsnos(:,:,filters),3);
    
    ELM_fsnos_tmp(isUS<1) = nan;
    MODSCAG_fsnos_tmp(isUS<1) = nan;
    SPIRES_fsnos_tmp(isUS<1) = nan;

        
    filters = ELM_swes_tmp<0 | UA_swes_tmp<0 | SNODAS_swes_tmp<0 | ...
        isnan(ELM_swes_tmp) | isnan(UA_swes_tmp)  | isnan(SNODAS_swes_tmp);
    
    ELM_swes_tmp(isUS<1 | filters) = nan;
    UA_swes_tmp(isUS<1 | filters) = nan;
    SNODAS_swes_tmp(isUS<1 | filters) = nan;
    
    
    
    
    if season_i==1
        snow_data_all_winter = [];
        
        for elevation_band = 1:7
            filters = elevations_bands>=0 & ELM_swes_tmp>0 & UA_swes_tmp>0 & SNODAS_swes_tmp>0 & ...
            ELM_fsnos_tmp>0;
            
            ELM_swes_tmp2 = ELM_swes_tmp(filters);
            UA_swes_tmp2 = UA_swes_tmp(filters);
            SNODAS_swes_tmp2 = SNODAS_swes_tmp(filters);
            ELM_swes_tmp2(elevations_bands(filters)~=elevation_band) = nan;
            UA_swes_tmp2(elevations_bands(filters)~=elevation_band) = nan;
            SNODAS_swes_tmp2(elevations_bands(filters)~=elevation_band) = nan;
            
            tmp = [ELM_swes_tmp2; UA_swes_tmp2; SNODAS_swes_tmp2];
            group_all_winter = [ones(size(tmp,1)/3, 1); 2*ones(size(tmp,1)/3, 1); 3*ones(size(tmp,1)/3, 1)];
            
            snow_data_all_winter =[snow_data_all_winter tmp];
        end
    else
        
        snow_data_all_spring = [];
        for elevation_band = 1:7
            filters = elevations_bands>=0 & ELM_swes_tmp>0 & UA_swes_tmp>0 & SNODAS_swes_tmp>0 & ...
            ELM_fsnos_tmp>0;
        
            ELM_swes_tmp2 = ELM_swes_tmp(filters);
            UA_swes_tmp2 = UA_swes_tmp(filters);
            SNODAS_swes_tmp2 = SNODAS_swes_tmp(filters);
            ELM_swes_tmp2(elevations_bands(filters)~=elevation_band) = nan;
            UA_swes_tmp2(elevations_bands(filters)~=elevation_band) = nan;
            SNODAS_swes_tmp2(elevations_bands(filters)~=elevation_band) = nan;
            
            tmp = [ELM_swes_tmp2; UA_swes_tmp2; SNODAS_swes_tmp2];
            
            
            %tmp(elevations_bands(filters)~=elevation_band) = nan;
            group_all_spring = [ones(size(tmp,1)/3, 1); 2*ones(size(tmp,1)/3, 1); 3*ones(size(tmp,1)/3, 1)];
            
            snow_data_all_spring =[snow_data_all_spring tmp];
        end
    end
end

%% forest cover



%% plot figure
set(gcf,'unit','normalized','position',[0.1,0.1,0.65,0.55]);
set(gca, 'Position', [0 0 1 1])

%% plot 1
subplot('position', [0.07 0.615 0.43 0.33]);
hold on

plot(2001:2019,ELM_swes_yearly(:,1),'color',colors(1,:),'linewidth',1.5)
plot(2001:2019,UA_swes_yearly(:,1),'color',colors(2,:),'linewidth',1.5)
plot(2004:2019,SNODAS_swes_yearly(4:end,1),'color',colors(3,:),'linewidth',1.5)


title('Winter')
box on
set(gca,'linewidth',1,'fontsize',10)
xlim([2001 2019])
ylim([0 100])
text(2001.5,100-100*0.06,'(a)','fontsize',12,'fontweight','bold')
ylabel('SWE (mm)','fontsize',12)
IAV1 = std(ELM_swes_yearly(:,1)-nanmean(ELM_swes_yearly(:,2)),1);
IAV2 = std(UA_swes_yearly(:,1)-nanmean(UA_swes_yearly(:,2)),1);
IAV3 = std(SNODAS_swes_yearly(4:end,1)-nanmean(SNODAS_swes_yearly(4:end,2)),1);
text(2001.7,100-100*0.18,['IAV\_ELM: ' num2str(IAV1,'%4.1f')],'fontsize',11,'fontweight','bold','color',colors(1,:))
text(2001.7,100-100*0.28,['IAV\_UA: ' num2str(IAV2,'%4.1f')],'fontsize',11,'fontweight','bold','color',colors(2,:))
text(2001.7,100-100*0.38,['IAV\_SNODAS: ' num2str(IAV3,'%4.1f')],'fontsize',11,'fontweight','bold','color',colors(3,:))

xlabel('Year')

sen_ps1 = nan(3,2);
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' ELM_swes_yearly(:,1)], 0.05, 0);
sen_ps1(1,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' UA_swes_yearly(:,1)], 0.05, 0);
sen_ps1(2,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' SNODAS_swes_yearly(:,1)], 0.05, 0);
sen_ps1(3,:) = [sen sig];

%% plot 2

subplot('position', [0.55 0.615 0.43 0.33]);
hold on


plot(2001:2019,ELM_swes_yearly(:,2),'color',colors(1,:),'linewidth',1.5)
plot(2001:2019,UA_swes_yearly(:,2),'color',colors(2,:),'linewidth',1.5)
plot(2004:2019,SNODAS_swes_yearly(4:end,2),'color',colors(3,:),'linewidth',1.5)
title('Spring')
xlim([2001 2019])
ylim([0 150])
box on
set(gca,'linewidth',1,'fontsize',10)
text(2001.5,150-150*0.06,'(b)','fontsize',12,'fontweight','bold')

IAV1 = std(ELM_swes_yearly(:,2)-nanmean(ELM_swes_yearly(:,2)),1);
IAV2 = std(UA_swes_yearly(:,2)-nanmean(UA_swes_yearly(:,2)),1);
IAV3 = std(SNODAS_swes_yearly(4:end,2)-nanmean(SNODAS_swes_yearly(4:end,2)),1);
text(2001.7,150-150*0.18,['IAV\_ELM: ' num2str(IAV1,'%4.1f')],'fontsize',11,'fontweight','bold','color',colors(1,:))
text(2001.7,150-150*0.28,['IAV\_UA: ' num2str(IAV2,'%4.1f')],'fontsize',11,'fontweight','bold','color',colors(2,:))
text(2001.7,150-150*0.38,['IAV\_SNODAS: ' num2str(IAV3,'%4.1f')],'fontsize',11,'fontweight','bold','color',colors(3,:))


xlabel('Year','fontsize',12)
sen_ps2 = nan(3,2);
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' ELM_swes_yearly(:,1)], 0.05, 0);
sen_ps2(1,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' UA_swes_yearly(:,1)], 0.05, 0);
sen_ps2(2,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' SNODAS_swes_yearly(:,1)], 0.05, 0);
sen_ps2(3,:) = [sen sig];


%% plot 3

subplot('position', [0.07 0.12 0.43 0.38]);
hold on

group_names = {'Humans', 'Dogs' , 'God'};
condition_names = {'0.5-1.0', '1.0-1.5', '1.5-2.0', '2.0-2.5' , '2.5-3.0', '>3.0'};


h = daboxplot(snow_data_all_winter(:,2:end),'groups',group_all_winter,'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',13,'mean',1, 'outliers',0,...
       'add_dis',-0.02,...
       'legend',{'ELM','UA','SNODAS'});

legend('boxoff')

box on
set(gca,'linewidth',1,'fontsize',10)
xlim([0.5 6.5])
ylim([0 600])
text(0.6,600-600*0.06,'(c)','fontsize',12,'fontweight','bold')
xlabel('Elevation (km)','fontsize',12)
ylabel('SWE (mm)','fontsize',12)

%% plot 4
subplot('position', [0.55 0.12 0.43 0.38]);
hold on

h = daboxplot(snow_data_all_spring(:,2:end),'groups',group_all_spring,'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',13,'mean',1, 'outliers',0);



box on
set(gca,'linewidth',1,'fontsize',10)
xlim([0.5 6.5])
ylim([0 600])
text(0.6,600-600*0.06,'(d)','fontsize',12,'fontweight','bold')


set(gca,'linewidth',1,'fontsize',10)
xlabel('Elevation (km)','fontsize',12)


print(gcf, '-dtiff', '-r300', '../../figure_all_tif/swe_trend_elevation_time.tif')

close all

