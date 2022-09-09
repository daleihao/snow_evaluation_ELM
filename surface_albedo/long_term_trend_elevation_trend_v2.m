
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

ELM_surfacealbedos_yearly = nan(19,2);
MCD43_surfacealbedos_yearly = nan(19,2);
ELM_surfacealbedos_yearly_std = nan(19,2);
MCD43_surfacealbedos_yearly_std = nan(19,2);


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
        
        ELM_surfacealbedos_tmp = nanmean(ELM_Albedos(:,:,filters),3);
        MCD43_surfacealbedo = MCD43_BSAs.* (1 - ELM_Skyls) + MCD43_WSAs .* ELM_Skyls;
        MCD43_surfacealbedos_tmp = nanmean(MCD43_surfacealbedo(:,:,filters),3)/1000;
         
        ELM_fsnos_tmp = nanmean(ELM_fsnos(:,:,filters),3);
        MODSCAG_fsnos_tmp = nanmean(MODSCAG_fsnos(:,:,filters),3);
        SPIRES_fsnos_tmp = nanmean(SPIRES_fsnos(:,:,filters),3);
        
        ELM_fsnos_tmp(isUS<1) = nan;
        MODSCAG_fsnos_tmp(isUS<1) = nan;
        SPIRES_fsnos_tmp(isUS<1) = nan;
        
        ELM_surfacealbedos_tmp(isUS<1) = nan;
        MCD43_surfacealbedos_tmp(isUS<1) = nan;

        
        
        filters = ELM_surfacealbedos_tmp>0 & MCD43_surfacealbedos_tmp>0 & ...
            ELM_fsnos_tmp>0;
        area_weight = Areas(filters)./sum(Areas(filters));
        
        ELM_surfacealbedos_yearly(year_i-2000,season_i)= nansum(ELM_surfacealbedos_tmp(filters).*area_weight);
        MCD43_surfacealbedos_yearly(year_i-2000,season_i) = nansum(MCD43_surfacealbedos_tmp(filters).*area_weight);
        
        ELM_surfacealbedos_yearly_std(year_i-2000,season_i)= std(ELM_surfacealbedos_tmp(filters),area_weight,1,'omitnan');
        MCD43_surfacealbedos_yearly_std(year_i-2000,season_i) = std(MCD43_surfacealbedos_tmp(filters),area_weight,1,'omitnan');
        
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
    ELM_surfacealbedos_tmp = nanmean(ELM_Albedos(:,:,filters),3);
    MCD43_surfacealbedo = MCD43_BSAs.* (1 - ELM_Skyls) + MCD43_WSAs .* ELM_Skyls;
    MCD43_surfacealbedos_tmp = nanmean(MCD43_surfacealbedo(:,:,filters),3)/1000;
    
      ELM_fsnos_tmp = nanmean(ELM_fsnos(:,:,filters),3);
        MODSCAG_fsnos_tmp = nanmean(MODSCAG_fsnos(:,:,filters),3);
        SPIRES_fsnos_tmp = nanmean(SPIRES_fsnos(:,:,filters),3);
        
        ELM_fsnos_tmp(isUS<1) = nan;
        MODSCAG_fsnos_tmp(isUS<1) = nan;
        SPIRES_fsnos_tmp(isUS<1) = nan;        
                
    filters = ELM_surfacealbedos_tmp<0 | MCD43_surfacealbedos_tmp<0  | ...
        isnan(ELM_surfacealbedos_tmp) | isnan(MCD43_surfacealbedos_tmp) |... 
        ELM_fsnos_tmp<=0 | MODSCAG_fsnos_tmp<=0 | SPIRES_fsnos_tmp<=0 | ...
        isnan(ELM_fsnos_tmp) | isnan(MODSCAG_fsnos_tmp)  | isnan(SPIRES_fsnos_tmp);
    
    ELM_surfacealbedos_tmp(isUS<1 | filters) = nan;
    MCD43_surfacealbedos_tmp(isUS<1 | filters) = nan;
    
    
    
    
    if season_i==1
        snow_data_all_winter = [];
        
        for elevation_band = 1:7
            filters = elevations_bands>=0 & ELM_surfacealbedos_tmp>0 & MCD43_surfacealbedos_tmp>0 & ...
               ELM_fsnos_tmp>0;
            
            ELM_surfacealbedos_tmp2 = ELM_surfacealbedos_tmp(filters);
            MCD43_surfacealbedos_tmp2 = MCD43_surfacealbedos_tmp(filters);

            ELM_surfacealbedos_tmp2(elevations_bands(filters)~=elevation_band) = nan;
            MCD43_surfacealbedos_tmp2(elevations_bands(filters)~=elevation_band) = nan;
            
            tmp = [ELM_surfacealbedos_tmp2; MCD43_surfacealbedos_tmp2];
            group_all_winter = [ones(size(tmp,1)/2, 1); 2*ones(size(tmp,1)/2, 1)];
            
            snow_data_all_winter =[snow_data_all_winter tmp];
        end
    else
        
        snow_data_all_spring = [];
        for elevation_band = 1:7
            filters = elevations_bands>=0 & ELM_surfacealbedos_tmp>0 & MCD43_surfacealbedos_tmp>0 & ...
               ELM_fsnos_tmp>0;
            ELM_surfacealbedos_tmp2 = ELM_surfacealbedos_tmp(filters);
            MCD43_surfacealbedos_tmp2 = MCD43_surfacealbedos_tmp(filters);

            ELM_surfacealbedos_tmp2(elevations_bands(filters)~=elevation_band) = nan;
            MCD43_surfacealbedos_tmp2(elevations_bands(filters)~=elevation_band) = nan;
            
            tmp = [ELM_surfacealbedos_tmp2; MCD43_surfacealbedos_tmp2];
            
            
            tmp(elevations_bands(filters)~=elevation_band) = nan;
            group_all_spring = [ones(size(tmp,1)/2, 1); 2*ones(size(tmp,1)/2, 1)];
            
            snow_data_all_spring =[snow_data_all_spring tmp];
        end
    end
end

%% snow cover plot

indexs  = repmat(0:18,3,1);

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
    ELM_surfacealbedos_tmp = nanmean(ELM_Albedos(:,:,filters),3);
    MCD43_surfacealbedo = MCD43_BSAs.* (1 - ELM_Skyls) + MCD43_WSAs .* ELM_Skyls;
    MCD43_surfacealbedos_tmp = nanmean(MCD43_surfacealbedo(:,:,filters),3)/1000;
    
      ELM_fsnos_tmp = nanmean(ELM_fsnos(:,:,filters),3);
        MODSCAG_fsnos_tmp = nanmean(MODSCAG_fsnos(:,:,filters),3);
        SPIRES_fsnos_tmp = nanmean(SPIRES_fsnos(:,:,filters),3);
        
        ELM_fsnos_tmp(isUS<1) = nan;
        MODSCAG_fsnos_tmp(isUS<1) = nan;
        SPIRES_fsnos_tmp(isUS<1) = nan;        
                
    filters = ELM_surfacealbedos_tmp<0 | MCD43_surfacealbedos_tmp<0  | ...
        isnan(ELM_surfacealbedos_tmp) | isnan(MCD43_surfacealbedos_tmp) |... 
        ELM_fsnos_tmp<=0 | MODSCAG_fsnos_tmp<=0 | SPIRES_fsnos_tmp<=0 | ...
        isnan(ELM_fsnos_tmp) | isnan(MODSCAG_fsnos_tmp)  | isnan(SPIRES_fsnos_tmp);
    
    ELM_surfacealbedos_tmp(isUS<1 | filters) = nan;
    MCD43_surfacealbedos_tmp(isUS<1 | filters) = nan;
    
    
    snowcover_bands = nan(size(ELM_fsnos_tmp));
for snowcover_band = 1:5
    filters = ELM_fsnos_tmp >=((snowcover_band-1)*0.2)  & ELM_fsnos_tmp<(snowcover_band*0.2)  & ...
               ELM_fsnos_tmp>0;
    if elevation_band == 5
        filters = ELM_fsnos_tmp >=((snowcover_band-1)*0.2) ;
    end
    snowcover_bands(filters) = snowcover_band;
end
    
    if season_i==1
        snow_data_all_winter_fsno = [];
        
        for snowcover_band = 1:5
            filters = snowcover_bands>=0 & ELM_surfacealbedos_tmp>=0 & MCD43_surfacealbedos_tmp>=0 & ...
               ELM_fsnos_tmp>0;
            
            ELM_surfacealbedos_tmp2 = ELM_surfacealbedos_tmp(filters);
            MCD43_surfacealbedos_tmp2 = MCD43_surfacealbedos_tmp(filters);

            ELM_surfacealbedos_tmp2(snowcover_bands(filters)~=snowcover_band) = nan;
            MCD43_surfacealbedos_tmp2(snowcover_bands(filters)~=snowcover_band) = nan;
            
            tmp = [ELM_surfacealbedos_tmp2; MCD43_surfacealbedos_tmp2];
            group_all_winter_fsno = [ones(size(tmp,1)/2, 1); 2*ones(size(tmp,1)/2, 1)];
            
            snow_data_all_winter_fsno =[snow_data_all_winter_fsno tmp];
        end
    else
        
        snow_data_all_spring_fsno = [];
        for snowcover_band = 1:5
            filters = snowcover_bands>=0 & ELM_surfacealbedos_tmp>0 & MCD43_surfacealbedos_tmp>0 & ...
               ELM_fsnos_tmp>0;
            ELM_surfacealbedos_tmp2 = ELM_surfacealbedos_tmp(filters);
            MCD43_surfacealbedos_tmp2 = MCD43_surfacealbedos_tmp(filters);

            ELM_surfacealbedos_tmp2(snowcover_bands(filters)~=snowcover_band) = nan;
            MCD43_surfacealbedos_tmp2(snowcover_bands(filters)~=snowcover_band) = nan;
            
            tmp = [ELM_surfacealbedos_tmp2; MCD43_surfacealbedos_tmp2];
            
            
            %tmp(snowcover_bands(filters)~=snowcover_band) = nan;
            group_all_spring_fsno = [ones(size(tmp,1)/2, 1); 2*ones(size(tmp,1)/2, 1)];
            
            snow_data_all_spring_fsno =[snow_data_all_spring_fsno tmp];
        end
    end
end
%% figure plot

figure
set(gcf,'unit','normalized','position',[0.1,0.1,0.65,0.78]);
set(gca, 'Position', [0 0 1 1])

%% plot 1
subplot('position', [0.07 0.72 0.43 0.24]);
hold on

plot(2001:2019,ELM_surfacealbedos_yearly(:,1),'color',colors(1,:),'linewidth',1.5)
plot(2001:2019,MCD43_surfacealbedos_yearly(:,1),'color',colors(2,:),'linewidth',1.5)


title('Winter')
box on
set(gca,'linewidth',1,'fontsize',10)
xlim([2001 2019])
ylim([0.2 0.5])
text(2001.2,0.5-0.3*0.06,'(a)','fontsize',12,'fontweight','bold')
ylabel('{\it \alpha}_{sur}','fontsize',12)
IAV1 = std(ELM_surfacealbedos_yearly(:,1)-nanmean(ELM_surfacealbedos_yearly(:,2)),1);
IAV2 = std(MCD43_surfacealbedos_yearly(:,1)-nanmean(MCD43_surfacealbedos_yearly(:,2)),1);
text(2001.7,0.5-0.3*0.18,['IAV\_ELM: ' num2str(IAV1,'%4.3f')],'fontsize',11,'fontweight','bold','color',colors(1,:))
text(2001.7,0.5-0.3*0.28,['IAV\_MCD43: ' num2str(IAV2,'%4.3f')],'fontsize',11,'fontweight','bold','color',colors(2,:))

xlabel('Year','fontsize',10)

sen_ps1 = nan(3,2);
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' ELM_surfacealbedos_yearly(:,1)], 0.05, 0);
sen_ps1(1,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' MCD43_surfacealbedos_yearly(:,1)], 0.05, 0);
sen_ps1(2,:) = [sen sig];

%% plot 2

subplot('position', [0.55 0.72 0.43 0.24]);
hold on


plot(2001:2019,ELM_surfacealbedos_yearly(:,2),'color',colors(1,:),'linewidth',1.5)
plot(2001:2019,MCD43_surfacealbedos_yearly(:,2),'color',colors(2,:),'linewidth',1.5)

title('Spring')
xlim([2001 2019])
ylim([0.1 0.30001])
box on
set(gca,'linewidth',1,'fontsize',10)
text(2001.2,0.3-0.2*0.06,'(b)','fontsize',12,'fontweight','bold')

IAV1 = std(ELM_surfacealbedos_yearly(:,2)-nanmean(ELM_surfacealbedos_yearly(:,2)),1);
IAV2 = std(MCD43_surfacealbedos_yearly(:,2)-nanmean(MCD43_surfacealbedos_yearly(:,2)),1);
text(2001.7,0.3-0.2*0.18,['IAV\_ELM: ' num2str(IAV1,'%4.3f')],'fontsize',11,'fontweight','bold','color',colors(1,:))
text(2001.7,0.3-0.2*0.28,['IAV\_MCD43: ' num2str(IAV2,'%4.3f')],'fontsize',11,'fontweight','bold','color',colors(2,:))


xlabel('Year','fontsize',10)
sen_ps2 = nan(3,2);
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' ELM_surfacealbedos_yearly(:,1)], 0.05, 0);
sen_ps2(1,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' MCD43_surfacealbedos_yearly(:,1)], 0.05, 0);
sen_ps2(2,:) = [sen sig];

%% plot 3

subplot('position', [0.07 0.40 0.43 0.24]);
hold on

group_names = {'Humans', 'Dogs' , 'God'};
condition_names = {'0.5-1.0', '1.0-1.5', '1.5-2.0', '2.0-2.5' , '2.5-3.0', '>3.0'};


h = daboxplot(snow_data_all_winter(:,2:end),'groups',group_all_winter,'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',13,'mean',1, 'outliers',0);



box on
set(gca,'linewidth',1,'fontsize',10)
xlim([0.5 6.5])
ylim([-0.1 0.9])
text(0.6,0.9-1*0.06,'(c)','fontsize',12,'fontweight','bold')
xlabel('Elevation (km)','fontsize',10)
ylabel('{\it \alpha}_{sur}','fontsize',12)

%% plot 4
subplot('position', [0.55 0.40 0.43 0.24]);
hold on

h = daboxplot(snow_data_all_spring(:,2:end),'groups',group_all_spring,'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',13,'mean',1, 'outliers',0);



box on
set(gca,'linewidth',1,'fontsize',10)
xlim([0.5 6.5])
ylim([-0.1 0.9])
text(0.6,0.9-1*0.06,'(d)','fontsize',12,'fontweight','bold')


set(gca,'linewidth',1,'fontsize',10)
xlabel('Elevation (km)','fontsize',10)



%% plot 5

subplot('position', [0.07 0.08 0.43 0.24]);
hold on

group_names = {'Humans', 'Dogs' , 'God'};
condition_names = {'0-0.2', '0.2-0.4', '0.4-0.6', '0.6-0.8', '0.8-1.0'};


h = daboxplot(snow_data_all_winter_fsno(:,:),'groups',group_all_winter_fsno,'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',13,'mean',1, 'outliers',0);



box on
set(gca,'linewidth',1,'fontsize',10)
xlim([0.5 5.5])
ylim([-0.1 0.9])
text(0.6,0.9-1*0.06,'(e)','fontsize',12,'fontweight','bold')
xlabel('{\it f}_{sno}','fontsize',10)
ylabel('{\it \alpha}_{sur}','fontsize',12)

%% plot 6
subplot('position', [0.55 0.08 0.43 0.24]);
hold on

h = daboxplot(snow_data_all_spring_fsno(:,:),'groups',group_all_spring_fsno,'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
   'add_dis',-0.02,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',13,'mean',1, 'outliers',0,'legend',{'ELM','MCD43'});



box on
set(gca,'linewidth',1,'fontsize',10)
xlim([0.5 5.5])
ylim([-0.1 0.9])
text(0.6,0.9-1*0.06,'(f)','fontsize',12,'fontweight','bold')


legend('boxoff')
set(gca,'linewidth',1,'fontsize',10)
xlabel('{\it f}_{sno}','fontsize',10)

print(gcf, '-dtiff', '-r300', '../../figure_all_tif/surfacealbedo_trend_elevation_time.tif')

close all

