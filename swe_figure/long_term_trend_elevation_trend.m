
clc;
clear all;
close all;
%% lat lon
load('../../../../data_processing/monthly_data_3000.mat');
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
        
        if(year_i>2003)
        filters = ELM_swes_tmp>=0 & UA_swes_tmp>=0 & SNODAS_swes_tmp>=0;
        else
        filters = ELM_swes_tmp>=0 & UA_swes_tmp>=0;

        end
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

ELM_swes_elevation = nan(7,2);
UA_swes_elevation = nan(7,2);
SNODAS_swes_elevation = nan(7,2);

ELM_swes_elevation_std = nan(7,2);
UA_swes_elevation_std = nan(7,2);
SNODAS_swes_elevation_std = nan(7,2);
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
    
    filters = ELM_swes_tmp<0 | UA_swes_tmp<0 | SNODAS_swes_tmp<0 | ...
        isnan(ELM_swes_tmp) | isnan(UA_swes_tmp)  | isnan(SNODAS_swes_tmp);
    
    ELM_swes_tmp(isUS<1 | filters) = nan;
    UA_swes_tmp(isUS<1 | filters) = nan;
    SNODAS_swes_tmp(isUS<1 | filters) = nan;

    for elevation_band = 1:7
        filters = elevations >=((elevation_band-1)*500)  & elevations<(elevation_band*500);
        if elevation_band == 7
            filters = elevations >=((elevation_band-1)*500) ;
        end
            area_weight = Areas(filters)./sum(Areas(filters));

        ELM_swes_elevation(elevation_band,season_i) = nansum(ELM_swes_tmp(filters).*area_weight);
        UA_swes_elevation(elevation_band,season_i) = nansum(UA_swes_tmp(filters).*area_weight);
        SNODAS_swes_elevation(elevation_band,season_i) = nansum(SNODAS_swes_tmp(filters).*area_weight);        
                
        ELM_swes_elevation_std(elevation_band,season_i) = std(ELM_swes_tmp(filters),area_weight,1,'omitnan');
        UA_swes_elevation_std(elevation_band,season_i) = std(UA_swes_tmp(filters),area_weight,1,'omitnan');
        SNODAS_swes_elevation_std(elevation_band,season_i) = std(SNODAS_swes_tmp(filters),area_weight,1,'omitnan');
    end
end


figure
set(gcf,'unit','normalized','position',[0.1,0.1,0.45,0.55]);
set(gca, 'Position', [0 0 1 1])

%% plot 1
subplot('position', [0.11 0.615 0.39 0.34]);
hold on
inBetween = [ELM_swes_yearly(:,1)+ELM_swes_yearly_std(:,1); flipud(ELM_swes_yearly(:,1)-ELM_swes_yearly_std(:,1))];
patch([[2001:2019]'; flipud([2001:2019]')], inBetween,'r','FaceColor', '#4DBEEE','LineStyle','none')
alpha(0.3)


inBetween = [UA_swes_yearly(:,1)+UA_swes_yearly_std(:,1); flipud(UA_swes_yearly(:,1)-UA_swes_yearly_std(:,1))];
patch([[2001:2019]'; flipud([2001:2019]')], inBetween,'r','FaceColor', '#EDB120','LineStyle','none')
alpha(0.3)

inBetween = [SNODAS_swes_yearly(4:end,1)+SNODAS_swes_yearly_std(4:end,1); flipud(SNODAS_swes_yearly(4:end,1)-SNODAS_swes_yearly_std(4:end,1))];
patch([[2004:2019]'; flipud([2004:2019]')], inBetween, 'm','LineStyle','none')
alpha(0.3)

plot(2001:2019,ELM_swes_yearly(:,1),'color','#4DBEEE','linewidth',1.5)
plot(2001:2019,UA_swes_yearly(:,1),'color','#EDB120','linewidth',1.5)
plot(2004:2019,SNODAS_swes_yearly(4:end,1),'m','linewidth',1.5)


title('Winter')
box on
set(gca,'linewidth',1,'fontsize',10)
xlim([2001 2019])
ylim([-100 200])
text(2001.5,200-300*0.06,'(a)','fontsize',12,'fontweight','bold')
ylabel('SWE (m)','fontsize',12)

xlabel('Year','fontsize',12)

sen_ps1 = nan(3,2);
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' ELM_swes_yearly(:,1)], 0.05, 0);
sen_ps1(1,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' UA_swes_yearly(:,1)], 0.05, 0);
sen_ps1(2,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2004:2019]' SNODAS_swes_yearly(4:end,1)], 0.05, 0);
sen_ps1(3,:) = [sen sig];

%% plot 2

subplot('position', [0.59 0.615 0.39 0.34]);
hold on
inBetween = [ELM_swes_yearly(:,2)+ELM_swes_yearly_std(:,2); flipud(ELM_swes_yearly(:,2)-ELM_swes_yearly_std(:,2))];
patch([[2001:2019]'; flipud([2001:2019]')], inBetween,'r','FaceColor', '#4DBEEE','LineStyle','none')
alpha(0.3)


inBetween = [UA_swes_yearly(:,2)+UA_swes_yearly_std(:,2); flipud(UA_swes_yearly(:,2)-UA_swes_yearly_std(:,2))];
patch([[2001:2019]'; flipud([2001:2019]')], inBetween,'r','FaceColor', '#EDB120','LineStyle','none')
alpha(0.3)


inBetween = [SNODAS_swes_yearly(4:end,2)+SNODAS_swes_yearly_std(4:end,2); flipud(SNODAS_swes_yearly(4:end,2)-SNODAS_swes_yearly_std(4:end,2))];
patch([[2004:2019]'; flipud([2004:2019]')], inBetween, 'm','LineStyle','none')
alpha(0.3)

plot(2001:2019,ELM_swes_yearly(:,2),'color','#4DBEEE','linewidth',1.5)
plot(2001:2019,UA_swes_yearly(:,2),'color','#EDB120','linewidth',1.5)
plot(2004:2019,SNODAS_swes_yearly(4:end,2),'m','linewidth',1.5)
title('Spring')
xlim([2001 2019])
ylim([-200 300])
box on
set(gca,'linewidth',1,'fontsize',10)
text(2001.5,300-500*0.06,'(b)','fontsize',12,'fontweight','bold')
xlabel('Year','fontsize',12)
sen_ps2 = nan(3,2);
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' ELM_swes_yearly(:,1)], 0.05, 0);
sen_ps2(1,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' UA_swes_yearly(:,1)], 0.05, 0);
sen_ps2(2,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2004:2019]' SNODAS_swes_yearly(4:end,1)], 0.05, 0);
sen_ps2(3,:) = [sen sig];


%% plot 3

subplot('position', [0.11 0.12 0.39 0.38]);
hold on

inBetween = [ELM_swes_elevation(:,1)+ELM_swes_elevation_std(:,1); flipud(ELM_swes_elevation(:,1)-ELM_swes_elevation_std(:,1))];
patch(inBetween,[[250:500:3500]'; flipud([250:500:3500]')],'r','FaceColor', '#4DBEEE','LineStyle','none')
alpha(0.3)

inBetween = [UA_swes_elevation(:,1)+UA_swes_elevation_std(:,1); flipud(UA_swes_elevation(:,1)-UA_swes_elevation_std(:,1))];
patch(inBetween,[[250:500:3500]'; flipud([250:500:3500]')],'r','FaceColor', '#EDB120','LineStyle','none')
alpha(0.3)

inBetween = [SNODAS_swes_elevation(:,1)+SNODAS_swes_elevation_std(:,1); flipud(SNODAS_swes_elevation(:,1)-SNODAS_swes_elevation_std(:,1))];
patch(inBetween,[[250:500:3500]'; flipud([250:500:3500]')],'r','FaceColor', 'm','LineStyle','none')
alpha(0.3)

plot(ELM_swes_elevation(:,1),[250:500:3500],'color','#4DBEEE','linewidth',1.5)
plot(UA_swes_elevation(:,1),[250:500:3500],'color','#EDB120','linewidth',1.5)
plot(SNODAS_swes_elevation(:,1),[250:500:3500],'m','linewidth',1.5)

box on
set(gca,'linewidth',1,'fontsize',10)
xlim([-100 420])
ylim([0 3500])
text(-90,3500-3500*0.06,'(c)','fontsize',12,'fontweight','bold')
ylabel('Elevation (m)','fontsize',12)
xlabel('SWE (m)','fontsize',12)

%% plot 4
subplot('position', [0.59 0.12 0.39 0.38]);
hold on

inBetween = [ELM_swes_elevation(:,2)+ELM_swes_elevation_std(:,2); flipud(ELM_swes_elevation(:,2)-ELM_swes_elevation_std(:,2))];
patch(inBetween,[[250:500:3500]'; flipud([250:500:3500]')],'r','FaceColor', '#4DBEEE','LineStyle','none')
alpha(0.3)

inBetween = [UA_swes_elevation(:,2)+UA_swes_elevation_std(:,2); flipud(UA_swes_elevation(:,2)-UA_swes_elevation_std(:,2))];
patch(inBetween,[[250:500:3500]'; flipud([250:500:3500]')],'r','FaceColor', '#EDB120','LineStyle','none')
alpha(0.3)

inBetween = [SNODAS_swes_elevation(:,2)+SNODAS_swes_elevation_std(:,2); flipud(SNODAS_swes_elevation(:,2)-SNODAS_swes_elevation_std(:,2))];
patch(inBetween,[[250:500:3500]'; flipud([250:500:3500]')],'r','FaceColor', 'm','LineStyle','none')
alpha(0.3)

plot(ELM_swes_elevation(:,2),[250:500:3500],'color','#4DBEEE','linewidth',1.5)
plot(UA_swes_elevation(:,2),[250:500:3500],'color','#EDB120','linewidth',1.5)
plot(SNODAS_swes_elevation(:,2),[250:500:3500],'m','linewidth',1.5)
box on
xlim([-100 420])
ylim([0 3500])
text(-90,3500-3500*0.06,'(d)','fontsize',12,'fontweight','bold')

legend('ELM','UA','SNODAS','location','southeast','fontsize',8)
legend('boxoff')
set(gca,'linewidth',1,'fontsize',10)
xlabel('SWE (m)','fontsize',12)


print(gcf, '-dtiff', '-r300', '../../figure_all_tif/swe_trend_elevation_time.tif')

close all

