
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


load('../../../../data_processing/snow_timing/snow_timing_alldata_v2.mat');
load('../../../../data_processing/snow_timing/snow_timing_num_alldata_v2.mat');



start_of_snow_mean = nan(19,3);
melt_of_snow_mean = nan(19,3);
midpoint_of_snow_mean = nan(19,3);
end_of_snow_mean = nan(19,3);
duration_of_snow_mean = nan(19,3);

start_of_snow_std = nan(19,3);
melt_of_snow_std = nan(19,3);
midpoint_of_snow_std = nan(19,3);
end_of_snow_std = nan(19,3);
duration_of_snow_std = nan(19,3);


for year_i = 2001:2019
    
    index_i = year_i - 2000;
    %% sos
    ELM_tmp = squeeze(ELM_start_of_snow_all(:,:,index_i))-122+365;
    MODSCAG_tmp = squeeze(MODSCAG_start_of_snow_all(:,:,index_i))-122+365;
    SPIRES_tmp = squeeze(SPIRES_start_of_snow_all(:,:,index_i))-122+365;
    
    filters = ELM_start_of_snow_num<18 | isUS<1 | MODSCAG_start_of_snow_num<18 | SPIRES_start_of_snow_num<18;
    ELM_tmp(filters) = nan;
    MODSCAG_tmp(filters) = nan;
    SPIRES_tmp(filters) = nan;
    
    filters = ELM_tmp>=0 & MODSCAG_tmp>=0 & SPIRES_tmp>=0;
    area_weight = Areas(filters)./sum(Areas(filters));
    
    start_of_snow_mean(year_i-2000,1)= nansum(ELM_tmp(filters).*area_weight);
    start_of_snow_mean(year_i-2000,2) = nansum(MODSCAG_tmp(filters).*area_weight);
    start_of_snow_mean(year_i-2000,3) = nansum(SPIRES_tmp(filters).*area_weight);
    
    start_of_snow_std(year_i-2000,1)= std(ELM_tmp(filters),area_weight,1,'omitnan');
    start_of_snow_std(year_i-2000,2) = std(MODSCAG_tmp(filters),area_weight,1,'omitnan');
    start_of_snow_std(year_i-2000,3) = std(SPIRES_tmp(filters),area_weight,1,'omitnan');
    %% snow melt
    ELM_tmp = squeeze(ELM_melt_of_snow_all(:,:,index_i));
    MODSCAG_tmp = squeeze(MODSCAG_melt_of_snow_all(:,:,index_i));
    SPIRES_tmp = squeeze(SPIRES_melt_of_snow_all(:,:,index_i));
    
    filters = ELM_melt_of_snow_num<18 | isUS<1 | MODSCAG_melt_of_snow_num<18 | SPIRES_melt_of_snow_num<18;
    ELM_tmp(filters) = nan;
    MODSCAG_tmp(filters) = nan;
    SPIRES_tmp(filters) = nan;
    
    filters = ELM_tmp>=0 & MODSCAG_tmp>=0 & SPIRES_tmp>=0;
    area_weight = Areas(filters)./sum(Areas(filters));
    
    melt_of_snow_mean(year_i-2000,1)= nansum(ELM_tmp(filters).*area_weight);
    melt_of_snow_mean(year_i-2000,2) = nansum(MODSCAG_tmp(filters).*area_weight);
    melt_of_snow_mean(year_i-2000,3) = nansum(SPIRES_tmp(filters).*area_weight);
    
    melt_of_snow_std(year_i-2000,1)= std(ELM_tmp(filters),area_weight,1,'omitnan');
    melt_of_snow_std(year_i-2000,2) = std(MODSCAG_tmp(filters),area_weight,1,'omitnan');
    melt_of_snow_std(year_i-2000,3) = std(SPIRES_tmp(filters),area_weight,1,'omitnan');
    
    %% midpoint
    ELM_tmp = squeeze(ELM_midpoint_of_snow_all(:,:,index_i));
    MODSCAG_tmp = squeeze(MODSCAG_midpoint_of_snow_all(:,:,index_i));
    SPIRES_tmp = squeeze(SPIRES_midpoint_of_snow_all(:,:,index_i));
    
    filters = ELM_midpoint_of_snow_num<18 | isUS<1 | MODSCAG_midpoint_of_snow_num<18 | SPIRES_midpoint_of_snow_num<18;
    ELM_tmp(filters) = nan;
    MODSCAG_tmp(filters) = nan;
    SPIRES_tmp(filters) = nan;
    
    filters = ELM_tmp>=0 & MODSCAG_tmp>=0 & SPIRES_tmp>=0;
    area_weight = Areas(filters)./sum(Areas(filters));
    
    midpoint_of_snow_mean(year_i-2000,1)= nansum(ELM_tmp(filters).*area_weight);
    midpoint_of_snow_mean(year_i-2000,2) = nansum(MODSCAG_tmp(filters).*area_weight);
    midpoint_of_snow_mean(year_i-2000,3) = nansum(SPIRES_tmp(filters).*area_weight);
    
    midpoint_of_snow_std(year_i-2000,1)= std(ELM_tmp(filters),area_weight,1,'omitnan');
    midpoint_of_snow_std(year_i-2000,2) = std(MODSCAG_tmp(filters),area_weight,1,'omitnan');
    midpoint_of_snow_std(year_i-2000,3) = std(SPIRES_tmp(filters),area_weight,1,'omitnan');
    %% end of snow
    ELM_tmp = squeeze(ELM_end_of_snow_all(:,:,index_i))-122;
    MODSCAG_tmp = squeeze(MODSCAG_end_of_snow_all(:,:,index_i))-122;
    SPIRES_tmp = squeeze(SPIRES_end_of_snow_all(:,:,index_i))-122;
    
    filters = ELM_end_of_snow_num<18 | isUS<1 | MODSCAG_end_of_snow_num<18 | SPIRES_end_of_snow_num<18;
    ELM_tmp(filters) = nan;
    MODSCAG_tmp(filters) = nan;
    SPIRES_tmp(filters) = nan;
    
    filters = ELM_tmp>=0 & MODSCAG_tmp>=0 & SPIRES_tmp>=0;
    area_weight = Areas(filters)./sum(Areas(filters));
    
    end_of_snow_mean(year_i-2000,1)= nansum(ELM_tmp(filters).*area_weight);
    end_of_snow_mean(year_i-2000,2) = nansum(MODSCAG_tmp(filters).*area_weight);
    end_of_snow_mean(year_i-2000,3) = nansum(SPIRES_tmp(filters).*area_weight);
    
    end_of_snow_std(year_i-2000,1)= std(ELM_tmp(filters),area_weight,1,'omitnan');
    end_of_snow_std(year_i-2000,2) = std(MODSCAG_tmp(filters),area_weight,1,'omitnan');
    end_of_snow_std(year_i-2000,3) = std(SPIRES_tmp(filters),area_weight,1,'omitnan');
    %% duration
    ELM_tmp = squeeze(ELM_duration_of_snow_all(:,:,index_i));
    MODSCAG_tmp = squeeze(MODSCAG_duration_of_snow_all(:,:,index_i));
    SPIRES_tmp = squeeze(SPIRES_duration_of_snow_all(:,:,index_i));
    
    filters = ELM_duration_of_snow_num<18 | isUS<1 | MODSCAG_duration_of_snow_num<18 | SPIRES_duration_of_snow_num<18;
    ELM_tmp(filters) = nan;
    MODSCAG_tmp(filters) = nan;
    SPIRES_tmp(filters) = nan;
    
    filters = ELM_tmp>=0 & MODSCAG_tmp>=0 & SPIRES_tmp>=0;
    area_weight = Areas(filters)./sum(Areas(filters));
    
    duration_of_snow_mean(year_i-2000,1)= nansum(ELM_tmp(filters).*area_weight);
    duration_of_snow_mean(year_i-2000,2) = nansum(MODSCAG_tmp(filters).*area_weight);
    duration_of_snow_mean(year_i-2000,3) = nansum(SPIRES_tmp(filters).*area_weight);
    
    duration_of_snow_std(year_i-2000,1)= std(ELM_tmp(filters),area_weight,1,'omitnan');
    duration_of_snow_std(year_i-2000,2) = std(MODSCAG_tmp(filters),area_weight,1,'omitnan');
    duration_of_snow_std(year_i-2000,3) = std(SPIRES_tmp(filters),area_weight,1,'omitnan');
end


%% elevation distribution
elevations = imread('../WUS_mean_of_Elevation.tif');
elevations = squeeze(elevations(:,:,1));

elevations_bands = nan(size(elevations));
for elevation_band = 1:7
    filters = elevations >=((elevation_band-1)*500)  & elevations<(elevation_band*500);
    if elevation_band == 7
        filters = elevations >=((elevation_band-1)*500) ;
    end
    elevations_bands(filters) = elevation_band;
end


ELM_melt_of_snow_mean = nanmean(ELM_melt_of_snow_all,3);
ELM_midpoint_of_snow_mean = nanmean(ELM_midpoint_of_snow_all,3);
ELM_start_of_snow_mean = nanmean(ELM_start_of_snow_all,3)-122+365;
ELM_end_of_snow_mean = nanmean(ELM_end_of_snow_all,3)-122;
ELM_duration_of_snow_mean = nanmean(ELM_duration_of_snow_all,3);

MODSCAG_melt_of_snow_mean = nanmean(MODSCAG_melt_of_snow_all,3);
MODSCAG_midpoint_of_snow_mean = nanmean(MODSCAG_midpoint_of_snow_all,3);
MODSCAG_start_of_snow_mean = nanmean(MODSCAG_start_of_snow_all,3)-122+365;
MODSCAG_end_of_snow_mean = nanmean(MODSCAG_end_of_snow_all,3)-122;
MODSCAG_duration_of_snow_mean = nanmean(MODSCAG_duration_of_snow_all,3);

SPIRES_melt_of_snow_mean = nanmean(SPIRES_melt_of_snow_all,3);
SPIRES_midpoint_of_snow_mean = nanmean(SPIRES_midpoint_of_snow_all,3);
SPIRES_start_of_snow_mean = nanmean(SPIRES_start_of_snow_all,3)-122+365;
SPIRES_end_of_snow_mean = nanmean(SPIRES_end_of_snow_all,3)-122;
SPIRES_duration_of_snow_mean = nanmean(SPIRES_duration_of_snow_all,3);

%% nan<-10 observation
ELM_melt_of_snow_mean(ELM_melt_of_snow_num<10 | isUS<1) = nan;
ELM_midpoint_of_snow_mean(ELM_midpoint_of_snow_num<10 | isUS<1)  = nan;
ELM_start_of_snow_mean(ELM_start_of_snow_num<10 | isUS<1)  = nan;
ELM_end_of_snow_mean(ELM_end_of_snow_num<10 | isUS<1)  = nan;
ELM_duration_of_snow_mean(ELM_duration_of_snow_num<10 | isUS<1)  = nan;

MODSCAG_melt_of_snow_mean(MODSCAG_melt_of_snow_num<10 | isUS<1) = nan;
MODSCAG_midpoint_of_snow_mean(MODSCAG_midpoint_of_snow_num<10 | isUS<1)  = nan;
MODSCAG_start_of_snow_mean(MODSCAG_start_of_snow_num<10 | isUS<1)  = nan;
MODSCAG_end_of_snow_mean(MODSCAG_end_of_snow_num<10 | isUS<1)  = nan;
MODSCAG_duration_of_snow_mean(MODSCAG_duration_of_snow_num<10 | isUS<1)  = nan;

SPIRES_melt_of_snow_mean(SPIRES_melt_of_snow_num<10 | isUS<1) = nan;
SPIRES_midpoint_of_snow_mean(SPIRES_midpoint_of_snow_num<10 | isUS<1)  = nan;
SPIRES_start_of_snow_mean(SPIRES_start_of_snow_num<10 | isUS<1)  = nan;
SPIRES_end_of_snow_mean(SPIRES_end_of_snow_num<10 | isUS<1)  = nan;
SPIRES_duration_of_snow_mean(SPIRES_duration_of_snow_num<10 | isUS<1)  = nan;


%% sos
ELM_tmp = ELM_start_of_snow_mean;
MODSCAG_tmp = MODSCAG_start_of_snow_mean;
SPIRES_tmp = SPIRES_start_of_snow_mean;


filters = ELM_tmp<0 | MODSCAG_tmp<0 | SPIRES_tmp<0 | ...
    isnan(ELM_tmp) | isnan(MODSCAG_tmp)  | isnan(SPIRES_tmp);

ELM_tmp(isUS<1 | filters) = nan;
MODSCAGs_tmp(isUS<1 | filters) = nan;
SPIRES_tmp(isUS<1 | filters) = nan;

snow_data_all_start = [];

for elevation_band = 1:7
    filters = elevations_bands>=0 & ELM_tmp>0 & MODSCAG_tmp>0 & SPIRES_tmp>0;
    
    ELM_tmp2 = ELM_tmp(filters);
    MODSCAG_tmp2 = MODSCAG_tmp(filters);
    SPIRES_tmp2 = SPIRES_tmp(filters);
    ELM_tmp2(elevations_bands(filters)~=elevation_band) = nan;
    MODSCAG_tmp2(elevations_bands(filters)~=elevation_band) = nan;
    SPIRES_tmp2(elevations_bands(filters)~=elevation_band) = nan;
    
    tmp = [ELM_tmp2; MODSCAG_tmp2; SPIRES_tmp2];
    group_all_start = [ones(size(tmp,1)/3, 1); 2*ones(size(tmp,1)/3, 1); 3*ones(size(tmp,1)/3, 1)];
    
    snow_data_all_start =[snow_data_all_start tmp];
end

%% sos
ELM_tmp = ELM_melt_of_snow_mean;
MODSCAG_tmp = MODSCAG_melt_of_snow_mean;
SPIRES_tmp = SPIRES_melt_of_snow_mean;


filters = ELM_tmp<0 | MODSCAG_tmp<0 | SPIRES_tmp<0 | ...
    isnan(ELM_tmp) | isnan(MODSCAG_tmp)  | isnan(SPIRES_tmp);

ELM_tmp(isUS<1 | filters) = nan;
MODSCAGs_tmp(isUS<1 | filters) = nan;
SPIRES_tmp(isUS<1 | filters) = nan;

snow_data_all_melt = [];

for elevation_band = 1:7
    filters = elevations_bands>=0 & ELM_tmp>0 & MODSCAG_tmp>0 & SPIRES_tmp>0;
    
    ELM_tmp2 = ELM_tmp(filters);
    MODSCAG_tmp2 = MODSCAG_tmp(filters);
    SPIRES_tmp2 = SPIRES_tmp(filters);
    ELM_tmp2(elevations_bands(filters)~=elevation_band) = nan;
    MODSCAG_tmp2(elevations_bands(filters)~=elevation_band) = nan;
    SPIRES_tmp2(elevations_bands(filters)~=elevation_band) = nan;
    
    tmp = [ELM_tmp2; MODSCAG_tmp2; SPIRES_tmp2];
    group_all_melt = [ones(size(tmp,1)/3, 1); 2*ones(size(tmp,1)/3, 1); 3*ones(size(tmp,1)/3, 1)];
    
    snow_data_all_melt =[snow_data_all_melt tmp];
end

%% sos
ELM_tmp = ELM_midpoint_of_snow_mean;
MODSCAG_tmp = MODSCAG_midpoint_of_snow_mean;
SPIRES_tmp = SPIRES_midpoint_of_snow_mean;


filters = ELM_tmp<0 | MODSCAG_tmp<0 | SPIRES_tmp<0 | ...
    isnan(ELM_tmp) | isnan(MODSCAG_tmp)  | isnan(SPIRES_tmp);

ELM_tmp(isUS<1 | filters) = nan;
MODSCAGs_tmp(isUS<1 | filters) = nan;
SPIRES_tmp(isUS<1 | filters) = nan;

snow_data_all_midpoint = [];

for elevation_band = 1:7
    filters = elevations_bands>=0 & ELM_tmp>0 & MODSCAG_tmp>0 & SPIRES_tmp>0;
    
    ELM_tmp2 = ELM_tmp(filters);
    MODSCAG_tmp2 = MODSCAG_tmp(filters);
    SPIRES_tmp2 = SPIRES_tmp(filters);
    ELM_tmp2(elevations_bands(filters)~=elevation_band) = nan;
    MODSCAG_tmp2(elevations_bands(filters)~=elevation_band) = nan;
    SPIRES_tmp2(elevations_bands(filters)~=elevation_band) = nan;
    
    tmp = [ELM_tmp2; MODSCAG_tmp2; SPIRES_tmp2];
    group_all_midpoint = [ones(size(tmp,1)/3, 1); 2*ones(size(tmp,1)/3, 1); 3*ones(size(tmp,1)/3, 1)];
    
    snow_data_all_midpoint =[snow_data_all_midpoint tmp];
end

%% sos
ELM_tmp = ELM_end_of_snow_mean;
MODSCAG_tmp = MODSCAG_end_of_snow_mean;
SPIRES_tmp = SPIRES_end_of_snow_mean;


filters = ELM_tmp<0 | MODSCAG_tmp<0 | SPIRES_tmp<0 | ...
    isnan(ELM_tmp) | isnan(MODSCAG_tmp)  | isnan(SPIRES_tmp);

ELM_tmp(isUS<1 | filters) = nan;
MODSCAGs_tmp(isUS<1 | filters) = nan;
SPIRES_tmp(isUS<1 | filters) = nan;

snow_data_all_end = [];

for elevation_band = 1:7
    filters = elevations_bands>=0 & ELM_tmp>0 & MODSCAG_tmp>0 & SPIRES_tmp>0;
    
    ELM_tmp2 = ELM_tmp(filters);
    MODSCAG_tmp2 = MODSCAG_tmp(filters);
    SPIRES_tmp2 = SPIRES_tmp(filters);
    ELM_tmp2(elevations_bands(filters)~=elevation_band) = nan;
    MODSCAG_tmp2(elevations_bands(filters)~=elevation_band) = nan;
    SPIRES_tmp2(elevations_bands(filters)~=elevation_band) = nan;
    
    tmp = [ELM_tmp2; MODSCAG_tmp2; SPIRES_tmp2];
    group_all_end = [ones(size(tmp,1)/3, 1); 2*ones(size(tmp,1)/3, 1); 3*ones(size(tmp,1)/3, 1)];
    
    snow_data_all_end =[snow_data_all_end tmp];
end

%% sos
ELM_tmp = ELM_duration_of_snow_mean;
MODSCAG_tmp = MODSCAG_duration_of_snow_mean;
SPIRES_tmp = SPIRES_duration_of_snow_mean;


filters = ELM_tmp<0 | MODSCAG_tmp<0 | SPIRES_tmp<0 | ...
    isnan(ELM_tmp) | isnan(MODSCAG_tmp)  | isnan(SPIRES_tmp);

ELM_tmp(isUS<1 | filters) = nan;
MODSCAGs_tmp(isUS<1 | filters) = nan;
SPIRES_tmp(isUS<1 | filters) = nan;

snow_data_all_duration = [];

for elevation_band = 1:7
    filters = elevations_bands>=0 & ELM_tmp>0 & MODSCAG_tmp>0 & SPIRES_tmp>0;
    
    ELM_tmp2 = ELM_tmp(filters);
    MODSCAG_tmp2 = MODSCAG_tmp(filters);
    SPIRES_tmp2 = SPIRES_tmp(filters);
    ELM_tmp2(elevations_bands(filters)~=elevation_band) = nan;
    MODSCAG_tmp2(elevations_bands(filters)~=elevation_band) = nan;
    SPIRES_tmp2(elevations_bands(filters)~=elevation_band) = nan;
    
    tmp = [ELM_tmp2; MODSCAG_tmp2; SPIRES_tmp2];
    group_all_duration = [ones(size(tmp,1)/3, 1); 2*ones(size(tmp,1)/3, 1); 3*ones(size(tmp,1)/3, 1)];
    
    snow_data_all_duration =[snow_data_all_duration tmp];
end

%% forest cover
%% elevation distribution
%% elevation distribution
forest_cover = imread('../forest_cover.tif');

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


ELM_melt_of_snow_mean = nanmean(ELM_melt_of_snow_all,3);
ELM_midpoint_of_snow_mean = nanmean(ELM_midpoint_of_snow_all,3);
ELM_start_of_snow_mean = nanmean(ELM_start_of_snow_all,3)-122+365;
ELM_end_of_snow_mean = nanmean(ELM_end_of_snow_all,3)-122;
ELM_duration_of_snow_mean = nanmean(ELM_duration_of_snow_all,3);

MODSCAG_melt_of_snow_mean = nanmean(MODSCAG_melt_of_snow_all,3);
MODSCAG_midpoint_of_snow_mean = nanmean(MODSCAG_midpoint_of_snow_all,3);
MODSCAG_start_of_snow_mean = nanmean(MODSCAG_start_of_snow_all,3)-122+365;
MODSCAG_end_of_snow_mean = nanmean(MODSCAG_end_of_snow_all,3)-122;
MODSCAG_duration_of_snow_mean = nanmean(MODSCAG_duration_of_snow_all,3);

SPIRES_melt_of_snow_mean = nanmean(SPIRES_melt_of_snow_all,3);
SPIRES_midpoint_of_snow_mean = nanmean(SPIRES_midpoint_of_snow_all,3);
SPIRES_start_of_snow_mean = nanmean(SPIRES_start_of_snow_all,3)-122+365;
SPIRES_end_of_snow_mean = nanmean(SPIRES_end_of_snow_all,3)-122;
SPIRES_duration_of_snow_mean = nanmean(SPIRES_duration_of_snow_all,3);

%% nan<-10 observation
ELM_melt_of_snow_mean(ELM_melt_of_snow_num<10 | isUS<1) = nan;
ELM_midpoint_of_snow_mean(ELM_midpoint_of_snow_num<10 | isUS<1)  = nan;
ELM_start_of_snow_mean(ELM_start_of_snow_num<10 | isUS<1)  = nan;
ELM_end_of_snow_mean(ELM_end_of_snow_num<10 | isUS<1)  = nan;
ELM_duration_of_snow_mean(ELM_duration_of_snow_num<10 | isUS<1)  = nan;

MODSCAG_melt_of_snow_mean(MODSCAG_melt_of_snow_num<10 | isUS<1) = nan;
MODSCAG_midpoint_of_snow_mean(MODSCAG_midpoint_of_snow_num<10 | isUS<1)  = nan;
MODSCAG_start_of_snow_mean(MODSCAG_start_of_snow_num<10 | isUS<1)  = nan;
MODSCAG_end_of_snow_mean(MODSCAG_end_of_snow_num<10 | isUS<1)  = nan;
MODSCAG_duration_of_snow_mean(MODSCAG_duration_of_snow_num<10 | isUS<1)  = nan;

SPIRES_melt_of_snow_mean(SPIRES_melt_of_snow_num<10 | isUS<1) = nan;
SPIRES_midpoint_of_snow_mean(SPIRES_midpoint_of_snow_num<10 | isUS<1)  = nan;
SPIRES_start_of_snow_mean(SPIRES_start_of_snow_num<10 | isUS<1)  = nan;
SPIRES_end_of_snow_mean(SPIRES_end_of_snow_num<10 | isUS<1)  = nan;
SPIRES_duration_of_snow_mean(SPIRES_duration_of_snow_num<10 | isUS<1)  = nan;


%% sos
ELM_tmp = ELM_start_of_snow_mean;
MODSCAG_tmp = MODSCAG_start_of_snow_mean;
SPIRES_tmp = SPIRES_start_of_snow_mean;


filters = ELM_tmp<0 | MODSCAG_tmp<0 | SPIRES_tmp<0 | ...
    isnan(ELM_tmp) | isnan(MODSCAG_tmp)  | isnan(SPIRES_tmp);

ELM_tmp(isUS<1 | filters) = nan;
MODSCAGs_tmp(isUS<1 | filters) = nan;
SPIRES_tmp(isUS<1 | filters) = nan;

snow_data_all2_start = [];

for forest_band = 1:5
    filters = fc_bands>=0 & ELM_tmp>0 & MODSCAG_tmp>0 & SPIRES_tmp>0;
    
    ELM_tmp2 = ELM_tmp(filters);
    MODSCAG_tmp2 = MODSCAG_tmp(filters);
    SPIRES_tmp2 = SPIRES_tmp(filters);
    ELM_tmp2(fc_bands(filters)~=forest_band) = nan;
    MODSCAG_tmp2(fc_bands(filters)~=forest_band) = nan;
    SPIRES_tmp2(fc_bands(filters)~=forest_band) = nan;
    
    tmp = [ELM_tmp2; MODSCAG_tmp2; SPIRES_tmp2];
    group_all2_start = [ones(size(tmp,1)/3, 1); 2*ones(size(tmp,1)/3, 1); 3*ones(size(tmp,1)/3, 1)];
    
    snow_data_all2_start =[snow_data_all2_start tmp];
end

%% sos
ELM_tmp = ELM_melt_of_snow_mean;
MODSCAG_tmp = MODSCAG_melt_of_snow_mean;
SPIRES_tmp = SPIRES_melt_of_snow_mean;


filters = ELM_tmp<0 | MODSCAG_tmp<0 | SPIRES_tmp<0 | ...
    isnan(ELM_tmp) | isnan(MODSCAG_tmp)  | isnan(SPIRES_tmp);

ELM_tmp(isUS<1 | filters) = nan;
MODSCAGs_tmp(isUS<1 | filters) = nan;
SPIRES_tmp(isUS<1 | filters) = nan;

snow_data_all2_melt = [];

for forest_band = 1:5
    filters = fc_bands>=0 & ELM_tmp>0 & MODSCAG_tmp>0 & SPIRES_tmp>0;
    
    ELM_tmp2 = ELM_tmp(filters);
    MODSCAG_tmp2 = MODSCAG_tmp(filters);
    SPIRES_tmp2 = SPIRES_tmp(filters);
    ELM_tmp2(fc_bands(filters)~=forest_band) = nan;
    MODSCAG_tmp2(fc_bands(filters)~=forest_band) = nan;
    SPIRES_tmp2(fc_bands(filters)~=forest_band) = nan;
    
    tmp = [ELM_tmp2; MODSCAG_tmp2; SPIRES_tmp2];
    group_all2_melt = [ones(size(tmp,1)/3, 1); 2*ones(size(tmp,1)/3, 1); 3*ones(size(tmp,1)/3, 1)];
    
    snow_data_all2_melt =[snow_data_all2_melt tmp];
end

%% sos
ELM_tmp = ELM_midpoint_of_snow_mean;
MODSCAG_tmp = MODSCAG_midpoint_of_snow_mean;
SPIRES_tmp = SPIRES_midpoint_of_snow_mean;


filters = ELM_tmp<0 | MODSCAG_tmp<0 | SPIRES_tmp<0 | ...
    isnan(ELM_tmp) | isnan(MODSCAG_tmp)  | isnan(SPIRES_tmp);

ELM_tmp(isUS<1 | filters) = nan;
MODSCAGs_tmp(isUS<1 | filters) = nan;
SPIRES_tmp(isUS<1 | filters) = nan;

snow_data_all2_midpoint = [];

for forest_band = 1:5
    filters = fc_bands>=0 & ELM_tmp>0 & MODSCAG_tmp>0 & SPIRES_tmp>0;
    
    ELM_tmp2 = ELM_tmp(filters);
    MODSCAG_tmp2 = MODSCAG_tmp(filters);
    SPIRES_tmp2 = SPIRES_tmp(filters);
    ELM_tmp2(fc_bands(filters)~=forest_band) = nan;
    MODSCAG_tmp2(fc_bands(filters)~=forest_band) = nan;
    SPIRES_tmp2(fc_bands(filters)~=forest_band) = nan;
    
    tmp = [ELM_tmp2; MODSCAG_tmp2; SPIRES_tmp2];
    group_all2_midpoint = [ones(size(tmp,1)/3, 1); 2*ones(size(tmp,1)/3, 1); 3*ones(size(tmp,1)/3, 1)];
    
    snow_data_all2_midpoint =[snow_data_all2_midpoint tmp];
end

%% sos
ELM_tmp = ELM_end_of_snow_mean;
MODSCAG_tmp = MODSCAG_end_of_snow_mean;
SPIRES_tmp = SPIRES_end_of_snow_mean;


filters = ELM_tmp<0 | MODSCAG_tmp<0 | SPIRES_tmp<0 | ...
    isnan(ELM_tmp) | isnan(MODSCAG_tmp)  | isnan(SPIRES_tmp);

ELM_tmp(isUS<1 | filters) = nan;
MODSCAGs_tmp(isUS<1 | filters) = nan;
SPIRES_tmp(isUS<1 | filters) = nan;

snow_data_all2_end = [];

for forest_band = 1:5
    filters = fc_bands>=0 & ELM_tmp>0 & MODSCAG_tmp>0 & SPIRES_tmp>0;
    
    ELM_tmp2 = ELM_tmp(filters);
    MODSCAG_tmp2 = MODSCAG_tmp(filters);
    SPIRES_tmp2 = SPIRES_tmp(filters);
    ELM_tmp2(fc_bands(filters)~=forest_band) = nan;
    MODSCAG_tmp2(fc_bands(filters)~=forest_band) = nan;
    SPIRES_tmp2(fc_bands(filters)~=forest_band) = nan;
    
    tmp = [ELM_tmp2; MODSCAG_tmp2; SPIRES_tmp2];
    group_all2_end = [ones(size(tmp,1)/3, 1); 2*ones(size(tmp,1)/3, 1); 3*ones(size(tmp,1)/3, 1)];
    
    snow_data_all2_end =[snow_data_all2_end tmp];
end

%% sos
ELM_tmp = ELM_duration_of_snow_mean;
MODSCAG_tmp = MODSCAG_duration_of_snow_mean;
SPIRES_tmp = SPIRES_duration_of_snow_mean;


filters = ELM_tmp<0 | MODSCAG_tmp<0 | SPIRES_tmp<0 | ...
    isnan(ELM_tmp) | isnan(MODSCAG_tmp)  | isnan(SPIRES_tmp);

ELM_tmp(isUS<1 | filters) = nan;
MODSCAGs_tmp(isUS<1 | filters) = nan;
SPIRES_tmp(isUS<1 | filters) = nan;

snow_data_all2_duration = [];

for forest_band = 1:5
    filters = fc_bands>=0 & ELM_tmp>0 & MODSCAG_tmp>0 & SPIRES_tmp>0;
    
    ELM_tmp2 = ELM_tmp(filters);
    MODSCAG_tmp2 = MODSCAG_tmp(filters);
    SPIRES_tmp2 = SPIRES_tmp(filters);
    ELM_tmp2(fc_bands(filters)~=forest_band) = nan;
    MODSCAG_tmp2(fc_bands(filters)~=forest_band) = nan;
    SPIRES_tmp2(fc_bands(filters)~=forest_band) = nan;
    
    tmp = [ELM_tmp2; MODSCAG_tmp2; SPIRES_tmp2];
    group_all2_duration = [ones(size(tmp,1)/3, 1); 2*ones(size(tmp,1)/3, 1); 3*ones(size(tmp,1)/3, 1)];
    
    snow_data_all2_duration =[snow_data_all2_duration tmp];
end
%% figure

figure
set(gcf,'unit','normalized','position',[0.1,0.1,0.8,0.9]);
set(gca, 'Position', [0 0 1 1])

%% plot 1
subplot('position', [0.04 0.81 0.3 0.17]);
hold on

plot(2002:2019,start_of_snow_mean(2:end,1),'color',colors(1,:),'linewidth',1.5)
plot(2002:2019,start_of_snow_mean(2:end,2),'color',colors(2,:),'linewidth',1.5)
plot(2002:2019,start_of_snow_mean(2:end,3),'color',colors(3,:),'linewidth',1.5)


%title('SOnDt')
box on
set(gca,'linewidth',1,'fontsize',8)
xlim([2002 2019])
ylim([280 360])
text(2002.5,360-80*0.06,'(a)','fontsize',11,'fontweight','bold')
ylabel('Accumulation\_onset\_date','fontsize',8,'fontweight','bold')
IAV1 = std(start_of_snow_mean(2:end,1)-nanmean(start_of_snow_mean(:,1)),1);
IAV2 = std(start_of_snow_mean(2:end,2)-nanmean(start_of_snow_mean(:,2)),1);
IAV3 = std(start_of_snow_mean(2:end,3)-nanmean(start_of_snow_mean(:,3)),1);
text(2002.7,360-80*0.18,['IAV\_ELM: ' num2str(IAV1,'%4.1f')],'fontsize',10,'fontweight','bold','color',colors(1,:))
text(2002.7,360-80*0.28,['IAV\_STC-MODSCAG: ' num2str(IAV2,'%4.1f')],'fontsize',10,'fontweight','bold','color',colors(2,:))
text(2002.7,360-80*0.38,['IAV\_SPIReS: ' num2str(IAV3,'%4.1f')],'fontsize',10,'fontweight','bold','color',colors(3,:))


sen_ps1 = nan(3,2);
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' start_of_snow_mean(:,1)], 0.05, 0);
sen_ps1(1,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' start_of_snow_mean(:,2)], 0.05, 0);
sen_ps1(2,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' start_of_snow_mean(:,3)], 0.05, 0);
sen_ps1(3,:) = [sen sig];

set(gca,'xticklabel',{})


subplot('position', [0.365 0.81 0.3 0.17]);
hold on

group_names = {'Humans', 'Dogs' , 'God'};
condition_names = {'0.5-1.0', '1.0-1.5', '1.5-2.0', '2.0-2.5' , '2.5-3.0', '>3.0'};


h = daboxplot(snow_data_all_start(:,2:end),'groups',group_all_start(:),'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',8,'mean',1, 'outliers',0);



box on
set(gca,'linewidth',1,'fontsize',8)
xlim([0.5 6.5])
ylim([240 380])
text(0.6,380-140*0.06,'(b)','fontsize',11,'fontweight','bold')
%xlabel('Elevation (km)','fontsize',11)
%ylabel('DOY','fontsize',10)
set(gca,'xticklabel',{})

subplot('position', [0.69 0.81 0.3 0.17]);
hold on

group_names = {'Humans', 'Dogs' , 'God'};
condition_names = {'0-10', '10-30', '30-50', '50-70' , '>70'};


h = daboxplot(snow_data_all2_start(:,1:end),'groups',group_all2_start(:),'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',8,'mean',1, 'outliers',0);



box on
set(gca,'linewidth',1,'fontsize',8)
xlim([0.5 5.5])
ylim([240 380])
text(0.6,380-140*0.06,'(c)','fontsize',11,'fontweight','bold')
%xlabel('Elevation (km)','fontsize',11)
%ylabel('DOY','fontsize',10)
set(gca,'xticklabel',{})
%% plot 2
subplot('position', [0.04 0.62 0.3 0.17]);
hold on

plot(2002:2019,melt_of_snow_mean(2:end,1),'color',colors(1,:),'linewidth',1.5)
plot(2002:2019,melt_of_snow_mean(2:end,2),'color',colors(2,:),'linewidth',1.5)
plot(2002:2019,melt_of_snow_mean(2:end,3),'color',colors(3,:),'linewidth',1.5)


%title('')
box on
set(gca,'linewidth',1,'fontsize',8)
xlim([2002 2019])
ylim([60 160])
ylabel('Depletion\_onset\_date','fontsize',8,'fontweight','bold')

text(2002.5,160-100*0.06,'(d)','fontsize',11,'fontweight','bold')
IAV1 = std(melt_of_snow_mean(2:end,1)-nanmean(melt_of_snow_mean(:,1)),1);
IAV2 = std(melt_of_snow_mean(2:end,2)-nanmean(melt_of_snow_mean(:,2)),1);
IAV3 = std(melt_of_snow_mean(2:end,3)-nanmean(melt_of_snow_mean(:,3)),1);
text(2002.7,160-100*0.18,['IAV\_ELM: ' num2str(IAV1,'%4.1f')],'fontsize',10,'fontweight','bold','color',colors(1,:))
text(2002.7,160-100*0.28,['IAV\_STC-MODSCAG: ' num2str(IAV2,'%4.1f')],'fontsize',10,'fontweight','bold','color',colors(2,:))
text(2002.7,160-100*0.38,['IAV\_SPIReS: ' num2str(IAV3,'%4.1f')],'fontsize',10,'fontweight','bold','color',colors(3,:))


sen_ps1 = nan(3,2);
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' melt_of_snow_mean(:,1)], 0.05, 0);
sen_ps1(1,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' melt_of_snow_mean(:,2)], 0.05, 0);
sen_ps1(2,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' melt_of_snow_mean(:,3)], 0.05, 0);
sen_ps1(3,:) = [sen sig];

set(gca,'xticklabel',{})

subplot('position', [0.365 0.62 0.3 0.17]);
hold on

group_names = {'Humans', 'Dogs' , 'God'};
condition_names = {'0.5-1.0', '1.0-1.5', '1.5-2.0', '2.0-2.5' , '2.5-3.0', '>3.0'};


h = daboxplot(snow_data_all_melt(:,2:end),'groups',group_all_melt(:),'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',8,'mean',1, 'outliers',0,'legend',{'ELM','MODSCAG','SPIReS'});



box on
set(gca,'linewidth',1,'fontsize',8)
xlim([0.5 6.5])
ylim([30 200])
text(0.6,200-170*0.06,'(e)','fontsize',11,'fontweight','bold')
%xlabel('Elevation (km)','fontsize',11)

set(gca,'xticklabel',{})


subplot('position', [0.69 0.62 0.3 0.17]);
hold on

group_names = {'Humans', 'Dogs' , 'God'};
condition_names = {'0-10', '10-30', '30-50', '50-70' , '>70'};


h = daboxplot(snow_data_all2_melt(:,1:end),'groups',group_all2_melt(:),'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',8,'mean',1, 'outliers',0);



box on
set(gca,'linewidth',1,'fontsize',8)
xlim([0.5 5.5])
ylim([30 200])
text(0.6,200-170*0.06,'(f)','fontsize',11,'fontweight','bold')
%xlabel('Elevation (km)','fontsize',11)

set(gca,'xticklabel',{})
%% plot 3
subplot('position', [0.04 0.43 0.3 0.17]);
hold on

plot(2002:2019,midpoint_of_snow_mean(2:end,1),'color',colors(1,:),'linewidth',1.5)
plot(2002:2019,midpoint_of_snow_mean(2:end,2),'color',colors(2,:),'linewidth',1.5)
plot(2002:2019,midpoint_of_snow_mean(2:end,3),'color',colors(3,:),'linewidth',1.5)


%title('')

ylabel('Midpoint\_date','fontsize',8,'fontweight','bold')

box on
set(gca,'linewidth',1,'fontsize',8)
xlim([2002 2019])
ylim([110 200])
text(2002.5,200-90*0.06,'(g)','fontsize',11,'fontweight','bold')
IAV1 = std(midpoint_of_snow_mean(2:end,1)-nanmean(midpoint_of_snow_mean(:,1)),1);
IAV2 = std(midpoint_of_snow_mean(2:end,2)-nanmean(midpoint_of_snow_mean(:,2)),1);
IAV3 = std(midpoint_of_snow_mean(2:end,3)-nanmean(midpoint_of_snow_mean(:,3)),1);
text(2002.7,200-90*0.18,['IAV\_ELM: ' num2str(IAV1,'%4.1f')],'fontsize',10,'fontweight','bold','color',colors(1,:))
text(2002.7,200-90*0.28,['IAV\_STC-MODSCAG: ' num2str(IAV2,'%4.1f')],'fontsize',10,'fontweight','bold','color',colors(2,:))
text(2002.7,200-90*0.38,['IAV\_SPIReS: ' num2str(IAV3,'%4.1f')],'fontsize',10,'fontweight','bold','color',colors(3,:))


sen_ps1 = nan(3,2);
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' midpoint_of_snow_mean(:,1)], 0.05, 0);
sen_ps1(1,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' midpoint_of_snow_mean(:,2)], 0.05, 0);
sen_ps1(2,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' midpoint_of_snow_mean(:,3)], 0.05, 0);
sen_ps1(3,:) = [sen sig];

set(gca,'xticklabel',{})


subplot('position', [0.365 0.43 0.3 0.17]);
hold on

group_names = {'Humans', 'Dogs' , 'God'};
condition_names = {'0.5-1.0', '1.0-1.5', '1.5-2.0', '2.0-2.5' , '2.5-3.0', '>3.0'};


h = daboxplot(snow_data_all_midpoint(:,2:end),'groups',group_all_midpoint(:),'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',8,'mean',1, 'outliers',0);



box on
set(gca,'linewidth',1,'fontsize',8)
xlim([0.5 6.5])
ylim([40 200])
text(0.6,200-160*0.06,'(h)','fontsize',11,'fontweight','bold')
%xlabel('Elevation (km)','fontsize',11)
set(gca,'xticklabel',{})

subplot('position', [0.69 0.43 0.3 0.17]);
hold on

group_names = {'Humans', 'Dogs' , 'God'};
condition_names = {'0-10', '10-30', '30-50', '50-70' , '>70'};


h = daboxplot(snow_data_all2_midpoint(:,1:end),'groups',group_all2_midpoint(:),'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',8,'mean',1, 'outliers',0);



box on
set(gca,'linewidth',1,'fontsize',8)
xlim([0.5 5.5])
ylim([40 200])
text(0.6,200-160*0.06,'(i)','fontsize',11,'fontweight','bold')
%xlabel('Elevation (km)','fontsize',11)

set(gca,'xticklabel',{})
%% plot 4
subplot('position', [0.04 0.24 0.3 0.17]);
hold on

plot(2002:2019,end_of_snow_mean(2:end,1),'color',colors(1,:),'linewidth',1.5)
plot(2002:2019,end_of_snow_mean(2:end,2),'color',colors(2,:),'linewidth',1.5)
plot(2002:2019,end_of_snow_mean(2:end,3),'color',colors(3,:),'linewidth',1.5)


%title('')
box on
set(gca,'linewidth',1,'fontsize',8)
xlim([2002 2019])
ylim([80 180])
text(2002.5,180-80*0.06,'(j)','fontsize',11,'fontweight','bold')

ylabel('End\_date','fontsize',8,'fontweight','bold')

IAV1 = std(end_of_snow_mean(2:end,1)-nanmean(end_of_snow_mean(:,1)),1);
IAV2 = std(end_of_snow_mean(2:end,2)-nanmean(end_of_snow_mean(:,2)),1);
IAV3 = std(end_of_snow_mean(2:end,3)-nanmean(end_of_snow_mean(:,3)),1);
text(2002.7,180-100*0.18,['IAV\_ELM: ' num2str(IAV1,'%4.1f')],'fontsize',10,'fontweight','bold','color',colors(1,:))
text(2002.7,180-100*0.28,['IAV\_STC-MODSCAG: ' num2str(IAV2,'%4.1f')],'fontsize',10,'fontweight','bold','color',colors(2,:))
text(2002.7,180-100*0.38,['IAV\_SPIReS: ' num2str(IAV3,'%4.1f')],'fontsize',10,'fontweight','bold','color',colors(3,:))


sen_ps1 = nan(3,2);
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' end_of_snow_mean(:,1)], 0.05, 0);
sen_ps1(1,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' end_of_snow_mean(:,2)], 0.05, 0);
sen_ps1(2,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' end_of_snow_mean(:,3)], 0.05, 0);
sen_ps1(3,:) = [sen sig];


set(gca,'xticklabel',{})

subplot('position', [0.365 0.24 0.3 0.17]);
hold on

group_names = {'Humans', 'Dogs' , 'God'};
condition_names = {'0.5-1.0', '1.0-1.5', '1.5-2.0', '2.0-2.5' , '2.5-3.0', '>3.0'};


h = daboxplot(snow_data_all_end(:,2:end),'groups',group_all_end(:),'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',8,'mean',1, 'outliers',0);



box on
set(gca,'linewidth',1,'fontsize',8)
xlim([0.5 6.5])
ylim([20 220])
text(0.6,220-200*0.06,'(k)','fontsize',11,'fontweight','bold')
%xlabel('Elevation (km)','fontsize',11)
set(gca,'xticklabel',{})

subplot('position', [0.69 0.24 0.3 0.17]);
hold on

group_names = {'Humans', 'Dogs' , 'God'};
condition_names = {'0-10', '10-30', '30-50', '50-70' , '>70'};


h = daboxplot(snow_data_all2_end(:,1:end),'groups',group_all2_end(:),'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',8,'mean',1, 'outliers',0);



box on
set(gca,'linewidth',1,'fontsize',8)
xlim([0.5 5.5])
ylim([20 220])
text(0.6,220-200*0.06,'(l)','fontsize',11,'fontweight','bold')
%xlabel('Elevation (km)','fontsize',11)
set(gca,'xticklabel',{})
%% plot 5
subplot('position', [0.04 0.05 0.3 0.17]);
hold on

plot(2002:2019,duration_of_snow_mean(2:end,1),'color',colors(1,:),'linewidth',1.5)
plot(2002:2019,duration_of_snow_mean(2:end,2),'color',colors(2,:),'linewidth',1.5)
plot(2002:2019,duration_of_snow_mean(2:end,3),'color',colors(3,:),'linewidth',1.5)


%title('')
box on
set(gca,'linewidth',1,'fontsize',8)
xlim([2002 2019])
ylim([130 270])
text(2002.5,270-130*0.06,'(m)','fontsize',11,'fontweight','bold')

ylabel('Duration','fontsize',8,'fontweight','bold')


IAV1 = std(duration_of_snow_mean(2:end,1)-nanmean(duration_of_snow_mean(:,1)),1);
IAV2 = std(duration_of_snow_mean(2:end,2)-nanmean(duration_of_snow_mean(:,2)),1);
IAV3 = std(duration_of_snow_mean(2:end,3)-nanmean(duration_of_snow_mean(:,3)),1);
text(2002.7,270-140*0.18,['IAV\_ELM: ' num2str(IAV1,'%4.1f')],'fontsize',10,'fontweight','bold','color',colors(1,:))
text(2002.7,270-140*0.28,['IAV\_STC-MODSCAG: ' num2str(IAV2,'%4.1f')],'fontsize',10,'fontweight','bold','color',colors(2,:))
text(2002.7,270-140*0.38,['IAV\_SPIReS: ' num2str(IAV3,'%4.1f')],'fontsize',10,'fontweight','bold','color',colors(3,:))


sen_ps1 = nan(3,2);
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' duration_of_snow_mean(:,1)], 0.05, 0);
sen_ps1(1,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' duration_of_snow_mean(:,2)], 0.05, 0);
sen_ps1(2,:) = [sen sig];
[taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub([[2001:2019]' duration_of_snow_mean(:,3)], 0.05, 0);
sen_ps1(3,:) = [sen sig];

xlabel('Year','fontsize',12,'fontweight','bold')


subplot('position', [0.365 0.05 0.3 0.17]);
hold on

group_names = {'Humans', 'Dogs' , 'God'};
condition_names = {'0.5-1.0', '1.0-1.5', '1.5-2.0', '2.0-2.5' , '2.5-3.0', '\geq3.0'};


h = daboxplot(snow_data_all_duration(:,2:end),'groups',group_all_duration(:),'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',8,'mean',1, 'outliers',0);



box on
set(gca,'linewidth',1,'fontsize',8)
xlim([0.5 6.5])
ylim([30 300])
text(0.6,300-270*0.06,'(n)','fontsize',11,'fontweight','bold')
xlabel('Elevation (km)','fontsize',11,'fontweight','bold')

subplot('position', [0.69 0.05 0.3 0.17]);
hold on

group_names = {'Humans', 'Dogs' , 'God'};
condition_names = {'0-10', '10-30', '30-50', '50-70' , '\geq70'};


h = daboxplot(snow_data_all2_duration(:,1:end),'groups',group_all2_duration(:),'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',8,'mean',1, 'outliers',0);



box on
set(gca,'linewidth',1,'fontsize',8)
xlim([0.5 5.5])
ylim([30 300])
text(0.6,300-270*0.06,'(o)','fontsize',11,'fontweight','bold')
xlabel('Forest cover (%)','fontsize',11,'fontweight','bold')


print(gcf, '-dtiff', '-r300', '../../figure_all_tif/SnowPhenology_trend_elevation_time.tif')

close all

