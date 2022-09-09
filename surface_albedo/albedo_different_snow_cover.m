
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

GridNames = {'Winter', 'Spring', 'Summer', 'Autumn'};


%% plot
colors_1 = flipud(brewermap(101, 'Spectral'));
colors_2 =  (brewermap(101, 'BrBG'));
colors_3 =  (brewermap(101, 'BuGn'));


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
    
    ELM_fsnos_tmp = nanmean(ELM_fsnos(:,:,filters),3);    
    ELM_fsnos_tmp(isUS<1) = nan;

filters = ~isnan(ELM_fsnos_tmp) & ~isnan(ELM_surfacealbedos_tmp) & ~isnan(MCD43_surfacealbedos_tmp);   
ELM_fsnos_all = ELM_fsnos_tmp(filters);
ELM_surfacealbedos_all = ELM_surfacealbedos_tmp(filters);
MCD43_surfacealbedos_all = MCD43_surfacealbedos_tmp(filters);

ELM_fsno_level = nan(size(ELM_fsnos_all));
ELM_fsno_level(ELM_fsnos_all>=0 & ELM_fsnos_all<0.2) = 1;
ELM_fsno_level(ELM_fsnos_all>=0.2 & ELM_fsnos_all<0.4) = 2;
ELM_fsno_level(ELM_fsnos_all>=0.4 & ELM_fsnos_all<0.6) = 3;
ELM_fsno_level(ELM_fsnos_all>=0.6 & ELM_fsnos_all<0.8) = 4;
ELM_fsno_level(ELM_fsnos_all>=0.8 & ELM_fsnos_all<=1.0) = 5;

AllData = [ELM_fsno_level ELM_fsnos_all ELM_surfacealbedos_all MCD43_surfacealbedos_all];
    
save([GridNames{season_i} '.mat'],'AllData');
end


