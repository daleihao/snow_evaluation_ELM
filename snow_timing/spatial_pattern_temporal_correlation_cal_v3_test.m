
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




%% plot
colors_1 = flipud(brewermap(10, 'Spectral'));
colors_2 =  flipud(brewermap(11, 'RdBu'));
colors_3 =  flipud((brewermap(10, 'Spectral')));

colors_2(6,:) = [0.85 0.85 0.85];


load('../../../../data_processing/snow_timing/snow_timing_alldata_v2_r1_0_05.mat');
load('../../../../data_processing/snow_timing/snow_timing_num_alldata_v2_r1_0_05.mat');

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


alloutputs = nan(10,5);
alloutputs(1,:) = calStat(ELM_start_of_snow_mean,MODSCAG_start_of_snow_mean, Areas);
alloutputs(2,:) = calStat(ELM_melt_of_snow_mean,MODSCAG_melt_of_snow_mean, Areas);
alloutputs(3,:) = calStat(ELM_midpoint_of_snow_mean,MODSCAG_midpoint_of_snow_mean, Areas);
alloutputs(4,:) = calStat(ELM_end_of_snow_mean,MODSCAG_end_of_snow_mean, Areas);
alloutputs(5,:) = calStat(ELM_duration_of_snow_mean,MODSCAG_duration_of_snow_mean, Areas);

alloutputs(6,:) = calStat(ELM_start_of_snow_mean,SPIRES_start_of_snow_mean, Areas);
alloutputs(7,:) = calStat(ELM_melt_of_snow_mean,SPIRES_melt_of_snow_mean, Areas);
alloutputs(8,:) = calStat(ELM_midpoint_of_snow_mean,SPIRES_midpoint_of_snow_mean, Areas);
alloutputs(9,:) = calStat(ELM_end_of_snow_mean,SPIRES_end_of_snow_mean, Areas);
alloutputs(10,:) = calStat(ELM_duration_of_snow_mean,SPIRES_duration_of_snow_mean, Areas);


min_values = [250 -60 -60];
max_values = [365 60 60];
%% snow start
%% plot 11
figure;
set(gcf,'unit','normalized','position',[0.1,0.1,0.46,0.9]);
set(gca, 'Position', [0 0 1 1])

Areas_tmp = Areas;
Areas_tmp(isnan(ELM_start_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));

ax1 = subplot('position', [0.06 + 0.31*(1-1) 0.78 0.31 0.18]);
colormap(ax1, colors_1);
hold on
plot_global_map_other(lats, lons, ELM_start_of_snow_mean, min_values(1), max_values(1),"(a)",1, 0);
m_text(-110,33.2,num2str(nansum(ELM_start_of_snow_mean(:).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');

set(gca,'fontsize',8,'fontname','time new roman')

ylabel({'Accumulation\_onset\_date'},'fontsize',8.5, 'fontweight', 'bold')

title("ELM",'fontsize',10, 'fontweight', 'bold')

hcb = colorbar;
x1=get(gca,'position');
x=get(hcb,'Position');
%x(1)=0.19*1;
set(hcb,'Position',x)

%% plot 22
Areas_tmp = Areas;
Areas_tmp(isnan(ELM_start_of_snow_mean) | isnan(MODSCAG_start_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));


ax2 = subplot('position', [0.06 + 0.31*(2-1) 0.78 0.31 0.18]);
colormap(ax2, colors_2);
hold on
    plot_global_map_other(lats, lons, ELM_start_of_snow_mean- MODSCAG_start_of_snow_mean, min_values(2), max_values(2),"(b)",0, 0);
    m_text(-110,33.2,num2str(nansum((ELM_start_of_snow_mean(:)- MODSCAG_start_of_snow_mean(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
    
set(gca,'fontsize',8,'fontname','time new roman')

  %  ylabel({''})
title("ELM - STC-MODSCAG",'fontsize',10, 'fontweight', 'bold')

    hcb = colorbar;
    x1=get(gca,'position');
    x=get(hcb,'Position');
   %x(1)=0.19*1;

    set(hcb,'Position',x)

    %% plot 33
Areas_tmp = Areas;
Areas_tmp(isnan(ELM_start_of_snow_mean) | isnan(SPIRES_start_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));


ax2 = subplot('position', [0.06 + 0.31*(3-1) 0.78 0.31 0.18]);
colormap(ax2, colors_2);
hold on
    plot_global_map_other(lats, lons, ELM_start_of_snow_mean- SPIRES_start_of_snow_mean, min_values(2), max_values(2),"(c)",0, 0);
    m_text(-110,33.2,num2str(nansum((ELM_start_of_snow_mean(:)- SPIRES_start_of_snow_mean(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
    
set(gca,'fontsize',8,'fontname','time new roman')

    %ylabel({''})
title("ELM - SPIReS",'fontsize',8, 'fontweight', 'bold')

    hcb = colorbar;
    x1=get(gca,'position');
    x=get(hcb,'Position');
 %  x(1)=0.19*1;

    set(hcb,'Position',x)
    
%% snow melt
min_values = [40 -60 -60];
max_values = [140 60 60];

%% plot 11

Areas_tmp = Areas;
Areas_tmp(isnan(ELM_melt_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));

ax1 = subplot('position', [0.06 + 0.31*(1-1) 0.59 0.31 0.18]);
colormap(ax1, colors_1);
hold on
plot_global_map_other(lats, lons, ELM_melt_of_snow_mean, min_values(1), max_values(1),"(d)",1, 0);
m_text(-110,33.2,num2str(nansum(ELM_melt_of_snow_mean(:).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');

set(gca,'fontsize',8,'fontname','time new roman')


%title("",'fontsize',10, 'fontweight', 'bold')
ylabel({'Depletion\_onset\_date'},'fontsize',8.5, 'fontweight', 'bold')

hcb = colorbar;
x1=get(gca,'position');
x=get(hcb,'Position');
%x(1)=0.19+0.185*1;

set(hcb,'Position',x)

%% plot 22
Areas_tmp = Areas;
Areas_tmp(isnan(ELM_melt_of_snow_mean) | isnan(MODSCAG_melt_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));


ax2 = subplot('position', [0.06 + 0.31*(2-1) 0.59 0.31 0.18]);
colormap(ax2, colors_2);
hold on
    plot_global_map_other(lats, lons, ELM_melt_of_snow_mean- MODSCAG_melt_of_snow_mean, min_values(2), max_values(2),"(e)",0, 0);
    m_text(-110,33.2,num2str(nansum((ELM_melt_of_snow_mean(:)- MODSCAG_melt_of_snow_mean(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
    
set(gca,'fontsize',8,'fontname','time new roman')


    hcb = colorbar;
    x1=get(gca,'position');
    x=get(hcb,'Position');
   % x(1)=0.19+0.185*1;

    set(hcb,'Position',x)

    %% plot 33
Areas_tmp = Areas;
Areas_tmp(isnan(ELM_melt_of_snow_mean) | isnan(SPIRES_melt_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));


ax2 = subplot('position', [0.06 + 0.31*(3-1) 0.59 0.31 0.18]);
colormap(ax2, colors_2);
hold on
    plot_global_map_other(lats, lons, ELM_melt_of_snow_mean- SPIRES_melt_of_snow_mean, min_values(2), max_values(2),"(f)",0, 0);
    m_text(-110,33.2,num2str(nansum((ELM_melt_of_snow_mean(:)- SPIRES_melt_of_snow_mean(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
    
set(gca,'fontsize',8,'fontname','time new roman')


    hcb = colorbar;
    x1=get(gca,'position');
    x=get(hcb,'Position');
   %x(1)=0.19+0.185*1;

    set(hcb,'Position',x)

%% snow midday
min_values = [60 -60 -60];
max_values = [180 60 60];

%% plot 11

Areas_tmp = Areas;
Areas_tmp(isnan(ELM_midpoint_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));

ax1 = subplot('position', [0.06 + 0.31*(1-1) 0.4 0.31 0.18]);
colormap(ax1, colors_1);
hold on
plot_global_map_other(lats, lons, ELM_midpoint_of_snow_mean, min_values(1), max_values(1),"(g)",1, 0);
m_text(-110,33.2,num2str(nansum(ELM_midpoint_of_snow_mean(:).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');

set(gca,'fontsize',8,'fontname','time new roman')

ylabel({'Midpoint\_date'},'fontsize',8.5, 'fontweight', 'bold')


hcb = colorbar;
x1=get(gca,'position');
x=get(hcb,'Position');
%x(1)=0.19+0.185*2;

set(hcb,'Position',x)

%% plot 22
Areas_tmp = Areas;
Areas_tmp(isnan(ELM_midpoint_of_snow_mean) | isnan(MODSCAG_midpoint_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));


ax2 = subplot('position', [0.06 + 0.31*(2-1) 0.4 0.31 0.18]);
colormap(ax2, colors_2);
hold on
    plot_global_map_other(lats, lons, ELM_midpoint_of_snow_mean- MODSCAG_midpoint_of_snow_mean, min_values(2), max_values(2),"(h)",0, 0);
    m_text(-110,33.2,num2str(nansum((ELM_midpoint_of_snow_mean(:)- MODSCAG_midpoint_of_snow_mean(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
    
set(gca,'fontsize',8,'fontname','time new roman')


    hcb = colorbar;
    x1=get(gca,'position');
    x=get(hcb,'Position');
 %  x(1)=0.19+0.185*2;

    set(hcb,'Position',x)

    %% plot 33
Areas_tmp = Areas;
Areas_tmp(isnan(ELM_midpoint_of_snow_mean) | isnan(SPIRES_midpoint_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));


ax2 = subplot('position', [0.06 + 0.31*(3-1) 0.4 0.31 0.18]);
colormap(ax2, colors_2);
hold on
    plot_global_map_other(lats, lons, ELM_midpoint_of_snow_mean- SPIRES_midpoint_of_snow_mean, min_values(2), max_values(2),"(i)",0, 0);
    m_text(-110,33.2,num2str(nansum((ELM_midpoint_of_snow_mean(:)- SPIRES_midpoint_of_snow_mean(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
    
set(gca,'fontsize',8,'fontname','time new roman')


    hcb = colorbar;
    x1=get(gca,'position');
    x=get(hcb,'Position');
  %  x(1)=0.19+0.185*2;

    set(hcb,'Position',x)

%% snow end 

min_values = [40 -60 -60];
max_values = [240 60 60];

%% plot 11

Areas_tmp = Areas;
Areas_tmp(isnan(ELM_end_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));

ax1 = subplot('position', [0.06 + 0.31*(1-1) 0.21 0.31 0.18]);
colormap(ax1, colors_1);
hold on
plot_global_map_other(lats, lons, ELM_end_of_snow_mean, min_values(1), max_values(1),"(j)",1, 0);
m_text(-110,33.2,num2str(nansum(ELM_end_of_snow_mean(:).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');

set(gca,'fontsize',8,'fontname','time new roman')

ylabel({'End\_date'},'fontsize',8.5, 'fontweight', 'bold')


hcb = colorbar;
x1=get(gca,'position');
x=get(hcb,'Position');
% x(1)=0.19+0.185*3;
set(hcb,'Position',x)

%% plot 22
Areas_tmp = Areas;
Areas_tmp(isnan(ELM_end_of_snow_mean) | isnan(MODSCAG_end_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));


ax2 = subplot('position', [0.06 + 0.31*(2-1) 0.21 0.31 0.18]);
colormap(ax2, colors_2);
hold on
    plot_global_map_other(lats, lons, ELM_end_of_snow_mean- MODSCAG_end_of_snow_mean, min_values(2), max_values(2),"(k)",0, 0);
    m_text(-110,33.2,num2str(nansum((ELM_end_of_snow_mean(:)- MODSCAG_end_of_snow_mean(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
    
set(gca,'fontsize',8,'fontname','time new roman')


    hcb = colorbar;
    x1=get(gca,'position');
    x=get(hcb,'Position');
 %    x(1)=0.19+0.185*3;
    set(hcb,'Position',x)

    %% plot 33
Areas_tmp = Areas;
Areas_tmp(isnan(ELM_end_of_snow_mean) | isnan(SPIRES_end_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));


ax2 = subplot('position', [0.06 + 0.31*(3-1) 0.21 0.31 0.18]);
colormap(ax2, colors_2);
hold on
    plot_global_map_other(lats, lons, ELM_end_of_snow_mean- SPIRES_end_of_snow_mean, min_values(2), max_values(2),"(l)",0, 0);
    m_text(-110,33.2,num2str(nansum((ELM_end_of_snow_mean(:)- SPIRES_end_of_snow_mean(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
    
set(gca,'fontsize',8,'fontname','time new roman')


    hcb = colorbar;
    x1=get(gca,'position');
    x=get(hcb,'Position');
  %   x(1)=0.19+0.185*3;
    set(hcb,'Position',x)
%% snow duration
min_values = [60 -60 -60];
max_values = [300 60 60];

%% plot 11

Areas_tmp = Areas;
Areas_tmp(isnan(ELM_duration_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));

ax1 = subplot('position', [0.06 + 0.31*(1-1) 0.02 0.31 0.18]);
colormap(ax1, colors_1);
hold on
plot_global_map_other(lats, lons, ELM_duration_of_snow_mean, min_values(1), max_values(1),"(m)",1, 1);
m_text(-110,33.2,num2str(nansum(ELM_duration_of_snow_mean(:).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');

set(gca,'fontsize',8,'fontname','time new roman')

ylabel({'Duration'},'fontsize',8.5, 'fontweight', 'bold')


hcb = colorbar;
x1=get(gca,'position');
x=get(hcb,'Position');
 %x(1)=0.19+0.185*4;
set(hcb,'Position',x)

%% plot 22
Areas_tmp = Areas;
Areas_tmp(isnan(ELM_duration_of_snow_mean) | isnan(MODSCAG_duration_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));


ax2 = subplot('position', [0.06 + 0.31*(2-1) 0.02 0.31 0.18]);
colormap(ax2, colors_2);
hold on
    plot_global_map_other(lats, lons, ELM_duration_of_snow_mean- MODSCAG_duration_of_snow_mean, min_values(2), max_values(2),"(n)",0, 1);
    m_text(-110,33.2,num2str(nansum((ELM_duration_of_snow_mean(:)- MODSCAG_duration_of_snow_mean(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
    
set(gca,'fontsize',8,'fontname','time new roman')


    hcb = colorbar;
    x1=get(gca,'position');
    x=get(hcb,'Position');
    % x(1)=0.19+0.185*4;
    set(hcb,'Position',x)

    %% plot 33
Areas_tmp = Areas;
Areas_tmp(isnan(ELM_duration_of_snow_mean) | isnan(SPIRES_duration_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));


ax2 = subplot('position', [0.06 + 0.31*(3-1) 0.02 0.31 0.18]);
colormap(ax2, colors_2);
hold on
    plot_global_map_other(lats, lons, ELM_duration_of_snow_mean- SPIRES_duration_of_snow_mean, min_values(2), max_values(2),"(o)",0, 1);
    m_text(-110,33.2,num2str(nansum((ELM_duration_of_snow_mean(:)- SPIRES_duration_of_snow_mean(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
   
mod = ELM_duration_of_snow_mean;
obs = SPIRES_duration_of_snow_mean;

    filters = mod >0 & obs >0;

    
mod = mod(filters);
obs = obs(filters);
areas = Areas(filters);

areas = areas./nansum(areas(:));

Bias = nansum((mod - obs).*areas(:));
rBias = Bias./nansum(obs.*areas(:))*100;



set(gca,'fontsize',8,'fontname','time new roman')


    hcb = colorbar;
    x1=get(gca,'position');
    x=get(hcb,'Position');
   % x(1)=0.19+0.185*4;
    set(hcb,'Position',x)

%print(gcf, '-dtiff', '-r300', '../../figure_all_tif/snowtiming_spaital_pattern_temporal_correlation.tif')

close all

