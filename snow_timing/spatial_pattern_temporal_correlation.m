
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
colors_1 = flipud(brewermap(101, 'Spectral'));
colors_2 =  flipud(brewermap(101, 'Spectral'));




load('../../../../data_processing/snow_timing/snow_timing_alldata.mat');
load('../../../../data_processing/snow_timing/snow_timing_num_alldata.mat');

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



min_values = [250 -40 -40];
max_values = [365 40 40];
%% snow start
%% plot 11
figure;
set(gcf,'unit','normalized','position',[0.1,0.1,0.74,0.58]);
set(gca, 'Position', [0 0 1 1])

Areas_tmp = Areas;
Areas_tmp(isnan(ELM_start_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));

ax1 = subplot('position', [0.03 + 0.185*(1-1) 0.65 0.17 0.28]);
colormap(ax1, colors_1);
hold on
plot_global_map_other(lats, lons, ELM_start_of_snow_mean, min_values(1), max_values(1),"(a)",1, 0);
m_text(-111,33.2,num2str(nansum(ELM_start_of_snow_mean(:).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');

set(gca,'fontsize',8,'fontname','time new roman')

ylabel({'ELM'})
title("Winter",'fontsize',10, 'fontweight', 'bold')

title("SOD",'fontsize',10, 'fontweight', 'bold')

hcb = colorbar;
x1=get(gca,'position');
x=get(hcb,'Position');
x(1)=0.19*1;
set(hcb,'Position',x)

%% plot 22
Areas_tmp = Areas;
Areas_tmp(isnan(ELM_start_of_snow_mean) | isnan(MODSCAG_start_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));


ax2 = subplot('position', [0.03 + 0.185*(1-1) 0.36 0.17 0.28]);
colormap(ax2, colors_2);
hold on
    plot_global_map_other(lats, lons, ELM_start_of_snow_mean- MODSCAG_start_of_snow_mean, min_values(2), max_values(2),"(f)",1, 0);
    m_text(-111,33.2,num2str(nansum((ELM_start_of_snow_mean(:)- MODSCAG_start_of_snow_mean(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
    
set(gca,'fontsize',8,'fontname','time new roman')

    ylabel({'ELM - MODSCAG'})

    hcb = colorbar;
    x1=get(gca,'position');
    x=get(hcb,'Position');
   x(1)=0.19*1;

    set(hcb,'Position',x)

    %% plot 33
Areas_tmp = Areas;
Areas_tmp(isnan(ELM_start_of_snow_mean) | isnan(SPIRES_start_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));


ax2 = subplot('position', [0.03 + 0.185*(1-1) 0.07 0.17 0.28]);
colormap(ax2, colors_2);
hold on
    plot_global_map_other(lats, lons, ELM_start_of_snow_mean- SPIRES_start_of_snow_mean, min_values(2), max_values(2),"(k)",1, 1);
    m_text(-111,33.2,num2str(nansum((ELM_start_of_snow_mean(:)- SPIRES_start_of_snow_mean(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
    
set(gca,'fontsize',8,'fontname','time new roman')

    ylabel({'ELM - SPIReS'})

    hcb = colorbar;
    x1=get(gca,'position');
    x=get(hcb,'Position');
   x(1)=0.19*1;

    set(hcb,'Position',x)
    
%% snow melt
min_values = [40 -40 -40];
max_values = [140 40 40];

%% plot 11

Areas_tmp = Areas;
Areas_tmp(isnan(ELM_melt_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));

ax1 = subplot('position', [0.03 + 0.185*(2-1) 0.65 0.17 0.28]);
colormap(ax1, colors_1);
hold on
plot_global_map_other(lats, lons, ELM_melt_of_snow_mean, min_values(1), max_values(1),"(b)",0, 0);
m_text(-111,33.2,num2str(nansum(ELM_melt_of_snow_mean(:).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');

set(gca,'fontsize',8,'fontname','time new roman')


title("SMOD",'fontsize',10, 'fontweight', 'bold')

hcb = colorbar;
x1=get(gca,'position');
x=get(hcb,'Position');
x(1)=0.19+0.185*1;

set(hcb,'Position',x)

%% plot 22
Areas_tmp = Areas;
Areas_tmp(isnan(ELM_melt_of_snow_mean) | isnan(MODSCAG_melt_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));


ax2 = subplot('position', [0.03 + 0.185*(2-1) 0.36 0.17 0.28]);
colormap(ax2, colors_2);
hold on
    plot_global_map_other(lats, lons, ELM_melt_of_snow_mean- MODSCAG_melt_of_snow_mean, min_values(2), max_values(2),"(g)",0, 0);
    m_text(-111,33.2,num2str(nansum((ELM_melt_of_snow_mean(:)- MODSCAG_melt_of_snow_mean(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
    
set(gca,'fontsize',8,'fontname','time new roman')


    hcb = colorbar;
    x1=get(gca,'position');
    x=get(hcb,'Position');
    x(1)=0.19+0.185*1;

    set(hcb,'Position',x)

    %% plot 33
Areas_tmp = Areas;
Areas_tmp(isnan(ELM_melt_of_snow_mean) | isnan(SPIRES_melt_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));


ax2 = subplot('position', [0.03 + 0.185*(2-1) 0.07 0.17 0.28]);
colormap(ax2, colors_2);
hold on
    plot_global_map_other(lats, lons, ELM_melt_of_snow_mean- SPIRES_melt_of_snow_mean, min_values(2), max_values(2),"(l)",0, 1);
    m_text(-111,33.2,num2str(nansum((ELM_melt_of_snow_mean(:)- SPIRES_melt_of_snow_mean(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
    
set(gca,'fontsize',8,'fontname','time new roman')


    hcb = colorbar;
    x1=get(gca,'position');
    x=get(hcb,'Position');
   x(1)=0.19+0.185*1;

    set(hcb,'Position',x)

%% snow midday
min_values = [60 -40 -40];
max_values = [180 40 40];

%% plot 11

Areas_tmp = Areas;
Areas_tmp(isnan(ELM_midpoint_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));

ax1 = subplot('position', [0.03 + 0.185*(3-1) 0.65 0.17 0.28]);
colormap(ax1, colors_1);
hold on
plot_global_map_other(lats, lons, ELM_midpoint_of_snow_mean, min_values(1), max_values(1),"(c)",0, 0);
m_text(-111,33.2,num2str(nansum(ELM_midpoint_of_snow_mean(:).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');

set(gca,'fontsize',8,'fontname','time new roman')

title("Winter",'fontsize',10, 'fontweight', 'bold')

title("SMMD",'fontsize',10, 'fontweight', 'bold')

hcb = colorbar;
x1=get(gca,'position');
x=get(hcb,'Position');
x(1)=0.19+0.185*2;

set(hcb,'Position',x)

%% plot 22
Areas_tmp = Areas;
Areas_tmp(isnan(ELM_midpoint_of_snow_mean) | isnan(MODSCAG_midpoint_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));


ax2 = subplot('position', [0.03 + 0.185*(3-1) 0.36 0.17 0.28]);
colormap(ax2, colors_2);
hold on
    plot_global_map_other(lats, lons, ELM_midpoint_of_snow_mean- MODSCAG_midpoint_of_snow_mean, min_values(2), max_values(2),"(h)",0, 0);
    m_text(-111,33.2,num2str(nansum((ELM_midpoint_of_snow_mean(:)- MODSCAG_midpoint_of_snow_mean(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
    
set(gca,'fontsize',8,'fontname','time new roman')


    hcb = colorbar;
    x1=get(gca,'position');
    x=get(hcb,'Position');
   x(1)=0.19+0.185*2;

    set(hcb,'Position',x)

    %% plot 33
Areas_tmp = Areas;
Areas_tmp(isnan(ELM_midpoint_of_snow_mean) | isnan(SPIRES_midpoint_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));


ax2 = subplot('position', [0.03 + 0.185*(3-1) 0.07 0.17 0.28]);
colormap(ax2, colors_2);
hold on
    plot_global_map_other(lats, lons, ELM_midpoint_of_snow_mean- SPIRES_midpoint_of_snow_mean, min_values(2), max_values(2),"(m)",0, 1);
    m_text(-111,33.2,num2str(nansum((ELM_midpoint_of_snow_mean(:)- SPIRES_midpoint_of_snow_mean(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
    
set(gca,'fontsize',8,'fontname','time new roman')


    hcb = colorbar;
    x1=get(gca,'position');
    x=get(hcb,'Position');
    x(1)=0.19+0.185*2;

    set(hcb,'Position',x)

%% snow end 

min_values = [40 -40 -40];
max_values = [240 40 40];

%% plot 11

Areas_tmp = Areas;
Areas_tmp(isnan(ELM_end_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));

ax1 = subplot('position', [0.03 + 0.185*(4-1) 0.65 0.17 0.28]);
colormap(ax1, colors_1);
hold on
plot_global_map_other(lats, lons, ELM_end_of_snow_mean, min_values(1), max_values(1),"(d)",0, 0);
m_text(-111,33.2,num2str(nansum(ELM_end_of_snow_mean(:).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');

set(gca,'fontsize',8,'fontname','time new roman')


title("SED",'fontsize',10, 'fontweight', 'bold')

hcb = colorbar;
x1=get(gca,'position');
x=get(hcb,'Position');
 x(1)=0.19+0.185*3;
set(hcb,'Position',x)

%% plot 22
Areas_tmp = Areas;
Areas_tmp(isnan(ELM_end_of_snow_mean) | isnan(MODSCAG_end_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));


ax2 = subplot('position', [0.03 + 0.185*(4-1) 0.36 0.17 0.28]);
colormap(ax2, colors_2);
hold on
    plot_global_map_other(lats, lons, ELM_end_of_snow_mean- MODSCAG_end_of_snow_mean, min_values(2), max_values(2),"(i)",0, 0);
    m_text(-111,33.2,num2str(nansum((ELM_end_of_snow_mean(:)- MODSCAG_end_of_snow_mean(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
    
set(gca,'fontsize',8,'fontname','time new roman')


    hcb = colorbar;
    x1=get(gca,'position');
    x=get(hcb,'Position');
     x(1)=0.19+0.185*3;
    set(hcb,'Position',x)

    %% plot 33
Areas_tmp = Areas;
Areas_tmp(isnan(ELM_end_of_snow_mean) | isnan(SPIRES_end_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));


ax2 = subplot('position', [0.03 + 0.185*(4-1) 0.07 0.17 0.28]);
colormap(ax2, colors_2);
hold on
    plot_global_map_other(lats, lons, ELM_end_of_snow_mean- SPIRES_end_of_snow_mean, min_values(2), max_values(2),"(n)",0, 1);
    m_text(-111,33.2,num2str(nansum((ELM_end_of_snow_mean(:)- SPIRES_end_of_snow_mean(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
    
set(gca,'fontsize',8,'fontname','time new roman')


    hcb = colorbar;
    x1=get(gca,'position');
    x=get(hcb,'Position');
     x(1)=0.19+0.185*3;
    set(hcb,'Position',x)
%% snow duration
min_values = [60 -40 -40];
max_values = [300 40 40];

%% plot 11

Areas_tmp = Areas;
Areas_tmp(isnan(ELM_duration_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));

ax1 = subplot('position', [0.03 + 0.185*(5-1) 0.65 0.17 0.28]);
colormap(ax1, colors_1);
hold on
plot_global_map_other(lats, lons, ELM_duration_of_snow_mean, min_values(1), max_values(1),"(e)",0, 0);
m_text(-111,33.2,num2str(nansum(ELM_duration_of_snow_mean(:).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');

set(gca,'fontsize',8,'fontname','time new roman')


title("SDD",'fontsize',10, 'fontweight', 'bold')

hcb = colorbar;
x1=get(gca,'position');
x=get(hcb,'Position');
 x(1)=0.19+0.185*4;
set(hcb,'Position',x)

%% plot 22
Areas_tmp = Areas;
Areas_tmp(isnan(ELM_duration_of_snow_mean) | isnan(MODSCAG_duration_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));


ax2 = subplot('position', [0.03 + 0.185*(5-1) 0.36 0.17 0.28]);
colormap(ax2, colors_2);
hold on
    plot_global_map_other(lats, lons, ELM_duration_of_snow_mean- MODSCAG_duration_of_snow_mean, min_values(2), max_values(2),"(j)",0, 0);
    m_text(-111,33.2,num2str(nansum((ELM_duration_of_snow_mean(:)- MODSCAG_duration_of_snow_mean(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
    
set(gca,'fontsize',8,'fontname','time new roman')


    hcb = colorbar;
    x1=get(gca,'position');
    x=get(hcb,'Position');
     x(1)=0.19+0.185*4;
    set(hcb,'Position',x)

    %% plot 33
Areas_tmp = Areas;
Areas_tmp(isnan(ELM_duration_of_snow_mean) | isnan(SPIRES_end_of_snow_mean)) = nan;
Areas_tmp = Areas_tmp./nansum(Areas_tmp(:));


ax2 = subplot('position', [0.03 + 0.185*(5-1) 0.07 0.17 0.28]);
colormap(ax2, colors_2);
hold on
    plot_global_map_other(lats, lons, ELM_duration_of_snow_mean- SPIRES_duration_of_snow_mean, min_values(2), max_values(2),"(o)",0, 1);
    m_text(-111,33.2,num2str(nansum((ELM_duration_of_snow_mean(:)- SPIRES_duration_of_snow_mean(:)).*Areas_tmp(:)),'%4.1f'),'fontsize',10, 'fontweight', 'bold');
    
set(gca,'fontsize',8,'fontname','time new roman')


    hcb = colorbar;
    x1=get(gca,'position');
    x=get(hcb,'Position');
    x(1)=0.19+0.185*4;
    set(hcb,'Position',x)

print(gcf, '-dtiff', '-r300', '../../figure_all_tif/snowtiming_spaital_pattern_temporal_correlation.tif')

close all

