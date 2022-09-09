
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



%% plot
colors = flipud(brewermap(101, 'Spectral'));
colors_2 =  flipud(brewermap(101, 'RdBu'));

figure;
set(gcf,'unit','normalized','position',[0.1,0.1,0.45,0.25]);
set(gca, 'Position', [0 0 1 1])

load('ELM_2004.mat');
midpoint_of_snow(isUS<1) = nan;




min_values = [30 30 30];
max_values = [180 180 180];

%% plot 1
ax1 = subplot('position', [0.02 + 0.3*(1-1) 0.07 0.29 0.8]);
colormap(ax1, colors);
hold on
plot_global_map_other(lats, lons, midpoint_of_snow, min_values(1), max_values(1),"ELM",1, 1);


set(gca,'fontsize',8,'fontname','time new roman')

load('MODSCAG_2004.mat');
midpoint_of_snow(isUS<1) = nan;

ax2 = subplot('position', [0.02 + 0.3*(2-1) 0.07 0.29 0.8]);
colormap(ax2, colors);
hold on
plot_global_map_other(lats, lons, midpoint_of_snow, min_values(2), max_values(2),"MODSCAG",0, 1);
set(gca,'fontsize',8,'fontname','time new roman')

load('SPIRES_2004.mat');
midpoint_of_snow(isUS<1) = nan;

ax3 = subplot('position', [0.02 + 0.3*(3-1) 0.07 0.29 0.8]);
colormap(ax3, colors);
hold on
plot_global_map_other(lats, lons, midpoint_of_snow, min_values(3), max_values(3),"SPIReS",0, 1);
set(gca,'fontsize',8,'fontname','time new roman')

hcb = colorbar;
hcb.Title.String = 'DOY';
hcb.Title.FontSize = 7;
hcb.Title.FontWeight = 'Bold';
x1=get(gca,'position');
x=get(hcb,'Position');
% x(3)=0.012;
x(4)=0.8;
x(1)=0.93;
x(2)=0.07;
set(hcb,'Position',x)
set(gca,'fontsize',8,'fontname','time new roman')


print(gcf, '-dtiff', '-r300', 'different_snowmidpoint_midpoint_date.tif')

close all

