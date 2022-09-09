
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

load('../../../../data_processing/daily_data_N_3000/ELM_daily_2004.mat');

ELM_fsnos = fsnos_daily;
load('../../../../data_processing/daily_data_N_3000/MODSCAG_daily_2004.mat');
MODSCAG_fsnos = fsnos_daily;

load('../../../../data_processing/daily_data_N_3000/SPIRES_daily_2004.mat');
SPIRES_fsnos = fsnos_daily;

Rs_1 = nan(144,168);
Rs_2 = nan(144,168);
Rs_3 = nan(144,168);

for row_i = 1:144
    for col_j = 1:168
        filter_day = [60:150];
        tmp1 = squeeze(ELM_fsnos(row_i, col_j,:));
        tmp2 = squeeze(MODSCAG_fsnos(row_i, col_j,:));
        tmp3 = squeeze(SPIRES_fsnos(row_i, col_j,:));
        
        Rs_1(row_i, col_j) = calculateR(tmp1, tmp2);
        Rs_2(row_i, col_j) = calculateR(tmp1, tmp3);
        Rs_3(row_i, col_j) = calculateR(tmp2, tmp3);
        
    end
end

min_values = [0 0 0];
max_values = [1 1 1];


%% plot
colors = flipud(brewermap(101, 'Spectral'));
colors =  (brewermap(101, 'YlGn'));

figure;
set(gcf,'unit','normalized','position',[0.1,0.1,0.8,0.55]);
set(gca, 'Position', [0 0 1 1])

%% plot 1
ax1 = subplot(1,3,1);
colormap(ax1, colors);
hold on
plot_global_map_other(lats, lons, Rs_1, min_values(1), max_values(1),'(a)',1, 0);

set(gca,'fontsize',8,'fontname','time new roman')

ax2 = subplot(1,3,2);
colormap(ax2, colors);
hold on
plot_global_map_other(lats, lons, Rs_2, min_values(1), max_values(1),'(b)',1, 0);

ax3 = subplot(1,3,3);
colormap(ax3, colors);
hold on
plot_global_map_other(lats, lons, Rs_3, min_values(1), max_values(1),'(c)',1, 0);

hcb = colorbar;
hcb.Title.String = 'Unitless';
hcb.Title.FontSize = 7;
hcb.Title.FontWeight = 'Bold';
x1=get(gca,'position');
x=get(hcb,'Position');

x(1)=0.945;
set(hcb,'Position',x)

set(gca,'fontsize',8,'fontname','time new roman')

print(gcf, '-dtiff', '-r300', '../../figure_all_tif/spaital_pattern_fsno.tif')

close all

