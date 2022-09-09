close all
clear all

load('isUS2.mat');
elevations = imread('TOP2_mean_of_Elevation.tif');

elevations(isUS2<1) = nan;

res_v = 18/2160;
res_h = 18/2160;
lon = (-125+res_h/2):res_h: (-104-res_h/2);
lat = (50-res_v/2):-res_v: (32 + res_v/2);
[lons,lats]=meshgrid(lon,lat);

res_v = 0.125;
res_h = 0.125;
lon = (-125+res_h/2):res_h: (-104-res_h/2);
lat = (50-res_v/2):-res_v: (32 + res_v/2);
[lons2,lats2]=meshgrid(lon,lat);


%% figue plot
figure;
set(gcf,'unit','normalized','position',[0.1,0.1,0.65,0.8]);
set(gca, 'Position', [0.08 0.05 0.9 0.9])
colors = flipud(brewermap(20, 'Spectral'));

a = subplot('Position',[0.05 0.45 0.45 0.5])
colormap(a, colors)

plot_global_map_elevation(lats, lons, elevations, 0, 4000, 'Elevation', 1, 1,'')
h = m_text(-123.5,42.5,'Cascade Range','color','k','fontweight','bold','fontsize',13);
set(h,'Rotation',80);

h = m_text(-119,36,'Sierra Nevada','color','k','fontweight','bold','fontsize',13);
set(h,'Rotation',120);

h = m_text(-110,45,'Rocky Mountains','color','k','fontweight','bold','fontsize',13);
set(h,'Rotation',-70);

hcb = colorbar;
hcb.Title.String = "m";

text(-0.17,0.6,'(a)','fontsize',14,'fontweight','bold')
set(gca,'fontsize',14)
%% forest cover
forest_cover = imread('forest_cover.tif');
load('isUS.mat');

forest_cover(isUS<1) = nan;

b = subplot('Position',[0.54 0.45 0.45 0.5])
colors = (brewermap(20, 'Greens'));
colormap(b, colors)
plot_global_map_forest_cover(lats2, lons2, forest_cover, 0, 80, 'Forest cover', 1, 0,'')
hcb = colorbar;
hcb.Title.String = "%";
text(-0.17,0.6,'(b)','fontsize',14,'fontweight','bold')
set(gca,'fontsize',14)

%% figure inset
elevations = imread('../WUS_mean_of_Elevation.tif');
load('../isUS.mat');

elevations = squeeze(elevations(:,:,1));

elevations_bands = nan(size(elevations));
for elevation_band = 1:7
    filters = elevations >=((elevation_band-1)*500)  & elevations<(elevation_band*500);
    if elevation_band == 7
        filters = elevations >=((elevation_band-1)*500) ;
    end
    elevations_bands(filters) = elevation_band;
end


elevations_bands(isUS<1) = nan;
t = [];
for i = 1:7
   t = [t; nansum(elevations_bands(:)==i)./nansum(elevations_bands(:)>=0)];
end

c = subplot('Position',[0.08 0.09 0.42 0.32])
box on
labels = {'0.0-0.5','0.5-1.0','1.0-1.5','1.5-2.0','2.0-2.5','2.5-3.0','\geq3.0'};
% p = pie(t)
% 
% 
% lgd = legend(labels,'location','east');
% lgd.Position(1) = lgd.Position(1)+0.11;
% legend('boxoff')
bar([1:7],t*100,'FaceColor',[0 .5 .5])
xlim([0.5 7.5])
ylabel('Area proportion (%)')

set(gca,'Fontsize',10,'xtick',1:7,'xticklabel',labels,'linewidth',1)  % THIS LINE ADDED
xlabel('Elevation (km)')
text(0.63,28,'(c)','fontsize',14,'fontweight','bold')
set(gca,'fontsize',14)

c.XAxis.FontSize = 11;

%% forest cover
forest_cover = imread('forest_cover.tif');
load('isUS.mat');

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


fc_bands(isUS<1) = nan;
tt = [];
for i = 1:5
   tt = [tt; nansum(fc_bands(:)==i)./nansum(fc_bands(:)>=0)];
end

d = subplot('Position',[0.55 0.09 0.42 0.32])
box on
labels = {'0-10','10-30','30-50','50-70','\geq70'};
% p = pie(t)
% 
% 
% lgd = legend(labels,'location','east');
% lgd.Position(1) = lgd.Position(1)+0.11;
% legend('boxoff')
bar([1:5],tt*100,'FaceColor',[0 .5 .5])
xlim([0.5 5.5])
%ylabel('Area proportion (%)')

set(gca,'Fontsize',10,'xtick',1:5,'xticklabel',labels,'linewidth',1)  % THIS LINE ADDED
xlabel('Forest cover (%)')
text(0.6,75,'(d)','fontsize',14,'fontweight','bold')
set(gca,'fontsize',14)
d.XAxis.FontSize = 11;
% figure;
print(gcf, '-dtiff', '-r300', ['../../figure_all_tif/figure_study_area.tif'])

close all
