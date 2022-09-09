clc
clear all
close all
load('ELM_UA_SNODAS_Data_v2.mat');
load('SNOTELData_v2.mat');

Rs = nan(829, 3);
RMSEs = nan(829, 3);
Biass = nan(829, 3);

sitenum = 829;

for site_i = 1:sitenum
    
    %% filter
    %  ELMData_monthl:,y<2000 & SNODASData_monthly<2000 & UAData_monthly <2000 & SNOTELData_monthly <2000;
    
    ELMData_filter = ELMData(:, site_i+3);
    SNODASData_filter = SNODASData(:, site_i+3);
    UAData_filter =  UAData(:, site_i+3);
    SNOTELData_filter = SNOTELData(:, site_i+3);
    
    filters = SNOTELData_filter >0 ;%& ...
    
    if(sum(filters)<10)
        continue;
    end
    ELMData_filter = ELMData_filter(filters);
    SNODASData_filter = SNODASData_filter(filters);
    UAData_filter =  UAData_filter(filters);
    SNOTELData_filter = SNOTELData_filter(filters);
    
    
    [R,RMSE,Bias] = calculateR2(SNOTELData_filter, ELMData_filter);%& ...
    Rs(site_i, 1) =R;
    RMSEs(site_i, 1) = RMSE;
    Biass(site_i, 1) = Bias;
    [R,RMSE,Bias] = calculateR2(SNOTELData_filter, UAData_filter);%& ...
    Rs(site_i, 2) =R;
    RMSEs(site_i, 2) = RMSE;
    Biass(site_i, 2) = Bias;
    
    [R,RMSE,Bias] = calculateR2(SNOTELData_filter, SNODASData_filter);%& ...
    Rs(site_i, 3) =R;
    RMSEs(site_i, 3) = RMSE;
    Biass(site_i, 3) = Bias;
    
end

%% lat lon
load('site_loc.mat');

load('../isUS.mat');


res_v = 0.125;
res_h = 0.125;
lon = (-125+res_h/2):res_h: (-104-res_h/2);
lat = (50-res_v/2):-res_v: (32 + res_v/2);
[lons,lats]=meshgrid(lon,lat);


indexs  = repmat(0:18,3,1);

GridNames = {'Winter', 'Spring', 'Summer', 'Autumn'};


%% plot
colors = flipud(brewermap(10, 'Spectral'));
colors_2 =  flipud(brewermap(11, 'RdBu'));

colors_2(6,:) = [0.85 0.85 0.85];

figure;
set(gcf,'unit','normalized','position',[0.1,0.1,0.8,0.42]);
set(gca, 'Position', [0 0 1 1])


%% plot 1
ax1 = subplot('position', [0.03  0.08 0.26 0.83]);
colormap(ax1, colors);
hold on
plot_global_map(lats, lons,0,1,'(a) R',1, 1,Rs(:,1));
set(gca,'fontsize',14,'fontname','time new roman')

hcb = colorbar;
hcb.Title.String = 'Unitless';
hcb.Title.FontSize = 10;
hcb.Title.FontWeight = 'Bold';
x1=get(gca,'position');
x=get(hcb,'Position');
x(3)=0.01;
x(4)=0.78;
x(1)=0.30;
x(2)=0.10;
set(hcb,'Position',x)


ax2 = subplot('position', [0.35  0.08 0.26 0.83]);
colormap(ax2, colors);
plot_global_map(lats, lons, 0,400,'(b) RMSE',0, 1,RMSEs(:,1));
set(gca,'fontsize',14,'fontname','time new roman')

hcb = colorbar;
hcb.Title.String = 'mm';
hcb.Title.FontSize = 10;
hcb.Title.FontWeight = 'Bold';
x1=get(gca,'position');
x=get(hcb,'Position');
x(3)=0.01;
x(4)=0.78;
x(1)=0.62;
x(2)=0.10;
set(hcb,'Position',x)


ax3 = subplot('position', [0.67  0.08 0.26 0.83]);
colormap(ax3, colors_2);
hold on
plot_global_map(lats, lons, -200,200,"(c) Bias",0, 1,Biass(:,1));
hcb = colorbar;
hcb.Title.String = 'mm';
hcb.Title.FontSize = 10;
hcb.Title.FontWeight = 'Bold';
x1=get(gca,'position');
x=get(hcb,'Position');
x(3)=0.01;
x(4)=0.78;
x(1)=0.94;
x(2)=0.10;
set(hcb,'Position',x)
set(gca,'fontsize',14,'fontname','time new roman')


print(gcf, '-dtiff', '-r300', '../../figure_all_tif/SNOTEL_SWE_figure.tif')

close all

