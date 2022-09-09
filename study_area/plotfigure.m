

close all
clear all

load('isUS.mat');

load('monthly_data.mat');
res_v = 0.125;
res_h = 0.125;
lon = (-125+res_h/2):res_h: (-104-res_h/2);
lat = (50-res_v/2):-res_v: (32 + res_v/2);
[lons,lats]=meshgrid(lon,lat);


%% figue plot
figure;
ha = tight_subplot(4,3,[.01 .03],[.1 .01],[.01 .01]);
for ii = 1:12;
    data_tmp = squeeze(ELM_swes(:,:,12*9+ii));

    axes(ha(ii));
    colors = flipud(brewermap(1000, 'Spectral'));
    colormap(colors)
    
    plot_global_map_other(lats, lons, data_tmp, 0, 1000, '', 1, 1,'')
    
    hcb = colorbar;
    hcb.Title.String = "m"; 
end

figure;
ha = tight_subplot(4,3,[.01 .03],[.1 .01],[.01 .01]);
for ii = 1:12;
    data_tmp = squeeze(UA_swes(:,:,12*9+ii));

    axes(ha(ii));
    colors = flipud(brewermap(1000, 'Spectral'));
    colormap(colors)
    
    plot_global_map_other(lats, lons, data_tmp, 0, 1000, '', 1, 1,'')
    
    hcb = colorbar;
    hcb.Title.String = "m"; 
end


figure;
ha = tight_subplot(4,3,[.01 .03],[.1 .01],[.01 .01]);
for ii = 1:12;
    data_tmp = squeeze(SNODAS_swes(:,:,12*9+ii));

    axes(ha(ii));
    colors = flipud(brewermap(1000, 'Spectral'));
    colormap(colors)
    
    plot_global_map_other(lats, lons, data_tmp, 0, 1000, '', 1, 1,'')
    
    hcb = colorbar;
    hcb.Title.String = "m"; 
end

% figure;
%print(gcf, '-dtiff', '-r300', ['figure_study_area.tif'])

%close all



