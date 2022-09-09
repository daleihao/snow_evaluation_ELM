function plot_global_map(lats, lons, min_clr, max_clr, seq_text, isyticklabel, isxticklabel, Stats)
axis equal
m_proj('miller','lat',[32 49.2],'lon',[-125 -104]); % robinson Mollweide
%m_coast('color','k','linewidth',1);
hold on


    yticks = [35:5:45];
    ytick_labels = ['35';'40';'45'];
    xticks = [-125:5:-110];
    xtick_labels = ['-125';'-120';'-115';'-110'];

if isyticklabel && isxticklabel
    m_grid('tickdir','out','linestyle','none','backcolor',[.9 .99 1], ...
        'fontsize',12,'tickstyle','dd',...
        'ytick',yticks,'yticklabel',ytick_labels,...
        'xtick',xticks,'xticklabel',xtick_labels...
        );
elseif isyticklabel && ~isxticklabel
    m_grid('tickdir','out','linestyle','none','backcolor',[.9 .99 1], 'xticklabels',[], ...
        'fontsize',12,'tickstyle','dd',...
        'ytick',yticks,'yticklabel',ytick_labels...
        );
elseif ~isyticklabel && isxticklabel
    m_grid('tickdir','out','linestyle','none','backcolor',[.9 .99 1], 'yticklabels',[], ...
        'fontsize',12,'tickstyle','dd',...
        'xtick',xticks,'xticklabel',xtick_labels,'yticklabel',[]...
        );
else
    m_grid('tickdir','out','linestyle','none','backcolor',[.9 .99 1], 'xticklabels',[], 'yticklabels',[], ...
        'fontsize',12,'tickstyle','dd','yticklabel',[]...
        );
end

%%
% load('elevations_TP.mat');
% sw_total(elevations<=1500 | isnan(elevations)) = nan;

%m_pcolor(lons,lats,sw_total);
% M=m_shaperead('../../data/TP_shp/ROTW_China');
% for k=1:length(M.ncst)
%     m_line(M.ncst{k}(:,1),M.ncst{k}(:,2),'color','k','linewidth',1);
% end

% M = shaperead('../tl_2017_us_state/tl_2017_us_state');
% for k=1:length(M)
%     m_line(M(k).X,M(k).Y,'color','k','linewidth',1);
% end
load('site_loc.mat');
% siteData = readmatrix('SNOTEL_summary.csv');
% filters = siteData(:,7)<-125 | siteData(:,7)>-104 | siteData(:,6)>50 | siteData(:,6)<32;
%
m_scatter(siteInfo(:,5),siteInfo(:,4),10,Stats,'filled')

res_v = 0.125;
res_h = 0.125;
load('boundary_WUS.mat');

for k=1:length(b)
    boundary = b{k};
    longs_k = boundary(:,2);
    lats_k = boundary(:,1);
    longs_k = (-125+res_h/2) + (longs_k-1)*res_h;
    lats_k = (50-res_v/2) - (lats_k-1)*res_v;
    m_line(longs_k,lats_k,'color','k','linewidth',1);
end

%% plot contour

% res_v = 0.125;
% res_h = 0.125;
% lon = (-180+res_h/2):res_h: (180-res_h/2);
% lat = (90-res_v/2):-res_v: (-90 + res_v/2);
% [lons2,lats2]=meshgrid(lon,lat);
%
% load("mean_DEM.mat");
% m_contour(lons2,lats2, mean_DEM_0125, [1500 4500], 'edgecolor',[0 0 0],'facecolor','none', 'linewidth', 1);
%
% load('IndusCHARISandHyrdoSheds.mat');
% m_plot( upperIndus.Vertices(:,1),upperIndus.Vertices(:,2),'color','r', 'linewidth', 1);
%

shading flat;
caxis([min_clr-1e-5,max_clr+1e-5])

%m_text(-124.5,33.5, seq_text,'fontsize',10);
%colormap(m_colmap('jet','step',10));
%m_text(66,24,num_label,'fontsize',12,'fontweight','bold');
%colorbar
%colormap(colors_single);

if seq_text ~= ""
   t = title(seq_text,'fontsize',12, 'fontweight', 'bold');
    set(t, 'horizontalAlignment', 'left');
    set(t, 'units', 'normalized');
    h1 = get(t, 'position');
    set(t, 'position', [0 h1(2) h1(3)]);
end

set(gca, 'FontName', 'Time New Roman');

view(0,90);
hold off