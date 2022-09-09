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

colors =  [0.45, 0.80, 0.69;...
    0.98, 0.40, 0.35;...
    0.55, 0.60, 0.79];
%% figure
figure
set(gcf,'unit','normalized','position',[0.1,0.1,0.65,0.35]);
set(gca, 'Position', [0 0 1 1])

subplot('position', [0.06 0.09 0.26 0.8]);
hold on

group_names = {'Humans', 'Dogs' , 'God'};
condition_names = {'ELM','UA','SNODAS'};
h = daboxplot(Rs,'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',13,'mean',1, 'outliers',0);


box on
set(gca,'linewidth',1,'fontsize',10)
xlim([0.5 3.5])
ylim([0 1.05])
text(0.6,1.05-1.05*0.06,'(a)','fontsize',12,'fontweight','bold')
h.bx(1).FaceColor = colors(1,:);
h.bx(2).FaceColor = colors(2,:);
h.bx(3).FaceColor = colors(3,:);

set(gca,'linewidth',1,'fontsize',10)
ylabel('R')
%% plot 2
subplot('position', [0.39 0.09 0.26 0.8]);
hold on

group_names = {'Humans', 'Dogs' , 'God'};
condition_names = {'ELM','UA','SNODAS'};
h = daboxplot(Biass,'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',13,'mean',1, 'outliers',0);


box on
set(gca,'linewidth',1,'fontsize',10)
xlim([0.5 3.5])
ylim([-400 200])
text(0.6,200-600*0.06,'(b)','fontsize',12,'fontweight','bold')
h.bx(1).FaceColor = colors(1,:);
h.bx(2).FaceColor = colors(2,:);
h.bx(3).FaceColor = colors(3,:);

set(gca,'linewidth',1,'fontsize',10)
ylabel('Bias (mm)')

%% 3
subplot('position', [0.72 0.09 0.26 0.8]);
hold on

group_names = {'Humans', 'Dogs' , 'God'};
condition_names = {'ELM','UA','SNODAS'};
h = daboxplot(RMSEs,'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,...
    'whiskers',1,'scatter',0,'jitter',0,'scattersize',13,'mean',1, 'outliers',0);


box on
set(gca,'linewidth',1,'fontsize',10)
xlim([0.5 3.5])
ylim([0 400])
text(0.6,400-400*0.06,'(c)','fontsize',12,'fontweight','bold')
h.bx(1).FaceColor = colors(1,:);
h.bx(2).FaceColor = colors(2,:);
h.bx(3).FaceColor = colors(3,:);

set(gca,'linewidth',1,'fontsize',10)

ylabel('RMSE (mm)')


print(gcf, '-dtiff', '-r300', '../../figure_all_tif/comparison_with_snotel.tif')

close all