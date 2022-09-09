clc
clear all
close all


load('../../../../data_processing/snowdepth_daily_3000.mat');
Rs_1 = nan(144,168,2);
Rs_2 = nan(144,168,2);
Rs_3 = nan(144,168,2);
for season_i = 1:2
    for row_i = 1:144
        for col_j = 1:168
            if season_i == 1
                filter_day = [1:59 335:366];
            else
                filter_day =  [60:150];
            end
            
            tmp = repmat(filter_day,1,19);
            for tmp_k = 1:19
                tmp((tmp_k -1)*91+(1:91)) = (tmp_k -1)*366 + filter_day;
            end
            filter_day = tmp;
            
            tmp1 = squeeze(ELM_snowdepths(row_i, col_j,filter_day))*1000;
            tmp2 = squeeze(UA_snowdepths(row_i, col_j,filter_day));
            tmp3 = squeeze(SNODAS_snowdepths(row_i, col_j,filter_day));
            
            Rs_1(row_i, col_j, season_i) = calculateR(tmp1, tmp2);
            Rs_2(row_i, col_j, season_i) = calculateR(tmp1, tmp3);
            Rs_3(row_i, col_j, season_i) = calculateR(tmp2, tmp3);
            
        end
    end
end

save('sd_temporal_R.mat','Rs_1','Rs_2','Rs_3');