clc
clear all
close all


load('../../../../data_processing/albedo_daily_3000.mat');
Rs_1 = nan(144,168,2);
Rs_2 = nan(144,168,2);
Rs_3 = nan(144,168,2);
for season_i = 1:2
    for row_i = 1:144
        for col_j = 1:168
            if season_i == 1
                filter_day = [1:59 335:366];
            else
                filter_day = [60:150];
                
                tmp = repmat(filter_day,1,19);
                for tmp_k = 1:19
                    tmp((tmp_k -1)*91+(1:91)) = (tmp_k -1)*366 + filter_day;
                end
                filter_day = tmp;
                
            end
            
    
            tmp1 = squeeze(ELM_Albedos(row_i, col_j,filter_day));
            MCD43_surfacealbedo = MCD43_BSAs(row_i, col_j,filter_day).* (1 - ELM_Skyls(row_i, col_j,filter_day)) + MCD43_WSAs(row_i, col_j,filter_day).* ELM_Skyls(row_i, col_j,filter_day);
            tmp2 = squeeze(MCD43_surfacealbedo)/1000;
            
            Rs_1(row_i, col_j, season_i) = calculateR(tmp1, tmp2);
            
        end
    end
end

save('Albedo_temporal_R','Rs_1');