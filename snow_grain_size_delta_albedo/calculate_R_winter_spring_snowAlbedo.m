clc
clear all
close all


load('../../../../data_processing/snowalbedo_daily_3000_v2.mat');
load('../../../../data_processing/fsno_daily_3000_v2.mat');
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
            tmp1 = squeeze(ELM_snowAlbedos(row_i, col_j,filter_day));
            tmp2 = squeeze(MODSCAG_snowAlbedos(row_i, col_j,filter_day));
            tmp3 = squeeze(SPIRES_snowAlbedos(row_i, col_j,filter_day));
            
            tmp11 = squeeze(ELM_fsnos(row_i, col_j,filter_day));
            tmp22 = squeeze(MODSCAG_fsnos(row_i, col_j,filter_day));
            tmp33 = squeeze(SPIRES_fsnos(row_i, col_j,filter_day));

            Rs_1(row_i, col_j, season_i) = calculateR2(tmp1, tmp2,tmp11,tmp22);
            Rs_2(row_i, col_j, season_i) = calculateR2(tmp1, tmp3,tmp11,tmp33);
            Rs_3(row_i, col_j, season_i) = calculateR2(tmp2, tmp3,tmp22,tmp33);
            
        end
    end
end

save('snowalbedo_temporal_R_v2.mat','Rs_1','Rs_2','Rs_3');