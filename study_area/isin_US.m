res_v = 0.125;
res_h = 0.125;
lon = (-125+res_h/2):res_h: (-104-res_h/2);
lat = (50-res_v/2):-res_v: (32 + res_v/2);
[lons,lats]=meshgrid(lon,lat);


M = shaperead('../tl_2017_us_state/tl_2017_us_state');
isins = zeros(size(lons(:),1), length(M));
for k=1:length(M)
isins(:,k) = inpolygon(lons(:),lats(:),M(k).X,M(k).Y);
end


isUS = reshape(sum(isins,2)>0, [144 168]);

%% get 1
res_v = 18/2160;
res_h = 18/2160;
lon = (-125+res_h/2):res_h: (-104-res_h/2);
lat = (50-res_v/2):-res_v: (32 + res_v/2);
[lons,lats]=meshgrid(lon,lat);

isUS2 = ones(size(lons));
for row = 1:size(lons,1)
    for col = 1:size(lons,2)
        lat_i = lats(row, col);
        lon_i = lons(row, col);
        row_1 = round(((50 - 0.125/2) - lat_i)/0.125)+1;
        col_1 = round((lon_i - (-125+0.125/2))/0.125)+1;
        isUS2(row, col) =isUS(row_1, col_1);
    end
end
