elevations = imread('../WUS_mean_of_Elevation.tif');
elevations = squeeze(elevations(:,:,1));
elevation_stds = imread('../WUS_std_of_Elevation.tif');

diff = ELM_swes_tmp - UA_swes_tmp;