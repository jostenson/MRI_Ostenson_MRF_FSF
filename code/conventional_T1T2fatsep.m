%% T1, T2, and fat-water separation using conventional methods

% clear, close all, clc;

addpath('./utils/')

BASEPATH = '../contrib/hernando/';
addpath([BASEPATH 'common/']);
addpath([BASEPATH 'graphcut/']);
addpath([BASEPATH 'descent/']);
addpath([BASEPATH 'mixed_fitting/']);
addpath([BASEPATH 'create_synthetic/']);
addpath([BASEPATH 'matlab_bgl/']);

%% parameters

disp_flag = 1;

% fat-water sep parameters
imDataParams.FieldStrength = 3.0;
imDataParams.PrecessionIsClockwise = 1;

ppm_ref = 4.65; % H2O
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = (ppm_ref - [0.9 1.3 1.3 2.1 2.1 2.75 4.2 5.3]);
algoParams.species(2).relAmps = [0.144 0.5*1.0 0.5*1.0 0.5*0.241 0.5*0.241 0.033 0.064 0.122];

% Algorithm-specific parameters
algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
algoParams.range_r2star = [0 100]; % Range of R2* values
algoParams.NUM_R2STARS = 11; % Numbre of R2* values for quantization
algoParams.range_fm = [-250 250]; % Range of field map values
algoParams.NUM_FMS = 301; % Number of field map values to discretize
algoParams.NUM_ITERS = 40; % Number of graph cut iterations
algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
algoParams.lambda = 0.05; % Regularization parameter
algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
algoParams.TRY_PERIODIC_RESIDUAL = 0;
algoParams.THRESHOLD = 0.01;

%% load data

load([dir_in fn_conv]);

n_row = size( spgr_data.img_stack_me_spgr, 1 );
n_col = n_row;

%% load/create ROI map

n_sl = size( spgr_data.img_stack_me_spgr, 4 );

img_conventional = zeros( n_row, n_col, n_sl );
for ii = 1:n_sl
    img_conventional(:,:,ii) = max( abs( squeeze( spgr_data.img_stack_me_spgr(:,:,:,ii) ) ), [],  3 );
end

img_mask = img_conventional(:,:,1);
img_mask( img_conventional(:,:,1) < multithresh(img_conventional(:,:,1),1) ) = 0;
img_mask(img_mask~=0) = 1;

if create_ROIs == 1
    
    ROI_stack = zeros( n_row, n_col, n_ROI );
    
    figure(1); clf;
    imagesc( img_conventional(:,:,1) ); axis image; colormap(gray);
    
    for ii = 1:n_ROI
        
        title( sprintf('draw ROI %d',ii) );
        my_roi = roipoly();
        ROI_stack(:,:,ii) = my_roi.* img_mask;
        
    end
    
    save( [dir_out fn_ROI], 'ROI_stack');
    
else
    
    load( [dir_out fn_ROI] )
    
end

img_mask = sum( ROI_stack, 3 );
img_mask( img_mask > 1 ) = 1;

figure(1); clf;
for ii = 1:n_sl
    imagesc( img_conventional(:,:,ii) ); axis image; colormap(gray);
    green = cat( 3, zeros(n_row,n_col), ones(n_row,n_col), zeros(n_row,n_col) );
    hold on;
    h = imshow(green);
    alpha_map = img_mask;
    set(h, 'AlphaData', alpha_map);
    hold off
    title('anatomical image with ROI overlays')
    pause(0.1)
end

conventional_maps.img_conventional = img_conventional;

%% T1 mapping

if do_T1 == 1
    
    lb = [0 -1 1];
    ub = [1.2*max(T1_data.img_stack_T1(:)) -0.2 5000];
    [T1_map_ms, RMSE_T1] = T1_IR_TD( T1_data.img_stack_T1, T1_data.TI_ms_v, T1_data.TD_ms, img_mask, lb, ub, disp_flag );
    
    figure(10); clf;
    imagesc( T1_map_ms ); axis image; colorbar();
    title('Conventional T1 map')
    drawnow;
    
    conventional_maps.T1_map_ms = T1_map_ms;
    conventional_maps.RMSE_T1 = RMSE_T1;
    
end

%% T2 mapping

if do_T2 == 1

    if n_T2_mod_params == 3
        lb = [0 0 0];
        ub = [5*max(T2_data.img_stack_T2(:)) 3000 0.5*max(T2_data.img_stack_T2(:))];
    elseif n_T2_mod_params == 2
        lb = [0 0];
        ub = [2*max(img_stack_T2(:)) 3000];
    end
    [T2_map_ms, RMSE_T2] = T2_MSE( T2_data.img_stack_T2, T2_data.TE_ms_v, img_mask, lb, ub, n_T2_mod_params, disp_flag );
    
    figure(20); clf;
    imagesc( T2_map_ms ); axis image; colorbar();
    title('Conventional T2 map')
    drawnow;
    
    conventional_maps.T2_map_ms = T2_map_ms;
    conventional_maps.RMSE_T2 = RMSE_T2;

end

%% graph cut fat-water separation

if do_fatsep == 1

    n_fatwater_sets = size( spgr_data.img_stack_me_spgr, 4 );
    fat_graphcut = zeros( n_row, n_col, n_fatwater_sets );
    water_graphcut = fat_graphcut;
    FSF_graphcut = fat_graphcut;
    B0_graphcut = fat_graphcut;
    R2star_graphcut = fat_graphcut;
    
    for ii = 1:n_fatwater_sets
        
        imDataParams.TE = spgr_data.TE_ms_v/1000; % sec
        n_TE = numel(spgr_data.TE_ms_v);
        imDataParams.images = reshape( spgr_data.img_stack_me_spgr(:,:,:,ii), [n_row, n_col, 1, 1, n_TE] );
        
        figure(29); clf;
        outParams = FattyRiot_fw_i2cm1i_3pluspoint_hernando_graphcut( imDataParams, algoParams );
        
        fat_graphcut(:,:,ii) = abs(outParams.species(2).amps);
        water_graphcut(:,:,ii) = abs(outParams.species(1).amps);
        FSF_graphcut(:,:,ii) = fat_graphcut(:,:,ii)./ ( fat_graphcut(:,:,ii) + water_graphcut(:,:,ii) );
        B0_graphcut(:,:,ii) = outParams.fieldmap;
        R2star_graphcut(:,:,ii) = outParams.r2starmap;
        
        figure(30); clf;
        imagesc( FSF_graphcut(:,:,ii) ); axis image; colorbar()
        title(sprintf ('graph cut FSF map series %d', ii ) );
        drawnow
        
    end

    conventional_maps.fat_graphcut = fat_graphcut;
    conventional_maps.water_graphcut = water_graphcut;
    conventional_maps.FSF_graphcut = FSF_graphcut;
    conventional_maps.B0_graphcut = B0_graphcut;
    conventional_maps.R2star_graphcut = R2star_graphcut;
    
end

%% save results

if exist('anatomical_img','var')
    conventional_maps.anatomical_img = anatomical_img;
end

save([dir_out fn_conv_maps], 'conventional_maps' );
