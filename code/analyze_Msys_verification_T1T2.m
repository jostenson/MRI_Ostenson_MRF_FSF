%%

clear, close all, clc;



%% parameters

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_MRF{1} = 'MRF_Msys_fixTE_proc_no_fat_sep.mat';
fn_MRF{2} = 'MRF_Msys_varTR_proc_no_fat_sep.mat';
fn_MRF{3} = 'MRF_Msys_varTE_proc_fat_sep.mat';
fn_conv = 'conventional_Msys.mat';
fn_ROI = 'conventional_Msys_ROIs.mat';
fn_save = 'Msys_T1T2_verification.mat';

create_ROIs = 0;
do_T1 = 1;
do_T2 = 1;
n_ROI = 14;
n_T2_mod_params = 2;
disp_flag = 1;

n_MRF = numel(fn_MRF);

%% load MRF data

T1_map = [];
T2_map = [];
M0_map = [];
for ii = 1:n_MRF
    
    load([dir_out fn_MRF{ii}]);
    try
        T1_map(:,:,ii) = output_MRF_match.T1_map;
        T2_map(:,:,ii) = output_MRF_match.T2_map;
        M0_map(:,:,ii) = output_MRF_match.M0_map;
        clear( 'output_MRF_match' );
    catch
        T1_map(:,:,ii) = output_MRFAT.T1_water_map;
        T2_map(:,:,ii) = output_MRFAT.T2_water_map;
        M0_map(:,:,ii) = output_MRFAT.M0_water_map;
        clear( 'output_MRFAT');
    end
    
end

%% load conventional

load([dir_in fn_conv])

img_conventional = sum( T2_data.SE.img_stack_T2, 3 );
n_row = size( img_conventional, 1 );
n_col = n_row;

%% load ROI map

copyfile([dir_in fn_ROI],[dir_out fn_ROI]);
load( [dir_out fn_ROI] );

if create_ROIs == 1
    
    ROI_stack = zeros( n_row, n_col, n_ROI );
    
    figure(1); clf;
    imagesc( img_conventional ); axis image; colormap(gray);
    
    for ii = 1:n_ROI
        
        title( sprintf('draw ROI %d',ii) );
        my_roi = roipoly();
%         ROI_stack(:,:,ii) = my_roi.* img_mask;
        ROI_stack(:,:,ii) = my_roi;
        
    end
    
    save( [dir_out fn_ROI], 'ROI_stack');
    
else
    
    load( [dir_out fn_ROI] )
    
end


img_mask = sum( ROI_stack, 3 );
img_mask( img_mask > 1 ) = 1;

figure(1); clf;
imagesc( img_conventional ); axis image; colormap(gray);
green = cat( 3, zeros(n_row,n_col), ones(n_row,n_col), zeros(n_row,n_col) );
hold on;
h = imshow(green);
alpha_map = img_mask;
set(h, 'AlphaData', alpha_map);
hold off
title('anatomical image with ROI overlays')
pause(0.1)


%% T1 mapping

if do_T1 == 1
    
    % IR
    
    lb = [0 -1 1];
    ub = [1.2*max(T1_data.SIR.img_stack_T1(:)) -0.2 5000];

    
    [T1_map_ms, RMSE_T1] = T1_IR_TD( T1_data.SIR.img_stack_T1, T1_data.SIR.TI_ms_v, T1_data.SIR.TD_ms, img_mask, lb, ub, disp_flag );
    
    figure(11); clf;
    imagesc( T1_map_ms ); axis image; colorbar();
    title('Conventional T1 map')
    drawnow;
    
    conventional_maps.SIR.T1_map_ms = T1_map_ms;
    conventional_maps.SIR.RMSE_T1 = RMSE_T1;
    
end

%% T2 mapping

if do_T2 == 1
        
        % SE
        
        if n_T2_mod_params == 3
            lb = [0 0 0];
            ub = [5*max(T2_data.SE.img_stack_T2(:)) 3000 0.5*max(T2_data.SE.img_stack_T2(:))];
        elseif n_T2_mod_params == 2
            lb = [0 0];
            ub = [5*max(T2_data.SE.img_stack_T2(:)) 3000];
        end
        [T2_map_ms, RMSE_T2] = T2_MSE( T2_data.SE.img_stack_T2, T2_data.SE.TE_ms_v, img_mask, lb, ub, n_T2_mod_params, disp_flag );
        
        figure(21); clf;
        imagesc( T2_map_ms ); axis image; colorbar();
        title('Conventional T2 map')
        drawnow;
        
        conventional_maps.SE.T2_map_ms = T2_map_ms;
        conventional_maps.SE.RMSE_T2 = RMSE_T2;
        
    
end

T1_ref = conventional_maps.SIR.T1_map_ms;
T2_ref = conventional_maps.SE.T2_map_ms;

%% get stats

T1_mean_meds = zeros(n_ROI,n_MRF+1);
T2_mean_meds = T1_mean_meds;
T1_std = T1_mean_meds;
T2_std = T1_mean_meds;

for ii = 1:n_ROI
    
    my_roi = logical( ROI_stack(:,:,ii) );
    for jj = 1:n_MRF + size( T2_ref, 3 )
        
        if jj <= n_MRF
            myT1 = T1_map(:,:,jj);
            myT2 = T2_map(:,:,jj);
            T1_mean_meds(ii,jj) = mean( myT1(my_roi) );
            T2_mean_meds(ii,jj) = mean( myT2(my_roi) );        
        elseif jj > n_MRF && jj <= n_MRF + 1
            myT1 = T1_ref;
            myT2 = T2_ref;
            T1_mean_meds(ii,jj) = median( myT1(my_roi) );
            T2_mean_meds(ii,jj) = median( myT2(my_roi) );      
        end
        
        T1_std(ii,jj) = std( myT1(my_roi) );
        T2_std(ii,jj) = std( myT2(my_roi) );
        
    end   
    
    
end

%% save results

save([dir_out fn_save],'T1_mean_meds','T2_mean_meds','T1_std','T2_std')




