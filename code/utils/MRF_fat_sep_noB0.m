function [output_MRF_fat_sep] = MRF_fat_sep_noB0( input_MRF_raw, input_MRF_img_recon, input_ksp_traj, input_MRF_dict, input_params )
%% apply ksp-img space SVD based MFI fat water separation

start_time = tic;

%% params

% MRF fat sep params
i_sign = input_params.i_sign;

% constants
gamma_bar_Hz_per_T = input_params.gamma_bar_Hz_per_T;
ppm_ref = input_params.ppm_ref;
B0_Tesla = input_params.B0_Tesla;
p_rel_ref_v = input_params.p_rel_ref_v;

fat_freq_Hz_v = p_rel_ref_v * gamma_bar_Hz_per_T * B0_Tesla;
fat_amps_v = input_params.fat_amps_v;
fat_amps_v = fat_amps_v./sum(fat_amps_v);
T1_fat_s_v = input_params.T1_fat_s_v;
T2_fat_s_v = input_params.T2_fat_s_v;
n_fat_pks = numel(fat_freq_Hz_v);

%% setup parallel pool

my_pool = parpool([1 128],'SpmdEnabled',false);

%% load data
tic;

MRF_img_data = input_MRF_img_recon.MRF_img_stack_coil_combined;
[n_row, n_col, n_TR ] = size( MRF_img_data );

disp(sprintf('MRF fat sep: completed coil-combined recon load in %.1f seconds', toc));

%% put into k-space
tic;

MRF_ksp_data = zeros( size( MRF_img_data )  );
for ii = 1:n_TR
    
    my_img = MRF_img_data(:,:,ii);
    my_ksp = fftshift( fftn( ifftshift( my_img ) ) );
    MRF_ksp_data(:,:,ii) = my_ksp;
    
end

clear MRF_img_data % save memory

% reshape data matrix
MRF_ksp_data = reshape( permute( MRF_ksp_data, [3 1 2] ), [ n_TR n_row*n_col ] );

disp(sprintf('MRF fat sep: completed coil-combined transform to k-space in %.1f seconds', toc));

%% get acq time map


%   INPUT: input.ksp_norm_1 = normalized k-space coords for dimension 1
%              ".ksp_norm_2 = normalized k-space coords for dimension 2
%              ".tacq = acquisition time of single trajectory interleaf
%              ".effMtx = pos integer is matrix size of k-space

%   OUTPUT: output.acq_time_map = acquisition time map

eval( sprintf( 'input_acq_time_map.ksp_norm_1 = input_ksp_traj.ksp_norm_%d;',input_params.permute_order(1) ) );
eval( sprintf( 'input_acq_time_map.ksp_norm_2 = input_ksp_traj.ksp_norm_%d;',input_params.permute_order(2) ) );
input_acq_time_map.tacq = input_MRF_raw.params.tacq_s;
input_acq_time_map.effMtx = n_row;

acq_time_map_struct = make_acq_time_map( input_acq_time_map );

acq_time_map_s = fftshift(acq_time_map_struct.acq_time_map);

output_MRF_fat_sep.acq_time_map_s = acq_time_map_s;


%% B1 constrained fitting

B1_unique_v = unique(input_MRF_dict.B1_v);
n_B1 = numel(B1_unique_v);
B1_map = input_MRF_img_recon.B1_map;

%% get water dictionary

tic;

r_water_SVD_space = size(input_MRF_dict.dict_compress,1);

Ur_water_dict_B1_stack = zeros( n_B1, n_TR, r_water_SVD_space );
for ii = 1:n_B1
    Ur_water_dict_B1_stack(ii,:,:) = input_MRF_dict.U_r(:,input_MRF_dict.B1_compress_dict_list_v == B1_unique_v(ii) ); % extract pefect B1
end


disp(sprintf('MRF fat sep: load of water dictionary in %.1f seconds', toc));

%% do k-space fat-water separation

tic;

TE_s_v = input_MRF_dict.TE_v/1000;

% get fat contrast by pks
TR_s_v = input_MRF_dict.TR_v/1000;
FA_deg_v = input_MRF_dict.FA_v;
delk = input_MRF_dict.delk;
szomega = input_MRF_dict.szomega;
phi_v = input_MRF_dict.phi_v;
TI_s = input_MRF_dict.TI/1000;
fat_contrast_by_pk = zeros( n_TR,n_fat_pks,n_B1 );
for jj = 1:n_B1
    for ii = 1:n_fat_pks
        
        T1_s = T1_fat_s_v(ii);
        T2_s = T2_fat_s_v(ii);
        my_FA_deg_v = FA_deg_v * B1_unique_v(jj);
        
        my_signal_v = EPG_MRF_SSFP( T1_s, T2_s, TE_s_v, TR_s_v, my_FA_deg_v, delk, n_TR, szomega, phi_v, TI_s );
        my_signal_v = my_signal_v./norm(my_signal_v);
        fat_contrast_by_pk(:,ii,jj) = fat_amps_v(ii).*my_signal_v(:);
        
    end
end

coeffs_measured_water = zeros( n_B1, n_row, n_col, r_water_SVD_space, 'single' );
coeffs_measured_fat = zeros(n_B1, n_row, n_col, 'single' );
parfor ii = 1:n_row
    disp( sprintf( 'MRF fat sep: beginning coefficient calc. in k-space for row %d of %d', ii, n_row ) );
    for jj = 1:n_col
        
        t1_v = TE_s_v + acq_time_map_s(ii,jj);
        
        % form signal basis
        my_idx = ( jj - 1 ) * n_col + ii;
        my_signal_v = MRF_ksp_data(:,my_idx);
        
        % form fat model
        my_phase_cor1 = exp( i_sign * 1i * 2 * pi * t1_v * fat_freq_Hz_v(:)' );
        
        for mm = 1:n_B1
            
            % form fat model continued
            my_fat_signal_model1_v = sum( fat_contrast_by_pk(:,:,mm).* my_phase_cor1, 2 );
            
            % normalize fat signal
            my_fat_signal_model1_v = my_fat_signal_model1_v./norm(my_fat_signal_model1_v);
            
            % form augmented design matrix
            my_water_basis = squeeze( Ur_water_dict_B1_stack( mm,:,: ) );
            A = [ my_water_basis my_fat_signal_model1_v ];
            
            % get water and fatq components
            beta_v = ( A' * A ) \ A' * my_signal_v;
            
            c_water_v = beta_v( 1:(end - 1) );
            c_fat_v = beta_v( end );
            
            % measured coefficients of k-space data
            my_water_coeffs = c_water_v;
            my_fat_coeffs = c_fat_v;
            coeffs_measured_water(mm,ii,jj,:) = my_water_coeffs;
            coeffs_measured_fat(mm,ii,jj) = my_fat_coeffs;
            
        end
    end
end

clear MRF_ksp_data % save memory

disp(sprintf('MRF fat sep: completed coefficient calc of k-space in %.1f seconds', toc));

%% if applicable, do IFT of coefficient k-spaces for all freq. mod. basis

tic;

coeffs_measured_maps_water = complex( zeros( size( coeffs_measured_water ), 'single' ) );
coeffs_measured_maps_fat = complex( zeros( size( coeffs_measured_fat ), 'single' ) );
for mm = 1:n_B1
    for ii = 1:r_water_SVD_space
        
        
        my_ksp = squeeze( coeffs_measured_water( mm,:,:,ii) );
        my_img = ifftshift( ifftn( fftshift( my_ksp ) ) );
        coeffs_measured_maps_water( mm,:,:,ii ) = my_img;
        
    end
    
    my_ksp = squeeze( coeffs_measured_fat (mm,:,:) );
    my_img = ifftshift( ifftn( fftshift( my_ksp ) ) );
    coeffs_measured_maps_fat( mm,:,:) = my_img;
    
end


clear coeffs_measured coeffs_measured_water coeffs_measured_fat % save memory

disp(sprintf('MRF fat sep: completed coefficient transformation to image domain in %.1f seconds', toc));


%% get T1 and T2 of water and M0 fat

tic;

B1_idx_map = zeros( n_row, n_col );
T1_fit_map = zeros( n_row, n_col );
T2_fit_map = zeros( n_row, n_col );
M0_fit_map = zeros( n_row, n_col );
M0_fat_fit_map = zeros( n_row, n_col );
for ii = 1:n_row
    for jj = 1:n_col
        
        % determine B1 constraint
        B1 = B1_map(ii,jj);
        [~, idx_B1_min] = min( abs( B1_unique_v(:) - B1 ) );
        B1_idx_map(ii,jj) = idx_B1_min;
        
        logical_idx_dict_v = input_MRF_dict.dict_list(:,3) == B1_unique_v(idx_B1_min);
        
        my_signal_v = squeeze( coeffs_measured_maps_water(idx_B1_min,ii,jj,:) );
        my_signal_norm_v = my_signal_v./ norm(my_signal_v);
        
        dict_compress = input_MRF_dict.dict_compress(:,logical_idx_dict_v);
        my_ip_v = dict_compress' * my_signal_norm_v(:); % determine basis coeff correlation
        
        [~,my_idx_max] = max( my_ip_v(:) );
        
        T1_fit_map(ii,jj) = input_MRF_dict.dict_list(my_idx_max,1);
        T2_fit_map(ii,jj) = input_MRF_dict.dict_list(my_idx_max,2);
        M0_fit_map(ii,jj) = dict_compress(:,my_idx_max)'*my_signal_v;
        M0_fat_fit_map(ii,jj) = coeffs_measured_maps_fat( idx_B1_min, ii,jj );
        
        
    end
end

FSF_fit_map = abs( M0_fat_fit_map )./ ( abs( M0_fat_fit_map ) + abs( M0_fit_map ) );

disp(sprintf('MRF fat sep: completed water T1 and T2 fits in %.1f seconds', toc));

%% save results

output_MRF_fat_sep.T1_water_map = T1_fit_map;
output_MRF_fat_sep.T2_water_map = T2_fit_map;
output_MRF_fat_sep.M0_water_map = M0_fit_map;
output_MRF_fat_sep.M0_fat_map = M0_fat_fit_map;
output_MRF_fat_sep.FSF_map = FSF_fit_map;

%% report total execution time

delete(my_pool);

disp(sprintf('MRF fat sep: TOTAL EXECUTION TIME: %.1f seconds', toc(start_time) ));

end