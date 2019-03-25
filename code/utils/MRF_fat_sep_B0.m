function [output_MRF_fat_sep] = MRF_fat_sep_B0( input_MRF_raw, input_MRF_img_recon, input_ksp_traj, input_MRF_dict, input_params )
%% apply ksp-img space SVD based MFI fat water separation
%  this code is not optimized

%  see calling scripts for input structure(s) specific details
%  input_MRF_raw is sruct that holds the raw data and parameters
%  input_MRF_img_recon is struct that holds image reconstructions
%  input_ksp_traj is struct that holds the k-space trajectory
%  input_MRF_dict is struct that holds the water dictionary
%  input_params is struct that holds fat-water fitting specifications (see
%  also below)

% output_MRF_fat_sep.img_water_cp is the B0-corrected water
% stack (save currently commented out)
% output_MRF_fat_sep.img_fat_cp is the B0-corrected fat stack
% (save currently commented out)
% output_MRF_fat_sep.B0_fit_map is the B0 map estimate
% output_MRF_fat_sep.T1_water_map is the fat-separated water T1 map estimate
% output_MRF_fat_sep.T2_water_map is the fat-separated water T2 map estimate
% output_MRF_fat_sep.M0_water_map is the (complex) magnetization density
% map estimate for water
% output_MRF_fat_sep.M0_fat_map is the (complex) magnetization density
% map estimate for fat
% output_MRF_fat_sep.FSF_map is the (real) fat signal fraction map estimate


start_time = tic;

%% params

% MRF fat sep params
i_sign = input_params.i_sign; % precession direction convention
n_B0_basis = input_params.n_B0_basis; % number of frequencies to use in MFI basis set
freq_MFI_max_Hz = input_params.freq_MFI_max_Hz; % maximum resonant frequency deviation from zero in Hz
freq_MFI_step_Hz = input_params.freq_MFI_step_Hz; % frequency step in B0 map
s_frac = input_params.s_frac; % fraction of singular value energy to keep when determining number of elements in basis derived from SVD
B0_range_Hz = input_params.B0_range_Hz; % maximum range of B0 (less than maximum frequency considered in basis)
try 
    sigma_filter = input_params.sigma_filter; % size of Gaussian blurring kernel to smooth B0 estimation
    sigma_filter_does_not_exist = 0;
    if sigma_filter == 0
        sigma_filter_does_not_exist = 1;
    end
catch
    sigma_filter = 0;
    sigma_filter_does_not_exist = 1;
end

% constants
gamma_bar_Hz_per_T = input_params.gamma_bar_Hz_per_T;
B0_Tesla = input_params.B0_Tesla;
p_rel_ref_v = input_params.p_rel_ref_v;

fat_freq_Hz_v = p_rel_ref_v * gamma_bar_Hz_per_T * B0_Tesla;
fat_amps_v = input_params.fat_amps_v;
fat_amps_v = fat_amps_v./sum(fat_amps_v);
T1_fat_s_v = input_params.T1_fat_s_v;
T2_fat_s_v = input_params.T2_fat_s_v; 
n_fat_pks = numel(fat_freq_Hz_v);

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
    my_ksp = ifftshift( fftn( fftshift( my_img ) ) );
    MRF_ksp_data(:,:,ii) = my_ksp;
    
end

clear MRF_img_data % save memory

% reshape data matrix
MRF_ksp_data = reshape( permute( MRF_ksp_data, [3 1 2] ), [ n_TR n_row*n_col ] );

disp(sprintf('MRF fat sep: completed coil-combined transform to k-space in %.1f seconds', toc));

%% get acq time map

eval( sprintf( 'input_acq_time_map.ksp_norm_1 = input_ksp_traj.ksp_norm_%d;',input_params.permute_order(1) ) );
eval( sprintf( 'input_acq_time_map.ksp_norm_2 = input_ksp_traj.ksp_norm_%d;',input_params.permute_order(2) ) );
input_acq_time_map.tacq = input_MRF_raw.params.tacq_s;
input_acq_time_map.effMtx = n_row;

acq_time_map_struct = make_acq_time_map( input_acq_time_map );

acq_time_map_s = fftshift(acq_time_map_struct.acq_time_map);

output_MRF_fat_sep.acq_time_map_s = acq_time_map_s;

%% determine freq domain via MFI

tic;

TE_s_v = input_MRF_dict.TE_v(:)/1000;
TE_base_s = min(TE_s_v);
acq_window_s = input_MRF_raw.params.tacq_s;
TE_max_s = max(TE_s_v);

params_MFI.method = 'DFT';
params_MFI.M = n_B0_basis;
params_MFI.t0 = TE_base_s;
params_MFI.t_end = acq_window_s + TE_max_s;
params_MFI.t1 = TE_base_s;
params_MFI.t2 = 1.2*( params_MFI.t_end ) - params_MFI.t_end;
params_MFI.N = input_ksp_traj.n_rd;
params_MFI.f_max = freq_MFI_max_Hz;
params_MFI.f_step = freq_MFI_step_Hz;

results_MFI_coeff = MFI_coeffs(params_MFI);

disp(sprintf('MRF fat sep: completed MFI coefficient calc. in %.1f seconds', toc));

%% B1 constrained fitting

B1_unique_v = unique(input_MRF_dict.B1_v);
n_B1 = numel(B1_unique_v);
B1_map = input_MRF_img_recon.B1_map;

%% create off-resonance dictionary

disp('MRF fat sep: beginning composite dictionary creation and compression');

tic;

for jj = 1:n_B1
    
    fprintf('MRF fat sep: computing composite dict and compression for B1 %.2f\n', B1_unique_v(jj) );
    
    my_dict_norm = input_MRF_dict.dict_norm( :,( input_MRF_dict.dict_list(:,3) == B1_unique_v(jj) ) ); % select B1

    n_dict = size(my_dict_norm,2);
    
    dict_composite = complex( zeros(n_TR, n_B0_basis*n_dict) );
    
    t_v = TE_s_v(:);
    for ii = 1:n_B0_basis
        
        my_basis_v = exp( i_sign * 1i * 2 * pi * t_v * results_MFI_coeff.delta_f_v(ii) );
        my_basis_mat = repmat( my_basis_v, [ 1 n_dict ] );
        my_index_v = ( ( ii - 1 )*n_dict + 1 ):( ii*n_dict );
        dict_composite(:,my_index_v) = my_basis_mat.* my_dict_norm;
        
    end
    
    % compress
    [U,S,~] = svd(dict_composite,'econ');
        
    if jj == 1 % first B1 determines rank of truncated SVD
        s_v = diag(S);
        s2_cum_rel_v = cumsum(s_v.^2)/sum(s_v.^2);
        r_SVD_space = sum(s2_cum_rel_v <= s_frac);
        Ur = complex( zeros( n_TR, r_SVD_space, n_B1 ) );
    end
    
    Ur(:,:,jj) = U(:,1:r_SVD_space);
    
%     error_mat = dict_composite - Ur(:,:,jj) * Ur(:,:,jj)' * dict_composite;
%     error_v = vecnorm( error_mat );
%     fprintf( 'MRF fat sep: max and median of residual norm for current B1 are %.4f and %.4f\n', max(error_v), median(error_v) )
    
    clear dict_composite
end

disp(sprintf('MRF fat sep: completed composite dictionary creation and compression in %.1f seconds', toc));


%% get water dictionary

tic;

Ur_water_dict_B1_stack = zeros( n_B1, n_TR, size(input_MRF_dict.dict_compress,1) );
for ii = 1:n_B1
    Ur_water_dict_B1_stack(ii,:,:) = input_MRF_dict.U_r(:,input_MRF_dict.B1_compress_dict_list_v == B1_unique_v(ii) ); % extract basis for given B1
end

disp(sprintf('MRF fat sep: load of water dictionary in %.1f seconds', toc));

%% setup parallel pool

my_pool = parpool([1 128],'SpmdEnabled',false);

%% get coeffs for all of k-space using basis from off-res dictionary for each freq. modulated basis
% and separate water and fat

tic;

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

delta_f_Hz_v = results_MFI_coeff.delta_f_v;
coeffs_measured_water = complex( zeros( n_B1, n_row, n_col, r_SVD_space, n_B0_basis, 'single' ) );
coeffs_measured_fat = complex( zeros(n_B1, n_row, n_col, r_SVD_space, n_B0_basis, 'single' ) );
coeffs_measured = complex( zeros( n_B1, n_row, n_col, r_SVD_space, n_B0_basis, 'single' ) );
cond_map = zeros( n_row, n_col, n_B1, 'single' );
for ii = 1:n_row
    disp( sprintf( 'MRF fat sep: beginning coefficient calc. in k-space for row %d of %d', ii, n_row ) );
    data_inds_v = sub2ind([n_row n_col],ii*ones(n_col,1),(1:n_col)');
    MRF_ksp_data_sub = MRF_ksp_data(:,data_inds_v);
    coeffs_measured_water_sub = complex( zeros( n_B1, n_col, r_SVD_space, n_B0_basis, 'single' ) );
    coeffs_measured_fat_sub = complex( zeros( n_B1, n_col, r_SVD_space, n_B0_basis, 'single' ) );
    coeffs_measured_sub = complex( zeros( n_B1, n_col, r_SVD_space, n_B0_basis, 'single' ) );
    cond_map_sub = zeros( n_col,n_B1, 'single' );
    parfor jj = 1:n_col
        
        t1_v = TE_s_v + acq_time_map_s(ii,jj);
        % start form fat model
        my_phase_cor1 = exp( i_sign * 1i * 2 * pi * t1_v * fat_freq_Hz_v(:)' );
        
        for mm = 1:n_B1
            
            % form fat model continued
            my_fat_signal_model1_v = sum( fat_contrast_by_pk(:,:,mm).* my_phase_cor1, 2 );
            
            % normalize fat signal
            my_fat_signal_model1_v = my_fat_signal_model1_v./norm(my_fat_signal_model1_v);
            
            % form augmented design matrix
            my_water_basis = squeeze( Ur_water_dict_B1_stack( mm,:,: ) );
            A = [ my_water_basis my_fat_signal_model1_v ];
            pA = ( A' * A ) \ A';
%             cond_map_sub(jj,mm) = cond(A);
            
            % form model without phase accrual due to acquisition timing
            t2_v = TE_s_v;
            my_phase_cor2 = exp( i_sign * 1i * 2 * pi * t2_v * fat_freq_Hz_v(:)' );
            my_fat_signal_model2_v = sum( fat_contrast_by_pk(:,:,mm).* my_phase_cor2, 2 );
            
            % normalize fat signal
            my_fat_signal_model2_v = my_fat_signal_model2_v./norm(my_fat_signal_model2_v);
            
            for kk = 1:n_B0_basis
                               
                % form phasor
                j_v = exp( -i_sign * 1i * 2 * pi * delta_f_Hz_v(kk) * t1_v(:) ); % conjugate i_sign
                
                % form signal basis
                my_signal_v = MRF_ksp_data_sub(:,jj);
                signal_basis_v = my_signal_v(:).* j_v;                             
                
                % get water and fat components
                % beta_v = pinv( A ) * signal_basis_v;
                beta_v = pA * signal_basis_v;
                
                c_water_v = beta_v( 1:(end - 1) );
                c_fat_v = beta_v( end );
                
                % get residual
                my_resid_v = signal_basis_v - A * beta_v;                             
                
                % measured coefficients of k-space data
                my_water_coeffs = Ur(:,:,mm)' * my_water_basis * c_water_v;
                my_fat_coeffs = Ur(:,:,mm)' * my_fat_signal_model2_v * c_fat_v;
                
                coeffs_measured_water_sub(mm,jj,:,kk) = my_water_coeffs;
                coeffs_measured_fat_sub(mm,jj,:,kk) = my_fat_coeffs;
                coeffs_measured_sub(mm,jj,:,kk) = ( Ur(:,:,mm)' * my_resid_v ) + my_water_coeffs + my_fat_coeffs;
                
            end
        end
        
    end % end parfor
    
    coeffs_measured_water(:,ii,:,:,:) = coeffs_measured_water_sub;
    coeffs_measured_fat(:,ii,:,:,:) = coeffs_measured_fat_sub;
    coeffs_measured(:,ii,:,:,:) = coeffs_measured_sub;
    cond_map(ii,:,:) = cond_map_sub;
    
end

clear MRF_ksp_data % save memory

disp(sprintf('MRF fat sep: completed coefficient calc of k-space in %.1f seconds', toc));

%% do IFT of coefficient k-spaces for all freq. mod. basis

tic;

coeffs_measured_maps_all_bases = complex( zeros( size( coeffs_measured ), 'single' ) );
coeffs_measured_maps_water = complex( zeros( size( coeffs_measured ), 'single' ) );
coeffs_measured_maps_fat = complex( zeros( size( coeffs_measured ), 'single' ) );
for ii = 1:r_SVD_space
    for jj = 1:n_B0_basis
        for mm = 1:n_B1
            
            my_ksp = squeeze( coeffs_measured(mm,:,:,ii,jj) );
            my_img = ifftshift( ifftn( fftshift( my_ksp ) ) );
            coeffs_measured_maps_all_bases( mm,:,:,ii,jj ) = my_img;
            
            my_ksp = squeeze( coeffs_measured_water( mm,:,:,ii,jj) );
            my_img = ifftshift( ifftn( fftshift( my_ksp ) ) );
            coeffs_measured_maps_water( mm,:,:,ii,jj ) = my_img;
            
            my_ksp = squeeze( coeffs_measured_fat (mm,:,:,ii,jj) );
            my_img = ifftshift( ifftn( fftshift( my_ksp ) ) );
            coeffs_measured_maps_fat( mm,:,:,ii,jj ) = my_img;
            
        end
    end
end

clear coeffs_measured coeffs_measured_water coeffs_measured_fat % save memory

disp(sprintf('MRF fat sep: completed coefficient transformation to image domain in %.1f seconds', toc));

%% get MFI coefficents over B0 range and relevant spacing

freq_fine_Hz_v = results_MFI_coeff.delta_f_fine_v;
coeffs_MFI = results_MFI_coeff.c_DFT;

%% model fat signal

fat_model_v = zeros(n_TR,n_B1);
for jj = 1:n_B1
    for ii = 1:n_fat_pks
        
        T1_s = T1_fat_s_v(ii);
        T2_s = T2_fat_s_v(ii);
        
        my_FA_deg_v = FA_deg_v * B1_unique_v(jj);
        my_signal_v = EPG_MRF_SSFP( T1_s, T2_s, TE_s_v, TR_s_v, my_FA_deg_v, delk, n_TR, szomega, phi_v, TI_s );
        my_signal_v = my_signal_v./norm(my_signal_v);
        fat_model_v(:,jj) = fat_model_v(:,jj) + fat_amps_v(ii).*my_signal_v(:).*exp( i_sign * 1i * 2 * pi * fat_freq_Hz_v(ii) * TE_s_v );
        
    end
    
    % normalize fat signal
    fat_model_v(:,jj) = fat_model_v(:,jj)./norm( fat_model_v(:,jj) );
end

%% form augmented water dictionary w/fat

r_water_SVD_space = size(input_MRF_dict.dict_compress,1);

dict_norm_augmented = zeros( n_TR, ( r_water_SVD_space + 1 ), n_B1 );
for ii = 1:n_B1
    
    dict_norm_augmented(:,:,ii) = [squeeze(Ur_water_dict_B1_stack(ii,:,:)) fat_model_v(:,ii)];

end

%% smooth coefficient maps
tic;

if sigma_filter_does_not_exist == 1
    
    coeffs_measured_maps_all_bases_smoothed = coeffs_measured_maps_all_bases;
    
else
    
    coeffs_measured_maps_all_bases_smoothed = zeros( size( coeffs_measured_maps_all_bases ), 'single' );
    PSF_coeff = fspecial( 'gaussian', round(sigma_filter * 5), sigma_filter );
    for ii = 1:r_SVD_space
        for jj = 1:n_B0_basis
            for mm = 1:n_B1
                
                my_coeff_map = squeeze( coeffs_measured_maps_all_bases( mm,:,:,ii,jj ) );
                my_coeff_map_smoothed = imfilter( my_coeff_map, PSF_coeff );
                
                coeffs_measured_maps_all_bases_smoothed( mm,:,:,ii,jj ) = my_coeff_map_smoothed;
                
            end
        end
    end
    
end

clear coeffs_measured_maps_all_bases; % save memory

disp(sprintf('MRF fat sep: completed coefficient filtering in image domain in %.1f seconds', toc));

%% for each data point (image domain) cycle through B0 to determine which B0 minimizes the objective

tic;

n_r_B0_fit = r_water_SVD_space + 1;

% determine # of coefficients for B0 fit
Ur_B0_fit = Ur(:,1:n_r_B0_fit,:);

% form operative VARPRO matrix
VARPRO_mat = zeros( n_r_B0_fit, n_TR, n_B1 );
for ii = 1:n_B1

    VARPRO_mat(:,:,ii) = Ur_B0_fit(:,:,ii)' - Ur_B0_fit(:,:,ii)' * dict_norm_augmented(:,:,ii) * pinv( dict_norm_augmented(:,:,ii) );
    
end

% restrict off-resonance range
idx_freq_restrict_v = find( freq_fine_Hz_v >= B0_range_Hz(1) & freq_fine_Hz_v <= B0_range_Hz(2) );
fine_freq_restrict_Hz_v = freq_fine_Hz_v( idx_freq_restrict_v ); 
n_freq_restrict = numel( fine_freq_restrict_Hz_v );
coeffs_MFI_restricted = coeffs_MFI(:,idx_freq_restrict_v);

% loop over test signals
B0_fit_map = zeros( n_row,n_col );
B1_idx_map = zeros( n_row,n_col );
idx_B0_solution_map = zeros( n_row,n_col );
for ii = 1:n_row
    disp( sprintf( 'MRF fat sep: beginning B0 fit for row %d of %d', ii, n_row ) );
    for jj = 1:n_col
        
        % determine B1 constraint
        B1 = B1_map(ii,jj);
        [~, idx_B1_min] = min( abs( B1_unique_v(:) - B1 ) );
        B1_idx_map(ii,jj) = idx_B1_min;
        
        my_coeffs = squeeze( coeffs_measured_maps_all_bases_smoothed(idx_B1_min,ii,jj,:,:) );
        my_residual_mat = VARPRO_mat(:,:,idx_B1_min) * Ur(:,:,idx_B1_min) * my_coeffs * coeffs_MFI_restricted;
        
        % loop over fine frequency discretization to residual magnitude
        r2_v = zeros( n_freq_restrict, 1);
        for kk = 1:n_freq_restrict;
            
            r2_v(kk) = my_residual_mat(:,kk)'*my_residual_mat(:,kk);
            
        end
        
        %         figure(10);
        %         plot(fine_freq_restrict_Hz_v,r2_v);
        %         xlabel('B0 (Hz)')
        %         ylabel('objective value')
        %         drawnow
        
        
        [r2_min, idx_min] = min( r2_v(:) );
        
        B0_fit_map(ii,jj) = -fine_freq_restrict_Hz_v(idx_min);
        idx_B0_solution_map(ii,jj) = idx_min;
        
               
    end
end

clear coeffs_measured_maps_all_bases % save memory

disp(sprintf('MRF fat sep: completed B0 fit in %.1f seconds', toc));

%% calculate water and fat data sets

tic;

img_water_cp = zeros( n_row, n_col, n_TR );
img_fat_cp = zeros( n_row, n_col, n_TR );
M0_fat_fit_map = zeros( n_row, n_col );
for ii = 1:n_row
    for jj = 1:n_col
        
        % determine B1 constraint
        idx_B1_min = B1_idx_map(ii,jj);
        
        my_coeffs_water = squeeze( coeffs_measured_maps_water(idx_B1_min,ii,jj,:,:) );
        my_coeffs_fat = squeeze( coeffs_measured_maps_fat(idx_B1_min,ii,jj,:,:) );
        
        my_coeffs_water = my_coeffs_water*coeffs_MFI_restricted(:,idx_B0_solution_map(ii,jj));
        my_coeffs_fat = my_coeffs_fat*coeffs_MFI_restricted(:,idx_B0_solution_map(ii,jj));
        
        signal_water_v = Ur(:,:,idx_B1_min)*my_coeffs_water;
        signal_fat_v = Ur(:,:,idx_B1_min)*my_coeffs_fat;
        
        img_water_cp(ii,jj,:) = signal_water_v;
        img_fat_cp(ii,jj,:) = signal_fat_v;
                
        M0_fat_fit_map(ii,jj) = fat_model_v(:,idx_B1_min)' * signal_fat_v(:);
        
    end
end

disp(sprintf('MRF fat sep: completed B0 corrected fat-water image stack formation in %.1f seconds', toc));

%% get T1 and T2 of water

tic;

T1_fit_map = zeros( n_row, n_col );
T2_fit_map = zeros( n_row, n_col );
M0_fit_map = zeros( n_row, n_col );
for ii = 1:n_row
    for jj = 1:n_col
        
        idx_B1_min = B1_idx_map(ii,jj);
        logical_idx_B1_col_v = input_MRF_dict.B1_compress_dict_list_v(:) == B1_unique_v(idx_B1_min);
    
        my_signal_v = squeeze( img_water_cp(ii,jj,:) );
        my_signal_norm_v = my_signal_v / norm( my_signal_v );
        
        my_water_basis = input_MRF_dict.U_r(:,logical_idx_B1_col_v);
        
        logical_idx_dict_v = input_MRF_dict.dict_list(:,3) == B1_unique_v(idx_B1_min);
        
        dict_compress = input_MRF_dict.dict_compress(:,logical_idx_dict_v);
        testR_v = my_water_basis' * my_signal_norm_v(:); % project onto SVD space
        my_ip_v = dict_compress' * testR_v; % determine basis coeff correlation
        
        [~,my_idx_max] = max( my_ip_v(:) );
        
        T1_fit_map(ii,jj) = input_MRF_dict.dict_list(my_idx_max,1);
        T2_fit_map(ii,jj) = input_MRF_dict.dict_list(my_idx_max,2);
        M0_fit_map(ii,jj) = my_signal_v'*input_MRF_dict.dict_norm(:,my_idx_max);
    
    end
end

FSF_fit_map = abs( M0_fat_fit_map )./ ( abs( M0_fat_fit_map ) + abs( M0_fit_map ) );

disp(sprintf('MRF fat sep: completed water T1 and T2 fits in %.1f seconds', toc));

%% save results

% output_MRF_fat_sep.img_water_cp = img_water_cp; % save storage
% output_MRF_fat_sep.img_fat_cp = img_fat_cp; % save storage
output_MRF_fat_sep.B0_fit_map = B0_fit_map;
output_MRF_fat_sep.T1_water_map = T1_fit_map;
output_MRF_fat_sep.T2_water_map = T2_fit_map;
output_MRF_fat_sep.M0_water_map = M0_fit_map;
output_MRF_fat_sep.M0_fat_map = M0_fat_fit_map;
output_MRF_fat_sep.FSF_map = FSF_fit_map;

%% report total execution time

delete(my_pool);

disp(sprintf('MRF fat sep: TOTAL EXECUTION TIME: %.1f seconds', toc(start_time) ));

end