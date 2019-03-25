%% create synthetic unbalanced SSFP variable TR MRF image
% intra-TR phase evolutions occur after TE

%% parameters specified in batch_simulations.m

%% confirm 1 B0

if B0_Hz_v(1) ~= B0_Hz_v(2)
    error( 'For MRF varTR, B0 must currently be uniform.' );
end

%% starting image

img_start = phantom( N );

img_start = round( img_start, 1 );

figure(1); clf;
imagesc( img_start ); axis image; colormap(gray);
title('starting image')
drawnow

%% get phantom image segmentation

mask = true( N, N, n_segs );
for ii = 1:n_segs
   
    my_mask = false( N );
    my_mask( img_start == segs_v(ii) ) = true;
    mask(:,:,ii) = my_mask;
    
end

img_segs = zeros( N, N, n_segs );
figure(10); clf;
for ii = 1:n_segs
    
    imagesc( mask(:,:,ii) ); axis image;
    title( sprintf( 'Phantom segment %d',ii ) )
    pause(0.5)
    
    img_segs(:,:,ii) = double( mask(:,:,ii) );
    
end

%% define fat components and weights

img_fat = zeros( N, N, n_segs );
img_water = img_fat;
for ii = 1:n_segs

    my_img_fat = img_segs(:,:,ii) * fsf_v(ii);
    my_img_water = img_segs(:,:,ii) * (1 - fsf_v(ii));
    
    figure(2); clf;
    subplot(121)
    imagesc( my_img_water ); axis image; colorbar(); caxis([0 max([my_img_water(:); my_img_fat(:)]) ] )
    title( sprintf( 'true water image segment %d', ii ) )
    subplot(122)
    imagesc( my_img_fat ); axis image; colorbar(); colormap(gray);caxis([0 max([my_img_water(:); my_img_fat(:)]) ] )
    title( sprintf( 'true fat image image segment %d', ii ) )
    pause(0.5)

    img_fat(:,:,ii) = my_img_fat;
    img_water(:,:,ii) = my_img_water;
    
end

img_fat = sum( img_fat, 3 );

%% replicate for total number of TRs

img_stack_water = single( repmat( img_water, [1 1 1 n_TR] ) );
img_stack_fat = single( repmat( img_fat, [1 1 n_TR n_fat_pks] ) );

img_stack_water = permute( img_stack_water, [1 2 4 3] );

%% define B1 map

[X, ~] = meshgrid( linspace( B1_v(1), B1_v(2), 2*N ) );
    
for ii = 1:2*N
    B1_grid(ii,:) = circshift( X(ii,:), -ii+1, 2 );
end

B1_true = B1_grid( 1:N,1:N );
B1_true = rot90(B1_true);

figure(8); clf;
imagesc( B1_true ); axis image; colorbar();
title('B1 map')
drawnow

output_B1_data.B1_map = B1_true;

%% define MRF contrast evolution due to T1/T2 for each water type

water_basis_signals = zeros( n_TR, n_segs );
for ii = 1:n_segs
    
    fprintf('Calculating MRF water contrast for segment %d\n',ii);
    
    B0_kHz = B0_Hz_v(1)/1000;
    
    my_seg = logical( img_segs(:,:,ii) );
    my_B1_seg = B1_true.*double(my_seg);
    my_B1_v = B1_true( my_seg );
    my_B1_unique_v = unique( my_B1_v );
    my_n_B1 = numel( my_B1_unique_v );
    for jj = 1:my_n_B1
        
        my_B1 = my_B1_unique_v(jj);
        my_B1_reps = sum( my_B1 == my_B1_v );
        my_idx_reps_v = find( my_B1 == my_B1_seg );  
        
        my_signal_v = EPG_MRF_SSFP_B0( T1_water_ms_v(ii), T2_water_ms_v(ii), TE_ms_v, TR_ms_v, my_B1*FA_deg_v, delk, n_TR, szomega, phi_v, TI_ms, B0_kHz, i_sign );
        my_signal_v = my_signal_v./norm(my_signal_v);

%       figure(10); clf;
%       plot(abs(my_signal_v));
%       title(sprintf( 'mag of true water signal for seg %d',ii ) )
%       drawnow
    
        [my_row_idx, my_col_idx] = ind2sub( [N N], my_idx_reps_v );
        for kk = 1:my_B1_reps
            
            my_row = my_row_idx(kk);
            my_col = my_col_idx(kk);  
            img_stack_water(my_row,my_col,:,ii) = my_signal_v(:).* squeeze( img_stack_water(my_row,my_col,:,ii) );
        
        end
        
    end
    
end

img_stack_water = sum( img_stack_water, 4 );

%% define acquisition time map

load(fn_ksp_traj)
input_ksp_traj = output_ksp_traj;

eval( sprintf( 'input_acq_time_map.ksp_norm_1 = input_ksp_traj.ksp_norm_%d;',permute_order(1) ) );
eval( sprintf( 'input_acq_time_map.ksp_norm_2 = input_ksp_traj.ksp_norm_%d;',permute_order(2) ) );
input_acq_time_map.tacq = tacq_ms;
input_acq_time_map.effMtx = N;

acq_time_map_struct = make_acq_time_map( input_acq_time_map );
acq_time_map_struct.acq_time_map = fftshift( acq_time_map_struct.acq_time_map );

acq_time_map_ms = acq_time_map_struct.acq_time_map;
acq_time_map_ms( isnan( acq_time_map_ms ) ) = 0;

figure(30); clf;
imagesc( acq_time_map_ms ); axis image; colorbar()
title('acquisition time map in ms')
drawnow

acq_time_map_s = acq_time_map_ms./1000;

%%

B1_unique_v = unique( B1_true(:) );
n_B1_unique = numel( B1_unique_v );
ksp_fat_comp = complex( zeros( N,N,n_TR, 'single' ) );
fat_mask = double( logical( img_fat ) );
B1_true_fat = B1_true.*fat_mask;
my_fat_stack_wghts = single( repmat( img_fat, [1 1 n_TR n_fat_pks] ) );

for ii = 1:n_B1_unique
    
    tic;
    fprintf('Calculating MRF fat contrast for unqiue B1 # %d\n',ii);
    
    my_B1 = B1_unique_v(ii);
    my_B1_mask = zeros(N);
    my_B1_mask( my_B1 == B1_true(:) ) = 1;
    
    if any( my_B1 == B1_true_fat(:) )
        
        my_B1_mask = repmat( my_B1_mask, [1 1 n_TR] );
        % image space
        
        my_fat_contrasts = complex( zeros( n_TR, n_fat_pks, 'single' ) );
        my_fat_stack = complex( zeros( N,N,n_TR,n_fat_pks, 'single' ) );
        fprintf('Calculating MRF fat contrast for unqiue B1 # %d, beginning MRF contrast...\n',ii);
        
        for jj = 1:n_fat_pks
            
            T1_ms = 1000*T1_fat_s_v(jj);
            T2_ms = 1000*T2_fat_s_v(jj);
            
            my_pk_B0_kHz = fat_freq_Hz_v(jj)/1000;
            my_signal_v = EPG_MRF_SSFP_B0( T1_ms, T2_ms, TE_ms_v, TR_ms_v, my_B1*FA_deg_v, delk, n_TR, szomega, phi_v, TI_ms, my_pk_B0_kHz, i_sign );
            my_signal_v = fat_amps_v(jj) * my_signal_v./ norm( my_signal_v );
            
            my_fat_contrasts(:,jj) = my_signal_v(:);
            
            my_contrast = repmat( my_signal_v(:), [1 N^2] );
            my_contrast = reshape( permute( my_contrast, [2 1] ), [N N n_TR] );
            
            my_fat_stack(:,:,:,jj) = my_contrast.* my_fat_stack_wghts(:,:,:,jj).* my_B1_mask;
            
        end
        fprintf('Calculating MRF fat contrast for unqiue B1 # %d, finished MRF contrast.\n',ii);       

        % determine fat norm map
        
        my_fat_norm_factors = zeros( N, N );
        for jj = 1:N
            for kk = 1:N
                
                acq_time_s = acq_time_map_s(jj,kk);
                
                my_signal_v = zeros(n_TR,1);
                for ll = 1:n_fat_pks
                    f_Hz = fat_freq_Hz_v(ll);
                    phasor_v = exp( i_sign * 1i * 2 * pi * f_Hz * acq_time_s ); % phase accrual up to TE already accounted for
                    fat_pk_v = my_fat_contrasts(:,ll).* phasor_v;
                    my_signal_v = my_signal_v + fat_pk_v;
                end
                
                my_fat_norm_factors(jj,kk) = 1./norm(my_signal_v);
                
            end
        end
        
        my_fat_norm_factors( my_fat_norm_factors == Inf ) = 0;
        
        % transform to k-space and modulate phase
        fprintf('Calculating MRF fat contrast for unqiue B1 # %d, beginning phase contrast and norm...\n',ii);
        for jj = 1:n_TR
            
            for kk = 1:n_fat_pks
                
                my_ksp = ifftshift( fftn( fftshift( my_fat_stack(:,:,jj,kk) ) ) );
                f_Hz = fat_freq_Hz_v(kk);
%                 phasor_map = exp( i_sign * 1i * 2 * pi * f_Hz * ( acq_time_map_s + TE_ms_v(jj)/1000 ) );
                phasor_map = exp( i_sign * 1i * 2 * pi * f_Hz * acq_time_map_s ); % phase accrual up to TE already accounted for
                my_ksp = my_ksp.* phasor_map.* my_fat_norm_factors;
                ksp_fat_comp(:,:,jj) = ksp_fat_comp(:,:,jj) + my_ksp;
                
            end
            
        end
        fprintf('Calculating MRF fat contrast for unqiue B1 # %d, finished contrast and norm.\n',ii);

    end
    toc;
    
end


%% transform water components to k-space

ksp_water = complex( zeros( size( img_stack_water ), 'single' ) );
for ii = 1:n_TR
    
    ksp_water(:,:,ii) = ifftshift( fftn( fftshift( img_stack_water(:,:,ii) ) ) );
    
end

figure(20); clf;
subplot(121)
imagesc( log( abs( ksp_water(:,:,1) ) ) ); axis image;
title('k-space of water TR 1');
subplot(122)
imagesc( log( abs( ksp_fat_comp(:,:,1) ) ) ); axis image;
title('k-space of fat TR 1');


%% transform to image domain

img_stack_fat_comp = complex( zeros( size( ksp_fat_comp ), 'single' ) );
for ii = 1:n_TR
    
    my_ksp_water = ksp_water(:,:,ii);
    my_ksp_fat = ksp_fat_comp(:,:,ii);
          
    img_stack_water(:,:,ii) = ifftshift( ifftn( fftshift( my_ksp_water ) ) );
    img_stack_fat_comp(:,:,ii) = ifftshift( ifftn( fftshift( my_ksp_fat ) ) );
    
end

figure(40); clf;
subplot(121);
imagesc( abs( sum( img_stack_water, 3 ) ) ); axis image; colormap(gray);
title('mag sum of water MRF stack')
subplot(122);
imagesc( abs( sum( img_stack_fat_comp, 3 ) ) ); axis image; colormap(gray);
title('mag sum of fat MRF stack')
drawnow

img_stack = img_stack_water + img_stack_fat_comp;

figure(41); clf;
imagesc( abs( sum( img_stack, 3 ) ) ); axis image; colormap(gray);
title('mag sum of fat-water MRF stack')
drawnow

%% B0 perturbations

if B0_pert == 1     

    %% angled B0 map
    
    [X, ~] = meshgrid( linspace( B0_Hz_v(1), B0_Hz_v(2), 2*N ) );
    
    for ii = 1:2*N
        B0_grid(ii,:) = circshift( X(ii,:), -ii+1, 2 );
    end
    
    B0_true_Hz = B0_grid( 1:N,1:N );
    
    B0_uniques_Hz_v = unique( B0_true_Hz(:) );
    n_B0 = numel( B0_uniques_Hz_v );
    
    ksp_stack_B0 = complex( zeros( size( img_stack ), 'single' ) );
    for ii = 1:n_B0
        fprintf( 'Solving for ksp for unique B0 value #%d\n',ii );
        
        idx_B0_v = find( B0_true_Hz == B0_uniques_Hz_v(ii) );
        for jj = 1:n_TR
            
            my_img = img_stack(:,:,jj);            
            my_img_B0 = zeros( N, 'single' );
            my_img_B0( idx_B0_v ) = my_img( idx_B0_v );
            my_ksp = ifftshift( fftn( fftshift( my_img_B0 ) ) );
            my_t_map_s = acq_time_map_s;
            my_ksp = my_ksp.* exp( my_t_map_s * ( i_sign * 1i * 2 * pi * B0_uniques_Hz_v(ii) ) );
            
            ksp_stack_B0(:,:,jj) = ksp_stack_B0(:,:,jj) + my_ksp;
            
        end
        
    end
    
    ksp_stack = ksp_stack_B0;
      
    
    % convert to image domain
    
    for ii = 1:n_TR
        
        img_stack(:,:,ii) = ifftshift( ifftn( fftshift( ksp_stack(:,:,ii) ) ) );
        
    end
    
    %%
    
    figure(42); clf;
    imagesc( abs( sum( img_stack, 3 ) ) ); axis image; colormap(gray);
    title('mag sum of fat-water MRF stack with B0')
    drawnow
    
end

%% perform forward transform on fully sampled image data

% get k-space trajs

Nrd = size( input_ksp_traj.ksp_norm_1, 1 );

k_1 = zeros( Nrd, Nint );
k_2 = k_1;

for ii = 1:Nint

    k_1(:,ii) = median( (N-1)*input_acq_time_map.ksp_norm_1(:,ii:Nint:end), 2 );
    k_2(:,ii) = median( (N-1)*input_acq_time_map.ksp_norm_2(:,ii:Nint:end), 2 ); 
    
end

for ii = 1:Nint
    figure(1); clf;
    plot( k_1(:,ii), k_2(:,ii) );
    title('Default spiral')
    drawnow
    
end

% setup NUFFTs

params_nufft.N = N;
params_nufft.Nrd = Nrd;
params_nufft.Nk = Nrd*Nint;
params_nufft.Nc = 1;
params_nufft.Nint = Nint;
FT = Gmri( [k_1(:) k_2(:)], true(N) );
params_nufft.FT = FT;


ksp_stack_nufft = complex( zeros( Nrd, n_TR ) );
ksp_stack_nufft_full = complex( zeros( Nrd*Nint, n_TR ) );
for ii = 1:n_TR
    my_interleaf = mod( ii - 1, Nint ) + 1;
    my_idx_v = ( my_interleaf - 1 ) * Nrd + 1:my_interleaf*Nrd;
    my_img = img_stack(:,:,ii);
    ksp_v = FT * my_img(:);
    ksp_stack_nufft(:,ii) = ksp_v(my_idx_v);
    ksp_stack_nufft_full(:,ii) = ksp_v;
%     figure(1); clf;
%     plot( abs( ksp_v ) );
%     pause();
end

%% add noise in k-space

ksp_stack_nufft = ksp_stack_nufft + sqrt( Nrd * Nint / n_TR ) * sigma * ( randn( size(ksp_stack_nufft) ) + 1i * randn( size(ksp_stack_nufft) ) );
ksp_stack_nufft_full = ksp_stack_nufft_full + sqrt( Nrd * Nint / n_TR ) * sigma * ( randn( size(ksp_stack_nufft_full) ) + 1i * randn( size(ksp_stack_nufft_full) ) );

%% for testing

% img_stack_nufft_full = complex( zeros( N, N, n_TR ) );
% parfor ii = 1:n_TR
%     [xs, info] = qpwls_pcg( zeros( N^2,1 ), FT, 1, ksp_stack_nufft_full(:,ii), 0, 0, 1, 20, 1);
%     my_img = reshape( xs(:,end), [N N] );
%     img_stack_nufft_full(:,:,ii) = my_img;
% end

%% package for processing by MRF_fat_sep_B0


output_MRF_raw.params.TE_s = TE_base_ms/1000;
output_MRF_raw.params.nomFlip_deg = nom_flip_ang_deg;
output_MRF_raw.params.TRbase_s = TR_base_ms/1000;
output_MRF_raw.params.tacq_s = tacq_ms/1000;
output_MRF_raw.data = reshape( permute( ksp_stack_nufft, [2 1] ), [1 n_TR Nrd] );
output_MRF_raw.data_full = reshape( permute( ksp_stack_nufft_full, [2 1] ),[1 n_TR Nrd*Nint] );

save(fn_out_raw,'output_MRF_raw');

T1_true_ms = zeros( N );
T2_true_ms = T1_true_ms;
FSF_true = T1_true_ms;
for ii = 1:n_segs
    
    T1_true_ms = T1_true_ms + mask(:,:,ii).*T1_water_ms_v(ii);
    T2_true_ms = T2_true_ms + mask(:,:,ii).*T2_water_ms_v(ii);
    FSF_true = FSF_true + mask(:,:,ii).*fsf_v(ii);
    
end

output_img_recon.MRF_img_stack_coil_combined = img_stack;
output_img_recon.B1_map = ones( N );

output_img_recon.T1_true_ms = T1_true_ms;
output_img_recon.T2_true_ms = T2_true_ms;
output_img_recon.FSF_true = FSF_true;
output_img_recon.B0_true_Hz = B0_true_Hz;
output_img_recon.mask = mask;

save(fn_out_img_cart,'output_img_recon');

save(fn_B1,'output_B1_data');
