%% process simulation for specified iterations

output_MC.n_iter = n_iter;
output_MC.sigma_v = sigma_v;
output_MC.SNR_dB_v = SNR_dB_v;
output_MC.last_iter = 0;
output_MC.T1_water_maps = zeros( input_img_recon.effMtx, input_img_recon.effMtx, n_iter );
output_MC.T2_water_maps = output_MC.T1_water_maps;
output_MC.M0_water_maps = output_MC.T1_water_maps;
output_MC.M0_fat_maps = output_MC.T1_water_maps;
output_MC.FSF_maps = output_MC.T1_water_maps;
output_MC.B0_maps = output_MC.T1_water_maps;

output_img_recon_cart = load( [dir_out fn_out_img_cart] );
output_img_recon_cart = output_img_recon_cart.output_img_recon;
mask_segs = output_img_recon_cart.mask;
n_segs = size( mask_segs, 3);

output_MC.mask_segs = mask_segs;
output_MC.n_segs = n_segs;
output_MC.T1_map_true_ms = output_img_recon_cart.T1_true_ms;
output_MC.T2_map_true_ms = output_img_recon_cart.T2_true_ms;
output_MC.FSF_map_true = output_img_recon_cart.FSF_true;
output_MC.B0_map_true_Hz = output_img_recon_cart.B0_true_Hz;

save( [dir_out fn_MC_output], 'output_MC' );

for gg = 1:n_iter
    
    sigma = sigma_v(gg);
    
    if gg == 1
        doMagdict = 1;
    else
        doMagdict = 0;
    end
    
    load( [dir_in fn_MRF_raw_nonoise] );
    Nint = input_img_recon.n_angles;
    n_TR = size( output_MRF_raw.data, 2 );
    Nrd = size( output_MRF_raw.data, 3 );
    output_MRF_raw.data = output_MRF_raw.data + sqrt( Nrd * Nint / n_TR ) * sigma * ( randn( size(output_MRF_raw.data) ) + 1i * randn( size(output_MRF_raw.data) ) );
    
    save( [dir_in fn_MRF_raw], 'output_MRF_raw' );
    
    MRF_master_proc_w_FSF

    %%
    
    load( [dir_out fn_MC_output] );
    
    if doTmaps ~= 1
        output_MC.T1_water_maps(:,:,gg) = output_MRFAT.T1_water_map;
        output_MC.T2_water_maps(:,:,gg) = output_MRFAT.T2_water_map;
        output_MC.M0_water_maps(:,:,gg) = output_MRFAT.M0_water_map;
        output_MC.M0_fat_maps(:,:,gg) = output_MRFAT.M0_fat_map;
        output_MC.FSF_maps(:,:,gg) = output_MRFAT.FSF_map;
        output_MC.B0_maps(:,:,gg) = output_MRFAT.B0_fit_map;
        if gg == 1
            output_MC.B1_map_true = output_img_recon.B1_map;
        end
        
        %% display restuls
        
        figure(10); clf;
        subplot(121);
        imagesc( output_img_recon_cart.T1_true_ms ); axis image; colorbar(); caxis( [0 1500] );
        title('true T1')
        subplot(122);
        imagesc( output_MRFAT.T1_water_map ); axis image; colorbar(); caxis( [0 1500] );
        title('fit T1')
        
        figure(20); clf;
        subplot(121);
        imagesc( output_img_recon_cart.T2_true_ms ); axis image; colorbar(); caxis( [0 150] );
        title('true T2')
        subplot(122);
        imagesc( output_MRFAT.T2_water_map ); axis image; colorbar(); caxis( [0 150] );
        title('fit T2')
        
        figure(30); clf;
        subplot(121);
        imagesc( output_img_recon_cart.B0_true_Hz ); axis image; colorbar(); caxis( [-250 250] );
        title('true B0')
        subplot(122);
        imagesc( -output_MRFAT.B0_fit_map ); axis image; colorbar(); caxis( [-250 250] );
        title('fit B0')
        
        figure(40); clf;
        subplot(121);
        imagesc( output_img_recon_cart.FSF_true ); axis image; colorbar(); caxis( [0 1] );
        title('true FSF')
        subplot(122);
        imagesc( output_MRFAT.FSF_map ); axis image; colorbar(); caxis( [0 1] );
        title('fit FSF')
        drawnow();
    
    else
        
        output_MC.T1_water_maps(:,:,gg) = output_MRF_match.T1_map;
        output_MC.T2_water_maps(:,:,gg) = output_MRF_match.T2_map;
        output_MC.M0_water_maps(:,:,gg) = output_MRF_match.M0_map;
        if gg == 1
            output_MC.B1_map_true = output_img_recon.B1_map;
        end
        
        figure(10); clf;
        subplot(121);
        imagesc( output_img_recon_cart.T1_true_ms ); axis image; colorbar(); caxis( [0 1500] );
        title('true T1')
        subplot(122);
        imagesc( output_MRF_match.T1_map ); axis image; colorbar(); caxis( [0 1500] );
        title('fit T1')
        
        figure(20); clf;
        subplot(121);
        imagesc( output_img_recon_cart.T2_true_ms ); axis image; colorbar(); caxis( [0 150] );
        title('true T2')
        subplot(122);
        imagesc( output_MRF_match.T2_map ); axis image; colorbar(); caxis( [0 150] );
        title('fit T2')
        drawnow();
        
    end
    %%

    output_MC.last_iter = gg;
    
    save( [dir_out fn_MC_output], 'output_MC' );
    
end