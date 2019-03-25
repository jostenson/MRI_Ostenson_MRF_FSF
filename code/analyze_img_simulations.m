%% analyze results of fat-water image simulations

clear, close all, clc;

%%

dir_in = '../data_out/';
dir_out = '../data_out/';

N = 240; % image size

fn_save = 'img_simulation_proc_stats.mat';

fn_simproc{1} = 'MRF_fatwatersep_simulation_MC_uniB1offB0_varTE.mat';
fn_simproc{2} = 'MRF_fatwater_simulation_MC_uniB1offB0_fixTE.mat';
fn_simproc{3} = 'MRF_fatwater_simulation_MC_uniB1offB0_varTR.mat';

fn_key = {'varTE_1B1_1B0','fixTE_1B1_1B0','varTR_1B1_1B0'};

n_metrics = 4; %T1,T2,(FSF),(B0)
nm_metrics = {'T1_water_maps','T2_water_maps','FSF_maps','B0_maps'};
nm_metrics_true = {'T1_map_true_ms','T2_map_true_ms','FSF_map_true','B0_map_true_Hz'};

eg_snr_pick = 2; % for example maps

%%

n_sims = numel( fn_simproc );

results = struct();
for ii = 1:n_sims
    load( [dir_in fn_simproc{ii}] );
    eval( sprintf( 'results.%s = output_MC;', fn_key{ii} ) );
    clear output_MC;
end

%%

n_noise = numel( eval( sprintf( 'results.%s.SNR_dB_v',fn_key{1} ) ) );
n_ROI = size( eval( sprintf( 'results.%s.mask_segs',fn_key{1} ) ), 3 );

simulation_means = zeros( n_ROI+1, n_metrics, n_noise, n_sims );
simulation_bias_means = simulation_means;
simulation_bias_std = simulation_means;
example_maps = zeros( N, N, n_metrics, n_sims, 2 );
example_snr_dB = eval( sprintf( 'results.%s.SNR_dB_v(%d);',fn_key{1},eg_snr_pick ) );

for ii = 1:n_sims
   
    eval( sprintf( 'my_result = results.%s;', fn_key{ii} ) );
    
    for jj = 1:n_noise
        
        for kk = 1:n_metrics
            
            eval( sprintf( 'my_maps = my_result.%s;', nm_metrics{kk} ) );
            eval( sprintf( 'my_map_true = my_result.%s;', nm_metrics_true{kk} ) );
            
            if strcmp( nm_metrics(kk),'B0_maps' ) % fix B0 sign convention
                my_maps = -my_maps;
            end
            
            example_maps(:,:,kk,ii,1) = my_maps(:,:,eg_snr_pick).*double( logical( sum(my_result.mask_segs,3) ) );
            example_maps(:,:,kk,ii,2) = my_map_true.*double( logical( sum(my_result.mask_segs,3) ) );
            
            for ll = 1:n_ROI+1
                
                if ll <= n_ROI
                    my_mask = my_result.mask_segs(:,:,ll);
                else
                    my_mask = logical( sum( my_result.mask_segs, 3 ) );
                end
                
                my_map = my_maps(:,:,jj);
                
                simulation_means(ll,kk,jj,ii) = mean( my_map(my_mask) );
                simulation_bias_means(ll,kk,jj,ii) = mean( my_map(my_mask) - my_map_true(my_mask) );
                simulation_bias_std(ll,kk,jj,ii) = std( my_map(my_mask) - my_map_true(my_mask) );
                
            end
            
        end
        
    end
    
end

%%

save( [dir_out fn_save], 'simulation_means','simulation_bias_means','simulation_bias_std','example_maps','example_snr_dB','fn_key','nm_metrics' );