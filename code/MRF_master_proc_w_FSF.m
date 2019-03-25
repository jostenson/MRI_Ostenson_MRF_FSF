%% process MRF data

%% dependencies

% should already be added in run_all.m


%% params

% WAT T1, T2 characterization from Hamilton et al., JMRI, 2011
input_fat_sep_params.p_rel_ref_v = (input_fat_sep_params.ppm_ref - [0.9 1.3 1.3 2.1 2.1 2.75 4.2 5.3]) * 1e-6; % parts
input_fat_sep_params.fat_amps_v = [0.144 0.5*1.0 0.5*1.0 0.5*0.241 0.5*0.241 0.033 0.064 0.122];
input_fat_sep_params.T1_fat_s_v = [543 280 240 249 202 284 154 421]/1000;
input_fat_sep_params.T2_fat_s_v = [80.1 54.7 54.7 51.9 51.9 46.2 50 44.1]/1000; % 7th entry not supplied, 50 ms is estimate

%% load k-space trajectories

load([dir_in fn_ksp])

%% load MRF raw data

load([dir_in fn_MRF_raw])

if exist('skiprecon','var') && skiprecon == 0
    %% reconstruct and coil combine raw data

    input_img_recon.img_data = output_MRF_raw.data;
    input_img_recon.output_ksp_traj = output_ksp_traj;
    input_img_recon.ecalib_threshold = 0.001; % sensitivity threshold

    output_img_recon = MRF_img_recon_bart( input_img_recon );

    %% load B1 map

    load([dir_in fn_B1])
    output_img_recon.B1_map = output_B1_data.B1_map;

    % output_img_recon.B1_map = ones(input_img_recon.effMtx); % test

    %% save image recon

    save([dir_out fn_MRF_img_recon],'output_img_recon')
    
else
    load([dir_out fn_MRF_img_recon])
end


%% make or load dictionary

TE_s = output_MRF_raw.params.TE_s;% s
nomFlip_deg = output_MRF_raw.params.nomFlip_deg; % nominal flip angle in degrees
TRbase_s = output_MRF_raw.params.TRbase_s; % s

if doMagdict == 1

    % load acquisition vectors
    data_struct = importdata([dir_in fn_MRF_seq_params]);

    % set dictionary construction parameters
    input_dict.FA_v = data_struct.data(:,1)*nomFlip_deg; % deg
    input_dict.phi_v = data_struct.data(:,2); % deg
    input_dict.TR_v = TRbase_s*1000 + data_struct.data(:,3); % ms
    input_dict.TE_v = TE_s*1000 + data_struct.data(:,4); % ms
    if ndims( output_MRF_raw.data ) == 2
        input_dict.nreps = size( output_MRF_raw.data, 1 );
    elseif ndims(  output_MRF_raw.data ) == 3
        input_dict.nreps = size( output_MRF_raw.data, 2 );
    end

    input_dict.delk = 1; % step between states equal to a full dephasing imparted by crusher gradient
    input_dict.szomega = 101; % number of factors of k to include in phase history

    % plot sequence parameters
    figure(1); clf;
    subplot(411)
    plot(input_dict.FA_v); ylabel('degrees'); title(['FA, TI is ' num2str(input_dict.TI) ' ms'])
    subplot(412)
    plot(input_dict.phi_v); ylabel('degrees'); title('phase')
    subplot(413)
    plot(input_dict.TR_v); ylabel('msec'); title('TR')
    subplot(414)
    plot(input_dict.TE_v); ylabel('msec'); title('TE')
    drawnow

    % do dictionary construction
    input_dict.reduce = 1;
    output_dict = MRF_dict_B1(input_dict);

    %% save result
    output_dict.fn = fn_MRF_water_dict;
    save([dir_out fn_MRF_water_dict],'output_dict')

else

    load([dir_out fn_MRF_water_dict]);

end

%% do a standard dictionary match if applicable

if doTmaps == 1

    reduce_flag = 1;

    output_MRF_match = MRF_dict_match_B1( output_img_recon, output_dict, reduce_flag );

    figure(1000)
    imagesc(output_MRF_match.T1_map)
    axis image
    colorbar()
    title('T1 map')

    figure(1001)
    imagesc(output_MRF_match.T2_map)
    axis image
    colorbar()
    title('T2 map')

    figure(1002)
    imagesc(abs(output_MRF_match.M0_map))
    axis image
    colorbar()
    title('Magnitude of M0')


    % save MRF match w/o fat sep

    save([dir_out fn_MRF_proc_no_fat_sep],'output_MRF_match')

end

%% do fat-water separated dictionary match if applicable

if doMRFAT == 1

    output_MRFAT = MRF_fat_sep_B0( output_MRF_raw, output_img_recon, output_ksp_traj, output_dict, input_fat_sep_params );

    figure(2000)
    imagesc(output_MRFAT.T1_water_map)
    axis image
    colorbar()
    title('T1 map')

    figure(2001)
    imagesc(output_MRFAT.T2_water_map)
    axis image
    colorbar()
    title('T2 map')

    figure(2002)
    imagesc(output_MRFAT.FSF_map)
    axis image
    colorbar()
    title('FSF map')

    figure(2003);
    imagesc(output_MRFAT.B0_fit_map)
    axis image; colormap(gray)
    colorbar()
    title('B0_map')


    % save MRF match w/ fat sep

    save([dir_out fn_MRF_proc_fat_sep],'output_MRFAT')

end
