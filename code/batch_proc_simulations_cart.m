%% batch process simulations


%% varTE uni-B1, off resonance B0

clear, close all, clc;

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_tra_varTE_brain.mat';
fn_MRF_raw = 'fatwater_simulation_rawparams_noisey.mat';
fn_MRF_raw_nonoise = 'fatwater_simulation_rawparams_uniB1offB0_varTE.mat';
fn_B1 = 'B1_240_sim0.mat';
fn_MRF_water_dict = 'MRF_dict_varTE.mat';
fn_MRF_seq_params = 'MRF103.csv';
fn_MRF_img_recon = 'fatwater_simulation_img_cart_uniB1offB0_varTE.mat';
fn_MRF_proc_no_fat_sep = 'MRF_simulation_varTE_proc_no_fat_sep_cart.mat';
fn_MRF_proc_fat_sep = 'MRF_simulation_varTE_proc_fat_sep_cart.mat';
fn_MC_output = 'MRF_fatwatersep_simulation_cart_MC_uniB1offB0_varTE.mat';
fn_out_img_cart = 'fatwater_simulation_img_cart_uniB1offB0_varTE.mat';

SNR_dB_v = [inf];
sigma_v = 10.^(-SNR_dB_v/20); %
n_iter = numel(SNR_dB_v);

skiprecon = 1;
doMagdict = 1; % if 1, create dictionary
doTmaps = 0; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 1; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [1 2 3];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
% input_dict.B1_v = 1;
input_dict.sfrac = 1 - 1e-4; % indirectly applied to fat sep algo

% MRF fat sep params
input_fat_sep_params.i_sign = -1;
input_fat_sep_params.n_B0_basis = 31;
input_fat_sep_params.freq_MFI_max_Hz = 700;
input_fat_sep_params.freq_MFI_step_Hz = 10;
input_fat_sep_params.s_frac = 1 - 1e-4; % fraction of dictionary energy to keep
input_fat_sep_params.B0_range_Hz = [-250 250];
input_fat_sep_params.gamma_bar_Hz_per_T = 42.5775e6;
input_fat_sep_params.ppm_ref = 4.65; % H2O
input_fat_sep_params.B0_Tesla = 3.0;
input_fat_sep_params.permute_order = input_img_recon.permute_order;
input_fat_sep_params.sigma_filter = 1.5;

% run processing simulation
proc_simulation


%% fixTE uni-B1, off resonance B0

clear, close all, clc;

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_tra_fixTE_brain.mat';
fn_MRF_raw = 'fatwater_simulation_rawparams_noisey.mat';
fn_MRF_raw_nonoise = 'fatwater_simulation_rawparams_uniB1offB0_fixTE.mat';
fn_B1 = 'B1_240_sim0.mat';
fn_MRF_water_dict = 'MRF_dict_fixTE.mat';
fn_MRF_seq_params = 'MRF104.csv';
fn_MRF_img_recon = 'fatwater_simulation_img_cart_uniB1offB0_fixTE.mat';
fn_MRF_proc_no_fat_sep = 'MRF_simulation_fixTE_proc_no_fat_sep_cart.mat';
fn_MC_output = 'MRF_fatwater_simulation_cart_MC_uniB1offB0_fixTE.mat';
fn_out_img_cart = 'fatwater_simulation_img_cart_uniB1offB0_fixTE.mat';

SNR_dB_v = [inf];
sigma_v = 10.^(-SNR_dB_v/20); %
n_iter = numel(SNR_dB_v);

skiprecon = 1;
doMagdict = 1; % if 1, create dictionary
doTmaps = 1; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 0; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [1 2 3];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
% input_dict.B1_v = 1;
input_dict.sfrac = 1 - 1e-4; % indirectly applied to fat sep algo

% MRF fat sep params
input_fat_sep_params.i_sign = -1;
input_fat_sep_params.n_B0_basis = 31;
input_fat_sep_params.freq_MFI_max_Hz = 700;
input_fat_sep_params.freq_MFI_step_Hz = 10;
input_fat_sep_params.s_frac = 1 - 1e-4; % fraction of dictionary energy to keep
input_fat_sep_params.B0_range_Hz = [-250 250];
input_fat_sep_params.gamma_bar_Hz_per_T = 42.5775e6;
input_fat_sep_params.ppm_ref = 4.65; % H2O
input_fat_sep_params.B0_Tesla = 3.0;
input_fat_sep_params.permute_order = input_img_recon.permute_order;
input_fat_sep_params.sigma_filter = 1.5;

% run processing simulation
proc_simulation

%% varTR uni-B1, off resonance B0

clear, close all, clc;

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_tra_varTR_brain.mat';
fn_MRF_raw = 'fatwater_simulation_rawparams_noisey.mat';
fn_MRF_raw_nonoise = 'fatwater_simulation_rawparams_uniB1offB0_varTR.mat';
fn_B1 = 'B1_240_sim0.mat';
fn_MRF_water_dict = 'MRF_dict_varTR.mat';
fn_MRF_seq_params = 'MRF001.csv';
fn_MRF_img_recon = 'fatwater_simulation_img_cart_uniB1offB0_varTR.mat';
fn_MRF_proc_no_fat_sep = 'MRF_simulation_varTR_proc_no_fat_sep_cart.mat';
fn_MC_output = 'MRF_fatwater_simulation_cart_MC_uniB1offB0_varTR.mat';
fn_out_img_cart = 'fatwater_simulation_img_cart_uniB1offB0_varTR.mat';

SNR_dB_v = [inf];
sigma_v = 10.^(-SNR_dB_v/20); %
n_iter = numel(SNR_dB_v);

skiprecon = 1;
doMagdict = 1; % if 1, create dictionary
doTmaps = 1; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 0; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [1 2 3];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
% input_dict.B1_v = 1;
input_dict.sfrac = 1 - 1e-4; % indirectly applied to fat sep algo

% MRF fat sep params
input_fat_sep_params.i_sign = -1;
input_fat_sep_params.n_B0_basis = 31;
input_fat_sep_params.freq_MFI_max_Hz = 700;
input_fat_sep_params.freq_MFI_step_Hz = 10;
input_fat_sep_params.s_frac = 1 - 1e-4; % fraction of dictionary energy to keep
input_fat_sep_params.B0_range_Hz = [-250 250];
input_fat_sep_params.gamma_bar_Hz_per_T = 42.5775e6;
input_fat_sep_params.ppm_ref = 4.65; % H2O
input_fat_sep_params.B0_Tesla = 3.0;
input_fat_sep_params.permute_order = input_img_recon.permute_order;
input_fat_sep_params.sigma_filter = 1.5;

% run processing simulation
proc_simulation

