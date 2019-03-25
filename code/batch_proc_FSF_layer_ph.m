%%

clear, close all, clc;

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_varTE.mat';
fn_MRF_raw = 'MRF_layer_varTE_1.mat';
fn_B1 = 'B1_layer_1.mat';
fn_MRF_water_dict = 'MRF_dict_layer_varTE.mat';
fn_MRF_seq_params = 'MRF103.csv';
fn_MRF_img_recon = 'MRF_layer_varTE_1_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_varTE_1_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = 'MRF_layer_varTE_1_proc_fat_sep.mat';

skiprecon = 0;
doMagdict = 1; % if 1, create dictionary
doTmaps = 0; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 1; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF

%%

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_varTE.mat';
fn_MRF_raw = 'MRF_layer_varTE_2.mat';
fn_B1 = 'B1_layer_2.mat';
fn_MRF_water_dict = 'MRF_dict_layer_varTE.mat';
fn_MRF_seq_params = 'MRF103.csv';
fn_MRF_img_recon = 'MRF_layer_varTE_2_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_varTE_2_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = 'MRF_layer_varTE_2_proc_fat_sep.mat';

skiprecon = 0;
doMagdict = 0; % if 1, create dictionary
doTmaps = 0; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 1; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF

%%

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_varTE.mat';
fn_MRF_raw = 'MRF_layer_varTE_3.mat';
fn_B1 = 'B1_layer_3.mat';
fn_MRF_water_dict = 'MRF_dict_layer_varTE.mat';
fn_MRF_seq_params = 'MRF103.csv';
fn_MRF_img_recon = 'MRF_layer_varTE_3_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_varTE_3_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = 'MRF_layer_varTE_3_proc_fat_sep.mat';

skiprecon = 0;
doMagdict = 0; % if 1, create dictionary
doTmaps = 0; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 1; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF

%%

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_varTE.mat';
fn_MRF_raw = 'MRF_layer_varTE_4.mat';
fn_B1 = 'B1_layer_4.mat';
fn_MRF_water_dict = 'MRF_dict_layer_varTE.mat';
fn_MRF_seq_params = 'MRF103.csv';
fn_MRF_img_recon = 'MRF_layer_varTE_4_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_varTE_4_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = 'MRF_layer_varTE_4_proc_fat_sep.mat';

skiprecon = 0;
doMagdict = 0; % if 1, create dictionary
doTmaps = 0; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 1; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF

%%

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_varTE.mat';
fn_MRF_raw = 'MRF_layer_varTE_5.mat';
fn_B1 = 'B1_layer_5.mat';
fn_MRF_water_dict = 'MRF_dict_layer_varTE.mat';
fn_MRF_seq_params = 'MRF103.csv';
fn_MRF_img_recon = 'MRF_layer_varTE_5_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_varTE_5_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = 'MRF_layer_varTE_5_proc_fat_sep.mat';

skiprecon = 0;
doMagdict = 0; % if 1, create dictionary
doTmaps = 0; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 1; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF

%%

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_varTE.mat';
fn_MRF_raw = 'MRF_layer_varTE_6.mat';
fn_B1 = 'B1_layer_6.mat';
fn_MRF_water_dict = 'MRF_dict_layer_varTE.mat';
fn_MRF_seq_params = 'MRF103.csv';
fn_MRF_img_recon = 'MRF_layer_varTE_6_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_varTE_6_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = 'MRF_layer_varTE_6_proc_fat_sep.mat';

skiprecon = 0;
doMagdict = 0; % if 1, create dictionary
doTmaps = 0; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 1; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF

%%

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_varTE.mat';
fn_MRF_raw = 'MRF_layer_varTE_7.mat';
fn_B1 = 'B1_layer_7.mat';
fn_MRF_water_dict = 'MRF_dict_layer_varTE.mat';
fn_MRF_seq_params = 'MRF103.csv';
fn_MRF_img_recon = 'MRF_layer_varTE_7_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_varTE_7_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = 'MRF_layer_varTE_7_proc_fat_sep.mat';

skiprecon = 0;
doMagdict = 0; % if 1, create dictionary
doTmaps = 0; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 1; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF


%%

clear, close all, clc;

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_varTR.mat';
fn_MRF_raw = 'MRF_layer_varTR_1.mat';
fn_B1 = 'B1_layer_1.mat';
fn_MRF_water_dict = 'MRF_dict_layer_varTR.mat';
fn_MRF_seq_params = 'MRF001.csv';
fn_MRF_img_recon = 'MRF_layer_varTR_1_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_varTR_1_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = '';

skiprecon = 0;
doMagdict = 1; % if 1, create dictionary
doTmaps = 1; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 0; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF

%%

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_varTR.mat';
fn_MRF_raw = 'MRF_layer_varTR_2.mat';
fn_B1 = 'B1_layer_2.mat';
fn_MRF_water_dict = 'MRF_dict_layer_varTR.mat';
fn_MRF_seq_params = 'MRF001.csv';
fn_MRF_img_recon = 'MRF_layer_varTR_2_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_varTR_2_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = '';

skiprecon = 0;
doMagdict = 0; % if 1, create dictionary
doTmaps = 1; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 0; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF

%%

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_varTR.mat';
fn_MRF_raw = 'MRF_layer_varTR_3.mat';
fn_B1 = 'B1_layer_3.mat';
fn_MRF_water_dict = 'MRF_dict_layer_varTR.mat';
fn_MRF_seq_params = 'MRF001.csv';
fn_MRF_img_recon = 'MRF_layer_varTR_3_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_varTR_3_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = '';

skiprecon = 0;
doMagdict = 0; % if 1, create dictionary
doTmaps = 1; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 0; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF

%%

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_varTR.mat';
fn_MRF_raw = 'MRF_layer_varTR_4.mat';
fn_B1 = 'B1_layer_4.mat';
fn_MRF_water_dict = 'MRF_dict_layer_varTR.mat';
fn_MRF_seq_params = 'MRF001.csv';
fn_MRF_img_recon = 'MRF_layer_varTR_4_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_varTR_4_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = '';

skiprecon = 0;
doMagdict = 0; % if 1, create dictionary
doTmaps = 1; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 0; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF

%%

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_varTR.mat';
fn_MRF_raw = 'MRF_layer_varTR_5.mat';
fn_B1 = 'B1_layer_5.mat';
fn_MRF_water_dict = 'MRF_dict_layer_varTR.mat';
fn_MRF_seq_params = 'MRF001.csv';
fn_MRF_img_recon = 'MRF_layer_varTR_5_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_varTR_5_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = '';

skiprecon = 0;
doMagdict = 0; % if 1, create dictionary
doTmaps = 1; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 0; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF

%%

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_varTR.mat';
fn_MRF_raw = 'MRF_layer_varTR_6.mat';
fn_B1 = 'B1_layer_6.mat';
fn_MRF_water_dict = 'MRF_dict_layer_varTR.mat';
fn_MRF_seq_params = 'MRF001.csv';
fn_MRF_img_recon = 'MRF_layer_varTR_6_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_varTR_6_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = '';

skiprecon = 0;
doMagdict = 0; % if 1, create dictionary
doTmaps = 1; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 0; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF

%%

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_varTR.mat';
fn_MRF_raw = 'MRF_layer_varTR_7.mat';
fn_B1 = 'B1_layer_7.mat';
fn_MRF_water_dict = 'MRF_dict_layer_varTR.mat';
fn_MRF_seq_params = 'MRF001.csv';
fn_MRF_img_recon = 'MRF_layer_varTR_7_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_varTR_7_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = '';

skiprecon = 0;
doMagdict = 0; % if 1, create dictionary
doTmaps = 1; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 0; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF

%%

clear, close all, clc;

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_fixTE.mat';
fn_MRF_raw = 'MRF_layer_fixTE_1.mat';
fn_B1 = 'B1_layer_1.mat';
fn_MRF_water_dict = 'MRF_dict_layer_fixTE.mat';
fn_MRF_seq_params = 'MRF104.csv';
fn_MRF_img_recon = 'MRF_layer_fixTE_1_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_fixTE_1_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = '';

skiprecon = 0;
doMagdict = 1; % if 1, create dictionary
doTmaps = 1; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 0; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF

%%

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_fixTE.mat';
fn_MRF_raw = 'MRF_layer_fixTE_2.mat';
fn_B1 = 'B1_layer_2.mat';
fn_MRF_water_dict = 'MRF_dict_layer_fixTE.mat';
fn_MRF_seq_params = 'MRF104.csv';
fn_MRF_img_recon = 'MRF_layer_fixTE_2_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_fixTE_2_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = '';

skiprecon = 0;
doMagdict = 0; % if 1, create dictionary
doTmaps = 1; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 0; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF

%%

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_fixTE.mat';
fn_MRF_raw = 'MRF_layer_fixTE_3.mat';
fn_B1 = 'B1_layer_3.mat';
fn_MRF_water_dict = 'MRF_dict_layer_fixTE.mat';
fn_MRF_seq_params = 'MRF104.csv';
fn_MRF_img_recon = 'MRF_layer_fixTE_3_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_fixTE_3_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = '';

skiprecon = 0;
doMagdict = 0; % if 1, create dictionary
doTmaps = 1; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 0; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF

%%

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_fixTE.mat';
fn_MRF_raw = 'MRF_layer_fixTE_4.mat';
fn_B1 = 'B1_layer_4.mat';
fn_MRF_water_dict = 'MRF_dict_layer_fixTE.mat';
fn_MRF_seq_params = 'MRF104.csv';
fn_MRF_img_recon = 'MRF_layer_fixTE_4_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_fixTE_4_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = '';

skiprecon = 0;
doMagdict = 0; % if 1, create dictionary
doTmaps = 1; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 0; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF

%%

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_fixTE.mat';
fn_MRF_raw = 'MRF_layer_fixTE_5.mat';
fn_B1 = 'B1_layer_5.mat';
fn_MRF_water_dict = 'MRF_dict_layer_fixTE.mat';
fn_MRF_seq_params = 'MRF104.csv';
fn_MRF_img_recon = 'MRF_layer_fixTE_5_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_fixTE_5_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = '';

skiprecon = 0;
doMagdict = 0; % if 1, create dictionary
doTmaps = 1; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 0; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF

%%

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_fixTE.mat';
fn_MRF_raw = 'MRF_layer_fixTE_6.mat';
fn_B1 = 'B1_layer_6.mat';
fn_MRF_water_dict = 'MRF_dict_layer_fixTE.mat';
fn_MRF_seq_params = 'MRF104.csv';
fn_MRF_img_recon = 'MRF_layer_fixTE_6_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_fixTE_6_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = '';

skiprecon = 0;
doMagdict = 0; % if 1, create dictionary
doTmaps = 1; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 0; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF

%%

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_ksp = 'MRF_ksp_traj_cor_fixTE.mat';
fn_MRF_raw = 'MRF_layer_fixTE_7.mat';
fn_B1 = 'B1_layer_7.mat';
fn_MRF_water_dict = 'MRF_dict_layer_fixTE.mat';
fn_MRF_seq_params = 'MRF104.csv';
fn_MRF_img_recon = 'MRF_layer_fixTE_7_img_recon.mat';
fn_MRF_proc_no_fat_sep = 'MRF_layer_fixTE_7_proc_no_fat_sep.mat';
fn_MRF_proc_fat_sep = '';

skiprecon = 0;
doMagdict = 0; % if 1, create dictionary
doTmaps = 1; % if 1, do dictionary matching without fat sep/B0 fitting
doMRFAT = 0; % if 1, do MRF with fat sep/B0 fitting

% img recon
input_img_recon.effMtx = 240;
input_img_recon.permute_order = [3 2 1];
input_img_recon.flip1 = 0;
input_img_recon.cc_factor = 2; % coil compression reduction factor
input_img_recon.use_median_traj = 1;
input_img_recon.n_angles = 32;

% dictionary params
input_dict.TI = 40.0; % ms
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500];
input_dict.B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
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


MRF_master_proc_w_FSF


%%

% exit;

