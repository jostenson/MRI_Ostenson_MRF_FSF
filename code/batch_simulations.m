%% batch generate simulation images

clear, close all, clc;

%% varTE uni-B1, off resonance B0

fn_B1 = '../data_in/B1_240_sim0.mat';
fn_mrf_csv = '../data_in/MRF103.csv';
fn_out_img_cart = '../data_out/fatwater_simulation_img_cart_uniB1offB0_varTE.mat';
fn_out_raw = '../data_in/fatwater_simulation_rawparams_uniB1offB0_varTE.mat';

TE_base_ms = 3.5;

B0_pert = 1;
B0_Hz_v = [150 150];
B1_v = [1 1];

% k-space trajectory from measurement
fn_ksp_traj = '../data_in/MRF_ksp_traj_tra_varTE_brain.mat';
permute_order = [1 2 3];
tacq_ms = 5.0;

sigma = 0; % noiseless

% simulation parameters
N = 240;
n_TR = 1500;
segs_v =        [ 0.1 0.2 0.3  0.40 1.00 ]; %intensities for segmentation; actual M0s = 1
fsf_v =         [ 0.6   0 0.5  0.25 0.80 ]; % for the different image segments
T1_water_ms_v = [ 2250 1200 1600 500 800 ];
T2_water_ms_v = [  100   50  100  30  30 ];
i_sign = -1; % precession direction convention
gamma_bar_Hz_per_T = 42.577e6;
B0_Tesla = 3;

us_flag = 1; % undersampling flag
Nint = 32; % number of interleaves

n_segs = numel( T1_water_ms_v );

% MRF seq
nom_flip_ang_deg = 60;
TR_base_ms = 16.0;
seq_struct = importdata(fn_mrf_csv);
seq_data = seq_struct.data;
FA_deg_v = nom_flip_ang_deg*seq_data(:,1);
phi_v = seq_data(:,2);
TR_ms_v = TR_base_ms + seq_data(:,3);
TE_ms_v = TE_base_ms + seq_data(:,4);
TI_ms = 40;
delk = 1;
szomega = 101;

% WAT T1, T2 characterization from Hamilton et al., JMRI, 2011
ppm_H2O = 4.65;
p_rel_H2O_v = (ppm_H2O - [0.9 1.3 1.3 2.1 2.1 2.75 4.2 5.3]) * 1e-6; % parts
fat_freq_Hz_v = p_rel_H2O_v * gamma_bar_Hz_per_T * B0_Tesla;
fat_amps_v = [0.144 0.5*1.0 0.5*1.0 0.5*0.241 0.5*0.241 0.033 0.064 0.122];
fat_amps_v = fat_amps_v./sum(fat_amps_v);
T1_fat_s_v = [543 280 240 249 202 284 154 421]/1000;
T2_fat_s_v = [80.1 54.7 54.7 51.9 51.9 46.2 50 44.1]/1000; % 7th entry not supplied, 50 ms is estimate
n_fat_pks = numel(fat_freq_Hz_v);

% run simulation
fat_water_MRF_image_simulation_nufft;

%% fixTE uni-B1, off resonance B0

clear, close all, clc;

fn_B1 = '../data_in/B1_240_sim0.mat';
fn_mrf_csv = '../data_in/MRF104.csv';
fn_out_img_cart = '../data_out/fatwater_simulation_img_cart_uniB1offB0_fixTE.mat';
fn_out_raw = '../data_in/fatwater_simulation_rawparams_uniB1offB0_fixTE.mat';

TE_base_ms = 4.65;

B0_pert = 1;
B0_Hz_v = [150 150];
B1_v = [1 1];

% k-space trajectory from measurement
fn_ksp_traj = '../data_in/MRF_ksp_traj_tra_fixTE_brain.mat';
permute_order = [1 2 3];
tacq_ms = 5.0;

sigma = 0; % noiseless

% simulation parameters
N = 240;
n_TR = 1500;
segs_v =        [ 0.1 0.2 0.3  0.40 1.00 ]; %intensities for segmentation; actual M0s = 1
fsf_v =         [ 0.6   0 0.5  0.25 0.80 ]; % for the different image segments
T1_water_ms_v = [ 2250 1200 1600 500 800 ];
T2_water_ms_v = [  100   50  100  30  30 ];
i_sign = -1; % precession direction convention
gamma_bar_Hz_per_T = 42.577e6;
B0_Tesla = 3;

us_flag = 1; % undersampling flag
Nint = 32; % number of interleaves

n_segs = numel( T1_water_ms_v );

% MRF seq
nom_flip_ang_deg = 60;
TR_base_ms = 16.0;
seq_struct = importdata(fn_mrf_csv);
seq_data = seq_struct.data;
FA_deg_v = nom_flip_ang_deg*seq_data(:,1);
phi_v = seq_data(:,2);
TR_ms_v = TR_base_ms + seq_data(:,3);
TE_ms_v = TE_base_ms + seq_data(:,4);
TI_ms = 40;
delk = 1;
szomega = 101;

% WAT T1, T2 characterization from Hamilton et al., JMRI, 2011
ppm_H2O = 4.65;
p_rel_H2O_v = (ppm_H2O - [0.9 1.3 1.3 2.1 2.1 2.75 4.2 5.3]) * 1e-6; % parts
fat_freq_Hz_v = p_rel_H2O_v * gamma_bar_Hz_per_T * B0_Tesla;
fat_amps_v = [0.144 0.5*1.0 0.5*1.0 0.5*0.241 0.5*0.241 0.033 0.064 0.122];
fat_amps_v = fat_amps_v./sum(fat_amps_v);
T1_fat_s_v = [543 280 240 249 202 284 154 421]/1000;
T2_fat_s_v = [80.1 54.7 54.7 51.9 51.9 46.2 50 44.1]/1000; % 7th entry not supplied, 50 ms is estimate
n_fat_pks = numel(fat_freq_Hz_v);

% run simulation
fat_water_MRF_image_simulation_nufft;

%% varTR uni-B1, off resonance B0

clear, close all, clc;

fn_B1 = '../data_in/B1_240_sim0.mat';
fn_mrf_csv = '../data_in/MRF001.csv';
fn_out_img_cart = '../data_out/fatwater_simulation_img_cart_uniB1offB0_varTR.mat';
fn_out_raw = '../data_in/fatwater_simulation_rawparams_uniB1offB0_varTR.mat';

TE_base_ms = 3.5;

B0_pert = 1;
B0_Hz_v = [150 150];
B1_v = [1 1];

% k-space trajectory from measurement
fn_ksp_traj = '../data_in/MRF_ksp_traj_tra_varTR_brain.mat';
permute_order = [1 2 3];
tacq_ms = 5.0;

sigma = 0; % noiseless

% simulation parameters
N = 240;
n_TR = 1000;
segs_v =        [ 0.1 0.2 0.3  0.40 1.00 ]; %intensities for segmentation; actual M0s = 1
fsf_v =         [ 0.6   0 0.5  0.25 0.80 ]; % for the different image segments
T1_water_ms_v = [ 2250 1200 1600 500 800 ];
T2_water_ms_v = [  100   50  100  30  30 ];
i_sign = -1; % precession direction convention
gamma_bar_Hz_per_T = 42.577e6;
B0_Tesla = 3;

us_flag = 1; % undersampling flag
Nint = 32; % number of interleaves

n_segs = numel( T1_water_ms_v );

% MRF seq
nom_flip_ang_deg = 60;
TR_base_ms = 16.0;
seq_struct = importdata(fn_mrf_csv);
seq_data = seq_struct.data;
FA_deg_v = nom_flip_ang_deg*seq_data(:,1);
phi_v = seq_data(:,2);
TR_ms_v = TR_base_ms + seq_data(:,3);
TE_ms_v = TE_base_ms + seq_data(:,4);
TI_ms = 40;
delk = 1;
szomega = 101;

% WAT T1, T2 characterization from Hamilton et al., JMRI, 2011
ppm_H2O = 4.65;
p_rel_H2O_v = (ppm_H2O - [0.9 1.3 1.3 2.1 2.1 2.75 4.2 5.3]) * 1e-6; % parts
fat_freq_Hz_v = p_rel_H2O_v * gamma_bar_Hz_per_T * B0_Tesla;
fat_amps_v = [0.144 0.5*1.0 0.5*1.0 0.5*0.241 0.5*0.241 0.033 0.064 0.122];
fat_amps_v = fat_amps_v./sum(fat_amps_v);
T1_fat_s_v = [543 280 240 249 202 284 154 421]/1000;
T2_fat_s_v = [80.1 54.7 54.7 51.9 51.9 46.2 50 44.1]/1000; % 7th entry not supplied, 50 ms is estimate
n_fat_pks = numel(fat_freq_Hz_v);

% run simulation
fat_water_MRF_image_simulation_nufft_varTR;
