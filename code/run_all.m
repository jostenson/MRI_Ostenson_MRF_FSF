%% top level batch script

% set(groot,'DefaultFigureVisible','off'); % option to turn off plot figures
clear, close all, clc;

% point to correct bart path
bart_path = '';
addpath([bart_path 'matlab']);
setenv('TOOLBOX_PATH', bart_path);
addpath('./utils/')
addpath('../contrib/sdc3_nrz_11aug/')
addpath('../contrib/grid3_dct_11aug/')
addpath('../contrib/irt'); % Michigan Image Reconstruction Toolbox
my_what = what('../contrib/irt');
irtdir = my_what.path;
setup_irt


%% process MRF data

T1T2_bias_with_fat
T1T2_bias_with_fat_varTR
T1T2_bias_with_fat_varTE_proposed
batch_proc_FSF_emulsion_ph
MRF_FSF_direct_dict_match
batch_proc_FSF_layer_ph
batch_proc_FSF_layer_ph_badshim
batch_proc_FSF_knee 
batch_proc_FSF_liver
batch_proc_FSF_brain
batch_proc_Msys
batch_simulations
batch_proc_simulations
batch_proc_simulations_cart

%% process conventional data and analyze results

analyze_direct_and_ksp_fat_sep %
analyze_oil_water_layer %
analyze_invivo_conventional %
analyze_Msys_verification_T1T2 %
analyze_img_simulations
analyze_img_simulations_cart

%% generate figures

% figure_generation % layout dependent on display settings
analyze_invivo_roi
