%% analyze fat-water separation using direct dictionary match k-space based method

clear, close all, clc;


%%

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_conv = 'emulsion_conventional.mat';
fn_conv_maps = 'emulsion_conventional_maps.mat';
fn_direct = 'MRF_direct_FSF_emulsion.mat';
fn_ksp = 'MRF_emulsion_varTE_proc_fat_sep.mat';
fn_ROI = 'emulsion_conventional_ROIs.mat';

fn_save = 'MRF_directksp_analysis.mat';

create_ROIs = 0;
n_ROI = 7;
do_T1 = 0;
do_T2 = 0;
n_T2_mod_params = 3;
do_fatsep = 1;

%% get MRF data

load([dir_out fn_direct])

FSF_direct = output_MRF_FSF_direct_match.FSF_map;

load([dir_out fn_ksp])

FSF_ksp = output_MRFAT.FSF_map;


%% get ROIs

copyfile([dir_in fn_ROI],[dir_out fn_ROI])
load( [dir_out fn_ROI] )
 
%% graph cut fat-water separation

conventional_T1T2fatsep

FSF_graphcut = conventional_maps.FSF_graphcut;

%% get stats

FSF_directksp_means = zeros( n_ROI, 3 );
FSF_directksp_std = zeros( n_ROI, 3 );

for ii = 1:n_ROI
    
    my_ROI = logical( ROI_stack(:,:,ii) );
    
    FSF_directksp_means(ii,1) = mean( FSF_graphcut( my_ROI ) );
    FSF_directksp_std(ii,1) = std( FSF_graphcut( my_ROI ) );
    FSF_directksp_means(ii,2) = mean( FSF_direct( my_ROI ) );
    FSF_directksp_std(ii,2) = std( FSF_direct( my_ROI ) );
    FSF_directksp_means(ii,3) = mean( FSF_ksp( my_ROI ) );
    FSF_directksp_std(ii,3) = std( FSF_ksp( my_ROI ) );
    
end

[CCC_direct, CCC_CI_direct, ~, ~]  = ccc( FSF_directksp_means(:,1), FSF_directksp_means(:,2) );
[CCC_ksp, CCC_CI_ksp, ~, ~]  = ccc( FSF_directksp_means(:,1), FSF_directksp_means(:,3) );

%%

save([dir_out fn_save],'FSF_graphcut','FSF_direct','FSF_ksp','CCC_direct','CCC_CI_direct',...
    'CCC_ksp','CCC_CI_ksp','img_conventional', 'FSF_directksp_means', 'FSF_directksp_std')
