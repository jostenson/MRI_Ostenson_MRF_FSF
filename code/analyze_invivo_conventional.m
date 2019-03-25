%% process conventional in vivo data

clear, close all, clc;

dir_in = '../data_in/';
dir_out = '../data_out/';

%% knee

fn_conv = 'knee_conventional.mat';
fn_conv_maps = 'knee_conventional_maps.mat';
fn_ROI = 'knee_conventional_ROIs.mat';

create_ROIs = 0;
n_ROI = 1;
do_T1 = 1;
do_T2 = 1;
n_T2_mod_params = 3;
do_fatsep = 1;

%% load ROI map

copyfile([dir_in fn_ROI],[dir_out fn_ROI]);
load( [dir_out fn_ROI] );


%%

conventional_T1T2fatsep


%% brain

fn_conv = 'brain_conventional.mat';
fn_conv_maps = 'brain_conventional_maps.mat';
fn_ROI = 'brain_conventional_ROIs.mat';

create_ROIs = 0;
n_ROI = 1;
do_T1 = 0;
do_T2 = 0;
n_T2_mod_params = 3;
do_fatsep = 1;

%% load ROI map

copyfile([dir_in fn_ROI],[dir_out fn_ROI]);
load( [dir_out fn_ROI] );


%%

conventional_T1T2fatsep

