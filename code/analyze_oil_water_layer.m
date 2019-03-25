%% analyze oil-water layer phantom experiment

clear, close all, clc;

%%

dir_in = '../data_in/';
dir_out = '../data_out/';

ls_MRF_varTE = dir([dir_out '*layer*varTE*proc_fat*.mat']);
ls_MRF_fixTE = dir([dir_out '*layer*fixTE*proc*.mat']);
ls_MRF_varTR = dir([dir_out '*layer*varTR*proc*.mat']);
fn_ROI = 'layer_conventional_ROIs.mat';
fn_conv = 'layer_conventional.mat';
fn_conv_maps = 'layer_conventional_maps.mat';

fn_save = 'layer_FSF_T1T2_analysis.mat';

create_ROIs = 0;
n_ROI = 9;
do_T1 = 0;
do_T2 = 0;
n_T2_mod_params = 3;
do_fatsep = 1;

%%

n_varTE = numel(ls_MRF_varTE);
n_fixTE = numel(ls_MRF_fixTE);
n_varTR = numel(ls_MRF_varTR);

for ii = 1:n_varTE
    
    load([dir_out ls_MRF_varTE(ii).name]);
    
    if ii == 1
        [n_r, n_c] = size( output_MRFAT.FSF_map );
        FSF_varTE_stack = zeros( n_r,n_c,n_varTE );
        T1_varTE_stack = FSF_varTE_stack;
        T2_varTE_stack = FSF_varTE_stack;
        B0_varTE_stack = FSF_varTE_stack;
    end
    
    FSF_varTE_stack(:,:,ii) = output_MRFAT.FSF_map;
    T1_varTE_stack(:,:,ii) = output_MRFAT.T1_water_map;
    T2_varTE_stack(:,:,ii) = output_MRFAT.T2_water_map;
    if isfield(output_MRFAT,'B0_fit_map')
        B0_varTE_stack(:,:,ii) = output_MRFAT.B0_fit_map;
    end
    
end

for ii = 1:n_fixTE
    
    load([dir_out ls_MRF_fixTE(ii).name]);
    
    if ii == 1
        
        T1_fixTE_stack = zeros( n_r,n_c,n_fixTE );
        T2_fixTE_stack = T1_fixTE_stack;
    end
    
    T1_fixTE_stack(:,:,ii) = output_MRF_match.T1_map;
    T2_fixTE_stack(:,:,ii) = output_MRF_match.T2_map;
    
end

for ii = 1:n_varTR
    
    load([dir_out ls_MRF_varTR(ii).name]);
    
    if ii == 1
        T1_varTR_stack = zeros( n_r,n_c,n_fixTE );
        T2_varTR_stack = T1_varTR_stack;
    end
    
    T1_varTR_stack(:,:,ii) = output_MRF_match.T1_map;
    T2_varTR_stack(:,:,ii) = output_MRF_match.T2_map;
    
end

%% get ROIs

copyfile([dir_in fn_ROI],[dir_out fn_ROI])
load( [dir_out fn_ROI] )

%% do fat-water sep

conventional_T1T2fatsep

%% get means for all metrics

T1_means_MRF = zeros( n_ROI, n_fixTE, 3 ); % n_ROI x n_layers x n_MRF
T1_std_MRF = T1_means_MRF;
T2_means_MRF = T1_means_MRF;
T2_std_MRF = T1_means_MRF;
FSF_means_MRF = zeros( n_ROI, n_fixTE );
FSF_std_MRF = FSF_means_MRF;
B0_means_MRF = FSF_means_MRF;
B0_std_MRF = FSF_means_MRF;
FSF_means_conv = FSF_means_MRF;
FSF_std_conv = FSF_means_MRF;
B0_means_conv = FSF_means_MRF;
B0_std_conv = FSF_means_MRF;
R2star_means_conv = FSF_means_MRF;
R2star_std_conv = FSF_means_MRF;

for ii = 1:n_ROI

    my_ROI = logical(ROI_stack(:,:,ii));
            
    for jj = 1:n_fixTE
    
        
        my_FSF = FSF_varTE_stack(:,:,jj);
        my_T1 = T1_varTE_stack(:,:,jj);
        my_T2 = T2_varTE_stack(:,:,jj);
        my_B0 = B0_varTE_stack(:,:,jj);
        
        FSF_means_MRF(ii,jj) = mean(my_FSF( my_ROI ));
        T1_means_MRF(ii,jj,1) = mean(my_T1( my_ROI ));
        T2_means_MRF(ii,jj,1) = mean(my_T2( my_ROI ));
        B0_means_MRF(ii,jj) = mean(my_B0( my_ROI ));
        
        FSF_std_MRF(ii,jj) = std(my_FSF( my_ROI ));
        T1_std_MRF(ii,jj,1) = std(my_T1( my_ROI ));
        T2_std_MRF(ii,jj,1) = std(my_T2( my_ROI ));
        B0_std_MRF(ii,jj) = std(my_B0( my_ROI ));
        
        my_T1 = T1_fixTE_stack(:,:,jj);
        my_T2 = T2_fixTE_stack(:,:,jj);
        
        T1_means_MRF(ii,jj,2) = mean(my_T1( my_ROI ));
        T2_means_MRF(ii,jj,2) = mean(my_T2( my_ROI ));
        
        T1_std_MRF(ii,jj,2) = std(my_T1( my_ROI ));
        T2_std_MRF(ii,jj,2) = std(my_T2( my_ROI ));
        
        my_T1 = T1_varTR_stack(:,:,jj);
        my_T2 = T2_varTR_stack(:,:,jj);
        
        T1_means_MRF(ii,jj,3) = mean(my_T1( my_ROI ));
        T2_means_MRF(ii,jj,3) = mean(my_T2( my_ROI ));
        
        T1_std_MRF(ii,jj,3) = std(my_T1( my_ROI ));
        T2_std_MRF(ii,jj,3) = std(my_T2( my_ROI ));
        
        my_FSF = conventional_maps.FSF_graphcut(:,:,jj);
        my_R2star = conventional_maps.R2star_graphcut(:,:,jj);
        my_B0 = conventional_maps.B0_graphcut(:,:,jj);
        
        FSF_means_conv(ii,jj) = mean(my_FSF( my_ROI ));
        R2star_means_conv(ii,jj) = mean(my_T2( my_ROI ));
        B0_means_conv(ii,jj) = mean(my_B0( my_ROI ));
        
        FSF_std_conv(ii,jj) = std(my_FSF( my_ROI ));
        R2star_std_conv(ii,jj) = std(my_T2( my_ROI ));
        B0_std_conv(ii,jj) = std(my_B0( my_ROI ));
    
    end
    
end


%% compute FSF concordance

[CCC_FSF, CCC_FSF_CI, ~, ~] = ccc( FSF_means_conv(:),FSF_means_MRF(:) );

%% compute consensus T1 and T2 for references

T1_consensus = mean( squeeze( T1_means_MRF(:,1,:) ), 2 );
T2_consensus = mean( squeeze( T2_means_MRF(:,1,:) ), 2 );


%% compute T1 and T2 bias for each ROI through fat fraction

T1_bias_MRF = T1_means_MRF - repmat( T1_consensus, [1 n_fixTE 3] );
T2_bias_MRF = T2_means_MRF - repmat( T2_consensus, [1 n_fixTE 3] );


%% save results

save([dir_out fn_save],'T1_means_MRF','T1_std_MRF','T2_means_MRF','T2_std_MRF', ...
    'FSF_means_MRF','FSF_std_MRF','FSF_means_conv','FSF_std_MRF','B0_means_MRF',...
    'B0_std_MRF','B0_means_conv','B0_std_conv','T1_consensus','T2_consensus',...
    'CCC_FSF','CCC_FSF_CI','T1_bias_MRF','T2_bias_MRF');