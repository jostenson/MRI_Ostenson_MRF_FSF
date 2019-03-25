%%

clear, close all, clc;

%% data files

dir_in = '../data_out/';

fn_liver{1} = 'MRF_liver_varTE_proc_fat_sep.mat';
fn_liver_roi = 'liver_ROIs_tissue.mat';


%%

copyfile(['../data_in/' fn_liver_roi],[dir_in fn_liver_roi])
load( [dir_in fn_liver_roi] );
roi_liver = roi_stack;
n_roi_liver = size( roi_liver, 3 );
n_sets_liver = numel( fn_liver );

%%

T1_liver_means = zeros( n_roi_liver, n_sets_liver );
T1_liver_std = T1_liver_means;
T2_liver_means = T1_liver_means;
T2_liver_std = T2_liver_means;

for ii = 1:n_sets_liver
    
    load( [dir_in fn_liver{ii}] );
    
    try
        T1_liver = output_MRF_match.T1_map;
        T2_liver = output_MRF_match.T2_map;
        clear('output_MRF_match');
    catch
        T1_liver = output_MRFAT.T1_water_map;
        T2_liver = output_MRFAT.T2_water_map;
        clear('output_MRFAT');
    end
    
    figure(1); clf;
    my_fn = fn_liver{ii};
    subplot(121);
    imagesc( T1_liver ); axis image; colorbar()
    title( sprintf( '%s : T1 map', my_fn ) );
    subplot(122);
    imagesc( T2_liver ); axis image; colorbar()
    title('T2 map')
    drawnow()
    
    for jj = 1:n_roi_liver
        my_liver_roi = logical( roi_liver(:,:,jj) );
        T1_liver_means(jj,ii) = mean( T1_liver( my_liver_roi & T1_liver~=0 ) );
        T1_liver_std(jj,ii) = std( T1_liver( my_liver_roi & T1_liver~=0 ) );
        T2_liver_means(jj,ii) = mean( T2_liver( my_liver_roi & T1_liver~=0 ) );
        T2_liver_std(jj,ii) = std( T2_liver( my_liver_roi & T1_liver~=0 ) );
    end
    
end
