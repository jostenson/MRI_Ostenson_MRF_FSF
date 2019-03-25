function [ output_args ] = MRF_img_recon_bart( input_args )
%% MRF recon with Bart


%% measure total execution time
start_time = tic;

%% params

effMtx = input_args.effMtx;
permute_order = input_args.permute_order;
try
    flip1 = input_args.flip1;
catch
end
try
    flip2 = input_args.flip2;
catch
end
cc_factor = input_args.cc_factor;
data = input_args.img_data;
output_ksp_traj = input_args.output_ksp_traj;
ecalib_threshold = input_args.ecalib_threshold;
try
    use_median_traj = input_args.use_median_traj;
    n_angles = input_args.n_angles;
catch
    use_median_traj = 0;
end

input_args = rmfield( input_args, 'img_data' );

%% determine MRF data format

try
    
    assert(doFullksp == 1)
    % single coil
    [n_dyn, n_ph, n_rd] = size(data);
    data_comb = permute(data,[3 1 2]);
    data_comb = reshape(data_comb, [n_rd*n_dyn n_ph]);
    n_Coil = 1;
    
catch
    
    
    if numel(size(data)) == 3 % multicoil
        [n_Coil,n_ph,n_rd] = size(data); %num channels, num acqs, num readout pnts
    elseif numel(size(data)) == 2 % single coil
        [n_ph,n_rd] = size(data);
        data = shiftdim( data, -1 );
        n_Coil = 1;
    else
        disp('Data dimensions unexpected')
        return
    end
    disp(['num channels ', num2str(n_Coil), ...
        ' :: n phase ' num2str(n_ph) ' :: n read ' num2str(n_rd)])
    
end

%% create trajectory to use with bart

for ii = 1:2
    txt = sprintf('k%d_coords_normalized = output_ksp_traj.ksp_norm_%d;',ii,permute_order(ii));
    eval(txt);
end

%% caclulate median trajectory
if use_median_traj == 1
    
    disp('Img recon: using median k-space trajectory')
    
    for ii = 1:n_ph
        idx_angle = mod(ii-1,n_angles) + 1;
        k1_coords_normalized(:,ii) = median( k1_coords_normalized( :, idx_angle:n_angles:end ), 2 );
        k2_coords_normalized(:,ii) = median( k2_coords_normalized( :, idx_angle:n_angles:end ), 2 );
    end
    
end

%% reformat traj for BART

traj_bart = zeros([3 n_ph n_rd]);
for ii = 1:n_ph
    traj_bart(1,ii,:) = ( (effMtx - 1) ) * k1_coords_normalized(:,ii);
    traj_bart(2,ii,:) = ( (effMtx - 1) ) * k2_coords_normalized(:,ii);
end

%% kx, ky (and kz=0) coordinates
k1_coords_normalized_transpose = k1_coords_normalized.';
k2_coords_normalized_transpose = k2_coords_normalized.';
crds = zeros(3, n_ph*n_rd);
crds(1,:) = k1_coords_normalized_transpose(:);
crds(2,:) = k2_coords_normalized_transpose(:);
crds(3,:) = 0;

clear  k1_coords_normalized_transpose k2_coords_normalized_transpose;

%% calculate density compensation
tic;
numIter = 25;
osf     = 2.1;
verbose = 1;
DCF = sdc3_MAT(crds,numIter,effMtx,verbose,osf);
DCF = reshape(DCF,[n_ph n_rd]);

disp(sprintf('Img recon: completed SDC in %.1f seconds', toc));

clear crds

%% create k-space to use with bart
ksp_bart = complex(zeros([1 n_ph n_rd n_Coil],'single'));
for ii = 1:n_Coil,
    ksp_bart(1,:,:,ii) = data(ii,:,:);
end

clear data

%% do coil compression

if cc_factor ~= 0 && cc_factor ~= 1
    
    n_cc = ceil( n_Coil / cc_factor ); % number of virtual coils
    
    bart_cmd_txt = sprintf( 'cc -p %d -S -A',n_cc );
    ksp_bart_cc = bart( bart_cmd_txt, ksp_bart);
    
    output_args.cc.n_cc = n_cc;
    output_args.cc.bart_cmd_txt = bart_cmd_txt;
    
else
    
    n_cc = n_Coil;
    ksp_bart_cc = ksp_bart;
    
end

clear ksp_bart

%% calculate sensitivity maps
tic;
bart_cmd_txt_1 = 'nufft -i -d30:30:1 -t';
lowres_img = bart(bart_cmd_txt_1, traj_bart, ksp_bart_cc);
bart_cmd_txt_2 = 'fft -u 7';
lowres_ksp = bart(bart_cmd_txt_2, lowres_img);
bart_cmd_txt_3 = sprintf('resize -c 0 %d 1 %d', effMtx, effMtx);
ksp_zerop = bart( bart_cmd_txt_3, lowres_ksp);
bart_cmd_txt_4 = sprintf('ecalib -I -t %.6f -m1',  ecalib_threshold);
sens = bart( bart_cmd_txt_4, ksp_zerop); % -I for inhomogeneity correction

output_args.sens.bart_cmd_txt_1 = bart_cmd_txt_1;
output_args.sens.bart_cmd_txt_2 = bart_cmd_txt_2;
output_args.sens.bart_cmd_txt_3 = bart_cmd_txt_3;
output_args.sens.bart_cmd_txt_4 = bart_cmd_txt_4;

disp(sprintf('Img recon: completed bart sensitivity map estimation in %.1f seconds', toc));

%% display sensitivities
for ii = 1:n_cc
    figure(2); clf;
    imagesc(squeeze(real(sens(:,:,1,ii)))); axis image
    title(sprintf('bart: estimated coil sensitivity maps coil %d', ii ) );
    pause(0.01)
end

%% perform undersampled reconstruction for each TR

delete(gcp('nocreate'))
my_pool = parpool([1 128],'SpmdEnabled',false);

bart_cmd_txt = sprintf('nufft -a -d%d:%d:1', effMtx, effMtx);
recon_each_TR = complex(zeros([effMtx effMtx n_ph],'single'));
sens_norm = 1./ ( sum( squeeze( conj( sens(:,:,1,:) ).* sens(:,:,1,:) ), 3 ) );
sens_norm( isinf( sens_norm ) ) = 0;
parfor TR = 1:n_ph,
    tic;
    recon_this_TR = bart( bart_cmd_txt, traj_bart(:,TR,:), ksp_bart_cc(1,TR,:,:) .* shiftdim(repmat(DCF(TR,:),[1 1 n_cc]),-1) );
    
    recon_this_TR_coil_combine = complex(zeros([effMtx effMtx],'single'));
    for coil = 1:n_cc,
        recon_this_TR_coil_combine = recon_this_TR_coil_combine + squeeze( recon_this_TR(:,:,1,coil) .* conj(sens(:,:,1,coil)) );
    end
    recon_this_TR_coil_combine = recon_this_TR_coil_combine.* sens_norm;
    
    recon_each_TR(:,:,TR) = recon_this_TR_coil_combine;
    
    disp(sprintf('Img recon: completed bart recon of TR %04d of %04d in %.1f seconds', TR, n_ph, toc));
end

output_args.recon_stack.bart_cmd_txt = bart_cmd_txt;

%% optional flipping
if exist('flip1','var')
    if flip1 == 1
        recon_each_TR = flip(recon_each_TR,1);
    end
end
if exist('flip2','var')
    if flip2 == 1
        recon_each_TR = flip(recon_each_TR,2);
    end
end


%% display all parts of reconstruction using complex sum of each TR reconstruction
recon_each_TR_complex_sum = sum(recon_each_TR,3);
hfig = figure(10); clf;
% set(hfig,'Position',[200 200 1200 1000]);

subplot(2,2,1); imagesc(abs(recon_each_TR_complex_sum)); axis image; colormap(gray); colorbar;
title('abs(recon_each_TR_complex_sum)','interpreter','none');

subplot(2,2,2); imagesc(angle(recon_each_TR_complex_sum)); axis image; colormap(gray); colorbar;
title('angle(recon_each_TR_complex_sum)','interpreter','none');

subplot(2,2,3); imagesc(real(recon_each_TR_complex_sum)); axis image; colormap(gray); colorbar;
title('real(recon_each_TR_complex_sum)','interpreter','none');

subplot(2,2,4); imagesc(imag(recon_each_TR_complex_sum)); axis image; colormap(gray); colorbar;
title('imag(recon_each_TR_complex_sum)','interpreter','none');

%% store results for return

output_args.MRF_img_stack_coil_combined = recon_each_TR;
output_args.ksp_traj_med = cat( 3, k1_coords_normalized, k2_coords_normalized, zeros( size( k1_coords_normalized ), 'single' ) );


%% report total execution time

delete(my_pool);

disp(sprintf('Img recon: TOTAL EXECUTION TIME: %.1f seconds', toc(start_time) ));

end