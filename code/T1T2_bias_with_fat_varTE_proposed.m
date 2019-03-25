%% MRF k-space - image-space simulation w/B0 effects

clear, close all, clc;

%% params

% dirs and files
dir_in = '../data_in/';
dir_out = '../data_out/';
fn_dict = 'MRF_dict_varTE.mat';
fn_MRF_params_csv = 'MRF103.csv';
fn_save = 'T1T2_bias_with_FSF_TEvar.mat';

i_sign = 1; % offresonance precession convention
flag_gen_dict = 1; % if 1 then generate dictionary, else load

% simulation and fit properties
n_B0_basis = 31;
sigma_noise = 0.; % stand. dev. of complex Gaussian noise
s_frac = 1 - 1e-4; % fraction of dictionary energy to maintain after compression of coefficients
s_frac_water = 1 - 1e-4; % fraction of water dictionary eneery to maintain after compression

% FSFs to simulate
fsf_q_v = [0:0.05:1.0]; % fat fraction range
n_fsf = numel(fsf_q_v);

% constants
gamma_bar_Hz_per_T = 42.5775e6;
ppm_H2O = 4.65;
B0_Tesla = 3.0;
delk = 1; % EPG steps
szomega = 101; % size of EPG state matrix

% MRF sequence params
n_TR = 1500; % length of MRF train
TE_base_s = 3.5e-3; % minimum TE
del_TE_s = 4.0e-3; % span of TE
TR_base_s = 16e-3; % minimum TR
FA_deg_nom = 60; % maximum flip angle
TI_s = 40e-3;

% get MRF modulation params from .csv file
MRF_data = csvread([dir_in fn_MRF_params_csv],1,0);
TE_s_v = MRF_data(:,4)*1e-3 + TE_base_s;
% TE_s_v = linspace( TE_base_s, TE_base_s + del_TE_s, n_TR )';
TR_s_v = TR_base_s * ones(n_TR,1); % fix TR
phi_v = zeros(n_TR,1);
FA_deg_v = FA_deg_nom * MRF_data(:,1);

% WAT T1, T2 characterization from Hamilton et al., JMRI, 2011
p_rel_H2O_v = (ppm_H2O - [0.9 1.3 1.3 2.1 2.1 2.75 4.2 5.3]) * 1e-6; % parts
fat_freq_Hz_v = p_rel_H2O_v * gamma_bar_Hz_per_T * B0_Tesla;
fat_amps_v = [0.144 0.5*1.0 0.5*1.0 0.5*0.241 0.5*0.241 0.033 0.064 0.122];
fat_amps_v = fat_amps_v./sum(fat_amps_v);
T1_fat_s_v = [543 280 240 249 202 284 154 421]/1000;
T2_fat_s_v = [80.1 54.7 54.7 51.9 51.9 46.2 50 44.1]/1000; % 7th entry not supplied, 50 ms is estimate
n_fat_pks = numel(fat_freq_Hz_v);


%% model fat signal

fat_signal_v = zeros(n_TR,1);
for ii = 1:n_fat_pks

    T1_s = T1_fat_s_v(ii);
    T2_s = T2_fat_s_v(ii);
    B0_Hz = fat_freq_Hz_v(ii);
    
    my_signal_v = EPG_MRF_SSFP_B0( T1_s, T2_s, TE_s_v, TR_s_v, FA_deg_v, delk, n_TR, szomega, phi_v, TI_s, B0_Hz, i_sign );
    my_signal_v = my_signal_v./norm(my_signal_v); % normalize each peak
    fat_signal_v = fat_signal_v + fat_amps_v(ii).*my_signal_v(:);
    
end

% normalize fat signal
fat_signal_v = fat_signal_v./norm(fat_signal_v);

%% create or load water dictionary

if flag_gen_dict == 1
    
    % dictionary params
    input_dict.TI = TI_s;
    input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000]/1000;
    input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500]/1000;
    input_dict.B1_v = [1.0];
    input_dict.FA_v = FA_deg_v;
    input_dict.phi_v = phi_v;
    input_dict.TR_v = TR_s_v;
    input_dict.TE_v = TE_s_v;
    input_dict.nreps = n_TR;
    
    input_dict.delk = delk; % step between states equal to a full dephasing imparted by crusher gradient
    input_dict.szomega = szomega; % number of factors of k to include in phase history
    
    
    % plot sequence parameters
    figure(1); clf;
    subplot(411)
    plot( input_dict.FA_v ); ylabel('degrees'); title(['FA, TI is ' num2str(input_dict.TI) ' ms'])
    subplot(412)
    plot( input_dict.phi_v*1e3 ); ylabel('degrees'); title('phase')
    subplot(413)
    plot( input_dict.TR_v*1e3 ); ylabel('msec'); title('TR')
    subplot(414)
    plot( input_dict.TE_v*1e3 ); ylabel('msec'); title('TE')
    drawnow
    
    % do dictionary construction
    input_dict.reduce = 1;
    input_dict.sfrac = s_frac_water;
    
    % construct
    output_dict = MRF_dict_B1(input_dict);
    
    save([dir_in fn_dict],'output_dict')
    
else
    
    load([dir_in fn_dict])
    
end

%% determine freq domain via MFI

%   INPUT: params.method = 'DFT';
%               ".M = pos integer is number of basis frequencie
%               ".params.t0 = start of acq time, TE for center out
%               ".t_end = end of acquisiiton time, TE + tacq for center out
%               ".t1 = time extension on the beginning of target vectors
%               ".t2 = time extension on the end of target vectors
%               ".N = pos integer is discretizations of time
%               ".f_max = maximum frequency to consider
%               ".f_step = frequency spacing in lookup table to be
%                   generated

params_MFI.method = 'DFT';
params_MFI.M = n_B0_basis;
params_MFI.t0 = TE_base_s;
params_MFI.t_end = 5e-3 + del_TE_s + TE_base_s;
params_MFI.t1 = TE_base_s;
params_MFI.t2 = 1.2*(params_MFI.t_end) - params_MFI.t_end;
params_MFI.N = 1600;
params_MFI.f_max = 700;
params_MFI.f_step = 10;

results_MFI_coeff = MFI_coeffs(params_MFI);

%% create off-resonance dictionary for compression

n_dict = size(output_dict.dict_norm,2);

dict_composite = complex( zeros(n_TR, n_B0_basis*n_dict) );

t_v = TE_s_v(:);
for ii = 1:n_B0_basis
    my_basis_v = exp( i_sign * 1i * 2 * pi * t_v * results_MFI_coeff.delta_f_v(ii) );
    my_basis_mat = repmat( my_basis_v, [ 1 n_dict ] );
    my_index_v = ( ( ii - 1 )*n_dict + 1 ):( ii*n_dict );
    dict_composite(:,my_index_v) = my_basis_mat.* output_dict.dict_norm;
end

%% get SVD of off-resonance dictionary

tic
[U,S,~] = svd(dict_composite,'econ');
toc

s_v = diag(S);
s2_cum_rel_v = cumsum(s_v.^2)/sum(s_v.^2);

r_SVD_space = sum(s2_cum_rel_v <= s_frac);

Ur = U(:,1:r_SVD_space);

%% generate synthetic data with off-resonance

signal_m = zeros( n_TR, n_dict, n_fsf );
for ii = 1:n_fsf

    % create on-resonance signals
    fat_m = repmat(fat_signal_v(:),[1 n_dict]);
    water_m = output_dict.dict_norm;
    max_signal = max([ max(abs(fat_m(:))); max(abs(water_m(:))) ] );
    noise_m = sigma_noise * max_signal * (randn(n_TR,n_dict) + 1i*randn(n_TR,n_dict));
    
    % create off-resonance
    del_B0_Hz = 0;
    j_v = exp( i_sign * 1i * 2 * pi * del_B0_Hz * TE_s_v(:) );
    j_mat = repmat( j_v, [ 1 n_dict ] );
    
    my_fsf = fsf_q_v(ii);
    
    % final composite synthetic signal
    signal_m(:,:,ii) = (my_fsf.* fat_m.* j_mat) + ((1 - my_fsf).* water_m.* j_mat) + noise_m;

end

%% form augmented water dictionary w/fat

[U_water_dict,S_water_dict, ~] = svd(output_dict.dict_norm,'econ');

s_v = diag(S_water_dict);
s2_cum_rel_v = cumsum(s_v.^2)/sum(s_v.^2);

r_water_SVD_space = sum(s2_cum_rel_v <= s_frac_water);

Ur_water_dict = U_water_dict(:,1:r_water_SVD_space);

dict_norm_augmented = [Ur_water_dict fat_signal_v(:)];
pinv_dict_norm_augmented = pinv(dict_norm_augmented);

%% get coeffs for all of all signals using basis from off-res dictionary for each freq. modulated basis

t_v = TE_s_v(:);

coeffs_measured = zeros(r_SVD_space, n_dict, n_fsf, n_B0_basis);
fatwater_coeffs = zeros(r_water_SVD_space+1, n_dict, n_fsf, n_B0_basis);

for jj = 1:n_fsf
    fprintf( 'T1T2 bias proposed simulation: calc coefficients for fsf %.2f\n', fsf_q_v(jj) );
    for ii = 1:n_B0_basis
        
        % form MFI basis
        j_v = exp( -i_sign * 1i * 2 * pi * results_MFI_coeff.delta_f_v(ii) * t_v(:) ); % note conjugate of i_sign
        j_mat = repmat( j_v, [ 1 n_dict ] );
        
        % form signal basis
        signal_basis_m = signal_m(:,:,jj).* j_mat;
        
        % single point simulation coefficients
        coeffs_measured(:,:,jj,ii) = Ur'*signal_basis_m;
        fatwater_coeffs(:,:,jj,ii) = pinv_dict_norm_augmented*signal_basis_m;
        
    end
end

%% get MFI coefficents over B0 range and relevant spacing

freq_fine_Hz_v = results_MFI_coeff.delta_f_fine_v;
n_freq_fine = numel(freq_fine_Hz_v);
coeffs_MFI = results_MFI_coeff.c_DFT;


%% for each signal cycle through B0 to determine which B0 minimizes the residual

% form operative VARPRO matrix
VARPRO_mat = Ur(:,1:r_water_SVD_space + 1)'*eye(n_TR) - Ur(:,1:r_water_SVD_space + 1)'*dict_norm_augmented*pinv_dict_norm_augmented;

% loop over test signals
B0_solution_Hz = zeros(n_dict,n_fsf);
idx_B0_solution = zeros(n_dict,n_fsf);

for kk = 1:n_fsf
    fprintf( 'T1T2 bias proposed simulation: B0 fit for fsf %.2f\n', fsf_q_v(kk) );
    for ii = 1:n_dict
        
        my_coeffs = squeeze( coeffs_measured(:,ii,kk,:) );
        
        my_test_signals = Ur*my_coeffs;
        my_obj_mat = VARPRO_mat*my_test_signals*coeffs_MFI;
        
        % loop over fine frequency discretization to residual magnitude
        r2_v = zeros( n_freq_fine, 1);
        for jj = 1:n_freq_fine;
            
            r2_v(jj) = norm( my_obj_mat(:,jj) )^2;
            
        end
        
        [r2_min, idx_min] = min( r2_v(:) );
        
%         figure(10);
%         plot(freq_fine_Hz_v,r2_v);
%         title( sprintf('min objective at %.1f Hz :: element %d',freq_fine_Hz_v(idx_min), ii ) );
%         xlabel('B0 (Hz)')
%         ylabel('objective value')
%         drawnow      
        
        B0_solution_Hz(ii,kk) = freq_fine_Hz_v(idx_min);
        idx_B0_solution(ii,kk) = idx_min;
        
    end
end


%% calculate water and fat components

water_fit = zeros(n_TR, n_dict, n_fsf);
water_fit_coeffs = zeros(r_water_SVD_space, n_dict, n_fsf);
fat_fit = zeros(n_TR, n_dict, n_fsf);
fsf_fit = zeros(n_dict,n_fsf);
for jj = 1:n_fsf
    fprintf( 'T1T2 bias proposed simulation: calc water-fat components for fsf %.2f\n', fsf_q_v(jj) );
    for ii = 1:n_dict
        
        my_coeffs_measured = squeeze( coeffs_measured(:,ii,jj,:) );
        
        signal_conj_phase_v = Ur*my_coeffs_measured*coeffs_MFI(:,idx_B0_solution(ii,jj));
        
        beta_v = pinv_dict_norm_augmented*signal_conj_phase_v;
        
        water_fit(:,ii,jj) = Ur_water_dict*beta_v( 1:(end-1) );
        water_fit_coeffs(:,ii,jj) = beta_v( 1:(end-1) );
        fat_fit(:,ii,jj) = fat_signal_v*beta_v(end);
        fsf_fit(ii,jj) = abs( beta_v(end) ) / ( abs( beta_v(end))  + norm( beta_v( 1:(end-1) ) ) );
        
    end
end

%% get T1 and T2 etc. of water

dict_red = Ur_water_dict' * output_dict.dict_norm;

T1_fit = zeros(n_dict,n_fsf);
T2_fit = T1_fit;
M0_fit = T1_fit;
for jj = 1:n_fsf
    fprintf( 'T1T2 bias proposed simulation: get water T1 and T2 estimates for fsf %.2f\n', fsf_q_v(jj) );
    for ii = 1:n_dict
        
        [r_max, my_idx_max] = max( dict_red' * water_fit_coeffs(:,ii,jj) );
        
        T1_fit(ii,jj) = output_dict.dict_list(my_idx_max,1);
        T2_fit(ii,jj) = output_dict.dict_list(my_idx_max,2);
        M0_fit(ii,jj) = water_fit(:,ii,jj)' * output_dict.dict_norm(:,my_idx_max);
        
    end
end

T1T2_fit_ms = 1000*reshape( [T1_fit; T2_fit], [n_dict 2 n_fsf] );

%% calculate absolute and relative bias for each estimate

T1T2_ref_ms = 1000*repmat( output_dict.dict_list(:,1:2), [1 1 n_fsf] );

T1T2_abs_bias_ms = T1T2_fit_ms - T1T2_ref_ms;
T1T2_rel_bias = abs(T1T2_abs_bias_ms)./T1T2_ref_ms;

save([dir_out fn_save],'T1T2_ref_ms','T1T2_abs_bias_ms','T1T2_rel_bias','fsf_q_v')