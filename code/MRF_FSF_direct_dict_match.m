%% direct dictionary match of FSF using previously processed MRF data

clear, close all, clc

%% parameters

dir_in = '../data_in/';
dir_out = '../data_out/';

fn_dict = 'MRF_dict_emulsion_varTE.mat';
fn_MRF_img = 'MRF_emulsion_varTE_img_recon.mat';
fn_MRF_params_csv = 'MRF103.csv';
fn_out = 'MRF_direct_FSF_emulsion.mat';
fn_B1 = 'B1_emulsion.mat';

flag_gen_dict = 1;

% FSFs to simulate
fsf_q_v = [0:0.05:1.0]; % fat fraction range
n_fsf = numel(fsf_q_v);

% constants
i_sign = -1; % off-resonance precession convention
gamma_bar_Hz_per_T = 42.5775e6;
ppm_H2O = 4.65;
B0_Tesla = 3.0;
delk = 1; % EPG steps
szomega = 101; % size of EPG state matrix

% MRF sequence params
n_TR = 1500; % length of MRF train
TE_base_s = 3.5e-3; % minimum TE
TR_base_s = 16e-3; % minimum TR
FA_deg_nom = 60; % maximum flip angle
TI_s = 40e-3;
B1_v = [0.5, 0.6, 0.7 0.75, 0.8, 0.825:0.025:1.2, 1.25 1.3, 1.4, 1.5];
n_B1 = numel( B1_v );

% get MRF modulation params from .csv file
MRF_data = csvread([dir_in fn_MRF_params_csv],1,0);
FA_deg_v = FA_deg_nom * MRF_data(:,1);
phi_v = MRF_data(:,2);
TR_s_v = TR_base_s + 1e-3 * MRF_data(:,3);
TE_s_v = TE_base_s + 1e-3 * MRF_data(:,4);
s_frac_water = 1 - 1e-4;

% WAT T1, T2 characterization from Hamilton et al., JMRI, 2011
p_rel_H2O_v = (ppm_H2O - [0.9 1.3 1.3 2.1 2.1 2.75 4.2 5.3]) * 1e-6; % parts
fat_freq_Hz_v = p_rel_H2O_v * gamma_bar_Hz_per_T * B0_Tesla;
fat_amps_v = [0.144 0.5*1.0 0.5*1.0 0.5*0.241 0.5*0.241 0.033 0.064 0.122];
fat_amps_v = fat_amps_v./sum(fat_amps_v);
T1_fat_s_v = [543 280 240 249 202 284 154 421]/1000;
T2_fat_s_v = [80.1 54.7 54.7 51.9 51.9 46.2 50 44.1]/1000; % 7th entry not supplied, 50 ms is estimate
n_fat_pks = numel(fat_freq_Hz_v);

%% model fat signal

fat_signal = zeros(n_TR,n_B1);
for jj = 1:n_B1
    
    fat_signal_v = zeros(n_TR,1);
    for ii = 1:n_fat_pks
    
        T1_s = T1_fat_s_v(ii);
        T2_s = T2_fat_s_v(ii);

        my_signal_v = EPG_MRF_SSFP( T1_s, T2_s, TE_s_v, TR_s_v, B1_v(jj)*FA_deg_v, delk, n_TR, szomega, phi_v, TI_s );
        my_signal_v = my_signal_v./ norm(my_signal_v);
        fat_signal_v = fat_signal_v + fat_amps_v(ii).*my_signal_v(:).*exp( i_sign * 1i * 2 * pi * fat_freq_Hz_v(ii) * TE_s_v );
    
    end
    % normalize fat signal
    fat_signal_v = fat_signal_v./norm(fat_signal_v);
    % figure(1); clf;
    % plot(abs(fat_signal_v))
    % xlabel('MRF TR'); ylabel('signal mag (au)')
    
    fat_signal(:,jj) = fat_signal_v;
    
end


%% create or load water dictionary

if flag_gen_dict == 1
    
    % dictionary params
    input_dict.TI = TI_s;
    input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000]/1000;
    input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800, 900:100:1500]/1000;
    input_dict.B1_v = B1_v;
    input_dict.FA_v = FA_deg_v;
    input_dict.phi_v = phi_v;
    input_dict.TR_v = TR_s_v;
    input_dict.TE_v = TE_s_v;
    input_dict.nreps = n_TR;
    
    input_dict.delk = delk; % step between states equal to a full dephasing imparted by crusher gradient
    input_dict.szomega = szomega; % number of factors of k to include in phase history
    
    
    % plot sequence parameters
    figure(2); clf;
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
    
    save([dir_out fn_dict],'output_dict')
    
else
    
    load([dir_out fn_dict])
    
end

%%

% init FSF dictionary
n_dict = size(output_dict.dict,2);
dict_FSF_norm = complex( zeros(n_TR, n_dict*n_fsf, 'single' ) );

% create FSF dictionary
dict_list_FSF = repmat(output_dict.dict_list,[n_fsf 1]);
dict_list_FSF(:,3) = 0;
dict_list_FSF(:,4) = 0;
for jj = 1:n_B1
    jj
    logical_B1_v = output_dict.dict_list(:,3) == B1_v(jj);
    n_dict_sub = sum( logical_B1_v );
    fat_signal_m = repmat(fat_signal(:,jj),[1 n_dict_sub]);
   
    for ii = 1:n_fsf

        my_idx_v = ((ii-1)*n_dict_sub + (jj-1)*n_dict_sub*n_fsf + 1):(ii*n_dict_sub + (jj-1)*n_dict_sub*n_fsf);
        dict_FSF_norm(:,my_idx_v) = fsf_q_v(ii)*fat_signal_m + (1 - fsf_q_v(ii))*output_dict.dict_norm(:,logical_B1_v);
        dict_list_FSF(my_idx_v,3) = B1_v(jj);
        dict_list_FSF(my_idx_v,4) = fsf_q_v(ii);
        
    end

end


% normalize
for ii = 1:size(dict_FSF_norm,2)
    dict_FSF_norm(:,ii) = dict_FSF_norm(:,ii)./norm(dict_FSF_norm(:,ii));
end


%% compress FSF dictionary

% calc SVD
Uk_FSF = [];
dict_FSF_svd = [];
dict_FSF_B1_svd_v = [];
for ii = 1:n_B1
    ii
    [U,S,~] = svd( dict_FSF_norm(:,dict_list_FSF(:,3) == B1_v(ii) ),'econ');
    s_v = diag(S);
    if ii == 1
        idx_svd_sig_v = cumsum(s_v.^2)./sum(s_v.^2) < s_frac_water*ones(size(s_v));
        max_idx = sum( idx_svd_sig_v ) + 1; % account for small variation in rank across B1_v
    end
    Uk = U(:,1:max_idx);

    % compress FSF dict
    Uk_FSF = [Uk_FSF Uk];
    dict_FSF_svd = [dict_FSF_svd Uk'*dict_FSF_norm(:,dict_list_FSF(:,3) == B1_v(ii)) ];
    dict_FSF_B1_svd_v = [dict_FSF_B1_svd_v B1_v(ii)*ones(1, max_idx)];
    
end

% % check error
% error_svd = dict_FSF - Uk*dict_FSF_svd;
% rnorm_v = zeros( 1, size( error_svd,2 ) );
% for ii = 1:size( error_svd, 2 );
%     rnorm_v(ii) = sqrt( error_svd(:,ii)'*error_svd(:,ii) );
% end

%% do dictionary match with FSF dictionary

load([dir_out fn_MRF_img])
load([dir_in fn_B1])
B1_map = output_B1_data.B1_map;

effMtx = size( output_img_recon.MRF_img_stack_coil_combined, 1 );

T1_map_FSF = zeros(effMtx);
T2_map_FSF = zeros(effMtx);
R_map_FSF = zeros(effMtx);
M0_map_FSF = zeros(effMtx);
FSF_map = zeros(effMtx);

for ii = 1:effMtx
    fprintf('MRF FSF direct match: Processing row %d\n',ii)
    for jj = 1:effMtx
        
        B1 = B1_map(ii,jj);
        [B1_min_diff, idx_min] = min( abs( B1_v(:) - B1 ) );
        logical_idx_B1_col_v = dict_FSF_B1_svd_v == B1_v(idx_min);
        my_Uk = Uk_FSF(:,logical_idx_B1_col_v);

        idx_dict_v = find( dict_list_FSF(:,3) == B1_v(idx_min) );
        my_dict_list = dict_list_FSF(idx_dict_v,:);
        
        test_v = squeeze(output_img_recon.MRF_img_stack_coil_combined(ii,jj,:));
        test_norm_v = test_v(:)./norm(test_v);
        test_compressed_v = my_Uk'*test_norm_v;
        
        my_dict = dict_FSF_svd(:,idx_dict_v);
        [max_ip,idx_ip] = max(my_dict'*test_compressed_v);
        
        T1_map_FSF(ii,jj) = my_dict_list(idx_ip,1);
        T2_map_FSF(ii,jj) = my_dict_list(idx_ip,2);
        B1_map(ii,jj) = my_dict_list(idx_ip,3);
        FSF_map(ii,jj) = my_dict_list(idx_ip,4);
        
        R_map_FSF(ii,jj) = max_ip;
        M0_map_FSF(ii,jj) = (my_Uk'*test_v)'*squeeze(my_dict(:,idx_ip));
         
    end
end

%% save results

% output_MRF_FSF_direct_match.dict_FSF_svd = dict_FSF_svd;
output_MRF_FSF_direct_match.dict_list_FSF = dict_list_FSF;
output_MRF_FSF_direct_match.f_MRF = fn_MRF_img;
output_MRF_FSF_direct_match.fat_amps_v = fat_amps_v;
output_MRF_FSF_direct_match.fat_freq_Hz_v = fat_freq_Hz_v;
output_MRF_FSF_direct_match.T1_fat_s_v = T1_fat_s_v;
output_MRF_FSF_direct_match.T2_fat_s_v = T2_fat_s_v;
output_MRF_FSF_direct_match.Uk = Uk;
output_MRF_FSF_direct_match.sfrac = s_frac_water;
output_MRF_FSF_direct_match.T1_map = T1_map_FSF;
output_MRF_FSF_direct_match.T2_map = T2_map_FSF;
output_MRF_FSF_direct_match.R_map = R_map_FSF;
output_MRF_FSF_direct_match.M0_map = M0_map_FSF;
output_MRF_FSF_direct_match.FSF_map = FSF_map;

save([dir_out fn_out],'output_MRF_FSF_direct_match')

