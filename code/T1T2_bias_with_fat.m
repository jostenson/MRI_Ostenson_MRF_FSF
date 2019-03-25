%% Script to simulate bias introduced by fat to MRF signal

clear all, close all, clc;


%% Parameters

dir_in = '../data_in/';
dir_out = '../data_out/';
fn_dict = 'water_dictionary_1500_TRfix_TEfix.mat';
fn_MRF_params_csv = 'MRF104.csv';
fn_save = 'T1T2_bias_with_FSF.mat';

flag_gen_dict = 1; % if 1 then generate dictionary, else load

% simulation and fit properties
n_iter = 1; % number of iterations
sigma_noise = 0.0; % stand. dev. of complex Gaussian noise
s_frac = 1 - 1e-4; % fraction of dictionary energy to maintain after compression
s_frac_water = 1 - 1e-4;

% FSFs to simulate
fsf_q_v = [0:0.05:1.0]; % fat fraction range
n_fsf = numel(fsf_q_v);

% constants
i_sign = 1; % off-resonance precession convention
gamma_bar_Hz_per_T = 42.5775e6;
ppm_H2O = 4.65;
B0_Tesla = 3.0;
delk = 1; % EPG steps
szomega = 100; % size of EPG state matrix

% MRF sequence params
n_TR = 1500; % length of MRF train
TE_base_s = 4.65e-3; % minimum TE
TR_base_s = 16e-3; % minimum TR
FA_deg_nom = 60; % maximum flip angle
TI_s = 40e-3;

% get MRF modulation params from .csv file
MRF_data = csvread([dir_in fn_MRF_params_csv],1,0);
FA_deg_v = FA_deg_nom * MRF_data(:,1);
phi_v = MRF_data(:,2);
TR_s_v = TR_base_s + 1e-3 * MRF_data(:,3);
TE_s_v = TE_base_s + 1e-3 * MRF_data(:,4);

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
    
    my_signal_v = EPG_MRF_SSFP( T1_s, T2_s, TE_s_v, TR_s_v, FA_deg_v, delk, n_TR, szomega, phi_v, TI_s );
    my_signal_v = my_signal_v./ norm(my_signal_v);
    fat_signal_v = fat_signal_v + fat_amps_v(ii).*my_signal_v(:).*exp( i_sign * 1i * 2 * pi * fat_freq_Hz_v(ii) * TE_s_v );
    
end

% normalize fat signal
fat_signal_v = fat_signal_v./norm(fat_signal_v);

% figure(1); clf;
% plot(abs(fat_signal_v))
% xlabel('MRF TR'); ylabel('signal mag (au)')

%% create or load water dictionary

if flag_gen_dict == 1
    
    % dictionary params
    input_dict.TI = TI_s;
    input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000]/1000;
    input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800 900:100:1500]/1000;
    input_dict.B1_v = [1.0];
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
    plot( input_dict.FA_v ); ylabel('degrees'); title(['FA, TI is ' num2str(input_dict.TI) ' s'])
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

%% for each dictionary element mix fat and water across all FSFs and estimate T1 and T2

n_dict = size( output_dict.dict_list, 1 );

T1T2_estimates_ms = zeros( n_dict, 2, n_fsf );
Ur = output_dict.U_r;
for ii = 1:n_fsf
    
    fprintf('T1/T2 bias prediction: calc T1/T2 estimates for FSF %.3f\n',fsf_q_v(ii));
    
    for jj = 1:n_dict
        
        fsf = fsf_q_v(ii);
        my_sig_v = (1 - fsf)*output_dict.dict_norm(:,jj) + fsf*fat_signal_v;
        
        my_sig_v = my_sig_v./norm(my_sig_v);
        
        r_v = (Ur'*my_sig_v)'*output_dict.dict_compress;
        
        [max_r,max_ind] = max(abs(r_v(:)));
        
        T1T2_estimates_ms(jj,:,ii) = output_dict.dict_list(max_ind,1:2)*1000;
        
    end
    
end

%% calculate absolute and relative bias for each estimate

T1T2_ref_ms = 1000*repmat( output_dict.dict_list(:,1:2), [1 1 n_fsf] );

T1T2_abs_bias_ms = T1T2_estimates_ms - T1T2_ref_ms;
T1T2_rel_bias = abs(T1T2_abs_bias_ms)./T1T2_ref_ms;

save([dir_out fn_save],'T1T2_ref_ms','T1T2_abs_bias_ms','T1T2_rel_bias','fsf_q_v')

%% plot results

for ii = 1:n_fsf
    figure(10); clf;
    imagesc( abs( T1T2_rel_bias(:,:,ii) ) ); colorbar(); caxis([0 0.5]);
    title(sprintf('Relative T1 and T2 bias for FSF %.3f',fsf_q_v(ii)))
    pause(.2)
end
