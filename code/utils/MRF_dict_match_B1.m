function output_match = MRF_dict_match_B1( output_recon, output_dict, reduce_flag )
% do MRF dictionary construction/match
% used in conjunction with MRF processing

%   INPUT: output_recon.MRF_img_stack_coil_combined = Nr x Nc x nSig matrix of MRF img data
%                 ".B1_map = provided B1 map
%          output_dict.dict_compress = SVD compressed dictionary
%                   ".U_r = significant left singular vectors
%                   ".dict_list = table of dictionary parameters
%                   ".dict_norm = normalized full dictionary
%                   ".B1_compress_dict_list_v = list of B1 for cols of
%                       compressed dictionary
%          reduce_flag = 1 then dictionary is reduced by SVD, else no
%   OUTPUT: output.T1_map = T1_map
%                ".T2_map = T2_map
%                ".B1_map = B1_map
%                ".M0_map = magnetization (proton) density
%                ".R_map = complex matrix is correlation map
%                ".dict_list = list of dictionary parameters that
%                   correspond to the entries in the dictionary
%                ".dict_fn = input dictionary filename
% assumes reduced dictionary space is already calculated
disp('Doing MRF dictionary match...');
tic;

img_stack = output_recon.MRF_img_stack_coil_combined;
[Nr, Nc, ~] = size(img_stack);
T1_map = zeros(Nr,Nc);
T2_map = T1_map;
R_map = T1_map;
M0_map = T1_map;

B1_map = output_recon.B1_map;

dict_list = output_dict.dict_list;
dict_norm = output_dict.dict_norm;
dict_norm_H = dict_norm';

B1_unique_v = unique(output_dict.B1_compress_dict_list_v);

for ii = 1:Nr
    disp(['MRF Tx map generating, row ' num2str(ii)])
    for jj = 1:Nc
        
        test_v = squeeze(img_stack(ii,jj,:));
        test_norm_v = test_v./sqrt(test_v'*test_v); % normalized
        
        % get closest discretized B1
        B1 = B1_map(ii,jj);
        [B1_min_diff, idx_min] = min( abs( B1_unique_v(:) - B1 ) );
        logical_idx_B1_col_v = output_dict.B1_compress_dict_list_v(:) == B1_unique_v(idx_min);        

        idx_dict_v = find( output_dict.dict_list(:,3) == B1_unique_v(idx_min) );   
        
        min_dict_idx = min( idx_dict_v(:) );
        
        if reduce_flag == 1
            U_r = output_dict.U_r(:,logical_idx_B1_col_v);
            dict_compress = output_dict.dict_compress(:,idx_dict_v);
            testR_v = U_r'*test_norm_v(:); % project onto SVD space
            ipDT_v = dict_compress'*testR_v; % determine basis coeff correlation
        else
            ipDT_v = dict_norm_H(idx_dict_v,:)*test_norm_v(:);
        end
        
        [maxIP, idxMax] = max(ipDT_v);
        
        idxMax = idxMax + min_dict_idx - 1;
        
        T1_map(ii,jj) = dict_list(idxMax,1);
        T2_map(ii,jj) = dict_list(idxMax,2);
        B1_map(ii,jj) = dict_list(idxMax,3);
        R_map(ii,jj) = maxIP;
        
        M0_map(ii,jj) = dict_norm_H(idxMax,:)*test_v(:);
        
    end
end

output_match.T1_map = T1_map;
output_match.T2_map = T2_map;
output_match.B1_map = B1_map;
output_match.R_map = R_map;
output_match.M0_map = M0_map;
output_match.dict_list = dict_list;
output_match.dict_fn = output_dict.fn;

t = toc;
disp(['Doing dictionary match complete. Elapsed time is ' num2str(t) ' s.']);

