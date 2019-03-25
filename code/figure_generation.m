
clear, close all, clc;

addpath('../contrib/arrow');
addpath('../contrib/export_fig-XXX');

dir_in = '../data_out/';
dir_out = '../figures/';

font_sz = 24;
line_w = 3;
marker_sz = 10;
dpi_img = 300;
dpi_img2 = 1200;

%% simulation of T1/T2 bias with FSF

n_buff = 40; % space between T1 and T2 row

T1_v = [ 500, 800, 1200, 1600, 2250];
T2_v = [  30,  30,   50,  100,  100];

fn_sim = 'T1T2_bias_with_FSF.mat';
fn_sim2 = 'T1T2_bias_with_FSF_TRvar.mat';
fn_sim3 = 'T1T2_bias_with_FSF_TEvar.mat';

load([dir_in fn_sim])

idx_logical_v = ismember(T1T2_ref_ms(:,:,1),[T1_v(:) T2_v(:)],'rows');
idx_T1T2_v = find(idx_logical_v);
n_idx = numel( idx_T1T2_v );

T1_bias = squeeze( T1T2_abs_bias_ms(idx_T1T2_v,1,:) )';
T2_bias = squeeze( T1T2_abs_bias_ms(idx_T1T2_v,2,:) )';

load([dir_in fn_sim2])

T1_bias2 = squeeze( T1T2_abs_bias_ms(idx_T1T2_v,1,:) )';
T2_bias2 = squeeze( T1T2_abs_bias_ms(idx_T1T2_v,2,:) )';

load([dir_in fn_sim3])

T1_bias3 = squeeze( T1T2_abs_bias_ms(idx_T1T2_v,1,:) )';
T2_bias3 = squeeze( T1T2_abs_bias_ms(idx_T1T2_v,2,:) )';

legend_txt = strsplit( sprintf( '%d/%d ', [T1_v(:) T2_v(:)]' ) );
% legend_txt = legend_txt(1:end-1);
legend_txt = [ legend_txt(1:end-1) ];


figure(1); clf;
set(gcf,'Position',[100 200 1000 900],'color','w')
plot( repmat( fsf_q_v(:) , [1 n_idx] ), T1_bias3, 's:', 'LineWidth', line_w)
ylim([-2250 3000])
ax = gca;
ax.ColorOrderIndex = 1;
set(gca,'FontSize',font_sz,'LineWidth',line_w)
xlabel('FSF');
ylabel('T_{1} bias (ms)')
title('MRF-varTE with Fat-Water Separation')
tix=get(gca,'xtick')';
set(gca,'xticklabel',num2str(tix,'%.1f'))
my_letter = text(0.05,2500,'a','FontSize',font_sz*2,'FontWeight','bold');
[my_leg,my_obj] = legend( legend_txt, 'Location', 'northeast' );
for ii = n_idx+1:2:n_idx+2*n_idx
    my_obj(ii).LineStyle = '-';
end
my_leg.Position(1) = 0.66; my_leg.Position(2) = 0.55;
my_txt = text(0,0,'T_{1}/T_{2} (ms)','FontSize',font_sz);
my_txt.Position = [0.70, 2200, 0];
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T1a = imread('tmp.png');
delete('tmp.png');

figure(2); clf;
set(gcf,'Position',[100 200 1000 900],'color','w')
plot( repmat( fsf_q_v(:) , [1 n_idx] ), T1_bias, 's--', 'LineWidth', line_w)
ylim([-2250 3000])
ax = gca;
ax.ColorOrderIndex = 1;
set(gca,'FontSize',font_sz,'LineWidth',line_w)
xlabel('FSF');
ylabel('T_{1} bias (ms)')
title('MRF-fixTE')
tix=get(gca,'xtick')';
set(gca,'xticklabel',num2str(tix,'%.1f'))
my_letter = text(0.05,2500,'b','FontSize',font_sz*2,'FontWeight','bold');
[my_leg,my_obj] = legend( legend_txt, 'Location', 'northeast' );
for ii = n_idx+1:2:n_idx+2*n_idx
    my_obj(ii).LineStyle = '-';
end
my_leg.Position(1) = 0.66; my_leg.Position(2) = 0.55;
my_txt = text(0,0,'T_{1}/T_{2} (ms)','FontSize',font_sz);
my_txt.Position = [0.70, 2200, 0];
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T1b = imread('tmp.png');
delete('tmp.png');

figure(3); clf;
set(gcf,'Position',[100 200 1000 900],'color','w')
plot( repmat( fsf_q_v(:) , [1 n_idx] ), T1_bias2, 's-.', 'LineWidth', line_w)
ylim([-2250 3000])
ax = gca;
ax.ColorOrderIndex = 1;
set(gca,'FontSize',font_sz,'LineWidth',line_w)
xlabel('FSF');
ylabel('T_{1} bias (ms)')
title('MRF-varTR')
tix=get(gca,'xtick')';
set(gca,'xticklabel',num2str(tix,'%.1f'))
my_letter = text(0.05,2500,'c','FontSize',font_sz*2,'FontWeight','bold');
[my_leg,my_obj] = legend( legend_txt, 'Location', 'northeast' );
for ii = n_idx+1:2:n_idx+2*n_idx
    my_obj(ii).LineStyle = '-';
end
my_leg.Position(1) = 0.66; my_leg.Position(2) = 0.55;
my_txt = text(0,0,'T_{1}/T_{2} (ms)','FontSize',font_sz);
my_txt.Position = [0.70, 2200, 0];
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T1c = imread('tmp.png');
delete('tmp.png');


figure(4); clf;
set(gcf,'Position',[100 200 1000 900],'color','w')
plot( repmat( fsf_q_v(:) , [1 n_idx] ), T2_bias3, 's:', 'LineWidth', line_w)
ylim([-100 100]) 
ax = gca;
ax.ColorOrderIndex = 1;
set(gca,'FontSize',font_sz,'LineWidth',line_w)
xlabel('FSF');
ylabel('T_{2} bias (ms)')
title('MRF-varTE with Fat-Water Separation')
tix=get(gca,'xtick')';
set(gca,'xticklabel',num2str(tix,'%.1f'))
my_letter = text(0.05,75,'d','FontSize',font_sz*2,'FontWeight','bold');
[my_leg,my_obj] = legend( legend_txt, 'Location', 'northeast' );
for ii = n_idx+1:2:n_idx+2*n_idx
    my_obj(ii).LineStyle = '-';
end
my_leg.Position(1) = 0.66; my_leg.Position(2) = 0.55;
my_txt = text(0,0,'T_{1}/T_{2} (ms)','FontSize',font_sz);
my_txt.Position = [0.70, 70, 0];
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T2d = imread('tmp.png');
delete('tmp.png');

figure(5); clf;
set(gcf,'Position',[100 200 1000 900],'color','w')
plot( repmat( fsf_q_v(:) , [1 n_idx] ), T2_bias, 's--', 'LineWidth', line_w)
ylim([-100 100]) 
ax = gca;
ax.ColorOrderIndex = 1;
set(gca,'FontSize',font_sz,'LineWidth',line_w)
xlabel('FSF');
ylabel('T_{2} bias (ms)')
title('MRF-fixTE')
tix=get(gca,'xtick')';
set(gca,'xticklabel',num2str(tix,'%.1f'))
my_letter = text(0.05,75,'e','FontSize',font_sz*2,'FontWeight','bold');
[my_leg,my_obj] = legend( legend_txt, 'Location', 'northeast' );
for ii = n_idx+1:2:n_idx+2*n_idx
    my_obj(ii).LineStyle = '-';
end
my_leg.Position(1) = 0.25; my_leg.Position(2) = 0.15;
my_txt = text(0,0,'T_{1}/T_{2} (ms)','FontSize',font_sz);
my_txt.Position = [0.15, -30, 0];
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T2e = imread('tmp.png');
delete('tmp.png');

figure(6); clf;
set(gcf,'Position',[100 200 1000 900],'color','w')
plot( repmat( fsf_q_v(:) , [1 n_idx] ), T2_bias2, 's-.', 'LineWidth', line_w)
ylim([-100 100]) 
ax = gca;
ax.ColorOrderIndex = 1;
set(gca,'FontSize',font_sz,'LineWidth',line_w)
xlabel('FSF');
ylabel('T_{2} bias (ms)')
title('MRF-varTR')
tix=get(gca,'xtick')';
set(gca,'xticklabel',num2str(tix,'%.1f'))
my_letter = text(0.05,75,'f','FontSize',font_sz*2,'FontWeight','bold');
[my_leg,my_obj] = legend( legend_txt, 'Location', 'northeast' );
for ii = n_idx+1:2:n_idx+2*n_idx
    my_obj(ii).LineStyle = '-';
end
my_leg.Position(1) = 0.25; my_leg.Position(2) = 0.15;
my_txt = text(0,0,'T_{1}/T_{2} (ms)','FontSize',font_sz);
my_txt.Position = [0.15, -30, 0];
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T2f = imread('tmp.png');
delete('tmp.png');

T1c = cat(1, T1c, 255*ones( size(T1a,1) - size(T1c,1), size(T1c,2), 3 ) );
T2f = cat(1, T2f, 255*ones( size(T2d,1) - size(T2f,1), size(T2f,2), 3 ) );

T1 = cat(2,T1a,T1b,T1c);
T1T2_sim_bias = cat( 1, T1, 255*ones(n_buff,size(T1,2),3), cat(2,T2d,T2e,T2f) );
figure(7); clf;
imshow( T1T2_sim_bias );
eval( sprintf( 'export_fig %sFigure2.png -r%d',dir_out,3*dpi_img) );

%% simulated image parameter maps

close all;

fn_sim_stats = 'img_simulation_proc_stats.mat';

load( [dir_in fn_sim_stats] );

T1_truth = squeeze( example_maps(:,:,1,1,2) );
T2_truth = squeeze( example_maps(:,:,2,1,2) );
FSF_truth = squeeze( example_maps(:,:,3,1,2) );
B0_truth = squeeze( example_maps(:,:,4,1,2) );

T1_maps = squeeze( example_maps(:,:,1,:,1) );
T2_maps = squeeze( example_maps(:,:,2,:,1) );
FSF_maps = squeeze( example_maps(:,:,3,:,1) );
B0_maps = squeeze( example_maps(:,:,4,:,1) );

N = size( T1_maps, 1 );

T1_maps = reshape( permute( T1_maps, [1 3 2] ), [N*size(T1_maps,3) N] );
T2_maps = reshape( permute( T2_maps, [1 3 2] ), [N*size(T2_maps,3) N] );
FSF_maps = reshape( permute( FSF_maps, [1 3 2] ), [N*size(FSF_maps,3) N] );
B0_maps = reshape( permute( B0_maps, [1 3 2] ), [N*size(B0_maps,3) N] );

T1_maps = [T1_truth; T1_maps; zeros(N/4,N)];
T2_maps = [T2_truth; T2_maps; zeros(N/4,N)];
FSF_maps = [FSF_truth; FSF_maps; zeros(N/4,N)];
B0_maps = [B0_truth; B0_maps; zeros(N/4,N)];

figure(1); clf;
set(gcf,'Position',[100 100 800 1000],'color','w');
imagesc( T1_maps ); axis image; axis off; caxis([0 2250])
c = colorbar('Location','south');
c.Color = 'w';
c.FontSize = font_sz/2;
c.LineWidth = line_w;
c.Position(2) = .125;
c.Label.String = '(ms)';
c.Label.Position(2) = -0.8;
colormap(hot);
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T1sim = imread('tmp.png');
delete('tmp.png');

figure(2); clf;
set(gcf,'Position',[100 100 800 1000],'color','w');
imagesc( T2_maps ); axis image; axis off; caxis([0 125])
c = colorbar('Location','south');
c.Color = 'w';
c.FontSize = font_sz/2;
c.LineWidth = line_w;
c.Position(2) = .125;
c.Label.String = '(ms)';
c.Label.Position(2) = -0.8;
colormap(hot);
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T2sim = imread('tmp.png');
delete('tmp.png');

figure(3); clf;
set(gcf,'Position',[100 100 800 1000],'color','w');
imagesc( FSF_maps ); axis image; axis off; caxis([0 1])
c = colorbar('Location','south');
c.Color = 'w';
c.FontSize = font_sz/2;
c.LineWidth = line_w;
c.Position(2) = .5;
c.Label.String = '(FSF)';
c.Label.Position(2) = -1.0;
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
FSFsim = imread('tmp.png');
delete('tmp.png');

figure(4); clf;
set(gcf,'Position',[100 100 800 1000],'color','w');
imagesc( B0_maps ); axis image; axis off; caxis([-250 250])
c = colorbar('Location','south');
c.Color = 'w';
c.FontSize = font_sz/2;
c.LineWidth = line_w;
c.Position(2) = .5;
c.Label.String = '(Hz)';
c.Label.Position(2) = -1.0;
colormap(gray)
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
B0sim = imread('tmp.png');
delete('tmp.png');

B0sim = repmat( B0sim, [1 1 3] );

n_buff = 400;
n_mini_buff = 200;

Tsims = [ T1sim T2sim ];
Tsims =  [ 255*ones(n_mini_buff,size(Tsims,2),3); Tsims; 255*ones(n_buff,size(Tsims,2),3) ];
Tsims = circshift( Tsims, round(3*n_buff/4), 1 );
FBsims = [ FSFsim B0sim ];
FBsims = [ 255*ones(n_mini_buff,size(FBsims,2),3); FBsims; 255*ones(n_buff,size(FBsims,2),3) ];
FBsims = circshift( FBsims, -round(n_buff/5), 1 );
sim_maps = [ Tsims 255*ones(size(Tsims,1),n_buff,3) FBsims 255*ones(size(Tsims,1),n_mini_buff,3) ];
sim_maps = [ 255*ones(size(sim_maps,1),n_mini_buff,3) sim_maps ];

my_scale = 0.85;

figure(5);
imshow(sim_maps);
txt_a = text(100,900,'True','FontSize',font_sz*my_scale,'Rotation',90);
txt_b = text(100,1650,'MRF-varTE','FontSize',font_sz*my_scale,'Rotation',90);
txt_c = text(100,2250,'MRF-fixTE','FontSize',font_sz*my_scale,'Rotation',90);
txt_d = text(100,2900,'MRF-varTR','FontSize',font_sz*my_scale,'Rotation',90);
txt_1 = text(450,400,'T1','FontSize',font_sz*my_scale);
txt_2 = text(1050,400,'T2','FontSize',font_sz*my_scale);
txt_3 = text(2000,50,'FSF','FontSize',font_sz*my_scale);
txt_4 = text(2650,50,'B0','FontSize',font_sz*my_scale);
txt_e = text(1700,600,'True','FontSize',font_sz*my_scale,'Rotation',90);
txt_f = text(1700,1400,'MRF-varTE','FontSize',font_sz*my_scale,'Rotation',90);
eval( sprintf( 'export_fig Fig_simmaps_temp.png -r%d',dpi_img) );


% img sim stats plots

close all

orange = [0.9100 0.4100 0.1700];

legend_txt = {'MRF-varTE w/ fat sep','MRF-fixTE w/o fat sep','MRF-varTR w/o fat sep'};
cats_snr = categorical( {'28','32','38','\infty'} );

% simulation_means [ n_ROI+1, n_metrics, n_noise, n_sims ]
T1sim_bias = squeeze( simulation_bias_means(end,1,:,:) );
T1sim_bias = circshift( T1sim_bias, -1,1 ); %rearrange lowest to highest SNR
T2sim_bias = squeeze( simulation_bias_means(end,2,:,:) );
T2sim_bias = circshift( T2sim_bias, -1,1 ); %rearrange lowest to highest SNR
FSFsim_bias = squeeze( simulation_bias_means(end,3,:,:) );
FSFsim_bias = circshift( FSFsim_bias, -1,1 ); %rearrange lowest to highest SNR
B0sim_bias = squeeze( simulation_bias_means(end,4,:,:) );
B0sim_bias = circshift( B0sim_bias, -1,1 ); %rearrange lowest to highest SNR

T1sim_bias_std = squeeze( simulation_bias_std(end,1,:,:) );
T1sim_bias_std = circshift( T1sim_bias_std, -1,1 ); %rearrange lowest to highest SNR
T2sim_bias_std = squeeze( simulation_bias_std(end,2,:,:) );
T2sim_bias_std = circshift( T2sim_bias_std, -1,1 ); %rearrange lowest to highest SNR
FSFsim_bias_std = squeeze( simulation_bias_std(end,3,:,:) );
FSFsim_bias_std = circshift( FSFsim_bias_std, -1,1 ); %rearrange lowest to highest SNR
B0sim_bias_std = squeeze( simulation_bias_std(end,4,:,:) );
B0sim_bias_std = circshift( B0sim_bias_std, -1,1 ); %rearrange lowest to highest SNR


figure(10); clf;
set(gcf,'color','w','Position',[100 200 800 800]);
my_hand = bar(cats_snr,T1sim_bias,'LineWidth',line_w/2);
my_hand(1).FaceColor = orange;
my_hand(2).FaceColor = 'b';
my_hand(3).FaceColor = 'm';
ylim([-175 50]);
legend(legend_txt)
xlabel('SNR (dB)')
ylabel('(ms)')
title('T1 mean bias')
set(gca,'FontSize',font_sz,'LineWidth',line_w)
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T1simbias_img = imread('tmp.png');
delete('tmp.png');


figure(11); clf;
set(gcf,'color','w','Position',[100 200 800 800]);
my_hand = bar(cats_snr,T2sim_bias,'LineWidth',line_w/2);
my_hand(1).FaceColor = orange;
my_hand(2).FaceColor = 'b';
my_hand(3).FaceColor = 'm';
ylim([-3 5]);
legend(legend_txt)
xlabel('SNR (dB)')
ylabel('(ms)')
title('T2 mean bias')
set(gca,'FontSize',font_sz,'LineWidth',line_w)
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T2simbias_img = imread('tmp.png');
delete('tmp.png');

figure(12); clf;
set(gcf,'color','w','Position',[100 200 800 800]);
my_hand = bar(cats_snr,T1sim_bias_std,'LineWidth',line_w/2);
my_hand(1).FaceColor = orange;
my_hand(2).FaceColor = 'b';
my_hand(3).FaceColor = 'm';
ylim([0 350]);
legend(legend_txt)
xlabel('SNR (dB)')
ylabel('(ms)')
title('T1 standard deviation of bias')
set(gca,'FontSize',font_sz,'LineWidth',line_w)
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T1simbiasstd_img = imread('tmp.png');
delete('tmp.png');

figure(13); clf;
set(gcf,'color','w','Position',[100 200 800 800]);
my_hand = bar(cats_snr,T2sim_bias_std,'LineWidth',line_w/2);
my_hand(1).FaceColor = orange;
my_hand(2).FaceColor = 'b';
my_hand(3).FaceColor = 'm';
ylim([0 12.5]);
legend(legend_txt)
xlabel('SNR (dB)')
ylabel('(ms)')
title('T2 standard deviation of bias')
set(gca,'FontSize',font_sz,'LineWidth',line_w)
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T2simbiasstd_img = imread('tmp.png');
delete('tmp.png');

n_mini_buff = 20;

Mean_bias_plots = [T1simbias_img 255*ones(size(T1simbias_img,1),n_mini_buff,3) T2simbias_img];
Std_bias_plots = [T1simbiasstd_img 255*ones(size(T1simbiasstd_img,1),n_mini_buff,3) T2simbiasstd_img];

my_buff = 255*ones( size(Mean_bias_plots,1), size(Std_bias_plots,2) - size(Mean_bias_plots,2), 3 );

Sim_bias_plots = [ [ my_buff(:,1:floor(end/2),:) Mean_bias_plots my_buff(:,floor(end/2)+1:end,:) ]; ...
    255*ones( n_mini_buff, size( Std_bias_plots, 2 ), 3 );
    Std_bias_plots];


Sim_bias_plots = [ 255*ones( 150,size(Sim_bias_plots,2),3 ); Sim_bias_plots ];
Sim_bias_plots = [ Sim_bias_plots 255*ones( size(Sim_bias_plots,1),150,3 ) ];

figure(20); clf;
imshow( Sim_bias_plots );
eval( sprintf( 'export_fig tmp.png -r%d',3*dpi_img) );
simbiasplot_rgb = imread('tmp.png');
delete('tmp.png');

simmaps_rgb = imread('Fig_simmaps_temp.png');

[my_row, my_col, ~] = size( simmaps_rgb );
my_row_rsz = my_row - 1190;
my_col_rsz = my_col - 1150;
simbiasplot_rgb = imresize( simbiasplot_rgb, [my_row_rsz my_col_rsz]);

sim_maps_and_plot = simmaps_rgb;
sim_maps_and_plot(end-my_row_rsz+1:end,end-my_col_rsz+1:end,:) = simbiasplot_rgb;

figure(22); clf;
imshow( sim_maps_and_plot );
text(100,280,'a','FontSize',1.5*font_sz,'Color','k')
text(1250,100,'b','FontSize',1.5*font_sz,'Color','k')
text(1125,1150,'c','FontSize',1.5*font_sz,'Color','k')

eval( sprintf( 'export_fig %sFigure3.png -r%d',dir_out,dpi_img) );


%% fat-water layer phantom data

close all;

fn_conv_proc = 'layer_conventional_maps.mat';
fn_layer_analysis = 'layer_FSF_T1T2_analysis.mat';
fn_MRF_varTE_ref = 'MRF_layer_varTE_1_proc_fat_sep.mat';
fn_MRF_fixTE_ref = 'MRF_layer_fixTE_1_proc_no_fat_sep.mat';
fn_MRF_varTR_ref = 'MRF_layer_varTR_1_proc_no_fat_sep.mat';
fn_MRF_varTE_eg = 'MRF_layer_varTE_6_proc_fat_sep.mat';
fn_MRF_fixTE_eg = 'MRF_layer_fixTE_6_proc_no_fat_sep.mat';
fn_MRF_varTR_eg = 'MRF_layer_varTR_6_proc_no_fat_sep.mat';

n_sl_eg = 6;

load([dir_in fn_conv_proc]);
FSF_conv = conventional_maps.FSF_graphcut(:,:,n_sl_eg);
B0_conv = conventional_maps.B0_graphcut(:,:,n_sl_eg);
load([dir_in fn_layer_analysis]);

load([dir_in fn_MRF_varTE_ref]);
T1_varTE_ref = output_MRFAT.T1_water_map;
T2_varTE_ref = output_MRFAT.T2_water_map;

load([dir_in fn_MRF_fixTE_ref])
T1_fixTE_ref = output_MRF_match.T1_map;
T2_fixTE_ref = output_MRF_match.T2_map;

load([dir_in fn_MRF_varTR_ref])
T1_varTR_ref = output_MRF_match.T1_map;
T2_varTR_ref = output_MRF_match.T2_map;

load([dir_in fn_MRF_varTE_eg]);
FSF_varTE = output_MRFAT.FSF_map;
B0_varTE = output_MRFAT.B0_fit_map;
T1_varTE = output_MRFAT.T1_water_map;
T2_varTE = output_MRFAT.T2_water_map;

load([dir_in fn_MRF_fixTE_eg])
T1_fixTE = output_MRF_match.T1_map;
T2_fixTE = output_MRF_match.T2_map;

load([dir_in fn_MRF_varTR_eg])
T1_varTR = output_MRF_match.T1_map;
T2_varTR = output_MRF_match.T2_map;

% mask example slices
mask1 = conventional_maps.img_conventional(:,:,1);
my_thresh = multithresh( mask1, 1 );
mask1( mask1 < my_thresh(1) ) = 0;
mask1( mask1 > 0 ) = 1;

mask2 = conventional_maps.img_conventional(:,:,n_sl_eg);
my_thresh = multithresh( mask2, 1 );
mask2( mask2 < my_thresh(1) ) = 0;
mask2( mask2 > 0 ) = 1;

T1_varTE_ref = T1_varTE_ref.*mask1;
T2_varTE_ref = T2_varTE_ref.*mask1;
T1_fixTE_ref = T1_fixTE_ref.*mask1;
T2_fixTE_ref = T2_fixTE_ref.*mask1;
T1_varTR_ref = T1_varTR_ref.*mask1;
T2_varTR_ref = T2_varTR_ref.*mask1;

FSF_conv = FSF_conv.*mask2;
B0_conv = B0_conv.*mask2;
FSF_varTE = FSF_varTE.*mask2;
B0_varTE = B0_varTE.*mask2;
T1_varTE = T1_varTE.*mask2;
T2_varTE = T2_varTE.*mask2;
T1_fixTE = T1_fixTE.*mask2;
T2_fixTE = T2_fixTE.*mask2;
T1_varTR = T1_varTR.*mask2;
T2_varTR = T2_varTR.*mask2;

T1_variation = squeeze( T1_means_MRF(:,1,:) ) - T1_consensus;
T2_variation = squeeze( T2_means_MRF(:,1,:) ) - T2_consensus;

fprintf( 'Max. abs MRF T1 variation from average %.1f ms\n', max( abs( T1_variation(:) ) ) );
fprintf( 'Max. abs MRF T2 variation from average %.1f ms\n', max( abs( T2_variation(:) ) ) );

% example T1 and T2 maps for all MRF sequences
figure(5); clf;
set(gcf,'Position',[50 200 1200 900],'color','w');
imagesc( [T1_varTE_ref T1_fixTE_ref T1_varTR_ref ; T1_varTE T1_fixTE T1_varTR] ); colormap(hot)
axis image; axis off; caxis([0 1750]);
c = colorbar('Location','south');
c.Color = 'w';
c.FontSize = font_sz*.8;
c.LineWidth = line_w;
c.Position(1) = 0.25;
c.Position(3) = 0.6;
c.Position(4) = 0.02;
c.Label.Position(1) = -150;
c.Label.Position(2) = -0.1;
c.Label.String = 'T1 (ms)';
my_letter = text(20,25,'a','FontSize',font_sz,'FontWeight','bold','Color','w');
my_letter = text(20+240,25,'b','FontSize',font_sz,'FontWeight','bold','Color','w');
my_letter = text(20+480,25,'c','FontSize',font_sz,'FontWeight','bold','Color','w');
my_letter = text(20,25+240,'d','FontSize',font_sz,'FontWeight','bold','Color','w');
my_letter = text(20+240,25+240,'e','FontSize',font_sz,'FontWeight','bold','Color','w');
my_letter = text(20+480,25+240,'f','FontSize',font_sz,'FontWeight','bold','Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T1_rgb = imread('tmp.png');

figure(6); clf;
set(gcf,'Position',[50 200 1200 900],'color','w');
imagesc( [T2_varTE_ref T2_fixTE_ref T2_varTR_ref ; T2_varTE T2_fixTE T2_varTR] ); colormap(hot)
axis image; axis off; caxis([0 200]);
c = colorbar('Location','south');
c.Color = 'w';
c.FontSize = font_sz*.8;
c.LineWidth = line_w;
c.Position(1) = 0.25;
c.Position(3) = 0.6;
c.Position(4) = 0.02;
c.Label.Position(1) = -18;
c.Label.Position(2) = -0.1;
c.Label.String = 'T2 (ms)';
my_letter = text(20,25,'g','FontSize',font_sz,'FontWeight','bold','Color','w');
my_letter = text(20+240,25,'h','FontSize',font_sz,'FontWeight','bold','Color','w');
my_letter = text(20+480,25,'i','FontSize',font_sz,'FontWeight','bold','Color','w');
my_letter = text(20,25+240,'j','FontSize',font_sz,'FontWeight','bold','Color','w');
my_letter = text(20+240,25+240,'k','FontSize',font_sz,'FontWeight','bold','Color','w');
my_letter = text(20+480,25+240,'l','FontSize',font_sz,'FontWeight','bold','Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T2_rgb = imread('tmp.png');

my_sp = 1000;
my_sp2 = 2000;
scale_font = 0.6;
figure(7); clf;
my_img = [T1_rgb; 255*ones( 40, size(T1_rgb,2), 3); T2_rgb];
[n_row,n_col,~] = size( my_img );
my_buff1 = 255*ones(n_row+250,250,3);
my_buff2 = 255*ones(250,n_col,3);
imshow( [my_buff1 [my_buff2; my_img]] );
txt1 = text(350,150,'MRF-varTE','FontSize',font_sz*scale_font,'FontWeight','bold','Color','k');
txt2 = text(350+my_sp,150,'MRF-fixTE','FontSize',font_sz*scale_font,'FontWeight','bold','Color','k');
txt3 = text(350+2*my_sp,150,'MRF-varTR','FontSize',font_sz*scale_font,'FontWeight','bold','Color','k');
txt4 = text(125,1000,'water only','FontSize',font_sz*scale_font,'FontWeight','bold','Color','k','Rotation',90);
txt5 = text(125,1000+my_sp2,'water only','FontSize',font_sz*scale_font,'FontWeight','bold','Color','k','Rotation',90);
txt6 = text(125,2000,'water + oil','FontSize',font_sz*scale_font,'FontWeight','bold','Color','k','Rotation',90);
txt7 = text(125,2000+my_sp2,'water + oil','FontSize',font_sz*scale_font,'FontWeight','bold','Color','k','Rotation',90);
eval( sprintf( 'export_fig %sFigure4.png -r%d',dir_out,dpi_img) );

% plot of T1 and T2 deviation from consensus at 0 FF vs conventional FSF
T1_bias_varTE = T1_bias_MRF(:,:,1);
T1_bias_fixTE = T1_bias_MRF(:,:,2);
T1_bias_varTR = T1_bias_MRF(:,:,3);
figure(8); clf;
set(gcf,'Position',[100 200 1000 800],'color','w');
orange = [0.9100 0.4100 0.1700];
plot( FSF_means_conv(:), T1_bias_varTE(:),'v', 'MarkerEdgeColor', orange, 'MarkerFaceColor', orange, 'LineWidth', line_w, 'MarkerSize', marker_sz );
hold on
plot( FSF_means_conv(:), T1_bias_fixTE(:),'bo','MarkerFaceColor', 'b', 'LineWidth', line_w, 'MarkerSize', marker_sz );
plot( FSF_means_conv(:), T1_bias_varTR(:),'ks','MarkerFaceColor', 'm', 'LineWidth', line_w, 'MarkerSize', marker_sz );
hold off
set(gca,'FontSize',font_sz,'LineWidth',line_w)
line([0 0.8],[0 0],'LineStyle','--','Color','r','LineWidth',line_w)
tix=get(gca,'xtick')';
set(gca,'xticklabel',num2str(tix,'%.1f'))
xlabel('Conventional FSF');
ylabel('Difference from water T1 (ms)')
legend('MRF-varTE w/ fat sep','MRF-fixTE w/o fat sep','MRF-varTR w/o fat sep', 'Location', 'northeast' )
my_letter = text(0.05,1250,'a','FontSize',font_sz*2,'FontWeight','bold','Color','k');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T1_bias_rgb = imread('tmp.png');

T2_bias_varTE = T2_bias_MRF(:,:,1);
T2_bias_fixTE = T2_bias_MRF(:,:,2);
T2_bias_varTR = T2_bias_MRF(:,:,3);
figure(9); clf;
set(gcf,'Position',[100 200 1000 800],'color','w');
plot( FSF_means_conv(:), T2_bias_varTE(:),'v', 'MarkerEdgeColor', orange, 'MarkerFaceColor', orange, 'LineWidth', line_w, 'MarkerSize', marker_sz );
hold on
plot( FSF_means_conv(:), T2_bias_fixTE(:),'bo','MarkerFaceColor', 'b', 'LineWidth', line_w, 'MarkerSize', marker_sz );
plot( FSF_means_conv(:), T2_bias_varTR(:),'ks','MarkerFaceColor', 'm', 'LineWidth', line_w, 'MarkerSize', marker_sz );
hold off
set(gca,'FontSize',font_sz,'LineWidth',line_w)
line([0 0.8],[0 0],'LineStyle','--','Color','r','LineWidth',line_w)
tix=get(gca,'xtick')';
set(gca,'xticklabel',num2str(tix,'%.1f'))
xlabel('Conventional FSF');
ylabel('Difference from water T2 (ms)')
legend('MRF-varTE w/ fat sep','MRF-fixTE w/o fat sep','MRF-varTR w/o fat sep', 'Location', 'northeast' );
my_letter = text(0.05,125,'b','FontSize',font_sz*2,'FontWeight','bold','Color','k');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T2_bias_rgb = imread('tmp.png');

my_buff = 255*ones( size(T2_bias_rgb,1), size(T1_bias_rgb,2) - size(T2_bias_rgb,2), 3 );
figure(10); clf;
imshow( [T1_bias_rgb; [T2_bias_rgb my_buff]] );
eval( sprintf( 'export_fig %sFigure5.png -r%d',dir_out,dpi_img) );

% example conventional and MRF FSF and B0
bg = zeros( size( [FSF_conv FSF_varTE] ) );
bg(~logical([mask2 mask2])) = 1;

figure(1); clf;
set(gcf,'Position',[100 200 1000 800],'color','w');
imagesc( [FSF_conv FSF_varTE] ); axis image; axis off; caxis([0 0.8])
c = colorbar('Location','south');
c.Color = 'w';
c.FontSize = font_sz;
c.LineWidth = line_w;
c.Position(1) = 0.30;
c.Position(3) = 0.45;
c.Position(4) = 0.03;
% c.Ticks = [0:0.2:1.0];
c.Label.String = 'FSF';
c.TickLabels{1} = num2str( str2num( c.TickLabels{1} ), '%.1f' );
black = cat(3, zeros(size(bg)), zeros(size(bg)), zeros(size(bg)));
hold on
bg_h = imshow(black);
hold off
set( bg_h, 'AlphaData', bg );
my_letter = text(20,25,'a','FontSize',font_sz*2,'FontWeight','bold','Color','w');
my_letter = text(20+240,25,'b','FontSize',font_sz*2,'FontWeight','bold','Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
FSF_fig_eg = imread('tmp.png');

figure(2); clf;
set(gcf,'Position',[100 200 1000 800],'color','w');
imagesc( [B0_conv B0_varTE] ); axis image; axis off; colormap(gray); caxis([-150 150])
c = colorbar('Location','south');
c.Color = 'w';
c.FontSize = font_sz;
c.LineWidth = line_w;
c.Position(1) = 0.30;
c.Position(3) = 0.45;
c.Position(4) = 0.03;
c.Ticks = [-150:50:150];
c.Label.String = 'B0 (Hz)';
my_letter = text(20,25,'c','FontSize',font_sz*2,'FontWeight','bold','Color','w');
my_letter = text(20+240,25,'d','FontSize',font_sz*2,'FontWeight','bold','Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
B0_fig_eg = imread('tmp.png');
B0_fig_eg = repmat( B0_fig_eg, [1 1 3] );

figure(3); clf;
set(gcf,'Position',[100 200 1000 800],'color','w');
plot( FSF_means_conv(:), FSF_means_MRF(:), 'v', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', orange, 'LineWidth', line_w, 'MarkerSize', marker_sz); axis square;
set(gca,'FontSize',font_sz,'LineWidth',line_w)
xlim([0 0.9]); ylim([0 0.9]);
line([0 0.9],[0 0.9],'LineStyle','--','Color','r','LineWidth',line_w)
tix=get(gca,'xtick')';
set(gca,'ytick',tix);
set(gca,'xticklabel',num2str(tix,'%.1f'),'yticklabel',num2str(tix,'%.1f'))
xlabel('Conventional FSF');
ylabel('MRF FSF')
my_txt = text( 0.06, 0.70, sprintf( 'CCC = %.3f',CCC_FSF ), 'FontSize', font_sz);
my_letter = text(0.03,0.85,'e','FontSize',font_sz*2,'FontWeight','bold','Color','k');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
CCC_layer = imread('tmp.png');
CCC_layer = imresize( CCC_layer, 2*size(B0_fig_eg,1)/size(CCC_layer,1) );

figure(4); clf;
imshow([ [FSF_fig_eg; 255*ones(20,size(FSF_fig_eg,2),3); B0_fig_eg] [255*ones(size(CCC_layer,1),20,3) CCC_layer; 255*ones(20,size(CCC_layer,2) + 20,3) ]])
eval( sprintf( 'export_fig %sFigure6.png -r%d',dir_out,dpi_img) );

%% MRF B0 experiment

close all;

fn_conv_proc = 'layer_conventional_maps.mat';
fn_layer_analysis = 'layer_FSF_T1T2_analysis.mat';
fn_MRF_b0 = 'MRF_layer_varTE_badshim_proc_fat_sep.mat';
fn_MRF_nob0 = 'MRF_layer_varTE_badshim_proc_fat_sep_noB0.mat';
fn_MRF_ref = 'MRF_layer_varTE_4_proc_fat_sep.mat';

n_sl_eg = 8;

load([dir_in fn_conv_proc]);
FSF_conv = conventional_maps.FSF_graphcut(:,:,n_sl_eg);
B0_conv = conventional_maps.B0_graphcut(:,:,n_sl_eg);
load([dir_in fn_layer_analysis]);

load([dir_in fn_MRF_b0]);
FSF_MRF_b0 = output_MRFAT.FSF_map;
B0_MRF_b0 = output_MRFAT.B0_fit_map;
T1_b0 = output_MRFAT.T1_water_map;
T2_b0 = output_MRFAT.T2_water_map;

load([dir_in fn_MRF_nob0])
FSF_MRF_nob0 = output_MRFAT.FSF_map;
T1_nob0 = output_MRFAT.T1_water_map;
T2_nob0 = output_MRFAT.T2_water_map;

load([dir_in fn_MRF_ref]);
T1_ref = output_MRFAT.T1_water_map;
T2_ref = output_MRFAT.T2_water_map;

% mask example slices
mask1 = conventional_maps.img_conventional(:,:,n_sl_eg);
my_thresh = multithresh( mask1, 1 );
mask1( mask1 < my_thresh(1) ) = 0;
mask1( mask1 > 0 ) = 1;

FSF_conv = FSF_conv.*mask1;
B0_conv = B0_conv.*mask1;
FSF_MRF_b0 = FSF_MRF_b0.*mask1;
B0_MRF_b0 = B0_MRF_b0.*mask1;
FSF_MRF_nob0 = FSF_MRF_nob0.*mask1;
T1_b0 = T1_b0.*mask1;
T2_b0 = T2_b0.*mask1;
T1_nob0 = T1_nob0.*mask1;
T2_nob0 = T2_nob0.*mask1;
T1_ref = T1_ref.*mask1;
T2_ref = T2_ref.*mask1;

% conventional and MRF FSF/B0
bg = zeros( size( [FSF_conv FSF_MRF_nob0 FSF_MRF_b0] ) );
bg(~logical( [mask1 mask1 mask1] ) ) = 1;

figure(1); clf;
set(gcf,'Position',[50 200 1450 1000],'color','w');
imagesc( [FSF_conv FSF_MRF_nob0 FSF_MRF_b0] ); axis image; axis off; caxis([0 0.4])
c = colorbar('Location','south');
c.Color = 'w';
c.FontSize = font_sz;
c.LineWidth = line_w;
c.Position(1) = 0.30;
c.Position(3) = 0.45;
c.Position(4) = 0.03;
% c.Ticks = [0:0.2:1.0];
c.TickLabels{1} = num2str( str2num( c.TickLabels{1} ), '%.1f');
c.Label.String = 'FSF';
c.Label.Position(1) = -0.03;
c.Label.Position(2) = -0.2;
black = cat(3, zeros(size(bg)), zeros(size(bg)), zeros(size(bg)));
hold on
bg_h = imshow(black);
hold off
set( bg_h, 'AlphaData', bg );
my_letter = text(20,25,'a','FontSize',font_sz*2,'FontWeight','bold','Color','w');
my_letter = text(20+240,25,'b','FontSize',font_sz*2,'FontWeight','bold','Color','w');
my_letter = text(20+480,25,'c','FontSize',font_sz*2,'FontWeight','bold','Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
FSF_rgb = imread('tmp.png');

figure(2); clf;
set(gcf,'Position',[50 200 1450 1000],'color','w');
imagesc( [B0_conv zeros(size(B0_MRF_b0)) B0_MRF_b0] ); axis image; axis off; colormap(gray); caxis([-250 250])
c = colorbar('Location','south');
c.Color = 'w';
c.FontSize = font_sz;
c.LineWidth = line_w;
c.Position(1) = 0.30;
c.Position(3) = 0.45;
c.Position(4) = 0.03;
% c.Ticks = [-150:50:150];
c.Label.String = 'B0 (Hz)';
my_letter = text(20,25,'d','FontSize',font_sz*2,'FontWeight','bold','Color','w');
my_letter = text(20+480,25,'e','FontSize',font_sz*2,'FontWeight','bold','Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
B0_rgb = imread('tmp.png');
B0_rgb = repmat( B0_rgb, [1 1 3] );

my_img = [ FSF_rgb; 255*ones(20,size(FSF_rgb,2),3); B0_rgb];
[n_row,n_col,~] = size( my_img );
my_buff1 = 255*ones(250,n_col,3);

figure(3); clf;
imshow([my_buff1; my_img])
txt1 = text(250,150,'SPGR w/ B0 Fit','FontSize',font_sz*scale_font*1.5,'FontWeight','bold','Color','k');
txt1 = text(1400,150,'MRF w/o B0 Fit','FontSize',font_sz*scale_font*1.5,'FontWeight','bold','Color','k');
txt1 = text(2550,150,'MRF w/ B0 Fit','FontSize',font_sz*scale_font*1.5,'FontWeight','bold','Color','k');
eval( sprintf( 'export_fig %sFigure7.png -r%d',dir_out,dpi_img) );

% MRF well shimmed T1, and poorly shimmed MRF w/o and w/B0
figure(4); clf;
set(gcf,'Position',[50 200 1450 1000],'color','w');
imagesc( [T1_ref T1_nob0 T1_b0] ); colormap(hot)
axis image; axis off; caxis([0 1750]);
c = colorbar('Location','south');
c.Color = 'w';
c.FontSize = font_sz;
c.LineWidth = line_w;
c.Position(1) = 0.3;
c.Position(3) = 0.45;
c.Position(4) = 0.03;
c.Label.Position(1) = -200;
c.Label.Position(2) = -0.1;
c.Label.String = 'T1 (ms)';
my_letter = text(20,25,'a','FontSize',font_sz*2,'FontWeight','bold','Color','w');
my_letter = text(20+240,25,'b','FontSize',font_sz*2,'FontWeight','bold','Color','w');
my_letter = text(20+480,25,'c','FontSize',font_sz*2,'FontWeight','bold','Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T1_rgb = imread('tmp.png');

% MRF well shimmed T2, and poorly shimmed MRF w/o and w/B0
figure(5); clf;
set(gcf,'Position',[50 200 1450 1000],'color','w');
imagesc( [T2_ref T2_nob0 T2_b0] ); colormap(hot)
axis image; axis off; caxis([0 250]);
c = colorbar('Location','south');
c.Color = 'w';
c.FontSize = font_sz;
c.LineWidth = line_w;
c.Position(1) = 0.3;
c.Position(3) = 0.45;
c.Position(4) = 0.03;
c.Label.Position(1) = -28;
c.Label.Position(2) = -0.1;
c.Label.String = 'T2 (ms)';
my_letter = text(20,25,'d','FontSize',font_sz*2,'FontWeight','bold','Color','w');
my_letter = text(20+240,25,'e','FontSize',font_sz*2,'FontWeight','bold','Color','w');
my_letter = text(20+480,25,'f','FontSize',font_sz*2,'FontWeight','bold','Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T2_rgb = imread('tmp.png');

my_buff = 255*ones( 20, size([ T1_rgb; T2_rgb],2), 3 );

my_img = [ T1_rgb; my_buff; T2_rgb ];
[n_row,n_col,~] = size( my_img );
my_buff1 = 255*ones(250,n_col,3);

figure(6); clf;
imshow([my_buff1; my_img])
txt1 = text(250,150,'MRF PB Shim','FontSize',font_sz*scale_font*1.5,'FontWeight','bold','Color','k');
txt2 = text(1400,90,'MRF Poor Shim','FontSize',font_sz*scale_font*1.5,'FontWeight','bold','Color','k');
txt2b = text(1500,210,'w/o B0 Fit','FontSize',font_sz*scale_font*1.5,'FontWeight','bold','Color','k');
txt3 = text(2550,90,'MRF Poor Shim','FontSize',font_sz*scale_font*1.5,'FontWeight','bold','Color','k');
txt3b = text(2650,210,'w/ B0 Fit','FontSize',font_sz*scale_font*1.5,'FontWeight','bold','Color','k');
eval( sprintf( 'export_fig %sFigure8.png -r%d',dir_out,dpi_img) );


%% knee data

close all;

trim_cols = 40;

fn_conv_proc = 'knee_conventional_maps.mat';
fn_MRF_varTE = 'MRF_knee_varTE_proc_fat_sep.mat';
fn_MRF_fixTE = 'MRF_knee_fixTE_proc_no_fat_sep.mat';
fn_MRF_varTR = 'MRF_knee_varTR_proc_no_fat_sep.mat';
fn_ROI = 'knee_conventional_ROIs.mat';

load([dir_in fn_ROI]);
img_mask = sum(ROI_stack,3);
img_mask( img_mask > 1 ) = 1;

load([dir_in fn_conv_proc]);
img_mask_water = img_mask;
water_graphcut = conventional_maps.water_graphcut;
img_mask_water( water_graphcut < multithresh(water_graphcut,1) ) = 0;
FSF_graphcut = conventional_maps.FSF_graphcut.*img_mask;
B0_graphcut = conventional_maps.B0_graphcut.*img_mask;
T1_map_ms = conventional_maps.T1_map_ms.*img_mask_water;
T2_map_ms = conventional_maps.T2_map_ms.*img_mask_water;

anatomical_img = conventional_maps.anatomical_img;

load([dir_in fn_MRF_varTE]);
T1_MRF_varTE = img_mask_water.*output_MRFAT.T1_water_map;
T2_MRF_varTE = img_mask_water.*output_MRFAT.T2_water_map;
FSF_MRF_varTE = img_mask.*output_MRFAT.FSF_map;
B0_MRF_varTE = img_mask.*output_MRFAT.B0_fit_map;

load([dir_in fn_MRF_fixTE]);
T1_MRF_fixTE = img_mask_water.*output_MRF_match.T1_map;
T2_MRF_fixTE = img_mask_water.*output_MRF_match.T2_map;

load([dir_in fn_MRF_varTR]);
T1_MRF_varTR = img_mask_water.*output_MRF_match.T1_map;
T2_MRF_varTR = img_mask_water.*output_MRF_match.T2_map;

% 

T1_row = [T1_map_ms(:,1:end-trim_cols) T1_MRF_varTR(:,1:end-trim_cols) T1_MRF_fixTE(:,1:end-trim_cols) T1_MRF_varTE];
T2_row = [T2_map_ms(:,1:end-trim_cols) T2_MRF_varTR(:,1:end-trim_cols) T2_MRF_fixTE(:,1:end-trim_cols) T2_MRF_varTE];

figure(5); clf;
set(gcf,'Position',[100 200 1400 800],'color','w')
imagesc( T1_row ); axis image; caxis([0 1750]); axis off; colormap(hot);
c = colorbar();
c.Color = 'w';
c.FontSize = font_sz;
c.LineWidth = line_w;
c.Position(1) = 0.83;
c.Position(2) = 0.35;
c.Position(4) = 0.35;
c.Label.String = 'T1 (ms)';
c.Label.Position(1) = -2.0;
my_offset = (240-trim_cols);
my_letter = text(20,25,'c','FontSize',font_sz*2,'FontWeight','bold','Color','w');
arrow([327,134],[312,136],'Color','w');
arrow([327-my_offset,134],[312-my_offset,136],'Color','w');
arrow([327+my_offset,134],[312+my_offset,136],'Color','w');
arrow([327+2*my_offset,134],[312+2*my_offset,136],'Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T1_knee_rgb = imread('tmp.png');
delete('tmp.png');

figure(6); clf;
set(gcf,'Position',[100 200 1400 800],'color','w'); colormap(hot);
imagesc( T2_row ); axis image; caxis([0 100]); axis off;
c = colorbar();
c.Color = 'w';
c.FontSize = font_sz;
c.LineWidth = line_w;
c.Position(1) = 0.84;
c.Position(2) = 0.35;
c.Position(4) = 0.34;
c.Label.String = 'T2 (ms)';
c.Label.Position(1) = -2.0;
my_letter = text(20,25,'d','FontSize',font_sz*2,'FontWeight','bold','Color','w');
arrow([327,134],[312,136],'Color','w');
arrow([327-my_offset,134],[312-my_offset,136],'Color','w');
arrow([327+my_offset,134],[312+my_offset,136],'Color','w');
arrow([327+2*my_offset,134],[312+2*my_offset,136],'Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T2_knee_rgb = imread('tmp.png');
delete('tmp.png');

%

fat_row = [FSF_graphcut(:,1:end-trim_cols) FSF_MRF_varTE];
bg = ~logical([img_mask(:,1:end-trim_cols) img_mask]);
n_col_fat = size(fat_row,2);
n_row_anat = size(anatomical_img,1);
n_col_xtra = n_col_fat - size(anatomical_img,2) + trim_cols;
anat_row = [zeros(n_row_anat,n_col_xtra/2) anatomical_img(:,1:end-trim_cols) zeros(n_row_anat,n_col_xtra/2)];

anat_row = [anat_row zeros(size(fat_row,1),size(T1_row,2)-size(fat_row,2))];
my_buff = [zeros(size(fat_row,1),size(T1_row,2)-size(fat_row,2))];
fat_row = [ my_buff fat_row ];
bg = [~logical(my_buff) bg];

figure(1); clf;
set(gcf,'Position',[100 200 1400 800],'color','w')
imagesc( anat_row ); axis image; axis off; colormap(gray);
my_letter = text(20,25,'a','FontSize',font_sz*2,'FontWeight','bold','Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
anat_rgb = imread('tmp.png');
anat_rgb = repmat( anat_rgb, [1 1 3] );

figure(2); clf;
set(gcf,'Position',[100 200 1400 800],'color','w')
imagesc( fat_row ); axis image; axis off; caxis([0 1]);
c = colorbar();
c.Color = 'w';
c.FontSize = font_sz;
c.LineWidth = line_w;
c.Position(1) = 0.84;
c.Position(2) = 0.36;
c.Position(4) = 0.30;
c.TickLabels{1} = num2str( str2num( c.TickLabels{1} ), '%.1f' );
c.TickLabels{end} = num2str( str2num( c.TickLabels{end} ), '%.1f' );
c.Label.String = 'FSF';
c.Label.Position(1) = -2.0;
black = cat(3, zeros(size(bg)), zeros(size(bg)), zeros(size(bg)));
hold on
bg_h = imshow(black);
hold off
set( bg_h, 'AlphaData', bg );
my_letter = text(20 + 350,25,'b','FontSize',font_sz*2,'FontWeight','bold','Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
fsf_rgb = imread('tmp.png');
delete('tmp.png');

n_b_r = 20;
n_row = size(T1_knee_rgb,1);
n_col = size(T1_knee_rgb,2);
my_buff1 = 255*ones(250,n_col,3);
my_buff2 = 255*ones(n_row,n_b_r,3);
my_buff3 = 255*ones(20,n_col,3);

col_div1 = round( n_col*2/5 - n_b_r/2 );
col_div2 = round( n_col*2/5 + n_b_r/2 );

top_row = [anat_rgb(1:size(fsf_rgb,1),1:col_div1,:) my_buff2 fsf_rgb(:,col_div2 + 1:end,:)];
my_img = [my_buff1; top_row; my_buff1; T1_knee_rgb; my_buff3; T2_knee_rgb];

figure(8); clf;
imshow(my_img)

my_sp = 800;
scale_font = 0.9;

txt1 = text(200+2*my_sp,150,'SPGR','FontSize',font_sz*scale_font,'FontWeight','bold','Color','k');
txt2 = text(60+3*my_sp,150,'MRF-varTE','FontSize',font_sz*scale_font,'FontWeight','bold','Color','k');

txt3 = text(60,1350,'IR/MSE','FontSize',font_sz*scale_font,'FontWeight','bold','Color','k');
txt4 = text(60+my_sp,1350,'MRF-fixTE','FontSize',font_sz*scale_font,'FontWeight','bold','Color','k');
txt5 = text(60+2*my_sp,1350,'MRF-varTR','FontSize',font_sz*scale_font,'FontWeight','bold','Color','k');
txt6 = text(60+3*my_sp,1350,'MRF-varTE','FontSize',font_sz*scale_font,'FontWeight','bold','Color','k');

eval( sprintf( 'export_fig %sFigure9.png -r%d',dir_out,dpi_img) );

%% brain

close all

fn_conv_proc = 'brain_conventional_maps.mat';
fn_ROI = 'brain_conventional_ROIs.mat';
fn_MRF_varTE = 'MRF_brain_varTE_proc_fat_sep.mat';
fn_MRF_fixTE = 'MRF_brain_fixTE_proc_no_fat_sep.mat';
fn_MRF_varTR = 'MRF_brain_varTR_proc_no_fat_sep.mat';

load([dir_in fn_ROI]);
img_mask = sum(ROI_stack,3);
img_mask( img_mask > 1 ) = 1;

load([dir_in fn_conv_proc]);
img_mask_water = img_mask;
water_graphcut = conventional_maps.water_graphcut;
img_mask_water( water_graphcut < multithresh(water_graphcut,1) ) = 0;
FSF_graphcut = conventional_maps.FSF_graphcut.*img_mask;
B0_graphcut = conventional_maps.B0_graphcut.*img_mask;

load([dir_in fn_MRF_varTE]);
T1_MRF_varTE = img_mask_water.*output_MRFAT.T1_water_map;
T2_MRF_varTE = img_mask_water.*output_MRFAT.T2_water_map;
FSF_MRF_varTE = img_mask.*output_MRFAT.FSF_map;
B0_MRF_varTE = img_mask.*output_MRFAT.B0_fit_map;

load([dir_in fn_MRF_fixTE]);
T1_MRF_fixTE = img_mask_water.*output_MRF_match.T1_map;
T2_MRF_fixTE = img_mask_water.*output_MRF_match.T2_map;

load([dir_in fn_MRF_varTR]);
T1_MRF_varTR = img_mask_water.*output_MRF_match.T1_map;
T2_MRF_varTR = img_mask_water.*output_MRF_match.T2_map;

n_row_buff = 20;
T1_row = [T1_MRF_fixTE T1_MRF_varTR T1_MRF_varTE];
T1_row = [T1_row; zeros( n_row_buff, size( T1_row, 2 ) ) ]; 
T2_row = [T2_MRF_fixTE T2_MRF_varTR T2_MRF_varTE];
T2_row = [T2_row; zeros( n_row_buff, size( T2_row, 2 ) ) ]; 


n_row = size( B0_graphcut, 1 );
n_buff1 = n_row/2;
my_buff = zeros( n_row, n_buff1 );

x0_b_arrow = 281;
figure(1); clf;
set(gcf,'Position',[100 200 1400 800],'color','w')
imagesc( [my_buff B0_graphcut B0_MRF_varTE my_buff] ); axis off; axis image; colormap( gray ); caxis([-250 250])
c = colorbar('Location','east');
c.Color = 'w';
c.FontSize = font_sz*0.9;
c.LineWidth = line_w;
c.Label.String = 'B0 (Hz)';
my_letter = text(25,40,'a','FontSize',font_sz*2,'FontWeight','bold','Color','w');
arrow([487,61],[496,63],'Color','m');
arrow([478,68],[469,64],'Color','m');
arrow([273,115],[x0_b_arrow,115],'Color','b');
arrow([273+240,115],[x0_b_arrow+240,115],'Color','b');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
B0_brain_rgb = imread('tmp.png');
delete('tmp.png');

fat_row = [my_buff FSF_graphcut FSF_MRF_varTE my_buff];
bg = ~logical([my_buff img_mask img_mask my_buff]);

figure(2); clf;
set(gcf,'Position',[100 200 1400 800],'color','w')
imagesc( fat_row ); axis off; axis image; caxis([0 1])
c = colorbar('Location','east');
c.Color = 'w';
c.FontSize = font_sz*0.9;
c.LineWidth = line_w;
c.TickLabels{1} = num2str( str2num( c.TickLabels{1} ), '%.1f' );
c.TickLabels{end} = num2str( str2num( c.TickLabels{end} ), '%.1f' );
c.Label.String = 'FSF';
black = cat(3, zeros(size(bg)), zeros(size(bg)), zeros(size(bg)));
hold on
bg_h = imshow(black);
hold off
set( bg_h, 'AlphaData', bg );
my_letter = text(25,40,'b','FontSize',font_sz*2,'FontWeight','bold','Color','w');
arrow([487,61],[496,63],'Color','m');
arrow([478,68],[469,64],'Color','m');
arrow([499,182],[490,184],'Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
fsf_brain_rgb = imread('tmp.png');
delete('tmp.png');

figure(3); clf;
set(gcf,'Position',[100 200 1400 800],'color','w')
imagesc( T1_row ); axis off; axis image; colormap( hot ); caxis([0 2500])
c = colorbar('Location','south');
c.Color = 'w';
c.FontSize = font_sz*0.9;
c.LineWidth = line_w;
c.Position(1) = 0.25;
c.Position(2) = 0.29;
c.Position(3) = 0.55;
c.Label.String = 'T1 (ms)';
c.Label.Position(1) = -200;
c.Label.Position(2) = 0;
my_letter = text(25,40,'c','FontSize',font_sz*2,'FontWeight','bold','Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T1_brain_rgb = imread('tmp.png');
delete('tmp.png');

figure(4); clf;
set(gcf,'Position',[100 200 1400 800],'color','w')
imagesc( T2_row ); axis off; axis image; colormap( hot ); caxis([0 150])
c = colorbar('Location','south');
c.Color = 'w';
c.FontSize = font_sz*0.9;
c.LineWidth = line_w;
c.Position(1) = 0.25;
c.Position(2) = 0.29;
c.Position(3) = 0.55;
c.Label.String = 'T2 (ms)';
c.Label.Position(1) = -15;
c.Label.Position(2) = 0;
my_letter = text(25,40,'d','FontSize',font_sz*2,'FontWeight','bold','Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T2_brain_rgb = imread('tmp.png');
delete('tmp.png');

n_b_r = 20;
n_row = size(T1_brain_rgb,1);
n_col = size(T1_brain_rgb,2);
my_buff1 = 255*ones(250,n_col,3);
my_buff2 = 255*ones(20,n_col,3);

col_div1 = round( n_col*2/5 - n_b_r/2 );
col_div2 = round( n_col*2/5 + n_b_r/2 );

my_img = [my_buff1; B0_brain_rgb; my_buff2; fsf_brain_rgb; my_buff1; T1_brain_rgb; my_buff2; T2_brain_rgb];

figure(8); clf;
imshow(my_img)

my_sp = 1150;
scale_font = 0.55;

txt1 = text(950,150,'SPGR','FontSize',font_sz*scale_font,'FontWeight','bold','Color','k');
txt2 = text(800 + my_sp,150,'MRF-varTE','FontSize',font_sz*scale_font,'FontWeight','bold','Color','k');

txt4 = text(240,2700,'MRF-fixTE','FontSize',font_sz*scale_font,'FontWeight','bold','Color','k');
txt5 = text(240+my_sp,2700,'MRF-varTR','FontSize',font_sz*scale_font,'FontWeight','bold','Color','k');
txt6 = text(240+2*my_sp,2700,'MRF-varTE','FontSize',font_sz*scale_font,'FontWeight','bold','Color','k');

eval( sprintf( 'export_fig %sFigure10.png -r%d',dir_out,dpi_img) );

%% Supplemental liver figure

close all

fn_liver = 'MRF_liver_varTE_proc_fat_sep.mat';

n_buff = 20;

load([dir_in fn_liver]);
mask = abs( abs( output_MRFAT.M0_fat_map ) + abs( output_MRFAT.M0_water_map ) );
my_thresh = multithresh( mask, 1);
mask( mask < my_thresh(1) ) = 0;
mask( mask > 0 ) = 1;

M0_H20_liver = abs( output_MRFAT.M0_water_map );
M0_fat_liver = abs( output_MRFAT.M0_fat_map );
FSF_liver = output_MRFAT.FSF_map.* mask;
B0_liver = output_MRFAT.B0_fit_map.* mask;
T1_liver = output_MRFAT.T1_water_map.* mask;
T2_liver = output_MRFAT.T2_water_map.* mask;

figure(1); clf;
set(gcf,'Position',[100 200 800 600],'color','w'); colormap(gray);
imagesc( M0_H20_liver ); axis image; axis off; caxis([0 max(M0_H20_liver(:))]);
my_letter = text(10,25,'a','FontSize',font_sz*2,'FontWeight','bold','Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
m0_h20_rgb = imread('tmp.png');
delete('tmp.png');

figure(2); clf;
set(gcf,'Position',[100 200 800 600],'color','w'); colormap(gray);
imagesc( M0_fat_liver ); axis image; axis off; caxis([0 max(M0_H20_liver(:))]);
my_letter = text(10,25,'b','FontSize',font_sz*2,'FontWeight','bold','Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
m0_fat_rgb = imread('tmp.png');
delete('tmp.png');

figure(3); clf;
set(gcf,'Position',[100 200 800 600],'color','w');
imagesc( FSF_liver ); axis image; axis off;
c = colorbar('Location','north');
c.Color = 'w';
c.FontSize = font_sz;
c.LineWidth = line_w;
c.Position(1) = 0.3;
c.Position(3) = 0.45;
c.Position(4) = 0.03;
c.Ticks = [0:0.2:1.0];
c.Label.String = 'FSF';
c.TickLabels{1} = num2str( str2num( c.TickLabels{1} ), '%.1f' );
black = cat(3, zeros(size(mask)), zeros(size(mask)), zeros(size(mask)));
hold on
bg_h = imshow(black);
hold off
set( bg_h, 'AlphaData', ~mask );
my_letter = text(10,25,'c','FontSize',font_sz*2,'FontWeight','bold','Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
FSF_liver_rgb = imread('tmp.png');
delete('tmp.png');

figure(4); clf;
set(gcf,'Position',[100 200 800 600],'color','w'); colormap(gray);
imagesc( B0_liver ); axis image; axis off;
c = colorbar('Location','north');
c.Color = 'w';
c.FontSize = font_sz;
c.LineWidth = line_w;
c.Position(1) = 0.3;
c.Position(3) = 0.45;
c.Position(4) = 0.03;
c.Label.String = 'B0 (Hz)';
my_letter = text(10,25,'d','FontSize',font_sz*2,'FontWeight','bold','Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
B0_liver_rgb = imread('tmp.png');
B0_liver_rgb = repmat( B0_liver_rgb, [1 1 3] );
delete('tmp.png');

figure(5); clf;
set(gcf,'Position',[100 200 800 600],'color','w'); colormap(hot);
imagesc( T1_liver ); axis image; axis off;  caxis([0 3000])
c = colorbar('Location','north');
c.Color = 'w';
c.FontSize = font_sz;
c.LineWidth = line_w;
c.Position(1) = 0.3;
c.Position(3) = 0.45;
c.Position(4) = 0.03;
c.Label.String = 'T1 (ms)';
my_letter = text(10,25,'e','FontSize',font_sz*2,'FontWeight','bold','Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T1_liver_rgb = imread('tmp.png');
delete('tmp.png');

figure(6); clf;
set(gcf,'Position',[100 200 800 600],'color','w'); colormap(hot);
imagesc( T2_liver ); axis image; axis off;  caxis([0 150])
c = colorbar('Location','north');
c.Color = 'w';
c.FontSize = font_sz;
c.LineWidth = line_w;
c.Position(1) = 0.3;
c.Position(3) = 0.45;
c.Position(4) = 0.03;
c.Label.String = 'T2 (ms)';
my_letter = text(10,25,'f','FontSize',font_sz*2,'FontWeight','bold','Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T2_liver_rgb = imread('tmp.png');
delete('tmp.png');

m0_h20_rgb = repmat( m0_h20_rgb, [1 1 3] );
m0_fat_rgb = repmat( m0_fat_rgb, [1 1 3] );

my_buff1 = 255*ones( size(m0_h20_rgb,1) , n_buff, 3 );
my_buff2 = 255*ones( n_buff, size([m0_h20_rgb m0_fat_rgb FSF_liver_rgb],2) + 2*n_buff, 3 );


figure(5); clf;
my_img = [m0_h20_rgb my_buff1 m0_fat_rgb my_buff1 FSF_liver_rgb; my_buff2; B0_liver_rgb my_buff1 T1_liver_rgb my_buff1 T2_liver_rgb];
imshow(my_img);
eval( sprintf( 'export_fig %sFigure11.png -r%d',dir_out,dpi_img) );

%% Supplemental sequence figures

close all;

% flip angle pattern for proposed MRF

fn_csv = 'MRF103.csv';
n_TR = 1500;

data = csvread(['../data_in/' fn_csv],1,0);
fa_deg_v = 60*data(:,1);

figure(1); clf;
set(gcf,'Position',[100 200 800 600],'color','w')
plot(1:n_TR,fa_deg_v,'LineWidth', line_w,'Color','k')
set(gca,'FontSize',font_sz,'LineWidth',line_w)
xlabel('TR #');
ylabel('flip angle (deg)')
eval( sprintf( 'export_fig %sFigureS1.png -r%d',dir_out,dpi_img) );


% flip angle TR pattern for variable TR

fn_csv = 'MRF001.csv';
n_TR = 1000;
TR_base_ms = 16;

data = csvread(['../data_in/' fn_csv],1,0);
fa_deg_v = 60*data(:,1);
TR_ms_v = TR_base_ms + data(:,3);

figure(2); clf;
set(gcf,'Position',[100 200 800 600],'color','w')
plot(1:n_TR,fa_deg_v,'LineWidth', line_w,'Color','k')
set(gca,'FontSize',font_sz,'LineWidth',line_w)
xlabel('TR #');
ylabel('flip angle (deg)')
my_letter = text(75,55,'a','FontSize',2*font_sz,'FontWeight','bold','Color','k');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
varTR_fa = imread('tmp.png');
delete('tmp.png');

figure(3); clf;
set(gcf,'Position',[100 200 800 600],'color','w')
plot(1:n_TR,TR_ms_v,'LineWidth', line_w,'Color','k')
set(gca,'FontSize',font_sz,'LineWidth',line_w)
xlabel('TR #');
ylabel('TR length (ms)')
my_letter = text(75,18.75,'b','FontSize',2*font_sz,'FontWeight','bold','Color','k');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
varTR_TR = imread('tmp.png');
delete('tmp.png');

varTR_fa = [varTR_fa 255*ones( size(varTR_fa,1), size(varTR_TR,2) - size(varTR_fa,2) )];

figure(4); clf;
imshow( [varTR_fa; varTR_TR] );
eval( sprintf( 'export_fig %sFigureS2.png -r%d',dir_out,dpi_img) );

%% Supplemental T1/T2 verification

close all;
marker_sz = 6;

T1_spec_ms_v = [2480, 2173, 1907, 1604, 1332, 1044, 801.7, 608.6, 458.4, 336.5, 244.2, 176.6, 126.9, 90.9];
T2_spec_ms_v = [581.3, 403.5, 278.1, 190.9, 133.3, 96.9, 64.1, 46.4, 32.0, 22.6, 15.8, 11.2 7.9, 5.6];

n_buff = 20;
orange = [0.9100 0.4100 0.1700];

fn_T1T2 = 'Msys_T1T2_verification.mat';
legend_txt = [{'MRF fixTE'},{'MRF varTR'},{'MRF varTE'},{'Conventional'}]; % order in data storage
n_MRF = numel(legend_txt);

load([dir_in fn_T1T2]);
% permute
col_idx = [3 2 1 4];
T1_mean_meds = T1_mean_meds(:,col_idx);
T2_mean_meds = T2_mean_meds(:,col_idx);
legend_txt = legend_txt(col_idx);


% CCC
CCC_T1_varTE = ccc( T1_mean_meds(:,1), T1_spec_ms_v(:) );
CCC_T1_varTR = ccc( T1_mean_meds(:,2), T1_spec_ms_v(:) );
CCC_T1_fixTE = ccc( T1_mean_meds(:,3), T1_spec_ms_v(:) );
CCC_T2_varTE = ccc( T2_mean_meds(:,1), T2_spec_ms_v(:) );
CCC_T2_varTR = ccc( T2_mean_meds(:,2), T2_spec_ms_v(:) );
CCC_T2_fixTE = ccc( T2_mean_meds(:,3), T2_spec_ms_v(:) );
CCC_T1_conv = ccc( T1_mean_meds(:,4) , T1_spec_ms_v(:) );
CCC_T2_conv = ccc( T2_mean_meds(:,4) , T2_spec_ms_v(:) );



% plot results
figure(1); clf;
set(gcf,'Position',[100 200 1000 800],'color','w') 
plot( T1_spec_ms_v(:), T1_mean_meds(:,1),'v', 'MarkerEdgeColor', orange, 'MarkerFaceColor', orange, 'LineWidth', line_w, 'MarkerSize', marker_sz );
hold on
plot( T1_spec_ms_v(:), T1_mean_meds(:,2),'bo','MarkerFaceColor', 'b', 'LineWidth', line_w, 'MarkerSize', marker_sz );
plot( T1_spec_ms_v(:), T1_mean_meds(:,3),'ks','MarkerFaceColor', 'm', 'LineWidth', line_w, 'MarkerSize', marker_sz );
plot( T1_spec_ms_v(:), T1_mean_meds(:,4),'ko','MarkerFaceColor', 'k', 'LineWidth', line_w, 'MarkerSize', marker_sz );
hold off
axis square; grid on;
xlim([0 1.05*max(T1_mean_meds(:))]); ylim([0 1.05*max(T1_mean_meds(:))]);
line( [ 0 max( T1_mean_meds(:) ) ], [ 0 max( T1_mean_meds(:) ) ],'LineStyle','--','LineWidth',line_w );
xlabel('T1 (ms)')
ylabel('MRF T1 (ms)')
% title('T1 means');
legend(legend_txt,'Location','northwest')
set(gca,'LineWidth',line_w,'FontSize',font_sz);
my_txt = text(1250, 600, {sprintf('CCC_{varTE} = %.3f',CCC_T1_varTE),sprintf('CCC_{fixTE} = %.3f',CCC_T1_fixTE),sprintf('CCC_{varTR} = %.3f',CCC_T1_varTR),sprintf('CCC_{conv} = %.3f',CCC_T1_conv)});
my_txt.FontSize = font_sz;
my_letter = text(250,1900,'a','FontSize',2*font_sz,'FontWeight','bold','Color','k');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T1_CCC = imread('tmp.png');
delete('tmp.png');

figure(2); clf;
set(gcf,'Position',[100 200 1000 800],'color','w')
plot( T2_spec_ms_v(:), T2_mean_meds(:,1),'v', 'MarkerEdgeColor', orange, 'MarkerFaceColor', orange, 'LineWidth', line_w, 'MarkerSize', marker_sz );
hold on
plot( T2_spec_ms_v(:), T2_mean_meds(:,2),'bo','MarkerFaceColor', 'b', 'LineWidth', line_w, 'MarkerSize', marker_sz );
plot( T2_spec_ms_v(:), T2_mean_meds(:,3),'ks','MarkerFaceColor', 'm', 'LineWidth', line_w, 'MarkerSize', marker_sz );
plot( T2_spec_ms_v(:), T2_mean_meds(:,4),'ko','MarkerFaceColor', 'k', 'LineWidth', line_w, 'MarkerSize', marker_sz );

hold off
axis square; grid on;
line( [ 0 max( T2_mean_meds(:) ) ], [ 0 max( T2_mean_meds(:) ) ],'LineStyle','--','LineWidth',line_w);
xlim([0 1.05*max(T2_mean_meds(:))]); ylim([0 1.05*max(T2_mean_meds(:))]);
xlabel('T2 (ms)')
ylabel('MRF T2 (ms)')
% title('T2 means');
legend(legend_txt,'Location','northwest')
set(gca,'LineWidth',line_w,'FontSize',font_sz);
my_txt = text(350, 200, {sprintf('CCC_{varTE} = %.3f',CCC_T2_varTE),sprintf('CCC_{fixTE} = %.3f',CCC_T2_fixTE),sprintf('CCC_{varTR} = %.3f',CCC_T2_varTR),sprintf('CCC_{conv} = %.3f',CCC_T2_conv)});
my_txt.FontSize = font_sz;
my_letter = text(60,460,'b','FontSize',2*font_sz,'FontWeight','bold','Color','k');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T2_CCC = imread('tmp.png');
delete('tmp.png');

T1_CCC = cat(1, T1_CCC, 255*ones( size(T2_CCC,1) - size(T1_CCC,1), size(T1_CCC,2), 3 ) );
my_buff = 255*ones( size(T1_CCC,1), n_buff, 3 );

figure(3); clf;
imshow([T1_CCC my_buff T2_CCC my_buff]);
eval( sprintf( 'export_fig %sFigureS5.png -r%d',dir_out,dpi_img) );

%% Supplemental label segments

close all

fn_sim_stats = 'img_simulation_proc_stats.mat';

load( [dir_in fn_sim_stats] );

T1_truth = squeeze( example_maps(:,:,1,1,2) );
T1_truth = T1_truth(:,:,1);

figure(100); clf;
set(gcf,'Position',[100 100 1000 1000],'color','w');
imagesc( T1_truth ); axis image; axis off; caxis([0 2250]); colormap(hot)

x = [73 54; 100 92; 115 98; 124 145; 128 128; 121 125; 122 125; 112 120; 121 124];
y = [34 31; 95 110; 133 125; 106 106; 192 182; 113 114; 132 124; 192 182; 192 182];
colors = {'w','w','w','w','w','w','w','w','w'};

n_lines = size(x,1);


for ii = 1:n_lines
    line( [ x(ii,1) x(ii,2) ], [ y(ii,1) y(ii,2) ], 'Color', colors{ii}, 'LineWidth', line_w);
end

numbers = {'5','1','3','4','3','3','2','2','2'};
colors_num = {'w','w','w','w','w','k','w','w','w','w','w'};

x_num = [ 42, 86, 124, 144, 120, 115, 115, 70, 158];
y_num = [ 30, 122, 118, 109, 174, 80, 34, 177, 177];
n_numbers = numel(numbers);

for ii = 1:n_numbers
    my_h = text( x_num(ii),y_num(ii), numbers{ii},'Color',colors_num{ii});
    my_h.FontSize = font_sz*2;
    my_h.FontWeight = 'bold';
end

eval( sprintf( 'export_fig %sFigureS3.png -r%d',dir_out,dpi_img) );

%% Supplemental fully-sampled simulated image parameter maps

close all;

fn_sim_stats = 'img_simulation_proc_stats_cart.mat';

load( [dir_in fn_sim_stats] );

T1_truth = squeeze( example_maps(:,:,1,1,2) );
T2_truth = squeeze( example_maps(:,:,2,1,2) );
FSF_truth = squeeze( example_maps(:,:,3,1,2) );
B0_truth = squeeze( example_maps(:,:,4,1,2) );

T1_maps = squeeze( example_maps(:,:,1,:,1) );
T2_maps = squeeze( example_maps(:,:,2,:,1) );
FSF_maps = squeeze( example_maps(:,:,3,:,1) );
B0_maps = squeeze( example_maps(:,:,4,:,1) );

[N,~,n_sets] = size( T1_maps );

T1_maps = reshape( permute( T1_maps, [1 3 2] ), [N*n_sets N] );
T2_maps = reshape( permute( T2_maps, [1 3 2] ), [N*n_sets N] );
FSF_maps = reshape( permute( FSF_maps, [1 3 2] ), [N*n_sets N] );
B0_maps = reshape( permute( B0_maps, [1 3 2] ), [N*n_sets N] );

T1_maps = [T1_truth; T1_maps; zeros(N/4,N)];
T2_maps = [T2_truth; T2_maps; zeros(N/4,N)];
FSF_maps = [FSF_truth; FSF_maps; zeros(N/4,N)];
B0_maps = [B0_truth; B0_maps; zeros(N/4,N)];

figure(1); clf;
set(gcf,'Position',[100 100 800 1000],'color','w');
imagesc( T1_maps ); axis image; axis off;  caxis([0 2250])
c = colorbar('Location','south');
c.Color = 'w';
c.FontSize = font_sz/2;
c.LineWidth = line_w;
c.Position(2) = .125;
c.Label.String = '(ms)';
c.Label.Position(2) = -0.8;
colormap(hot);
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T1sim = imread('tmp.png');
delete('tmp.png');

figure(2); clf;
set(gcf,'Position',[100 100 800 1000],'color','w'); 
imagesc( T2_maps ); axis image; axis off; caxis([0 125])
c = colorbar('Location','south');
c.Color = 'w';
c.FontSize = font_sz/2;
c.LineWidth = line_w;
c.Position(2) = .125;
c.Label.String = '(ms)';
c.Label.Position(2) = -0.8;
colormap(hot);
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T2sim = imread('tmp.png');
delete('tmp.png');

figure(3); clf;
set(gcf,'Position',[100 100 800 1000],'color','w');
imagesc( FSF_maps ); axis image; axis off; caxis([0 1])
c = colorbar('Location','south');
c.Color = 'w';
c.FontSize = font_sz/2;
c.LineWidth = line_w;
c.Position(2) = .5;
c.Label.String = '(FSF)';
c.Label.Position(2) = -1.0;
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
FSFsim = imread('tmp.png');
delete('tmp.png');

figure(4); clf;
set(gcf,'Position',[100 100 800 1000],'color','w');
imagesc( B0_maps ); axis image; axis off; caxis([-250 250])
c = colorbar('Location','south');
c.Color = 'w';
c.FontSize = font_sz/2;
c.LineWidth = line_w;
c.Position(2) = .5;
c.Label.String = '(Hz)';
c.Label.Position(2) = -1.0;
colormap(gray)
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
B0sim = imread('tmp.png');
delete('tmp.png');

B0sim = repmat( B0sim, [1 1 3] );

n_buff = 200;
n_mini_buff = 10;

sim_maps = [ T1sim T2sim 255*ones(size(T1sim,1),n_mini_buff,3) FSFsim B0sim ];
sim_maps = [ 255*ones(n_buff,size(sim_maps,2),3); sim_maps ];
sim_maps = [ 255*ones(size(sim_maps,1),n_buff,3) sim_maps ];

figure(5);
imshow(sim_maps);
txt_a = text(100,600,'True','FontSize',font_sz,'Rotation',90);
txt_b = text(100,1350,'MRF-varTE','FontSize',font_sz,'Rotation',90);
txt_c = text(100,1950,'MRF-fixTE','FontSize',font_sz,'Rotation',90);
txt_d = text(100,2550,'MRF-varTR','FontSize',font_sz,'Rotation',90);
txt_1 = text(450,100,'T1','FontSize',font_sz);
txt_2 = text(1050,100,'T2','FontSize',font_sz);
txt_3 = text(1625,100,'FSF','FontSize',font_sz);
txt_4 = text(2250,100,'B0','FontSize',font_sz);
eval( sprintf( 'export_fig %sFigureS6.png -r%d',dir_out,dpi_img) );


%% Supplemental fully-sampled img sim stats plots

close all

fn_sim_stats = 'img_simulation_proc_stats_cart.mat';

load( [dir_in fn_sim_stats] );

orange = [0.9100 0.4100 0.1700];

legend_txt = {'MRF-varTE w/ fat sep','MRF-fixTE w/o fat sep','MRF-varTR w/o fat sep'};
cats_snr = categorical( {'1','2','3','4','5','all'} );

T1sim_bias = squeeze( simulation_bias_means(:,1,:,:) );
T2sim_bias = squeeze( simulation_bias_means(:,2,:,:) );
FSFsim_bias = squeeze( simulation_bias_means(:,3,:,:) );
B0sim_bias = squeeze( simulation_bias_means(:,4,:,:) );

T1sim_bias_std = squeeze( simulation_bias_std(:,1,:,:) );
T2sim_bias_std = squeeze( simulation_bias_std(:,2,:,:) );
FSFsim_bias_std = squeeze( simulation_bias_std(:,3,:,:) );
B0sim_bias_std = squeeze( simulation_bias_std(:,4,:,:) );


figure(10); clf;
set(gcf,'color','w','Position',[100 200 800 800]);
my_hand = bar(cats_snr,T1sim_bias,'LineWidth',line_w/2);
my_hand(1).FaceColor = orange;
my_hand(2).FaceColor = 'b';
my_hand(3).FaceColor = 'm';
ylim([-1250 200]);
legend(legend_txt,'Location','southeast')
xlabel('Segment')
ylabel('(ms)')
title('T1 mean bias')
set(gca,'FontSize',font_sz,'LineWidth',line_w)
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T1simbias_img = imread('tmp.png');
delete('tmp.png');


figure(11); clf;
set(gcf,'color','w','Position',[100 200 800 800]);
my_hand = bar(cats_snr,T2sim_bias,'LineWidth',line_w/2);
my_hand(1).FaceColor = orange;
my_hand(2).FaceColor = 'b';
my_hand(3).FaceColor = 'm';
ylim([-20 30]);
legend(legend_txt)
xlabel('Segment')
ylabel('(ms)')
title('T2 mean bias')
set(gca,'FontSize',font_sz,'LineWidth',line_w)
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T2simbias_img = imread('tmp.png');
delete('tmp.png');

figure(12); clf;
set(gcf,'color','w','Position',[100 200 800 800]);
my_hand = bar(cats_snr,T1sim_bias_std,'LineWidth',line_w/2);
my_hand(1).FaceColor = orange;
my_hand(2).FaceColor = 'b';
my_hand(3).FaceColor = 'm';
ylim([0 400]);
legend(legend_txt,'Location','northwest')
xlabel('Segment')
ylabel('(ms)')
title('T1 standard deviation of bias')
set(gca,'FontSize',font_sz,'LineWidth',line_w)
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T1simbiasstd_img = imread('tmp.png');
delete('tmp.png');

figure(13); clf;
set(gcf,'color','w','Position',[100 200 800 800]);
my_hand = bar(cats_snr,T2sim_bias_std,'LineWidth',line_w/2);
my_hand(1).FaceColor = orange;
my_hand(2).FaceColor = 'b';
my_hand(3).FaceColor = 'm';
ylim([0 20]);
legend(legend_txt)
xlabel('Segment')
ylabel('(ms)')
title('T2 standard deviation of bias')
set(gca,'FontSize',font_sz,'LineWidth',line_w)
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
T2simbiasstd_img = imread('tmp.png');
delete('tmp.png');

n_mini_buff = 20;

Mean_bias_plots = [T1simbias_img 255*ones(size(T1simbias_img,1),n_mini_buff,3) T2simbias_img];
Std_bias_plots = [T1simbiasstd_img 255*ones(size(T1simbiasstd_img,1),n_mini_buff,3) T2simbiasstd_img];

my_buff = 255*ones( size(Std_bias_plots,1),size(Mean_bias_plots,2) -  size(Std_bias_plots,2), 3 );

Sim_bias_plots = [  Mean_bias_plots; ...
    255*ones( n_mini_buff, size( Mean_bias_plots, 2 ), 3 ); ...
    [my_buff(:,(1:round(size(my_buff,2)/2)),:) Std_bias_plots my_buff(:,round(size(my_buff,2)/2)+1:end,:) ] ];

figure(20); clf;
imshow( Sim_bias_plots );
eval( sprintf( 'export_fig %sFigureS7.png -r%d',dir_out,dpi_img) );

%% Supplemental fat-water direct match vs. k-space fitting

close all;

my_buff = 5;
my_buff2 = 40;

fn_directksp_analysis = 'MRF_directksp_analysis.mat';

load([dir_in fn_directksp_analysis])

my_thresh = multithresh( img_conventional, 1 );
FSF_graphcut( img_conventional < my_thresh ) = 0;
FSF_graphcut = FSF_graphcut(:,1:end-my_buff2);
FSF_direct( img_conventional < my_thresh ) = 0;
FSF_direct = FSF_direct(:,1:end-my_buff2);
FSF_ksp( img_conventional < my_thresh ) = 0;
FSF_ksp = FSF_ksp(:,1:end-my_buff2);
bg = ones( size( FSF_graphcut ) ); bg( img_conventional >= my_thresh ) = 0;
bg = [ bg bg; bg ones(size(bg)) ];
[row_bg, col_bg] = size( bg );
bg_cross = zeros( size( bg ) );
bg_cross( (row_bg/2 - my_buff/2):(row_bg/2 + my_buff/2), : ) = 1;
bg_cross( :, (col_bg/2 - my_buff/2):(col_bg/2 + my_buff/2) ) = 1;
bg_cross( row_bg/2:end, col_bg/2:end) = 1;

figure(1); clf;
set(gcf,'Position',[100 200 1000 800],'color','w');
imagesc( [FSF_graphcut FSF_direct; FSF_ksp zeros(size(FSF_ksp))] ); axis image; axis off;
c = colorbar('Location','westoutside');
c.Color = 'k';
c.FontSize = font_sz;
c.LineWidth = line_w;
c.Position(1) = 0.18;
c.Position(2) = 0.27;
c.Position(4) = 0.5;
c.Ticks = [0:0.2:1.0];
c.TickLabels{1} = num2str( str2num(c.TickLabels{1}), '%.1f');
c.TickLabels{end} = num2str( str2num(c.TickLabels{end}), '%.1f');
c.Label.String = 'FSF';
c.Label.Position(1) = 2.7;
black = cat(3, zeros(size(bg)), zeros(size(bg)), zeros(size(bg)));
white = 255*cat(3, ones(size(bg)), ones(size(bg)), ones(size(bg)));
hold on
bg_h = imshow(black);
hold off
set( bg_h, 'AlphaData', bg );
hold on
bg_cross_h = imshow(white);
hold off
set( bg_cross_h, 'AlphaData', bg_cross );
my_letter = text(20,25,'a','FontSize',font_sz,'FontWeight','bold','Color','w');
my_letter = text(20+240-my_buff2,25,'b','FontSize',font_sz,'FontWeight','bold','Color','w');
my_letter = text(20,25+240,'c','FontSize',font_sz,'FontWeight','bold','Color','w');
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
FSF_directVksp = imread('tmp.png');
delete('tmp.png');


figure(2); clf;
set(gcf,'Position',[100 200 1000 800],'color','w');
plot(FSF_directksp_means(:,1), FSF_directksp_means(:,2), 'bo','MarkerFaceColor','b','MarkerSize',marker_sz,'LineWidth',line_w); axis square;
hold on;
orange = [0.9100 0.4100 0.1700];
plot(FSF_directksp_means(:,1), FSF_directksp_means(:,3), 'v','MarkerFaceColor',orange,'MarkerEdgeColor','k','MarkerSize',marker_sz,'LineWidth',line_w); axis square;
hold off;
line([0 1],[0 1],'LineStyle','--','Color','r','LineWidth',line_w);
set(gca,'FontSize',font_sz,'LineWidth',line_w)
xlabel('Conventional FSF');
ylabel('MRF FSF')
tix=get(gca,'xtick')';
set(gca,'ytick',tix);
set(gca,'xticklabel',num2str(tix,'%.1f'),'yticklabel',num2str(tix,'%.1f'))
my_leg = legend('dictionary fat sep.','k-space fat sep.','Location','southeast');
my_leg.Position(1) = 0.52;
my_leg.Position(2) = 0.20;
my_txt = text(0.56, 0.30, {sprintf('CCC_{dict} = %.3f',CCC_direct),sprintf('CCC_{ksp} = %.3f',CCC_ksp)});
my_txt.FontSize = font_sz;
eval( sprintf( 'export_fig tmp.png -r%d',dpi_img) );
CCC_FSF_directvksp = imread('tmp.png');
delete('tmp.png');

CCC_FSF_directvksp_sc = imresize( CCC_FSF_directvksp, 0.35 );

[n_r_d, n_c_d, ~] = size( CCC_FSF_directvksp_sc );
FSF_directVksp( (end - n_r_d + 1 - my_buff2):(end-my_buff2), (end - n_c_d + 1):end, :) = CCC_FSF_directvksp_sc;
figure(3); clf;
imshow( FSF_directVksp )
my_letter = text(1300,1120,'d','FontSize',font_sz,'FontWeight','bold','Color','k');
eval( sprintf( 'export_fig %sFigureS8.png -r%d',dir_out,dpi_img) );