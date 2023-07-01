function LF_MS_correlation_for_Fig3bc_of_Mobbing
% addpath('../utility')
% 
% load('../data_after_processing/Jackdaw_mob_01_all_frame_new_present_bird_0.25.mat');
% 
% all_time = unique(tracks_filt(:,5));
% all_bird = unique(tracks_filt(:,1));
% 
% 
% Frame_index = 2; 
% Frame = which_Frame(Frame_index);
% Frame_matrix = all_Frame_matrix{Frame_index};
% Frame_time = all_Frame_time{Frame_index};
% 
% anis_factor = 1;
% findRou = 'max';
% tau_threshold = 0.25;
% Corr_Type = 'Spearman';
% correlation_threshold = 0.8;
% time_slice = ceil(linspace(1,size(Frame_matrix,2),51));
% tau_slice  = time_slice(9:2:end);
% time_slice = time_slice(10:2:end);
% 
% [temporal_Mij, temporal_Delay,temporal_Delay_neg,corr_M_LF] =...
%     Bird_Motion_salience_vs_LF_Anisotropy_without_nan_L2F(anis_factor,Corr_Type,tau_slice,time_slice,Frame_matrix,tracks_filt,findRou,tau_threshold,correlation_threshold);
% save(['One_frame_Temporal_results_anis=' num2str(anis_factor) '_frame' num2str(Frame) '_' Corr_Type '_L2Fnet.mat'],'-v7.3')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_anis0 = 'One_frame_Temporal_results_anis=0_frame2673_Spearman_L2Fnet.mat';
file_anis1 = 'One_frame_Temporal_results_anis=1_frame2673_Spearman_L2Fnet.mat';

Plot_LF_MS_anis01(file_anis0,file_anis1);

end

function Plot_LF_MS_anis01(file_anis0,file_anis1)

load(file_anis0)

time_slice = ceil(linspace(1,size(Frame_matrix,2),51));
tau_slice  = time_slice(9:2:end);
time_slice = time_slice(10:2:end); 



%%%%%%%%%%%%%%%%%%% LFMS matrix for anis = 0
figure;
set(gcf,'position',[1958 148 490 1449])

blue = hex2rgb('#0046F6');
red = hex2rgb('#F74461');
w2r= [linspace(1,red(1),50)' linspace(1,red(2),50)' linspace(1,red(3),50)'];
w2r = w2r(length(w2r):-1:1,:);
b2w = [linspace(blue(1),1,50)' linspace(blue(2),1,50)' linspace(blue(3),1,50)'];
b2w = b2w(length(b2w):-1:1,:);
colormapp = [b2w(length(b2w):-1:1,:);w2r(length(w2r):-1:1,:)];

h1 = subplot('position',[0.12 0.34 0.8 0.18]);
imagesc(Frame_time(time_slice),Frame_time(tau_slice)-Frame_time(1),corr_M_LF)
xlabel('t (s)');ylabel('\tau  (s)');
colormap(h1,[hex2rgb('CCCCCC'); colormapp] );
set(gca,'YDir','reverse')
caxis([-1 1])
title([num2str(Frame) ', ' Corr_Type ', \alpha = ' num2str(anis_factor)])
set(gca,'fontsize',14)

clearvars -except file_anis1 
%% for anis = 1
load(file_anis1)

time_slice = ceil(linspace(1,size(Frame_matrix,2),51));
tau_slice  = time_slice(9:2:end);
time_slice = time_slice(10:2:end); 

%%%%%%%%%%%%%%%%%%% LFMS matrix for anis = 1
blue = hex2rgb('#0046F6');
red = hex2rgb('#F74461');
w2r= [linspace(1,red(1),50)' linspace(1,red(2),50)' linspace(1,red(3),50)'];
w2r = w2r(length(w2r):-1:1,:);
b2w = [linspace(blue(1),1,50)' linspace(blue(2),1,50)' linspace(blue(3),1,50)'];
b2w = b2w(length(b2w):-1:1,:);
colormapp = [b2w(length(b2w):-1:1,:);w2r(length(w2r):-1:1,:)];

h2 = subplot('position',[0.12 0.11 0.8 0.18]);
imagesc(Frame_time(time_slice),Frame_time(tau_slice)-Frame_time(1),corr_M_LF)
xlabel('t (s)');ylabel('\tau  (s)');
colormap(h2,[hex2rgb('CCCCCC'); colormapp] );
set(gca,'YDir','reverse')
caxis([-1 1])
title([num2str(Frame) ', ' Corr_Type ', \alpha = ' num2str(anis_factor)])
set(gca,'fontsize',14)


end








