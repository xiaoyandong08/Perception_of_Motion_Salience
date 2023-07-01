function Plot_Mobbing_of_Fig1a
addpath('../utility')

load('../data_after_processing/Jackdaw_mob_01_all_frame_new_present_bird_0.25.mat');

all_time = unique(tracks_filt(:,5));
all_bird = unique(tracks_filt(:,1));


Frame_index = 2; 
Frame = which_Frame(Frame_index);
Frame_matrix = all_Frame_matrix{Frame_index};
Frame_time = all_Frame_time{Frame_index};
Plot_order_trajectory_for_One_Frame(all_Frame_matrix{Frame_index},all_Frame_time{Frame_index},tracks_filt,Frame)


end