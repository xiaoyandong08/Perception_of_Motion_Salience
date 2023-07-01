function Generate_Modularity_of_Trajectory_Hover_for_anyFrame
addpath(genpath('../utility'))

% load('../../../../XYD_Data/Chimney_swift/Swift_flocking_A-Xiao_all_frame_new_present_bird_0.25_persist_time_150.mat');
% 
% Frame_index = find(which_Frame == 6221);
% Frame = which_Frame(Frame_index);
% Frame_matrix = all_Frame_matrix{Frame_index};
% Frame_time = tracks_filt(Frame_matrix(1,:),5)';
% tracks_filt = tracks_filt(Frame_matrix(:),:);
% Frame_matrix = reshape([1:1:size(tracks_filt,1)],[size(Frame_matrix,1),size(Frame_matrix,2)]);
% 
% save(['../data_after_processing/Circling_traj_of_Frame' num2str(Frame) '_for_Modularity' '.mat'],'Frame','Frame_matrix','tracks_filt',...
%     'Frame_time','-v7.3')

%%

% load ../data_after_processing/Circling_traj_of_Frame6221_for_Modularity.mat
% 
% findRou = {'max'};
% tau_threshold = [0.25];
% correlation_threshold = 0.8;
% 
% plot_LFT = 0;
% [Delay,Delay_neg] = Mapping_Leader_follow_network_anis_factor_simplified(plot_LFT,Frame_matrix,tracks_filt,findRou,tau_threshold,correlation_threshold);
% 
% save('../data_after_processing/LFT_of_Circling_of_Frame6221.mat','-v7.3')

load ../data_after_processing/LFT_of_Circling_of_Frame6221.mat
One_LFnet = Delay_neg;
%%
flock_tag = 'Circling';
json_type = 'onlyTraj';


rng(123)
bp = Bipartite((abs(One_LFnet)));
bp.community.Detect();
plotFormat = PlotFormat();
plotFormat.back_color = [41,38,78]/255;
plotFormat.cell_color = 'white';
plotFormat.use_labels = true;
plotFormat.font_size = 8;
% figure;
% bp.plotter.SetPlotFormat(plotFormat);
% bp.plotter.PlotModularMatrix();
% figure;
% subplot(121);PlotWebs.PLOT_MATRIX((abs(One_LFnet)), plotFormat);
% subplot(122);PlotWebs.PLOT_NESTED_MATRIX((abs(One_LFnet)), plotFormat);


row_label_in_plot = bp.community.index_rows;%row label in plot
col_label_in_plot = bp.community.index_cols;
row_module_index = bp.community.row_modules;% module index for each bird
col_module_index = bp.community.col_modules;

% modularity trajectory
for i = 1 : max(row_module_index)^2
    [index_row,index_col] = ind2sub([max(row_module_index) max(row_module_index)],i);
    index1 = find(row_module_index(row_label_in_plot)==index_row);
    a1 = row_label_in_plot(index1);
    index1 = find(col_module_index(col_label_in_plot)==max(row_module_index)-index_col+1);
    a2 = col_label_in_plot(index1);
    
    module_One_LFnet{index_row,index_col} = One_LFnet(a1,a2);
    module_row{index_row,index_col} =  a1;
    module_col{index_row,index_col} =  a2;
    
    module_intersect_bird{index_row,index_col} = intersect(module_row{index_row,index_col},module_col{index_row,index_col});
    module_intersect_One_LFnet{index_row,index_col} = One_LFnet(module_intersect_bird{index_row,index_col},module_intersect_bird{index_row,index_col});

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(json_type,'onlyTraj')
    modularity_index = zeros(1,size(Frame_matrix,1));
    count = 1;
    for i = 1 : size(module_intersect_bird,1)
        for j = 1 : size(module_intersect_bird,2)
            if length(module_intersect_bird{i,j})>0
                modularity_index(module_intersect_bird{i,j}) = count;
                count = count + 1;
            end
        end
    end
    fun = @(m)sRGB_to_OSAUCS(m,true,true);
    modularity_color = maxdistcolor(max(modularity_index),fun);
end


figure;
subplot('position',[0.1 0.68 0.8 0.30])
Color = jet(size(Frame_matrix,2));
xlabel('x');ylabel('y');zlabel('z')
set(gca,'FontSize',14,'TickLength',[0.03, 0.01],...
    'XMinorTick','on','YMinorTick','on','boxstyle','full'); 
view([-154 62])
xlim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2)))])
ylim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3)))])
zlim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4)))])
hold on;box on
for i = 1 : size(Frame_matrix,2)
    Id = Frame_matrix(find(Frame_matrix(:,i)>0),i);
    xyz=tracks_filt(Id,2:4);
    v_xyz = tracks_filt(Id,6:8);
    scatter3(xyz(:,1),xyz(:,2),xyz(:,3),5,modularity_color(modularity_index,:),'filled');
    if i == size(Frame_matrix,2)
        scatter3(xyz(:,1),xyz(:,2),xyz(:,3),65,[0.2 0.2 0.2],'o','filled');
    end
%     im(i) = getframe;
end
if exist('im')
    a=VideoWriter('0.5_0.45','MPEG-4');
    a.FrameRate = 10;
    a.Quality = 100;
    open(a);
    writeVideo(a,im);
    close(a)
end

% modulairty of LFT 
subplot('position',[0.1 0.08 0.8 0.5]);hold on;box on;
imagesc(One_LFnet(row_label_in_plot,col_label_in_plot));
colormap([parula;[1 1 1]])
title(['Frame = ' num2str(Frame)])
set(gca,'fontsize',14)
ylabel('Modularitys row index')
xlabel('Modularitys column index')
index1 = 0;index2 = size(One_LFnet,1);
for i = 1 : max(row_module_index)-1
    index1 = index1 + sum(row_module_index==i);
    plot([0.5 size(One_LFnet,1)],[index1+0.5 index1+0.5],'k-')
    
    index2 = index2 - sum(col_module_index==i);
    plot([index2+0.5 index2+0.5],[0.5 size(One_LFnet,1)],'k-')
end
xlim([0.3 size(One_LFnet,1)]);
ylim([0.3 size(One_LFnet,1)]);
set(gca,'Ydir','reverse')
set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gcf,'position',[205 192 1086 1131])
set(gcf,'position',[102 481 800 1131])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlims = [floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2)))];
ylims = [floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3)))];
zlims = [floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4)))];

fig = figure;
module_ave_order = nan(max(row_module_index),max(row_module_index));
module_nestedness = nan(max(row_module_index),max(row_module_index));
for i = 1 : max(row_module_index)^2
    [index_row,index_col] = ind2sub([max(row_module_index) max(row_module_index)],i);
    if length( module_intersect_bird{index_row,index_col})>0
        count = unique(modularity_index(module_intersect_bird{index_row,index_col}));
       [module_ave_order(index_row,index_col),module_nestedness(index_row,index_col)] = Plot_order_trajectory_for_all_module(fig,max(row_module_index),index_row,index_col,module_intersect_One_LFnet{index_row,index_col},...
            Frame_matrix(module_intersect_bird{index_row,index_col},:),tracks_filt,modularity_color(count,:),xlims,ylims,zlims);
    else
        %subplot(max(row_module_index),max(row_module_index),(index_row-1)*max(row_module_index)+index_col)
        subaxis(max(row_module_index),max(row_module_index),(index_row-1)*max(row_module_index)+index_col,'SpacingVertical',0.02,'SpacingHorizontal',0.02)
        view([-162 78])
        %view([-154 62])
        set(gca,'FontSize',12,'TickLength',[0.03, 0.01],...
            'XMinorTick','on','YMinorTick','on','boxstyle','full');
        box on    
    end
    set(gca,'XTick',[],'YTick',[],'ZTick',[])
    %xlabel('x');ylabel('y');zlabel('z')
end
set(gcf,'position',[102 481 1318 1276])



end

function [ave_order,One_LFnet_nestedness] = Plot_order_trajectory_for_all_module(fig,num_module,index_row,index_col,One_LFnet,Frame_matrix,tracks_filt,modularity_color,xlims,ylims,zlims)

order = nan(1,size(Frame_matrix,2));
for i = 1 : size(Frame_matrix,2)
    Id = Frame_matrix(find(Frame_matrix(:,i)>0),i);
    if length(Id)>0
        v_xyz = tracks_filt(Id,6:8);
        v_xyz = v_xyz./vecnorm(v_xyz')';
        order(i) = norm(sum(v_xyz,1))/size(v_xyz,1);
        
        p_xyz = tracks_filt(Id,2:4);
        mean_p_xyz = mean(p_xyz,1);
        p_xyz = (p_xyz - mean_p_xyz);
        p_xyz = p_xyz./vecnorm(p_xyz')';
        moment(i) = norm(sum(cross(p_xyz,v_xyz),1))/length(Id);

    end
end
ave_order = mean(order);

Color = jet(size(Frame_matrix,2));


%subplot(num_module,num_module,(index_row-1)*num_module+index_col)
subaxis(num_module,num_module,(index_row-1)*num_module+index_col,'SpacingVertical',0.02,'SpacingHorizontal',0.02)
% subplot('position',[0.1 0.45 0.8 0.36])
% xlabel('x');ylabel('y');zlabel('z')
set(gca,'FontSize',12,'TickLength',[0.03, 0.01],...
    'XMinorTick','on','YMinorTick','on','boxstyle','full'); 
view([-162 78])
% view([-154 62])
% plot3(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2),tracks_filt(Frame_matrix(Frame_matrix(:)>0),3),tracks_filt(Frame_matrix(Frame_matrix(:)>0),4),'.')

% xlim(xlims)
% ylim(ylims)
% zlim(zlims)
xlim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2)))])
ylim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3)))])
zlim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4)))])


hold on;box on
for i = 1 : size(Frame_matrix,2)
    
    Id = Frame_matrix(find(Frame_matrix(:,i)>0),i);
    xyz=tracks_filt(Id,2:4);
    v_xyz = tracks_filt(Id,6:8);
    %plot3(xyz(:,1),xyz(:,2),xyz(:,3),'.','color',Color(i,:));
    scatter3(xyz(:,1),xyz(:,2),xyz(:,3),5,modularity_color,'filled');
    
%     if i==1 & length(unique(sum(logical(Frame_matrix),1)))==1
%         text(xyz(:,1),xyz(:,2),xyz(:,3),num2str([1:size(Frame_matrix,1)]'))
%     end
%     im(i) = getframe;
end
for i = 1 : size(Frame_matrix,1)
    Id = Frame_matrix(i,find(Frame_matrix(i,:)>0));
    if length(Id)>0
        xyz=tracks_filt(Id(end),2:4);
        scatter3(xyz(:,1),xyz(:,2),xyz(:,3),65,[0.2 0.2 0.2],'o','filled');
    end
end
set(gca,'fontsize',14)

Delay2 = One_LFnet;
% Delay2(1:size(Delay2,1)+1:size(Delay2,1)^2) = 1;
One_LFnet_nestedness = Nestedness.NODF(abs(sign(Delay2)));
One_LFnet_nestedness = One_LFnet_nestedness.N;

end


