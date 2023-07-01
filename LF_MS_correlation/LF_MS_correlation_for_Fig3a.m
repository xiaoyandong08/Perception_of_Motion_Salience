function LF_MS_correlation_for_Fig3a
% addpath('../utility')

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
% Plot_order_trajectory_for_One_Frame(all_Frame_matrix{Frame_index},all_Frame_time{Frame_index},tracks_filt,Frame)
% 
% 
% 
% anis_factor = 0;
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

%%

load One_frame_Temporal_results_anis=0_frame2673_Spearman_L2Fnet.mat

given_tau = tau_slice(14);
time = time_slice(17);

tau = Frame_time(given_tau)-Frame_time(1);
now_time = Frame_time(time);

One_LFnet = temporal_Delay_neg{find(tau_slice==given_tau),find(time_slice==time)};

One_Mij = temporal_Mij{find(tau_slice==given_tau),find(time_slice==time)};
mean_MS = nanmean(One_Mij,1);

%%%%%%%%%%%%%%%%%%% Highlight a period of flock traj for Fig.3a
Color = jet(size(Frame_matrix,2));
figure
%xlabel('x');ylabel('y');zlabel('z')
set(gca,'FontSize',12,'TickLength',[0.03, 0.01],...
    'XMinorTick','on','YMinorTick','on');
grid on;box off
view([-67 73])
title({['Frames = ' num2str(Frame)];['T = ' num2str(round(now_time,4)) 's, \tau = ' num2str(round(tau,4)) 's']},'fontweight','normal')
xlim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2)))])
ylim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3)))])
zlim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4)))])
hold on;box on
for i = 1 : size(Frame_matrix,2)
    
    Id = Frame_matrix(find(Frame_matrix(:,i)>0),i);
    xyz=tracks_filt(Id,2:4);
    v_xyz = tracks_filt(Id,6:8);
    if i<=time && i>=time- given_tau
        plot3(xyz(:,1),xyz(:,2),xyz(:,3),'.','color',hex2rgb('F74421'));
    else
        plot3(xyz(:,1),xyz(:,2),xyz(:,3),'.','color',hex2rgb('A6A8AB'));
    end
%     if i==1 & length(unique(sum(logical(Frame_matrix),1)))==1
%         text(xyz(:,1),xyz(:,2),xyz(:,3),num2str([1:size(Frame_matrix,1)]'))
%     end
    if i == size(Frame_matrix,2)
        scatter3(xyz(:,1),xyz(:,2),xyz(:,3),65,[0.2 0.2 0.2],'o','filled');
    end
%     im(i) = getframe;
end
set(gca,'fontsize',14)
set(gcf,'position',[102 756 479 404])

%%%%%%%%%%%%%%%%%%% LF network for Fig.3a
[Gr,Cr] = global_reaching_centrality(sign(abs(One_LFnet))');
[~,sort_column] = sort(Cr,'descend');
[~,sort_row] = sort(sum(One_LFnet~=0,2),'descend');

G = digraph(abs(One_LFnet)');
figure
LWidths = 7*G.Edges.Weight/max(G.Edges.Weight);
p = plot(G,'layout','layered','NodeLabel',{},'LineWidth',LWidths,'ArrowSize',8,'MarkerSize',24,'NodeColor',hex2rgb('F7941D'),'EdgeColor',hex2rgb('262626'));axis off
text(p.XData,p.YData,num2str([1:size(Frame_matrix,1)]'),'Color','w','HorizontalAlignment','center','FontSize',18);
set(gca,'fontsize',14)
set(gcf,'position',[95 395 288 275])
title(['T = ' num2str(round(now_time,4)) 's, \tau = ' num2str(round(tau,4)) 's'],'fontweight','normal')

%%%%%%%%%%%%%%%%%%% Leading tier for each individual
figure
radius_Mij = Cr';
for i = 1 : size(radius_Mij,1)
    for j = 1 : size(radius_Mij,2)
        if ~isnan(radius_Mij(i,j))
            r_agent = 0.1+0.4*radius_Mij(i,j);
            arrayfun(@(x,y) rectangle('Position', [x-r_agent, y-r_agent, r_agent*2, r_agent*2], 'Curvature', [1 1], 'EdgeColor', hex2rgb('262626'), 'FaceColor', hex2rgb('ffa500')), i,j)
            text(i,j,num2str(j),'fontsize',18,'HorizontalAlignment','center')
            text(i*2,j,num2str(round(radius_Mij(i,j),2)),'fontsize',18,'HorizontalAlignment','center')
        end
    end
end
axis equal
xlim([0.5 size(radius_Mij,1)+0.5]);ylim([0.5 size(radius_Mij,2)+0.5])
axis off
set(gcf,'position',[57 412 408 290])
set(gca,'YDir','reverse')
title('Leading tier','fontsize',18)

%%%%%%%%%%%%%%%%%%% MS matrix for Fig.3a
radius_Mij = One_Mij;
figure;
subplot('position',[0.1 0.35 0.6 0.6]);box on
for i = 1 : size(radius_Mij,1)
    for j = 1 : size(radius_Mij,2)
        if ~isnan(radius_Mij(i,j))
            r_agent = radius_Mij(i,j);
            arrayfun(@(x,y) rectangle('Position', [x-r_agent, y-r_agent, r_agent*2, r_agent*2], 'Curvature', [1 1], 'EdgeColor', hex2rgb('262626'), 'FaceColor', hex2rgb('ffa500')), j,i)
            %text(j+r_agent+0.03,i,num2str(round(radius_Mij(i,j),2)),'fontsize',10,'HorizontalAlignment','center','VerticalAlignment','baseline','Rotation',270)
        end
    end
end
set(gca,'YDir','reverse')
axis equal
xlim([0.5 size(radius_Mij,2)+0.5]);ylim([0.5 size(radius_Mij,1)+0.5])
set(gca,'fontsize',14,'xtick',[1:1:size(radius_Mij,1)])
xlabel('Bird index');ylabel('Bird index');
title('MS of a period','fontsize',18)

radius_Mij = mean_MS;
subplot('position',[0.1 0.05 0.6 0.25]);box on
for i = 1 : size(radius_Mij,1)
    for j = 1 : size(radius_Mij,2)
        if ~isnan(radius_Mij(i,j))
            r_agent = radius_Mij(i,j);
            arrayfun(@(x,y) rectangle('Position', [x-r_agent, y-r_agent, r_agent*2, r_agent*2], 'Curvature', [1 1], 'EdgeColor', hex2rgb('262626'), 'FaceColor', hex2rgb('ffa500')), j,i)
            text(j,i,num2str(round(radius_Mij(i,j),2)),'fontsize',10,'HorizontalAlignment','center')
        end
    end
end
axis equal
xlim([0.5 size(radius_Mij,2)+0.5]);ylim([0.5 size(radius_Mij,1)+0.5])
axis off
set(gcf,'position',[57 412 408 484])

%%%%%%%%%%%%%%%%%%% average MS for each individual
figure
radius_Mij = mean_MS;
for i = 1 : size(radius_Mij,1)
    for j = 1 : size(radius_Mij,2)
        if ~isnan(radius_Mij(i,j))
            r_agent = 1.3*radius_Mij(i,j);
            arrayfun(@(x,y) rectangle('Position', [x-r_agent, y-r_agent, r_agent*2, r_agent*2], 'Curvature', [1 1], 'EdgeColor', hex2rgb('262626'), 'FaceColor', hex2rgb('ffa500')), i,j)
            text(i,j,num2str(j),'fontsize',18,'HorizontalAlignment','center')
            text(i*2,j,num2str(round(radius_Mij(i,j),2)),'fontsize',18,'HorizontalAlignment','center')
        end
    end
end
axis equal
xlim([0.5 size(radius_Mij,1)+0.5]);ylim([0.5 size(radius_Mij,2)+0.5])
axis off
set(gcf,'position',[57 412 408 358])
set(gca,'YDir','reverse')
title('average MS','fontsize',18)


%%%%%%%%%%%%%%%%%%% LF MS scatter for Fig.3a

Mj = nanmean(One_Mij)';
[R,pvalue] = corr(Cr,Mj,'Type','Spearman');
figure;hold on;box on;
[a,b] = sort(Mj,'ascend');
stem(a,'MarkerSize',14,'MarkerFaceColor',hex2rgb('0072BD'),'MarkerEdgeColor',hex2rgb('262626'))
stem(Cr(b),'MarkerSize',14,'MarkerFaceColor',hex2rgb('F74461'),'MarkerEdgeColor',hex2rgb('262626'))
xlim([0.5 length(Mj)+0.5])
legend('M_i(t,\tau)','L_i(t,\tau)','location','best')
xlabel('Bird index sorted by ascending M_i')
ylabel('M_i or L_i')
set(gca,'fontsize',18,'xtick',[1:1:7],'XtickLabel',b);
% title(['Spearman Correlation = ' num2str(round(R,4)) ', p-value = ' num2str(pvalue)],'fontweight','normal')
set(gcf,'position',[166 263 413 336])
end










