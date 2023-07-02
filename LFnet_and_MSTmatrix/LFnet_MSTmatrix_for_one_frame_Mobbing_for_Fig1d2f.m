function LFnet_MSTmatrix_for_one_frame_Mobbing_for_Fig1d2f
addpath('../utility')

load('../data_after_processing/Jackdaw_mob_01_all_frame_new_present_bird_0.25.mat');

all_time = unique(tracks_filt(:,5));
all_bird = unique(tracks_filt(:,1));


Frame_index = 2; % frame 2 is 2095;frame 8 is 3815;frame 10 is 4202;
Frame = which_Frame(Frame_index);
Frame_matrix = all_Frame_matrix{Frame_index};

Plot_order_trajectory_for_One_Frame(all_Frame_matrix{Frame_index},all_Frame_time{Frame_index},tracks_filt,Frame)

anis_factor = 0;
findRou = {'max'};
tau_threshold = [0.25];
correlation_threshold = 0.8;

% load /Users/yandongxiao/XYD_Data/nchoosek_2to1000.mat
for i = 2 : 1000
    time_pair{i} = nchoosek(1:i,2);
end
[Delay,Delay_neg,max_motion_delay,rou,max_rou] = Mapping_Leader_follow_network_anis_factor_with_MST(Frame_matrix,tracks_filt,anis_factor,time_pair,findRou,tau_threshold,correlation_threshold);
Delay_neg(isnan(Delay_neg))=0;

% clearvars time_pair
% save(['One_frame_results_anis=' num2str(anis_factor) '_frame' num2str(Frame) '_Mobbing.mat'],'Frame','Frame_matrix','tracks_filt',...
%     'anis_factor','Delay','Delay_neg','max_motion_delay','rou','max_rou','-v7.3')


% load One_frame_results_anis=0_frame2673_Mobbing.mat



One_LFnet = Delay_neg;
One_anis_factor = anis_factor;
One_delay = Delay;
One_MS = max_motion_delay;
Frame_time = tracks_filt(Frame_matrix(find(sum(Frame_matrix~=0,2)==size(Frame_matrix,2),1),:),5)';
One_rou = rou;
One_max_rou = max_rou;


[Gr,Cr] = global_reaching_centrality(sign(abs(One_LFnet')));
[~,sort_column] = sort(Cr,'descend');
[~,sort_row] = sort(sum(One_LFnet~=0,2),'descend');

%% LFT matrix and graph | Fig.1d inset, Fig.1d

G = digraph(abs(One_LFnet)');
figure
LWidths = 7*G.Edges.Weight/max(G.Edges.Weight);
p = plot(G,'layout','layered','NodeLabel',{},'LineWidth',LWidths,'ArrowSize',8,'MarkerSize',16,'NodeColor',hex2rgb('F7941D'),'EdgeColor',hex2rgb('262626'));axis off
text(p.XData,p.YData,num2str([1:size(Frame_matrix,1)]'),'Color','w','HorizontalAlignment','center','FontSize',14);
set(gca,'fontsize',14)
set(gcf,'position',[95 395 244 266])


% togml(sign(abs(One_LFnet))',ones(1,size(One_LFnet,1)),'LFnet_mobbing.gml')

%%%%%% LFT matrix
abs_LFnet = abs(One_LFnet);
Final_M0 = abs_LFnet(sort_row,sort_column);
figure;
h = imagesc(sign(Final_M0));
caxis([-1 1]);colormap(hex2rgb('262626'))
set(h,'AlphaData',Final_M0/max(Final_M0(:)))
for i = 1:size(One_LFnet,1)
    for j = 1:size(One_LFnet,1)
        if Final_M0(i,j)~=0
            t = text(j,i,sprintf('%.2f',-1*Final_M0(i,j)),'horizontalAlignment','center');
            %t = text(j,i,num2str(round(-1*Final_M0(i,j),2)),'horizontalAlignment','center');
            t.FontSize = 9;
            line([j-0.5 j-0.5],[i-0.5 i+0.5],'Color',hex2rgb('FFFFFF'),'LineWidth',5)
            line([j+0.5 j+0.5],[i-0.5 i+0.5],'Color',hex2rgb('FFFFFF'),'LineWidth',5)
            line([j-0.5 j+0.5],[i-0.5 i-0.5],'Color',hex2rgb('FFFFFF'),'LineWidth',7)
            line([j-0.5 j+0.5],[i+0.5 i+0.5],'Color',hex2rgb('FFFFFF'),'LineWidth',7)
            
            line([j-0.4 j-0.4],[i-0.4 i+0.4],'Color',hex2rgb('262626'),'LineWidth',1)
            line([j+0.4 j+0.4],[i-0.4 i+0.4],'Color',hex2rgb('262626'),'LineWidth',1)
            line([j-0.4 j+0.4],[i-0.4 i-0.4],'Color',hex2rgb('262626'),'LineWidth',1)
            line([j-0.4 j+0.4],[i+0.4 i+0.4],'Color',hex2rgb('262626'),'LineWidth',1)
        end
        
    end
end
set(gca,'fontsize',14,'xtick',[1:1:size(One_LFnet,1)],'xticklabel',sort_column,'ytick',[1:1:size(One_LFnet,1)],'yticklabel',sort_row);
set(gcf,'position',[113 113 562 515])
xlabel('Bird index');
ylabel('Bird index')
set(gcf,'position',[173 124 450 450])
set(gcf,'position',[173 200 450 374])
set(gcf,'position',[173 200 285 374])

aa = reshape(zeros(1,size(Frame_matrix,1)^2),size(Frame_matrix,1),size(Frame_matrix,1));
aa(1:size(Frame_matrix,1)+1:size(Frame_matrix,1)^2) = 1;

nestedness = Nestedness.NODF(sign(abs(One_LFnet)));
nestedness = nestedness.N;

nestedness1 = Nestedness.NODF(sign(abs(One_LFnet')));
nestedness1 = nestedness1.N;

figure;
subplot(121);PlotWebs.PLOT_MATRIX(sign(abs(One_LFnet)));
subplot(122);PlotWebs.PLOT_NESTED_MATRIX(sign(abs(One_LFnet)));


%%
%%%%%% MST matrix

Final_M0 = abs(One_MS);
figure;
h = imagesc(sign(Final_M0));
caxis([0 1]);
colormap(hex2rgb('262626'))
set(h,'AlphaData',Final_M0/max(Final_M0(:)))
% set(h,'AlphaData',Final_M0)
set(gca,'fontsize',14,'xtick',[1:1:size(Final_M0,1)],'ytick',[1:1:size(Final_M0,1)]);
xlabel('Bird index');ylabel('Bird index')
set(gcf,'position',[44 166 462 340])
for i = 1:size(Final_M0,1)
    for j = 1:size(Final_M0,1)
        if Final_M0(i,j)~=0
            if Final_M0(i,j)>1.5
                t = text(j,i,sprintf('%.2f',1*Final_M0(i,j)),'horizontalAlignment','center','color','w');
            else
                t = text(j,i,sprintf('%.2f',1*Final_M0(i,j)),'horizontalAlignment','center');                
            end
            %t = text(j,i,num2str(round(-1*Final_M0(i,j),2)),'horizontalAlignment','center');
            t.FontSize = 10;
            line([j-0.5 j-0.5],[i-0.5 i+0.5],'Color',hex2rgb('FFFFFF'),'LineWidth',8)
            line([j+0.5 j+0.5],[i-0.5 i+0.5],'Color',hex2rgb('FFFFFF'),'LineWidth',8)
            line([j-0.5 j+0.5],[i-0.5 i-0.5],'Color',hex2rgb('FFFFFF'),'LineWidth',6)
            line([j-0.5 j+0.5],[i+0.5 i+0.5],'Color',hex2rgb('FFFFFF'),'LineWidth',6)
            
            line([j-0.4 j-0.4],[i-0.4 i+0.4],'Color',hex2rgb('262626'),'LineWidth',1)
            line([j+0.4 j+0.4],[i-0.4 i+0.4],'Color',hex2rgb('262626'),'LineWidth',1)
            line([j-0.4 j+0.4],[i-0.4 i-0.4],'Color',hex2rgb('262626'),'LineWidth',1)
            line([j-0.4 j+0.4],[i+0.4 i+0.4],'Color',hex2rgb('262626'),'LineWidth',1)
        end
        
    end
end


end
