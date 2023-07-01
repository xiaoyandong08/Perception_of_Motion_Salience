function [temporal_Mij, Delay,Delay_neg,corr_M_LF] ...
    = Bird_Motion_salience_vs_LF_Anisotropy_without_nan_L2F(Anisotropy,Corr_Type,tau_slice,time_slice,Frame_matrix,tracks_filt,findRou,tau_threshold,correlation_threshold)

temporal_Mij = cell(length(tau_slice),length(time_slice));
Delay = cell(length(tau_slice),length(time_slice));
Delay_neg = cell(length(tau_slice),length(time_slice));

corr_M_LF  = nan(length(tau_slice),length(time_slice));
temporal_ave_leading_time  = nan(length(tau_slice),length(time_slice));
temporal_ave_MS_time = nan(length(tau_slice),length(time_slice));

for i = 1 : length(tau_slice)
    i;
    for j = 1 : length(time_slice)
        if tau_slice(i) <  time_slice(j)
            [temporal_Mij{i,j}] = Calculate_Mij(Anisotropy,Frame_matrix(:,[time_slice(j)-tau_slice(i)+1 time_slice(j)]),tracks_filt);
            [a,b] = sort(nansum(temporal_Mij{i,j}));

            %[Delay{i,j},Delay_neg{i,j},max_motion_delay{i,j},~,~] = Mapping_Leader_follow_network_anis_factor(Frame_matrix(:,time_slice(j)-tau_slice(i)+1:time_slice(j)),tracks_filt,Anisotropy,time_pair);
            [Delay{i,j},Delay_neg{i,j}] = Mapping_Leader_follow_network_anis_factor_simplified(Frame_matrix(:,time_slice(j)-tau_slice(i)+1:time_slice(j)),tracks_filt,findRou,tau_threshold,correlation_threshold);

            
            [Gr,Cr] = global_reaching_centrality(sign(abs(Delay_neg{i,j}))');
            
%             [~, index_rows] = sort(sum(abs(sign(Delay_neg{i,j})),2),'descend');
%             [~, index_cols] = sort(sum(abs(sign(Delay_neg{i,j})),1),'descend');
%             nested_matrix = Delay_neg{i,j}(index_rows,index_cols);
%             figure;imagesc(sign(nested_matrix))
            
            if strcmp(Corr_Type,'Spearman')
                corr_M_LF(i,j) = corr(nanmean(temporal_Mij{i,j},1)',Cr,'Type','Spearman');
            elseif strcmp(Corr_Type,'Pearson')
                corr_M_LF(i,j) = corr(nanmean(temporal_Mij{i,j},1)',Cr);
            end


            
%             figure
%             node_size = nanmean(temporal_Mij{i,j},1)';
%             node_size = node_size/max(node_size)*20;
%             G = digraph(abs(Delay_neg{i,j}));
%             LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
%             p = plot(G,'layout','layered','NodeLabel',{},'LineWidth',LWidths,'ArrowSize',12);axis off
%             p.MarkerSize = 2*indegree(G) + 14;
%             p.MarkerSize = node_size;
%             text(p.XData,p.YData,num2str([1:size(Frame_matrix,1)]'),'Color','w','HorizontalAlignment','center','FontSize',14);
%             title(['Corr-Cr-M = '  num2str(round( corr_M_LF(i,j),3))])
%             set(gca,'fontsize',14)
%             set(gcf,'position',[216 324 481 383])
%             
%             Nestedness.NODF(abs(sign(Delay1{i,j})))
%             PlotWebs.font_size = 14;
%             figure;subplot(121);PlotWebs.PLOT_MATRIX(abs(sign(Delay1{i,j})));
%             subplot(122);PlotWebs.PLOT_NESTED_MATRIX(abs(sign(Delay1{i,j})));
            
           
        end
    end
end

end



function [Mij] = Calculate_Mij(Anisotropy,Frame_matrix,tracks_filt)

Mij = nan(size(Frame_matrix,1),size(Frame_matrix,1));


for i = 1 : size(Frame_matrix,1)
    for j = 1 : size(Frame_matrix,1)
        if i ~= j
            temp_Frame_matrix = Frame_matrix([i j],:);
            index = find(sum(temp_Frame_matrix~=0,1) == 2);
            temp_Frame_matrix = temp_Frame_matrix(:,index);
            
            x1_pre = tracks_filt(temp_Frame_matrix(1,1),2:4);
            x2_pre = tracks_filt(temp_Frame_matrix(2,1),2:4);
            
            x1_now = tracks_filt(temp_Frame_matrix(1,2),2:4);
            x2_now = tracks_filt(temp_Frame_matrix(2,2),2:4);
            
            a = (x2_pre - x1_pre)'./vecnorm((x2_pre - x1_pre)');
            b = (x2_now - x1_now)'./vecnorm((x2_now - x1_now)');
            
            vi = tracks_filt(temp_Frame_matrix(1,2),6:8)';
            vj = tracks_filt(temp_Frame_matrix(2,2),6:8)';
            vi = vi./vecnorm(vi);
            vj = vj./vecnorm(vj);
            
            
            %Mij(i,j) = acos(a'*b) / (tracks_filt(temp_Frame_matrix(1,2),5)-tracks_filt(temp_Frame_matrix(1,1),5));

%             v1_pre = tracks_filt(temp_Frame_matrix(1,1),6:8);
%             v1_pre = v1_pre'./vecnorm(v1_pre');
%             A1 = acos(a'*b) / (tracks_filt(temp_Frame_matrix(1,2),5)-tracks_filt(temp_Frame_matrix(1,1),5));
%             A2 = (1+sum(v1_pre.*a))/2;
%             Mij(i,j) = A1.*(A2.^Anisotropy);

            v1_pre = tracks_filt(temp_Frame_matrix(1,1),6:8);
            v1_pre = v1_pre'./vecnorm(v1_pre');
            v1_now = tracks_filt(temp_Frame_matrix(1,2),6:8);
            v1_now = v1_now'./vecnorm(v1_now');
            A1 = acos(a'*b) / (tracks_filt(temp_Frame_matrix(1,2),5)-tracks_filt(temp_Frame_matrix(1,1),5));
            A2 = (1+sum(v1_pre.*a))/2;
            A3 = (1+sum(v1_now.*b))/2;
            Mij(i,j) = A1.*(A2.^Anisotropy)*(A3.^Anisotropy);
        end
    end
end
   
end
