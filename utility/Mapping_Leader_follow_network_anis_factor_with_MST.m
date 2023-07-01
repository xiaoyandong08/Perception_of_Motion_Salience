function [Delay,Delay1,max_motion_delay,rou,max_rou] = Mapping_Leader_follow_network_anis_factor_with_MST(Frame_matrix,tracks_filt,anis_factor,time_pair,findRou,tau_threshold,correlation_threshold)

Delay = zeros(size(Frame_matrix,1),size(Frame_matrix,1));
max_motion_delay = zeros(size(Frame_matrix,1),size(Frame_matrix,1));
max_rou = zeros(size(Frame_matrix,1),size(Frame_matrix,1));

for i = 1 : size(Frame_matrix,1)
    for j = 1 : size(Frame_matrix,1)
        if i ~= j
            temp_Frame_matrix = Frame_matrix([i j],:);
            index = find(sum(temp_Frame_matrix~=0,1) == 2);
            temp_Frame_matrix = temp_Frame_matrix(:,index);
            if size(temp_Frame_matrix,2)>0 % >=20
                time_delay_slice = [-(size(temp_Frame_matrix,2)-1):size(temp_Frame_matrix,2)-1];
                if size(temp_Frame_matrix,2)>1000
                    one_time_pair = nchoosek(1:size(temp_Frame_matrix,2),2);
                else
                    one_time_pair = time_pair{size(temp_Frame_matrix,2)};
                end
                %[rou,Delay(i,j),max_motion_delay(i,j)] = Cal_Rou_Perception_v3_anis_factor_withoutNaN(temp_Frame_matrix,time_delay_slice,one_time_pair,tracks_filt,anis_factor);
%                 profile on;
                [rou{i,j},max_rou(i,j),Delay(i,j),max_motion_delay(i,j)] = Cal_Rou_Perception_Anisotropy_2part_without_nan_with_MST(anis_factor,temp_Frame_matrix,time_delay_slice,one_time_pair,tracks_filt,findRou,tau_threshold,correlation_threshold);
%                 profile viewer;
            end
            %Plot_trajectory(rou,time_delay_slice,temp_Frame_matrix,tracks_filt,all_time(index))

        end
    end
end

Delay1 = Delay;


Delay1(abs(Delay1)<0.0) = 0;
Delay1(Delay1>0) = 0;

Delay2 = Delay1;
Delay2(1:size(Delay2,1)+1:size(Delay2,1)^2) = 1;

% nestedness = Nestedness.NODF(abs(sign(Delay2)));
% nestedness = nestedness.N;
% 
% 
% isSymtry = sum(sum(triu(Delay1)' + tril(Delay1)))==0;
% 
% 
% Gr = global_reaching_centrality(sign(abs(Delay1)));
    
% for i = size(Frame_matrix,2) : size(Frame_matrix,2)
%     Id = Frame_matrix(find(Frame_matrix(:,i)>0),i);
%     if length(Id)>0
%         v_xyz = tracks_filt(Id,6:8);
%         v_xyz = v_xyz./vecnorm(v_xyz')';
%         order(i) = norm(sum(v_xyz,1))/size(v_xyz,1);
%         
%         p_xyz = tracks_filt(Id,2:4);
%         dist = squareform(pdist(p_xyz));
%         density(i) = 6*size(p_xyz,1)/(pi*mean(max(dist,[],2))^3);
%         
%         order_each = nan(1,length(Id));
%         for j = 1 : length(Id)
%             part_v_xyz = v_xyz(setdiff(1 : length(Id),j),:);
%             part_order = norm(sum(part_v_xyz,1))/size(part_v_xyz,1);
%             order_each(1,j) = part_order - order(i);
%         end
% %         [a,b] = min(order_each);
% %         v_xyz(setdiff(1 : length(Id),b),:)
% %         v_xyz(b,:) 
%     end
%     
%     all_order_each{i} = order_each;
% end
    
end