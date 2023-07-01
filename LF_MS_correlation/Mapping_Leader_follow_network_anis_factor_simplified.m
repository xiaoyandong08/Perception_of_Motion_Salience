function [Delay,Delay1,rou,max_rou] = Mapping_Leader_follow_network_anis_factor_simplified(Frame_matrix,tracks_filt,findRou,tau_threshold,correlation_threshold)

Delay = zeros(size(Frame_matrix,1),size(Frame_matrix,1));
max_rou = zeros(size(Frame_matrix,1),size(Frame_matrix,1));

for i = 1 : size(Frame_matrix,1)
    for j = 1 : size(Frame_matrix,1)
        if i ~= j
            temp_Frame_matrix = Frame_matrix([i j],:);
            index = find(sum(temp_Frame_matrix~=0,1) == 2);
            temp_Frame_matrix = temp_Frame_matrix(:,index);
            if size(temp_Frame_matrix,2)>0 % >=20
                time_delay_slice = [-(size(temp_Frame_matrix,2)-1):size(temp_Frame_matrix,2)-1];
                %[rou,Delay(i,j),max_motion_delay(i,j)] = Cal_Rou_Perception_v3_anis_factor_withoutNaN(temp_Frame_matrix,time_delay_slice,one_time_pair,tracks_filt,anis_factor);
%                 profile on;
                [rou{i,j},max_rou(i,j),Delay(i,j)] = Cal_Rou_Perception_Anisotropy_2part_without_nan_simplified(temp_Frame_matrix,time_delay_slice,tracks_filt,findRou,tau_threshold,correlation_threshold);
%                 profile viewer;
            end

        end
    end
end

Delay1 = Delay;
Delay1(isnan(Delay1)) = 0;
Delay1(abs(Delay1)<0.0) = 0;
Delay1(Delay1>0) = 0;

Delay2 = Delay1;
Delay2(1:size(Delay2,1)+1:size(Delay2,1)^2) = 1;

% minn = floor(length(time_delay_slice)*tau_threshold);
% maxx = floor(length(time_delay_slice)*(1-tau_threshold));
% figure;
% for i = 1 : size(Delay1,1)
%     for j = 1 : size(Delay1,2)
%         if j ~= i
%             %subplot(size(Delay1,1),size(Delay1,2),(i-1)*size(Delay1,2)+j)
%             subaxis(size(Delay1,1),size(Delay1,2),(i-1)*size(Delay1,2)+j,'SpacingVertical',0.04,'SpacingHorizontal',0.02,'MarginLeft',.02,'MarginRight',.02,'MarginTop',.04,'MarginBottom',.04)
%             
%             hold on
%             plot(rou{i,j}(2,:),rou{i,j}(1,:),'-')
%             plot(Delay(i,j),max_rou(i,j),'ro')
%             xlim([rou{1,2}(2,1), rou{1,2}(2,end)])
%             ax = gca;
%             line([rou{i,j}(2,minn) rou{i,j}(2,minn)],ax.YLim,'color','r','linestyle','-')
%             line([rou{i,j}(2,maxx) rou{i,j}(2,maxx)],ax.YLim,'color','r','linestyle','-')
%         end
%     end
% end

end