function [rou,max_rou,Delay,max_motion_delay,time_delay,isLeader,Leader_follow_matrix] = Cal_Rou_Perception_Anisotropy_2part_without_nan(Anisotropy,Frame_matrix,time_delay_slice,time_pair,tracks_filt,findRou,tau_threshold,correlation_threshold)

time_delay = zeros(1,length(time_delay_slice));
rou = zeros(2,length(time_delay_slice));
for i = 1 : length(time_delay_slice)
    temp = Frame_matrix(2,:);
    if time_delay_slice(i)>0
        temp(1:time_delay_slice(i)) = [];
        index_i = Frame_matrix(1,1:length(temp));
        index_j = temp;
        
        index = find(time_pair(:,2)-time_pair(:,1)==time_delay_slice(i));
        index = time_pair(index,:);
        
        x1_pre = tracks_filt(Frame_matrix(1,index(:,1)),2:4);
        x2_pre = tracks_filt(Frame_matrix(2,index(:,1)),2:4);
        
        x1_now = tracks_filt(Frame_matrix(1,index(:,2)),2:4);
        x2_now = tracks_filt(Frame_matrix(2,index(:,2)),2:4);
        
        a = (x2_pre - x1_pre)'./vecnorm((x2_pre - x1_pre)');
        b = (x2_now - x1_now)'./vecnorm((x2_now - x1_now)');
        
        %rou(2,i) = mean(acos(sum(a.*b)) / time_delay_slice(i));
        %rou(2,i) = mean(acos(sum(a.*b)) / (tracks_filt(index_j(1),5)-tracks_filt(index_i(1),5)));
        
%         v1_pre = tracks_filt(Frame_matrix(1,index(:,1)),6:8);
%         v1_pre = v1_pre'./vecnorm(v1_pre');
%         A1 = acos(sum(a.*b)) / (tracks_filt(index_j(1),5)-tracks_filt(index_i(1),5));
%         A2 = (1+sum(v1_pre.*a))/2;
%         rou(2,i) = mean(A1.*(A2.^Anisotropy));
        
        
        v1_pre = tracks_filt(Frame_matrix(1,index(:,1)),6:8);
        v1_pre = v1_pre'./vecnorm(v1_pre');
        v1_now = tracks_filt(Frame_matrix(1,index(:,2)),6:8);
        v1_now = v1_now'./vecnorm(v1_now');
        A1 = acos(sum(a.*b)) / (tracks_filt(index_j(1),5)-tracks_filt(index_i(1),5));
        A2 = (1+sum(v1_pre.*a))/2;
        A3 = (1+sum(v1_now.*b))/2;
        rou(2,i) = mean(A1.*(A2.^Anisotropy).*(A3.^Anisotropy));
        %rou(2,i) = mean(A1.*(A2.^Anisotropy));    
    else
        temp(end+time_delay_slice(i)+1:end) = [];
        index_i = Frame_matrix(1,end-length(temp)+1:end);
        index_j = temp;
    end
    time_delay(i) = tracks_filt(index_j(1),5)-tracks_filt(index_i(1),5);
    
    vi = tracks_filt(index_i,6:8)';
    vj = tracks_filt(index_j,6:8)';
    vi = vi./vecnorm(vi);
    vj = vj./vecnorm(vj);
    rou(1,i) =  mean(sum(vi.*vj));%sum(sum(vi.*vj)./(vecnorm(vi).*vecnorm(vj)))/size(temp,2);
end
index0 = find(time_delay_slice==0);
rou(2,1:index0-1) = rou(2,end:-1:index0+1);
rou(2,time_delay_slice==0) = nan;
rou(3,:) = time_delay;
rou(4,:) = zeros(1,length(time_delay));
rou(5,:) = zeros(1,length(time_delay));

%% Delay time
 
% follow_correlation = rou(1,:);
% [follow_peak,peak_index] = findpeaks(follow_correlation);
% index = peak_index>floor(length(time_delay_slice)/4)&peak_index<floor(3*length(time_delay_slice)/4);
% peak_index = peak_index(index);
% follow_peak = follow_peak(index);
% [max_peak,index] = max(follow_peak);
% max_delay_index = peak_index(index);
% rou(4,max_delay_index) = 1;

if strcmp(findRou,'peaks') | strcmp(findRou,'findpeaks')
    follow_correlation = rou(1,:);
    [follow_peak,peak_index] = findpeaks(follow_correlation);
    index = peak_index>floor(length(time_delay_slice)*tau_threshold)&peak_index<floor(length(time_delay_slice)*(1-tau_threshold));
    peak_index = peak_index(index);
    follow_peak = follow_peak(index);
    [max_peak,index] = max(follow_peak);
    max_delay_index = peak_index(index);
    rou(4,max_delay_index) = 1;

elseif strcmp(findRou,'max')
    % without nan
    follow_correlation = rou(1,:);
    temp_follow_correlation = follow_correlation;
    temp_follow_correlation(temp_follow_correlation<correlation_threshold)=-2;
    
    [follow_peak,peak_index] = max(temp_follow_correlation);
    max_delay_index = peak_index(1);
    index = max_delay_index>floor(length(time_delay_slice)*tau_threshold)&max_delay_index<floor(length(time_delay_slice)*(1-tau_threshold));
    max_delay_index = max_delay_index(index);
    rou(4,max_delay_index) = 1;

else
    error('Error!')
end

% without nan
% follow_correlation = rou(1,:);
% [follow_peak,peak_index] = max(follow_correlation);
% max_delay_index = peak_index(end);
% index = max_delay_index>floor(length(time_delay_slice)/4)&max_delay_index<floor(3*length(time_delay_slice)/4);
% max_delay_index = max_delay_index(index);
% rou(4,max_delay_index) = 1;

isLeader = zeros(2,1);
if length(max_delay_index)>0
    if time_delay_slice(max_delay_index)>0
        isLeader(1) = 1;
        
        temp2 = Frame_matrix(2,:);
        temp2(1:time_delay_slice(max_delay_index)) = [];
        temp1 = Frame_matrix(1,1:length(temp2));
        
    else
        isLeader(2) = 1;
        
        temp2 = Frame_matrix(2,:);
        temp2(end+time_delay_slice(max_delay_index)+1:end) = [];
        temp1 = Frame_matrix(1,end-length(temp2)+1:end);
    end
    Leader_follow_matrix = [temp1;temp2];
    Delay = time_delay(max_delay_index);%tracks_filt(Leader_follow_matrix(2,1),5) - tracks_filt(Leader_follow_matrix(1,1),5);
    max_rou = rou(1,max_delay_index);
else
    Delay = nan;
    max_rou = nan;
end
%% Visual Motion

visual_motion = rou(2,:);
max_visual_index = find(visual_motion==max(visual_motion));
if time_delay_slice(max_delay_index)>=0
    index = find(time_delay_slice(max_visual_index)>=0);
    max_visual_index = max_visual_index(index);
    
else
    index = find(time_delay_slice(max_visual_index)<0);
    max_visual_index = max_visual_index(index);
end
max_motion_delay = time_delay(max_visual_index);
rou(5,max_visual_index) = 1;

% visual_motion = rou(2,:);
% slope = visual_motion(2:end) - visual_motion(1:end-1);
% [~,max_visual_index] = find(abs(slope) == min(abs(slope)));
% if time_delay_slice(max_delay_index)>=0
%     index = find(time_delay_slice(max_visual_index)>=0);
%     max_visual_index = max_visual_index(index);
%     
% else
%     index = find(time_delay_slice(max_visual_index)<0);
%     max_visual_index = max_visual_index(index);
% end
% max_motion_delay = time_delay(max_visual_index);
% rou(5,max_visual_index) = 1;

% X = rou(1,:);
% Y = rou(2,:);
% t = rou(3,:);
% figure;hold on;
% yyaxis left
% plot(t,X,'-')
% plot(t(rou(4,:)==1),X(rou(4,:)==1),'o')
% yyaxis right
% plot(t,Y,'-')
% plot(t(rou(5,:)==1),Y(rou(5,:)==1),'o')

end
