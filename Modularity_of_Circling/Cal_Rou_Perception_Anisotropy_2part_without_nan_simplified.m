function [rou,max_rou,Delay,time_delay,isLeader,Leader_follow_matrix] = Cal_Rou_Perception_Anisotropy_2part_without_nan_simplified(Frame_matrix,time_delay_slice,tracks_filt,findRou,tau_threshold,correlation_threshold)

time_delay = zeros(1,length(time_delay_slice));
rou = zeros(2,length(time_delay_slice));
for i = 1 : length(time_delay_slice)
    temp = Frame_matrix(2,:);
    if time_delay_slice(i)>0
        temp(1:time_delay_slice(i)) = [];
        index_i = Frame_matrix(1,1:length(temp));
        index_j = temp;
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
rou(2,:) = time_delay;

%% Delay time

if strcmp(findRou,'peaks') | strcmp(findRou,'findpeaks')
    follow_correlation = rou(1,:);
    [follow_peak,peak_index] = findpeaks(follow_correlation);
    index = peak_index>floor(length(time_delay_slice)*tau_threshold)&peak_index<floor(length(time_delay_slice)*(1-tau_threshold));
    peak_index = peak_index(index);
    follow_peak = follow_peak(index);
    [max_peak,index] = max(follow_peak);
    max_delay_index = peak_index(index);
elseif strcmp(findRou,'max')
    % without nan
    follow_correlation = rou(1,:);
    temp_follow_correlation = follow_correlation;
    temp_follow_correlation(temp_follow_correlation<correlation_threshold)=-2;
    
    [follow_peak,peak_index] = max(temp_follow_correlation);
    max_delay_index = peak_index(1);
    index = max_delay_index>floor(length(time_delay_slice)*tau_threshold)&max_delay_index<floor(length(time_delay_slice)*(1-tau_threshold));
    max_delay_index = max_delay_index(index);
else
    error('Error!')
end

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




end
