function Generate_tracks_filt_with_block_border_Jackdaw
addpath('../utility')

file_tag = 'mob_01';

load(['../data_original/' file_tag '.mat']);

all_time = unique(tracks_filt(:,5));
all_bird = unique(tracks_filt(:,1));

bird_time_matrix = zeros(length(all_bird),length(all_time));
for i = 1 : length(all_bird)
    for j = 1 : length(all_time)
        temp = intersect(find(all_bird(i)==tracks_filt(:,1)), find(all_time(j)==tracks_filt(:,5)) );
        if ~isempty(temp)
            bird_time_matrix(i,j) = temp;
        end
    end
end

u=tracks_filt(:,6:8);
U=mean(u,1);   
theta_temp=atan(U(2)/U(1));
R = [cos(theta_temp) -sin(theta_temp); sin(theta_temp) cos(theta_temp)];
UR=U(1:2)*R;UR(3)=U(3);
for j=1:size(tracks_filt,1)
    if UR(1)>0
        tracks_filt(j,2:3) = tracks_filt(j,2:3)*R; 
        tracks_filt(j,6:7)=tracks_filt(j,6:7)*R;
    elseif UR(1)<0
        tracks_filt(j,2:3) = tracks_filt(j,2:3)*R*(-1); 
        tracks_filt(j,6:7)=tracks_filt(j,6:7)*R*(-1);  
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fraction_new_present_bird = 0.25;

previous_birds_in_frame = [];
count = 1;
for i = 1 : size(bird_time_matrix,2)
    Frame = i;
    birds_in_frame = find(bird_time_matrix(:,Frame));
    if length(birds_in_frame)>10 & length(setdiff(birds_in_frame,previous_birds_in_frame))>length(birds_in_frame)*fraction_new_present_bird
        
        block_row = birds_in_frame;
        block_col = [1:size(bird_time_matrix,2)];

        Frame_matrix = bird_time_matrix(birds_in_frame,:);
        index = find(sum(Frame_matrix~=0,1)>size(Frame_matrix,1)/2);
        Frame_matrix = Frame_matrix(:,index);
        temp_all_time = all_time(index);
        
        block_col = block_col(index);
        block_row(sum(Frame_matrix~=0,2)<mean(sum(Frame_matrix~=0,2)),:) = [];
        Frame_matrix(sum(Frame_matrix~=0,2)<mean(sum(Frame_matrix~=0,2)),:) = [];
        
        index = sum(logical(Frame_matrix),1) == size(Frame_matrix,1);
        Frame_matrix = Frame_matrix(:,index);
        block_col = block_col(index);
        if size(Frame_matrix,1)>0
            all_Frame_matrix{count} = Frame_matrix;
            all_Frame_time{count} = temp_all_time(index);
            which_Frame(count) = i;
            all_block_col{count} = block_col;
            all_block_row{count} = block_row;
            count = count + 1;
        end
        previous_birds_in_frame = birds_in_frame;
    end
end

clear temp
for i = 1 : length(all_Frame_matrix)
    temp(i) = sum(all_Frame_matrix{i}(:));
end
[a,b,c] = unique(temp,'stable');
all_Frame_matrix = all_Frame_matrix(b);
all_Frame_time = all_Frame_time(b);
which_Frame = which_Frame(b);
all_block_col = all_block_col(b);
all_block_row = all_block_row(b);
sprintf('Delete %d frames',length(c)-length(b))

for k = 1 : length(all_Frame_matrix)
    
    Frame_matrix = all_Frame_matrix{k};
    order = nan(1,size(Frame_matrix,2));
    density = nan(1,size(Frame_matrix,2));
    for i = 1 : size(Frame_matrix,2)
        Id = Frame_matrix(find(Frame_matrix(:,i)>0),i);
        if length(Id)>0
            v_xyz = tracks_filt(Id,6:8);
            v_xyz = v_xyz./vecnorm(v_xyz')';
            order(i) = norm(sum(v_xyz,1))/size(v_xyz,1);
            
            p_xyz = tracks_filt(Id,2:4);
            dist = squareform(pdist(p_xyz));
            density(i) = 6*size(p_xyz,1)/(pi*mean(max(dist,[],2))^3);
            all_dist{k}{i} = dist;
        end
    end
    
    all_order{k} = order;
    all_density{k} = density;
end


save(['../data_original/' 'Jackdaw_' file_tag '_all_frame_new_present_bird_' num2str(fraction_new_present_bird) '.mat'],'tracks_filt','bird_time_matrix',...
    'all_dist','all_Frame_matrix','all_Frame_time','which_Frame','all_order','all_density','fraction_new_present_bird','all_block_col','all_block_row','-v7.3')


cmap = [[1 1 1];hex2rgb('262626')];
figure;imagesc(all_time,[1:1:size(bird_time_matrix,1)],logical(bird_time_matrix))
colormap(cmap)
xlabel('Time (s)');ylabel('Bird index')
hold on;
for i = 1 : length(all_Frame_matrix)
    for j = 1 : size(all_Frame_matrix{i},1)
        plot([all_time(all_block_col{i}(1)) all_time(all_block_col{i}(end))],[all_block_row{i}(j) all_block_row{i}(j)],'r-','LineWidth',0.5)
    end
end
for i = 1 : length(all_Frame_matrix)
    plot([all_time(all_block_col{i}(1)) all_time(all_block_col{i}(end))],[all_block_row{i}(1) all_block_row{i}(1)],'r-','LineWidth',1)
    plot([all_time(all_block_col{i}(1)) all_time(all_block_col{i}(1))],[all_block_row{i}(1) all_block_row{i}(end)],'r-','LineWidth',1)
    plot([all_time(all_block_col{i}(1)) all_time(all_block_col{i}(end))],[all_block_row{i}(end) all_block_row{i}(end)],'r-','LineWidth',1)
    plot([all_time(all_block_col{i}(end)) all_time(all_block_col{i}(end))],[all_block_row{i}(1) all_block_row{i}(end)],'r-','LineWidth',1)

end
set(gca,'fontsize',14)

figure;
for k = 1 : length(all_Frame_matrix)
    
    subplot(ceil(sqrt(length(all_Frame_matrix))),ceil(sqrt(length(all_Frame_matrix))),k)
    set(gca,'FontSize',12,'TickLength',[0.03, 0.01],'XMinorTick','on','YMinorTick','on','boxstyle','full'); 
    view([-162 78]);
    
    Frame_matrix = all_Frame_matrix{k};
    Color = jet(size(Frame_matrix,2));
    xlim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2)))])
    ylim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3)))])
    zlim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4)))])
    hold on;box on
    for i = 1 : size(Frame_matrix,2)
        Id = Frame_matrix(find(Frame_matrix(:,i)>0),i);
        xyz=tracks_filt(Id,2:4);
        plot3(xyz(:,1),xyz(:,2),xyz(:,3),'.','color',Color(i,:));
    end
    title(['Frame time = ' num2str(round(all_time(which_Frame(k)),4),'%4.4f') ' s'],'fontweight','normal','fontsize',12)
end
set(gcf,'position',[65 672 1804 849])

end