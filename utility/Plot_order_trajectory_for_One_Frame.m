function Plot_order_trajectory_for_One_Frame(Frame_matrix,all_time,tracks_filt,Frame)

order = nan(1,size(Frame_matrix,2));
moment = nan(1,size(Frame_matrix,2));
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

h = 1;
for id = 1 : size(Frame_matrix,1)
    id_xyz = tracks_filt(Frame_matrix(id,:),2:4);
    r = id_xyz';
    r1 = gradient(r)./h;
    r2 = gradient(r1)./h;
    r1 = r1';r2 = r2';
    v = cross(r1,r2,2);
    c = vecnorm(v');
    d = vecnorm(r1');
    k = (c./(d.^3));
    curvature(id,:) = k;
end
all_curvature = nanmean(curvature,1);

Color = jet(size(Frame_matrix,2));
figure('units','inches','position',[5 2 10 8]);
box on
subplot('position',[0.1 0.85 0.8 0.12]);
hold on;box on;
% for i = 1 : length(order)
%     scatter(all_time(i),order(i),8,Color(i,:),'o','filled')
% end
p1 = plot(all_time,order,'LineWidth',2);
cd = [uint8(Color*255) uint8(ones(size(Color,1),1))]';
drawnow
set(p1.Edge,'ColorBinding','interpolated','ColorData',cd);
%xlabel('Time');
ylabel('Order')
legend(['Ave order = ' num2str(round(mean(order),3))],'location','best');legend boxoff
set(gca,'fontsize',14)

subplot('position',[0.1 0.68 0.8 0.12])
p1 = plot(all_time,moment,'LineWidth',2);
cd = [uint8(Color*255) uint8(ones(size(Color,1),1))]';
drawnow
set(p1.Edge,'ColorBinding','interpolated','ColorData',cd);
% xlabel('Time');
ylabel('Momentum')
legend(['Ave moment = ' num2str(round(mean(moment),3))],'location','best');legend boxoff
set(gca,'fontsize',14)

subplot('position',[0.1 0.51 0.8 0.12])
p1 = plot(all_time,all_curvature,'LineWidth',2);
cd = [uint8(Color*255) uint8(ones(size(Color,1),1))]';
drawnow
set(p1.Edge,'ColorBinding','interpolated','ColorData',cd);
xlabel('Time (s)');
ylabel('Curvature')
legend(['Ave curvature = ' num2str(round(mean(all_curvature),3))],'location','best');legend boxoff
set(gca,'fontsize',14)

subplot('position',[0.1 0.04 0.8 0.36])
% xlabel('x');ylabel('y');zlabel('z')
set(gca,'FontSize',12,'TickLength',[0.03, 0.01],...
    'XMinorTick','on','YMinorTick','on','boxstyle','full'); 
view([-162 78])
title(['Frames = ' num2str(Frame)])
% plot3(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2),tracks_filt(Frame_matrix(Frame_matrix(:)>0),3),tracks_filt(Frame_matrix(Frame_matrix(:)>0),4),'.')
xlim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2)))])
ylim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3)))])
zlim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4)))])
hold on;box on
for i = 1 : size(Frame_matrix,2)
    
    Id = Frame_matrix(find(Frame_matrix(:,i)>0),i);
    xyz=tracks_filt(Id,2:4);
    v_xyz = tracks_filt(Id,6:8);
    plot3(xyz(:,1),xyz(:,2),xyz(:,3),'.','color',Color(i,:));
    
%     if i==1 & length(unique(sum(logical(Frame_matrix),1)))==1
%         text(xyz(:,1),xyz(:,2),xyz(:,3),num2str([1:size(Frame_matrix,1)]'))
%     end
    if i == size(Frame_matrix,2)
        scatter3(xyz(:,1),xyz(:,2),xyz(:,3),65,[0.2 0.2 0.2],'o','filled');
    end
%     im(i) = getframe;
end
% for i = 1 : size(Frame_matrix,1)
%     Id = Frame_matrix(i,find(Frame_matrix(i,:)>0));
%     if length(Id)>0
%         xyz=tracks_filt(Id(end),2:4);
%         scatter3(xyz(:,1),xyz(:,2),xyz(:,3),65,[0.2 0.2 0.2],'o','filled');
%     end
% end
set(gca,'fontsize',14)

if exist('im')
    a=VideoWriter('0.5_0.45','MPEG-4');
    a.FrameRate = 10;
    a.Quality = 100;
    open(a);
    writeVideo(a,im);
    close(a)
end

end