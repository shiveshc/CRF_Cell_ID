%%%% make figures for position consistency
addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\violin'))

xtick = [1:size(uniq_annotated_neurons,1)];
ytick = [1:size(uniq_annotated_neurons,1)];
figure,imagesc(true_PA_matrix)
set(gca,'ytick',ytick,'yticklabel',uniq_annotated_neurons(sort_index_X(:,1),:),'FontSize',5)
% ytickangle(45)
set(gca,'xtick',xtick,'xticklabel',uniq_annotated_neurons(sort_index_X(:,1),:),'FontSize',5)
xtickangle(90)
colorbar

figure,imagesc(mean_PA_consistency_matrix)
set(gca,'ytick',ytick,'yticklabel',uniq_annotated_neurons(sort_index_X(:,1),:),'FontSize',5)
% ytickangle(45)
set(gca,'xtick',xtick,'xticklabel',uniq_annotated_neurons(sort_index_X(:,1),:),'FontSize',5)
xtickangle(90)
colorbar

figure,imagesc(true_LR_matrix)
set(gca,'ytick',ytick,'yticklabel',uniq_annotated_neurons(sort_index_Y(:,1),:),'FontSize',5)
% ytickangle(45)
set(gca,'xtick',xtick,'xticklabel',uniq_annotated_neurons(sort_index_Y(:,1),:),'FontSize',5)
xtickangle(90)
colorbar

figure,imagesc(mean_LR_consistency_matrix)
set(gca,'ytick',ytick,'yticklabel',uniq_annotated_neurons(sort_index_Y(:,1),:),'FontSize',5)
% ytickangle(45)
set(gca,'xtick',xtick,'xticklabel',uniq_annotated_neurons(sort_index_Y(:,1),:),'FontSize',5)
xtickangle(90)
colorbar

figure,imagesc(true_DV_matrix)
set(gca,'ytick',ytick,'yticklabel',uniq_annotated_neurons(sort_index_Z(:,1),:),'FontSize',5)
% ytickangle(45)
set(gca,'xtick',xtick,'xticklabel',uniq_annotated_neurons(sort_index_Z(:,1),:),'FontSize',5)
xtickangle(90)
colorbar

figure,imagesc(mean_DV_consistency_matrix)
set(gca,'ytick',ytick,'yticklabel',uniq_annotated_neurons(sort_index_Z(:,1),:),'FontSize',5)
% ytickangle(45)
set(gca,'xtick',xtick,'xticklabel',uniq_annotated_neurons(sort_index_Z(:,1),:),'FontSize',5)
xtickangle(90)
colorbar

figure,imagesc(true_angle_matrix)
set(gca,'ytick',ytick,'yticklabel',uniq_annotated_neurons(sort_index_Z(:,1),:),'FontSize',5)
% ytickangle(45)
set(gca,'xtick',xtick,'xticklabel',uniq_annotated_neurons(sort_index_Z(:,1),:),'FontSize',5)
xtickangle(90)
colorbar

figure,imagesc(mean_angle_matrix)
set(gca,'ytick',ytick,'yticklabel',uniq_annotated_neurons(sort_index_Z(:,1),:),'FontSize',5)
% ytickangle(45)
set(gca,'xtick',xtick,'xticklabel',uniq_annotated_neurons(sort_index_Z(:,1),:),'FontSize',5)
xtickangle(90)
colorbar

figure,boxplot(data_PA_relationship,true_PA_relationship)
xlabel('true PA relationship')
ylabel('observed PA relationship')
title('3,655 relationships across 4 worms')
ylim([-0.5,1.5])
index_0 = find(true_PA_relationship == 0);
index_1 = find(true_PA_relationship == 1);
figure,violin({data_PA_relationship(index_0,:),data_PA_relationship(index_1,:)},'xlabel',{'0','1'},'facecolor',[0,0,1],'edgecolor',[1,1,1],'facealpha','1','bw',0.04)
ylim([-0.5,1.5])
xlabel('true PA relationship')
ylabel('observed PA relationship')
title('3,655 relationships across 4 worms')
% figure,plotSpread({data_PA_relationship(index_0,:),data_PA_relationship(index_1,:)})
    
figure,boxplot(data_LR_relationship,true_LR_relationship)
xlabel('true LR relationship')
ylabel('observed LR relationship')
title('3,655 relationships across 4 worms')
ylim([-0.5,1.5])
index_0 = find(true_LR_relationship == 0);
index_1 = find(true_LR_relationship == 1);
figure,violin({data_LR_relationship(index_0,:),data_LR_relationship(index_1,:)},'xlabel',{'0','1'},'facecolor',[0,0,1],'edgecolor',[1,1,1],'facealpha','1','bw',0.04)
ylim([-0.5,1.5])
xlabel('true LR relationship')
ylabel('observed LR relationship')
title('3,655 relationships across 4 worms')

figure,boxplot(data_DV_relationship,true_DV_relationship)
xlabel('true DV relationship')
ylabel('observed DV relationship')
title('3,655 relationships across 4 worms')
ylim([-0.5,1.5])
hold on
index_0 = find(true_DV_relationship == 0);
index_1 = find(true_DV_relationship == 1);
figure,violin({data_DV_relationship(index_0,:),data_DV_relationship(index_1,:)},'xlabel',{'0','1'},'facecolor',[0,0,1],'edgecolor',[1,1,1],'facealpha','1','bw',0.04)
ylim([-0.5,1.5])
xlabel('true DV relationship')
ylabel('observed DV relationship')
title('3,655 relationships across 4 worms')
% scatter(ones(size(index_0,1),1)+(rand(size(index_0,1),1) - 0.5)/10,data_DV_relationship(index_0,:),'r','filled')
% scatter(2*ones(size(index_1,1),1)+(rand(size(index_0,1),1) - 0.5)/10,data_DV_relationship(index_1,:),'.r')


% histogram2(true_angle_relationship,data_angle_relationship,'NumBins',[40,40])
% xlabel('true angle')
% ylabel('observed angles')
% colorbar
% view(2)

figure,scatter(true_angle_relationship,data_angle_relationship,'.r')
xlabel('true angle')
ylabel('observed angles')
hold on
plot([1:1:90],[1:1:90],'--b','LineWidth',3)

figure,histogram(true_angle_relationship-data_angle_relationship,'Normalization','probability')
xlabel('angle deviation from true angle (deg)')
ylabel('fraction')

figure,scatter(true_angle_relationship,abs(true_angle_relationship - data_angle_relationship),'.r')