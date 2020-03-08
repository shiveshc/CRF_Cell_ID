%%%% function to plot potential of all neurons being assinged to a
%%%% particular state given one neuron's state is fixed
%%%% inputs - fixed neuron n1 and state1
%%%%        - state 2
% function plot_node_pot(n1,state1,state2)

n1 = 117;
state1 = find(strcmp('SIAVL',Neuron_head));
state2 = find(strcmp('AINL',Neuron_head));

scaled_r = imadjustn(mat2gray(max(thisimage_r,[],3)),[0,0.4],[0,1],1);
temp_img = -max(scaled_r,[],3);
temp_img = (temp_img - min(temp_img(:)))/(max(temp_img(:))-min(temp_img(:)));
score_img = zeros(size(scaled_r));

n1_at_end1 = find(edgeStruct.edgeEnds(:,1) == n1);
scores_n1_at_end1 = edge_pot(state1,state2,n1_at_end1);
scores_n1_at_end1 = scores_n1_at_end1(:);
n1_at_end2 = find(edgeStruct.edgeEnds(:,2) == n1);
scores_n1_at_end2 = edge_pot(state2,state1,n1_at_end2);
scores_n1_at_end2 = scores_n1_at_end2(:);
other_end = edgeStruct.edgeEnds(n1_at_end1,2);
other_end = [other_end;edgeStruct.edgeEnds(n1_at_end2,1)];
[n1_edges,sort_index] = sort([n1_at_end1;n1_at_end2]);
n1_other_end = other_end(sort_index,:);
scores = [scores_n1_at_end1;scores_n1_at_end2];
scores = scores(sort_index,:);
scores = (scores - min(scores))/(max(scores) - min(scores));

cmap = hot(100);
col1 = [0,0,0];
col2 = cmap(1,:);
cmap_app = [linspace(col1(1,1),col2(1,1),25);linspace(col1(1,2),col2(1,2),25);linspace(col1(1,3),col2(1,3),25)]';
cmap_app = [cmap_app;cmap];
intervals = linspace(0,1,size(cmap_app,1));

for i = 1:size(scores,1)
    for n = 1:size(intervals,2)-1
        if scores(i,1) > intervals(1,n) && scores(i,1) <= intervals(1,n+1)
            scores(i,2) = n;
        end
    end
end
allCenters_wo_n1 = allCenters;
allCenters_wo_n1{1,1}(n1) = [];
for i = 1:size(allCenters_wo_n1{1,1},2)
    curr_PixelList = allCenters_wo_n1{1,1}(i).PixelList;
    for n = 1:size(curr_PixelList,1)
        score_img(curr_PixelList(n,1),curr_PixelList(n,2),curr_PixelList(n,3)) = scores(i,2);
    end
end

B = labeloverlay(temp_img,max(score_img,[],3),'Colormap',cmap_app,'Transparency',0);
figure,imshow(B,[])
hold on
text(30,30,Neuron_head{state2,1},'Color',[0,1,1],'FontSize',14)

figure,imshow(max(mat2gray(score_img),[],3),[])
colormap(cmap_app)
colorbar