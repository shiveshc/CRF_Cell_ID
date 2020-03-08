%%%% function to plot potential of all neurons being assinged to a
%%%% particular state
% function plot_node_pot(state)

state = find(strcmp('AWBL',Neuron_head));
scaled_r = imadjustn(mat2gray(max(thisimage_r,[],3)),[0,0.4],[0,1],1);
temp_img = -max(scaled_r,[],3);
temp_img = (temp_img - min(temp_img(:)))/(max(temp_img(:))-min(temp_img(:)));
score_img = zeros(size(temp_img));
scores = node_pot(:,state);

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
for i = 1:size(allCenters{1,1},2)
    curr_PixelList = allCenters{1,1}(i).PixelList;
    for n = 1:size(curr_PixelList,1)
        score_img(curr_PixelList(n,1),curr_PixelList(n,2),curr_PixelList(n,3)) = scores(i,2);
    end
end

B = labeloverlay(temp_img,max(score_img,[],3),'Colormap',cmap_app,'Transparency',0);
figure,imshow(B,[])
hold on
text(30,30,Neuron_head{state,1},'Color',[0,1,1],'FontSize',14)

figure,imshow(max(mat2gray(score_img),[],3),[])
colormap(cmap_app)
colorbar