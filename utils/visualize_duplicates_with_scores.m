function visualize_duplicates_with_scores(thisimage_r,duplicates,scores,allCenters,labeled_img_r,segNorm,name)

scaled_r = imadjustn(mat2gray(max(thisimage_r,[],3)),[0,0.4],[0,1],1);
temp_img = -max(scaled_r,[],3);
temp_img = (temp_img - min(temp_img(:)))/(max(temp_img(:))-min(temp_img(:)));
score_img = zeros(size(scaled_r));
lim1 = max(0,min(scores)-0.03);
lim2 = min(max(scores)+0.01,1);

cmap = hot(100);
col1 = [0,0,0];
col2 = cmap(1,:);
cmap_app = [linspace(col1(1,1),col2(1,1),25);linspace(col1(1,2),col2(1,2),25);linspace(col1(1,3),col2(1,3),25)]';
cmap_app = [cmap_app;cmap];
% intervals = linspace(min(floor(scores))-0.1,max(ceil(scores))+0.1,size(cmap_app,1));
intervals = linspace(lim1,lim2,size(cmap_app,1));
for i = 1:size(scores,1)
    for n = 1:size(intervals,2)-1
        if scores(i,1) > intervals(1,n) && scores(i,1) <= intervals(1,n+1)
            scores(i,2) = n;
        end
    end
end

fixed_val = 0.94;
for n = 1:size(intervals,2)-1
    if fixed_val > intervals(1,n) && fixed_val <= intervals(1,n+1)
        fixed_val = n;
    end
end
% for i = 1:size(allCenters{1,1},2)
%     curr_PixelList = allCenters{1,1}(i).PixelList;
%     for n = 1:size(curr_PixelList,1)
%         score_img(curr_PixelList(n,1),curr_PixelList(n,2),curr_PixelList(n,3)) = fixed_val;
%     end
% end
for i = 1:size(duplicates,1)
    curr_PixelList = allCenters{1,1}(duplicates(i)).PixelList;
    for n = 1:size(curr_PixelList,1)
        score_img(curr_PixelList(n,1),curr_PixelList(n,2),curr_PixelList(n,3)) = scores(i,2);
    end
end



figure,imshow(max(score_img,[],3),[])
colormap(cmap_app)
caxis([lim1,lim2])
colorbar

B = labeloverlay(temp_img,max(score_img,[],3),'Colormap',cmap_app,'Transparency',0);
% C = imfuse(scaled_r,max(score_img,[],3),'blend');
% B = labeloverlay(scaled_r,max(score_img,[],3));
figure,imshow(B,[])
% figure,imshow(C,[])
hold on
text(30,30,name,'Color',[0,1,1],'FontSize',14)