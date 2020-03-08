%%%% function to create merged 3D image of NeuroPAL strain from R, G, and B
%%%% channels. This 3D merged image will be used to create annotation in
%%%% Vaa3D (neuron names and neuron positions).

%%%% Note use the same normalization as the one  used for creating the rgb
%%%% image. Make r_norm, g_norm, b_norm using rgb_img_all_files_in_foder.m

% merge_img = zeros(size(thisimage_g));
% for i = 1:size(thisimage_g,3)
%     merge_img(:,:,i) = max(cat(3,r_norm(:,:,i),g_norm(:,:,i),b_norm(:,:,i)),[],3);
% end
% indexed_img_to_tiff(merge_img(:,:,:,1),[],['C:\Users\Shivesh\Dropbox
% (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\21\merge_img.tif'])
scaled_r = imadjustn(mat2gray(thisimage_r),[0,0.6],[0,1],0.8);
scaled_g = imadjustn(mat2gray(thisimage_g),[0,0.1],[0,1],1);
scaled_b = imadjustn(mat2gray(thisimage_b),[0,0.2],[0,1],1);
index = 1:1:size(thisimage_g,1)*size(thisimage_g,2);
index_img = reshape(index,size(thisimage_g,1),size(thisimage_g,2));

for i = 1:size(thisimage_g,3)
%     g_norm = imadjust(scaled_g(:,:,i),[0;0.7],[0;1],0.8);
%     b_norm = imadjust(scaled_b(:,:,i),[0;0.3],[0;1],1);
%     r_norm = imadjust(scaled_r(:,:,i),[0;0.3],[0;1],1);
    g_norm = thisimage_g(:,:,i);
    b_norm = thisimage_b(:,:,i);
    r_norm = thisimage_r(:,:,i);
    rgb_img = cat(3,r_norm,g_norm,b_norm);
    
%     rgb_img = im2uint8(rgb_img);    
%     [rgb_ind,cmap] = rgb2ind(rgb_img,256,'nodither');
    
%     rgb_img = imadjust(rgb_img,[],[],[0.9,0.9,0.9]);
    [rgb_ind,cmap] = rgb2ind(rgb_img,65536);
%     cmap = cat(2,r_norm(:),g_norm(:),b_norm(:));
    indexed_img_to_tiff(rgb_ind,cmap,['C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\13\index_images\',num2str(i),'.tif'])
end