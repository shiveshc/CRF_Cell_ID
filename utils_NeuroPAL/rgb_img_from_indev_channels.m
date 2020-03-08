%%% create rgb image of NeuroPAL strain
%%% Input - thisimage_g, thisimage_r, thisimage_b (read using
%%% 'readtif_Brucker_v2.m' or 'readCyOFP_SD.m')

g_norm = imadjust(max(mat2gray(thisimage_g),[],3),[0,0.2],[0;1]);
b_norm = imadjust(max(mat2gray(thisimage_b),[],3),[0,0.2],[0;1]);
r_norm = imadjust(max(mat2gray(thisimage_r),[],3),[0,0.2],[0;1]);
figure,imshow(g_norm)
figure,imshow(b_norm)
figure,imshow(r_norm)
rgb_img = cat(3,r_norm,g_norm,b_norm);
% rgb_img(30:31,330:345,:) = 1;
figure,imshow(rgb_img)
saveas(gcf,['C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190702_OH15500\3\','rgb_img.png'])