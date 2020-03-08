%%%% Code to import CyOFP data for landmarks acuired at Spinning Disk 
%%%% Note here reading order of zplanes is flipped similar to
%%%% readtif_Brucker_v2.m for labelling consistency.
%%%% Inputs - 

% Create a single file by combining volocity exported individual files
% each timepoint in a single folder in the directory

% direc = 'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190702_OH15500\3\3_mNeptune.tif';
direc = input_1{1,1};

z_plane_list = dir(direc);
z_plane_names = {z_plane_list(:).name};
rem_index = [find(strcmp(z_plane_names(:),'.')),find(strcmp(z_plane_names(:),'..'))];
z_plane_names(rem_index) = [];
num_z_planes = numel(z_plane_names);

img_size = size(imread([direc,'\',z_plane_names{1}]));

full_img = zeros(img_size(1),img_size(2),num_z_planes,1,'uint16');

for j = 1:num_z_planes
%     full_img(:,:,num_z_planes - j + 1,1) = imread([direc,'\',z_plane_names{j}]);
    full_img(:,:,j,1) = imread([direc,'\',z_plane_names{j}]);
end

%%%%% DONT forget to turn on these options%%%%%%
% addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\GMM_segmetation'))
% indexed_img_to_tiff(full_img(:,:,:,1),[],['C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190702_OH15500\3\thisimage_r.tif'])