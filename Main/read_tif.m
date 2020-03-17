%%%% Function to read an image stack
%%%% Inputs - directory where stack files are saved. each channel must be
%%%% saved separately


function img = read_tif(style,input)
if nargin < 2
    ["Please specify both input and style arguments"]
else
    if strcmp(style,'one')
        info = imfinfo(input);
        img = zeros(info(1).Height,info(1).Width,size(info,1),'uint16');
        for z = 1:size(img,3)
            img(:,:,z) = imread(input,z);
        end
    elseif strcmp(style,'separate')
        z_plane_list = dir(input);
        z_plane_names = {z_plane_list(:).name};
        rem_index = [find(strcmp(z_plane_names(:),'.')),find(strcmp(z_plane_names(:),'..'))];
        z_plane_names(rem_index) = [];
        num_z_planes = numel(z_plane_names);

        img_size = size(imread([input,'\',z_plane_names{1}]));

        img = zeros(img_size(1),img_size(2),num_z_planes,1,'uint16');

        for z = 1:num_z_planes
    %         img(:,:,num_z_planes - j + 1,1) = imread([direc,'\',z_plane_names{j}]);
            img(:,:,z,1) = imread([input,'\',z_plane_names{z}]);
        end
    else
        ["style option should be 'one' or 'separate' only"]
    end
end

%%%%% DONT forget to turn on these options%%%%%%
% addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\GMM_segmetation'))
% indexed_img_to_tiff(full_img(:,:,:,1),[],['C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190702_OH15500\3\thisimage_r.tif'])