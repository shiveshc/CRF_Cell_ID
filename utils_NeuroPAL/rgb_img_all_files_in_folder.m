%%%%% make rgb image of NeuroPAL strain
%%%%% This code uses 'readCyOFP_SD.m' and 'rgb_img_from_indev_channels.m'

%%%%% NOTE - the folders must be empty except individual image folders. And
%%%%% they must contain exactly 4 folders (RFP, BFP, CyOFP, and mNeptune
%%%%% channels)

addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\GMM_segmetation'))

input_dir = 'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\';
fileList = dir(input_dir);
short_list = [13]; %%%% re-run for the shortlisted images. These folders must be empty

if ~isempty(short_list)
    for i = 1:size(short_list,2)
        name = num2str(short_list(1,i));
        curr_files = dir([input_dir,'\',name]);
        for n = 1:size(curr_files,1)
            curr_name = curr_files(n).name;
            if ~(strcmp('.',curr_name)) && ~(strcmp('..',curr_name)) && isfolder([input_dir,'\',name,'\',curr_name]) && isempty(regexp(curr_name,'index','match'))
                z_plane_list = dir([input_dir,'\',name,'\',curr_name]);
                z_plane_names = {z_plane_list(:).name};
                rem_index = [find(strcmp(z_plane_names(:),'.')),find(strcmp(z_plane_names(:),'..'))];
                z_plane_names(rem_index) = [];
                num_z_planes = numel(z_plane_names);

                img_size = size(imread([input_dir,'\',name,'\',curr_name,'\',z_plane_names{1}]));

                full_img = zeros(img_size(1),img_size(2),num_z_planes,1,'uint16');

                for j = 1:num_z_planes
                    full_img(:,:,num_z_planes - j + 1,1) = imread([input_dir,'\',name,'\',curr_name,'\',z_plane_names{j}]);
%                     full_img(:,:,j,1) = imread([input_dir,'\',name,'\',curr_name,'\',z_plane_names{j}]);
                end
            
                if ~isempty(regexp(curr_name,'[rR]','match'))
                    indexed_img_to_tiff(full_img(:,:,:,1),[],[input_dir,'\',name,'\thisimage_rfp.tif'])
                elseif ~isempty(regexp(curr_name,'[bB]','match')) 
                    indexed_img_to_tiff(full_img(:,:,:,1),[],[input_dir,'\',name,'\thisimage_b.tif'])
                    thisimage_b = full_img;
                elseif ~isempty(regexp(curr_name,'[cC]','match')) 
                    indexed_img_to_tiff(full_img(:,:,:,1),[],[input_dir,'\',name,'\thisimage_c.tif'])
                    thisimage_g = full_img;
                elseif ~isempty(regexp(curr_name,'[mM]','match'))
                    indexed_img_to_tiff(full_img(:,:,:,1),[],[input_dir,'\',name,'\thisimage_r.tif'])
                    thisimage_r = full_img;
                end
            end
        end
        g_norm = imadjust(max(mat2gray(thisimage_g),[],3),[0;0.7],[0;1],0.8);
        b_norm = imadjust(max(mat2gray(thisimage_b),[],3),[0;0.3],[0;1],0.8);
        r_norm = imadjust(max(mat2gray(thisimage_r),[],3),[0;0.4],[0;1],0.8);
        rgb_img = cat(3,r_norm,g_norm,b_norm);
        % rgb_img(30:31,330:345,:) = 1;
        figure,imshow(rgb_img)
        saveas(gcf,[input_dir,'\',name,'\rgb_img.png'])
    end
else
    for i = 1:size(fileList,1)
        name = fileList(i).name;
        if ~(strcmp('.',name)) && ~(strcmp('..',name))
            curr_files = dir([input_dir,'\',name]);
            for n = 1:size(curr_files,1)
                curr_name = curr_files(n).name;
                if ~(strcmp('.',curr_name)) && ~(strcmp('..',curr_name)) && isfolder([input_dir,'\',name,'\',curr_name]) && isempty(regexp(curr_name,'index','match'))
                    z_plane_list = dir([input_dir,'\',name,'\',curr_name]);
                    z_plane_names = {z_plane_list(:).name};
                    rem_index = [find(strcmp(z_plane_names(:),'.')),find(strcmp(z_plane_names(:),'..'))];
                    z_plane_names(rem_index) = [];
                    num_z_planes = numel(z_plane_names);

                    img_size = size(imread([input_dir,'\',name,'\',curr_name,'\',z_plane_names{1}]));

                    full_img = zeros(img_size(1),img_size(2),num_z_planes,1,'uint16');

                    for j = 1:num_z_planes
                        full_img(:,:,num_z_planes - j + 1,1) = imread([input_dir,'\',name,'\',curr_name,'\',z_plane_names{j}]);
                    end

                    if ~isempty(regexp(curr_name,'[rR]','match'))
                        indexed_img_to_tiff(full_img(:,:,:,1),[],[input_dir,'\',name,'\thisimage_rfp.tif'])
                        thisimage_rfp = full_img;
                    elseif ~isempty(regexp(curr_name,'[bB]','match')) 
                        indexed_img_to_tiff(full_img(:,:,:,1),[],[input_dir,'\',name,'\thisimage_b.tif'])
                        thisimage_b = full_img;
                    elseif ~isempty(regexp(curr_name,'[cC]','match')) 
                        indexed_img_to_tiff(full_img(:,:,:,1),[],[input_dir,'\',name,'\thisimage_c.tif'])
                        thisimage_g = full_img;
                    elseif ~isempty(regexp(curr_name,'[mM]','match'))
                        indexed_img_to_tiff(full_img(:,:,:,1),[],[input_dir,'\',name,'\thisimage_r.tif'])
                        thisimage_r = full_img;
                    end
                end
            end
            g_norm = imadjust(max(mat2gray(thisimage_g),[],3),[0;0.8],[0;1],0.8);
            b_norm = imadjust(max(mat2gray(thisimage_b),[],3),[0;0.4],[0;1],0.8);
            r_norm = imadjust(max(mat2gray(thisimage_r),[],3),[0;0.5],[0;1],0.8);
            rgb_img = cat(3,r_norm,g_norm,b_norm);
            % rgb_img(30:31,330:345,:) = 1;
            figure,imshow(rgb_img)
            saveas(gcf,[input_dir,'\',name,'\rgb_img.png'])
        end
    end
end
                
            