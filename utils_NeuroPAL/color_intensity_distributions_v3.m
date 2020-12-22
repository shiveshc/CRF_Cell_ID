%%%% function to get RGB intensity distributions of annotated in NeuroPAL
%%%% strain. These distributions are used to assign node-potentials in
%%%% annotation method.
%%%% color constancy normalization before leave-one-out atlas building

% in_direc = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\5',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\21',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\2',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\5',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\8',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\10',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\22',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\24',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\27'};

src_img = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\27'};

in_direc = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\2',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\5',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\8',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\10',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\22',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\24'};

%%% make list of all annotated neurons
all_annotated_neurons = {};
for i = 1:size(in_direc,2)
    [X,Y,Z,marker_name,marker_index] = read_marker_files(in_direc{1,i});
    all_annotated_neurons = cat(1,all_annotated_neurons,marker_name);
end
uniq_annotated_neurons = unique(all_annotated_neurons);
    
%%% read  src image
curr_files = dir(src_img{1,1});
for n = 1:size(curr_files,1)
    curr_name = curr_files(n).name;
    if ~(strcmp('.',curr_name)) && ~(strcmp('..',curr_name)) && isfolder([src_img{1,1},'\',curr_name]) && isempty(regexp(curr_name,'index','match'))
        z_plane_list = dir([src_img{1,1},'\',curr_name]);
        z_plane_names = {z_plane_list(:).name};
        rem_index = [find(strcmp(z_plane_names(:),'.')),find(strcmp(z_plane_names(:),'..'))];
        z_plane_names(rem_index) = [];
        num_z_planes = numel(z_plane_names);

        img_size = size(imread([src_img{1,1},'\',curr_name,'\',z_plane_names{1}]));

        full_img = zeros(img_size(1),img_size(2),num_z_planes,1,'uint16');

        for j = 1:num_z_planes
%                 full_img(:,:,num_z_planes - j + 1,1) = imread([in_direc{1,i},'\',curr_name,'\',z_plane_names{j}]);
            full_img(:,:,j,1) = imread([src_img{1,1},'\',curr_name,'\',z_plane_names{j}]);
        end

        if ~isempty(regexp(curr_name,'[rR]','match'))
        elseif ~isempty(regexp(curr_name,'[bB]','match')) 
            thisimage_b = full_img;
        elseif ~isempty(regexp(curr_name,'[cC]','match')) 
            thisimage_g = full_img;
        elseif ~isempty(regexp(curr_name,'[mM]','match'))
            thisimage_r = full_img;
        end
    end
end
src_thisimage_r = double(thisimage_r);
src_thisimage_g = double(thisimage_g);
src_thisimage_b = double(thisimage_b);
norm_mat = [src_thisimage_r(:),src_thisimage_g(:),src_thisimage_b(:)];
for n = 1:10
    norm_mat = norm_mat./repmat(sum(norm_mat),size(norm_mat,1),1);
    norm_mat = norm_mat./repmat(sum(norm_mat,2),1,size(norm_mat,2));
end
src_thisimage_r_norm = reshape(norm_mat(:,1),size(thisimage_r,1),size(thisimage_r,2),size(thisimage_r,3));
src_thisimage_g_norm = reshape(norm_mat(:,2),size(thisimage_g,1),size(thisimage_g,2),size(thisimage_g,3));
src_thisimage_b_norm = reshape(norm_mat(:,3),size(thisimage_b,1),size(thisimage_b,2),size(thisimage_b,3));
% r_norm = max(mat2gray(thisimage_r),[],3);
% g_norm = max(mat2gray(thisimage_g),[],3);
% b_norm = max(mat2gray(thisimage_b),[],3);
% rgb_img = cat(3,r_norm,g_norm,b_norm);
% figure,imshow(rgb_img)
% r_norm = max(mat2gray(src_thisimage_r_norm),[],3);
% g_norm = max(mat2gray(src_thisimage_g_norm),[],3);
% b_norm = max(mat2gray(src_thisimage_b_norm),[],3);
% rgb_img = cat(3,r_norm,g_norm,b_norm);
% figure,imshow(rgb_img)

%%%% create RGB distributions for annotated neurons 
all_col_int = struct();
for i = 1:size(uniq_annotated_neurons,1)
    all_col_int(i).name = uniq_annotated_neurons{i,1};
    all_col_int(i).R_int = [];
    all_col_int(i).G_int = [];
    all_col_int(i).B_int = [];
end
for i = 1:size(in_direc,2)
    %%%% read images
    curr_files = dir([in_direc{1,i}]);
    for n = 1:size(curr_files,1)
        curr_name = curr_files(n).name;
        if ~(strcmp('.',curr_name)) && ~(strcmp('..',curr_name)) && isfolder([in_direc{1,i},'\',curr_name]) && isempty(regexp(curr_name,'index','match'))
            z_plane_list = dir([in_direc{1,i},'\',curr_name]);
            z_plane_names = {z_plane_list(:).name};
            rem_index = [find(strcmp(z_plane_names(:),'.')),find(strcmp(z_plane_names(:),'..'))];
            z_plane_names(rem_index) = [];
            num_z_planes = numel(z_plane_names);

            img_size = size(imread([in_direc{1,i},'\',curr_name,'\',z_plane_names{1}]));

            full_img = zeros(img_size(1),img_size(2),num_z_planes,1,'uint16');

            for j = 1:num_z_planes
%                 full_img(:,:,num_z_planes - j + 1,1) = imread([in_direc{1,i},'\',curr_name,'\',z_plane_names{j}]);
                full_img(:,:,j,1) = imread([in_direc{1,i},'\',curr_name,'\',z_plane_names{j}]);
            end

            if ~isempty(regexp(curr_name,'[rR]','match'))
            elseif ~isempty(regexp(curr_name,'[bB]','match')) 
                thisimage_b = full_img;
            elseif ~isempty(regexp(curr_name,'[cC]','match')) 
                thisimage_g = full_img;
            elseif ~isempty(regexp(curr_name,'[mM]','match'))
                thisimage_r = full_img;
            end
        end
    end
%     thisimage_r = mat2gray(thisimage_r);
%     thisimage_g = mat2gray(thisimage_g);
%     thisimage_b = mat2gray(thisimage_b);
%     g_norm = imadjust(max(mat2gray(thisimage_g),[],3),[0;0.8],[0;1],0.8);
%     b_norm = imadjust(max(mat2gray(thisimage_b),[],3),[0;0.4],[0;1],0.8);
%     r_norm = imadjust(max(mat2gray(thisimage_r),[],3),[0;0.5],[0;1],0.8);
%     rgb_img = cat(3,r_norm,g_norm,b_norm);
%     % rgb_img(30:31,330:345,:) = 1;
%     figure,imshow(rgb_img)
    
    %%% color cont img
    thisimage_r = double(thisimage_r);
    thisimage_g = double(thisimage_g);
    thisimage_b = double(thisimage_b);
    norm_mat = [thisimage_r(:),thisimage_g(:),thisimage_b(:)];
    for n = 1:10
        norm_mat = norm_mat./repmat(sum(norm_mat),size(norm_mat,1),1);
        norm_mat = norm_mat./repmat(sum(norm_mat,2),1,size(norm_mat,2));
    end
    thisimage_r_norm = reshape(norm_mat(:,1),size(thisimage_r,1),size(thisimage_r,2),size(thisimage_r,3));
    thisimage_g_norm = reshape(norm_mat(:,2),size(thisimage_g,1),size(thisimage_g,2),size(thisimage_g,3));
    thisimage_b_norm = reshape(norm_mat(:,3),size(thisimage_b,1),size(thisimage_b,2),size(thisimage_b,3));
%     r_norm = max(mat2gray(thisimage_r),[],3);
%     g_norm = max(mat2gray(thisimage_g),[],3);
%     b_norm = max(mat2gray(thisimage_b),[],3);
%     rgb_img = cat(3,r_norm,g_norm,b_norm);
%     figure,imshow(rgb_img)
%     r_norm = max(mat2gray(thisimage_r_norm),[],3);
%     g_norm = max(mat2gray(thisimage_g_norm),[],3);
%     b_norm = max(mat2gray(thisimage_b_norm),[],3);
%     rgb_img = cat(3,r_norm,g_norm,b_norm);
%     figure,imshow(rgb_img)
    
    %%% hist match RGB images
    thisimage_r_histmatch = imhistmatchn(mat2gray(thisimage_r_norm),mat2gray(src_thisimage_r_norm));
    thisimage_g_histmatch = imhistmatchn(mat2gray(thisimage_g_norm),mat2gray(src_thisimage_g_norm));
    thisimage_b_histmatch = imhistmatchn(mat2gray(thisimage_b_norm),mat2gray(src_thisimage_b_norm));
    
    %%% read curr marker file
    [X,Y,Z,marker_name,marker_index] = read_marker_files(in_direc{1,i});
    temp_img = zeros(size(thisimage_r));
    for an = 1:size(marker_index,1)
        X_range = max(1,round(X(marker_index(an,1),1)) - 2):1:min(size(thisimage_r,2),round(X(marker_index(an,1),1))+2);
        Y_range = max(1,round(Y(marker_index(an,1),1)) - 2):1:min(size(thisimage_r,1),round(Y(marker_index(an,1),1))+2);
        Z_range = max(1,round(Z(marker_index(an,1),1)) - 1):1:min(size(thisimage_r,3),round(Z(marker_index(an,1),1))+1);
        [x,y,z] = meshgrid(Y_range,X_range,Z_range);
        pixels = sub2ind(size(thisimage_r),x(:),y(:),z(:));
        temp_img(pixels) = 1;
        R_int = thisimage_r_histmatch(pixels);
        G_int = thisimage_g_histmatch(pixels);
        B_int = thisimage_b_histmatch(pixels);
        
        index_annotated_neurons = find(strcmp(marker_name{an,1},uniq_annotated_neurons));
        all_col_int(index_annotated_neurons).R_int = cat(1,all_col_int(index_annotated_neurons).R_int,R_int);
        all_col_int(index_annotated_neurons).G_int = cat(1,all_col_int(index_annotated_neurons).G_int,G_int);
        all_col_int(index_annotated_neurons).B_int = cat(1,all_col_int(index_annotated_neurons).B_int,B_int);
    end
end
save('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos_Col\all_col_int_ex9_colconst_samesess','all_col_int')        