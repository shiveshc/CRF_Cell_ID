%%%% function to get RGB intensity distributions of annotated in NeuroPAL
%%%% strain. These distributions are used to assign node-potentials in
%%%% annotation method

% in_direc = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\5',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\21',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\2',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\5',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\8',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\10',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\22',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\24',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\27'};

in_direc = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\5',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\21',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\2',...
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
    g_norm = imadjustn(mat2gray(thisimage_g),[0;0.8],[0;1],0.8);
    b_norm = imadjustn(mat2gray(thisimage_b),[0;0.4],[0;1],0.8);
    r_norm = imadjustn(mat2gray(thisimage_r),[0;0.5],[0;1],0.8);
%     rgb_img = cat(3,r_norm,g_norm,b_norm);
%     % rgb_img(30:31,330:345,:) = 1;
%     figure,imshow(rgb_img)
    
    
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
        R_int = r_norm(pixels);
        G_int = g_norm(pixels);
        B_int = b_norm(pixels);
        
        index_annotated_neurons = find(strcmp(marker_name{an,1},uniq_annotated_neurons));
        all_col_int(index_annotated_neurons).R_int = cat(1,all_col_int(index_annotated_neurons).R_int,R_int);
        all_col_int(index_annotated_neurons).G_int = cat(1,all_col_int(index_annotated_neurons).G_int,G_int);
        all_col_int(index_annotated_neurons).B_int = cat(1,all_col_int(index_annotated_neurons).B_int,B_int);
    end
end
save('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos_Col\all_col_int_ex9_imadjust','all_col_int')        