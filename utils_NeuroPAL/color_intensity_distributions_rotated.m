%%%% function to get RGB intensity distributions of annotated in NeuroPAL
%%%% strain. These distributions are used to assign node-potentials in
%%%% annotation method.
%%%% histogram matching done before leave-one-out atlas building

% in_direc = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\5',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\21',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\2',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\5',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\8',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\10',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\22',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\24',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\27'};

src_img = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\13'};

in_direc = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\5',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\21',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\2',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\5',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\8',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\10',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\22',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\24',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\27',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190705_OH15495_array\27',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\16',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\20',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\21',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\1',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\4',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\8',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\9',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\10',...
            'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\11'};

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
% src_thisimage_g = imadjustn(mat2gray(thisimage_g),[0;0.8],[0;1],0.8);
% src_thisimage_b = imadjustn(mat2gray(thisimage_b),[0;0.4],[0;1],0.8);
% src_thisimage_r = imadjustn(mat2gray(thisimage_r),[0;0.5],[0;1],0.8);
    
% % col constancy
% norm_mat = [double(thisimage_r(:)),double(thisimage_g(:)),double(thisimage_b(:))];
% for n = 1:10
%     norm_mat = norm_mat./repmat(sum(norm_mat),size(norm_mat,1),1);
%     norm_mat = norm_mat./repmat(sum(norm_mat,2),1,size(norm_mat,2));
% end
% src_thisimage_r_norm = reshape(norm_mat(:,1),size(thisimage_r,1),size(thisimage_r,2),size(thisimage_r,3));
% src_thisimage_g_norm = reshape(norm_mat(:,2),size(thisimage_g,1),size(thisimage_g,2),size(thisimage_g,3));
% src_thisimage_b_norm = reshape(norm_mat(:,3),size(thisimage_b,1),size(thisimage_b,2),size(thisimage_b,3));

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
    
%     thisimage_g = imadjustn(mat2gray(thisimage_g),[0;0.8],[0;1],0.8);
%     thisimage_b = imadjustn(mat2gray(thisimage_b),[0;0.4],[0;1],0.8);
%     thisimage_r = imadjustn(mat2gray(thisimage_r),[0;0.5],[0;1],0.8);
    
%     % col constancy
%     norm_mat = [double(thisimage_r(:)),double(thisimage_g(:)),double(thisimage_b(:))];
%     for n = 1:10
%         norm_mat = norm_mat./repmat(sum(norm_mat),size(norm_mat,1),1);
%         norm_mat = norm_mat./repmat(sum(norm_mat,2),1,size(norm_mat,2));
%     end
%     thisimage_r_norm = reshape(norm_mat(:,1),size(thisimage_r,1),size(thisimage_r,2),size(thisimage_r,3));
%     thisimage_g_norm = reshape(norm_mat(:,2),size(thisimage_g,1),size(thisimage_g,2),size(thisimage_g,3));
%     thisimage_b_norm = reshape(norm_mat(:,3),size(thisimage_b,1),size(thisimage_b,2),size(thisimage_b,3));
% 
%     %%% hist match RGB images
%     thisimage_r_histmatch = imhistmatchn(mat2gray(thisimage_r),mat2gray(src_thisimage_r));
%     thisimage_g_histmatch = imhistmatchn(mat2gray(thisimage_g),mat2gray(src_thisimage_g));
%     thisimage_b_histmatch = imhistmatchn(mat2gray(thisimage_b),mat2gray(src_thisimage_b));
    
    
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
        R_int = thisimage_r(pixels);
        G_int = thisimage_g(pixels);
        B_int = thisimage_b(pixels);
        
        index_annotated_neurons = find(strcmp(marker_name{an,1},uniq_annotated_neurons));
        all_col_int(index_annotated_neurons).R_int = cat(1,all_col_int(index_annotated_neurons).R_int,R_int);
        all_col_int(index_annotated_neurons).G_int = cat(1,all_col_int(index_annotated_neurons).G_int,G_int);
        all_col_int(index_annotated_neurons).B_int = cat(1,all_col_int(index_annotated_neurons).B_int,B_int);
    end
end
save('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\all_col_int_rot_ex11','all_col_int')        