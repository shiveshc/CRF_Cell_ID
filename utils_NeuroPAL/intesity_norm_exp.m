%%%%% check intensity values of '20190720_OH15496_array' imaging
%%%%% intensity values should be normalized so that color values match
%%%%% across imaging sessions

input_dir = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\16'};
for f = 1:size(input_dir,2)
    curr_files = dir([input_dir{1,f}]);
    for n = 1:size(curr_files,1)
        curr_name = curr_files(n).name;
        if ~(strcmp('.',curr_name)) && ~(strcmp('..',curr_name)) && isfolder([input_dir{1,f},'\',curr_name]) && isempty(regexp(curr_name,'index','match'))
            z_plane_list = dir([input_dir{1,f},'\',curr_name]);
            z_plane_names = {z_plane_list(:).name};
            rem_index = [find(strcmp(z_plane_names(:),'.')),find(strcmp(z_plane_names(:),'..'))];
            z_plane_names(rem_index) = [];
            num_z_planes = numel(z_plane_names);

            img_size = size(imread([input_dir{1,f},'\',curr_name,'\',z_plane_names{1}]));

            full_img = zeros(img_size(1),img_size(2),num_z_planes,1,'uint16');

            for j = 1:num_z_planes
                full_img(:,:,num_z_planes - j + 1,1) = imread([input_dir{1,f},'\',curr_name,'\',z_plane_names{j}]);
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

    [X,Y,Z] = read_marker_files_wo_marker_name(input_dir{1,f});
    R_int = [];
    G_int = [];
    B_int = [];
    for an = 1:size(X,1)
        X_range = max(1,round(X(an,1)) - 2):1:min(size(thisimage_r,1),round(X(an,1))+2);
        Y_range = max(1,round(Y(an,1)) - 2):1:min(size(thisimage_r,2),round(Y(an,1))+2);
        Z_range = max(1,round(Z(an,1)) - 1):1:min(size(thisimage_r,3),round(Z(an,1))+1);
        [x,y,z] = meshgrid(Y_range,X_range,Z_range);
        pixels = sub2ind(size(thisimage_r),x(:),y(:),z(:));
        R_int = [R_int;thisimage_r(pixels)];
        G_int = [G_int;thisimage_g(pixels)];
        B_int = [B_int;thisimage_b(pixels)];
    end
end

%%% compare intenisty distributions of annotated and new data
load('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\all_col_int.mat')
all_R_int = [];
all_G_int = [];
all_B_int = [];
for i = 1:size(all_col_int,2)
    all_R_int = [all_R_int;double(all_col_int(i).R_int)];
    all_G_int = [all_G_int;double(all_col_int(i).G_int)];
    all_B_int = [all_B_int;double(all_col_int(i).B_int)];
end

[~,centers] = hist(all_R_int,50);
h1 = hist(all_R_int,50);
h1_p = hist(double(R_int),centers);
figure,bar(h1)
hold on
bar(h1_p)
alpha(0.7)

[~,centers] = hist(all_G_int,50);
h2 = hist(all_G_int,50);
h2_p = hist(double(G_int),centers);
figure,bar(h2)
hold on
bar(h2_p)
alpha(0.7)

[~,centers] = hist(all_B_int,50);
h3 = hist(all_B_int,50);
h3_p = hist(double(B_int),centers);
figure,bar(h3)
hold on
bar(h3_p)
alpha(0.7)