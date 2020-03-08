%%%% function to create rgb image of annotated neuroPAL data
%%%% Inputs - 1. input_dir (directory name containing B, C, M channel folders
%%%%          2. X, Y read from 'markers' file in the input_dir using
%%%%             'read_marker_files.m'
%%%%          3. 'marker_name' and 'marker_index' variable read using
%%%%             'read_marker_files.m'

function rgb_img_annotated_data(input_dir,X,Y,marker_name,marker_index)  
    
    curr_files = dir([input_dir]);
    for n = 1:size(curr_files,1)
        curr_name = curr_files(n).name;
        if ~(strcmp('.',curr_name)) && ~(strcmp('..',curr_name)) && isfolder([input_dir,'\',curr_name]) && isempty(regexp(curr_name,'index','match'))
            z_plane_list = dir([input_dir,'\',curr_name]);
            z_plane_names = {z_plane_list(:).name};
            rem_index = [find(strcmp(z_plane_names(:),'.')),find(strcmp(z_plane_names(:),'..'))];
            z_plane_names(rem_index) = [];
            num_z_planes = numel(z_plane_names);

            img_size = size(imread([input_dir,'\',curr_name,'\',z_plane_names{1}]));

            full_img = zeros(img_size(1),img_size(2),num_z_planes,1,'uint16');

            for j = 1:num_z_planes
                full_img(:,:,num_z_planes - j + 1,1) = imread([input_dir,'\',curr_name,'\',z_plane_names{j}]);
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
    g_norm = imadjust(max(mat2gray(thisimage_g),[],3),[0;0.9],[0;1],0.8);
    b_norm = imadjust(max(mat2gray(thisimage_b),[],3),[0;0.3],[0;1],0.8);
    r_norm = imadjust(max(mat2gray(thisimage_r),[],3),[0;0.4],[0;1],0.8);
    rgb_img = cat(3,r_norm,g_norm,b_norm);
    % rgb_img(30:31,330:345,:) = 1;
    figure,imshow(rgb_img)
    hold on
    for n = 1:size(X,1)
        annot = find(n == marker_index(:,1));
        if ~isempty(annot)
            text(X(n,1) ,Y(n,1),marker_name{annot,1},'Color',[1,1,1],'FontSize',5)
        else
            text(X(n,1) ,Y(n,1),num2str(n),'Color',[1,1,1],'FontSize',5)
        end
    end
    fig = gcf;
    fig.InvertHardcopy = 'off';
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf,[input_dir,'\rgb_img_annotated.png'])
end