%%%% function to create rgb image of neuroPAL data based on annotation
%%%% results

function plot_labels_on_rgb_img()  
    results_dir = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\Results_RGB\non_pca_based_axes\Results_multi_3D_data_20190705_27_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\Results_RGB\non_pca_based_axes\Results_multi_3D_data_20190706_16_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\Results_RGB\non_pca_based_axes\Results_multi_3D_data_20190706_20_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\Results_RGB\non_pca_based_axes\Results_multi_3D_data_20190710_21_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\Results_RGB\non_pca_based_axes\Results_multi_3D_data_20190720_1_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\Results_RGB\non_pca_based_axes\Results_multi_3D_data_20190720_4_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\Results_RGB\non_pca_based_axes\Results_multi_3D_data_20190720_8_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\Results_RGB\non_pca_based_axes\Results_multi_3D_data_20190720_9_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\Results_RGB\non_pca_based_axes\Results_multi_3D_data_20190720_10_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\Results_RGB\non_pca_based_axes\Results_multi_3D_data_20190720_11_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\Results_RGB\non_pca_based_axes\Results_multi_3D_data_20190720_13_numLand0\'};
    
    data = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190705_27.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190706_16.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190706_20.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190710_21.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190720_1.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190720_4.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190720_8.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190720_9.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190720_10.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190720_11.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190720_13.mat'};

    input_dir = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190705_OH15495_array\27',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\16',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\20',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\21',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\1',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\4',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\8',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\9',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\10',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\11',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\13'};
    
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
        g_norm = imadjust(max(mat2gray(thisimage_g),[],3),[0;0.9],[0;1],0.8);
        b_norm = imadjust(max(mat2gray(thisimage_b),[],3),[0;0.3],[0;1],0.8);
        r_norm = imadjust(max(mat2gray(thisimage_r),[],3),[0;0.4],[0;1],0.8);
        rgb_img = cat(3,r_norm,g_norm,b_norm);
        % rgb_img(30:31,330:345,:) = 1;
        figure,imshow(rgb_img),hold on
    
    
        fileList = dir(results_dir{1,f});        
        labels = {};
        cnt = 1;
        for i = 1:size(fileList,1)
            name = fileList(i).name;
            if ~strcmp(name,'.') && ~strcmp(name,'..')
                load([results_dir{1,f},name])

                for k = 1:size(experiments(1).node_label,1)
                    if experiments(1).node_label(k,1) == 0
                        labels{k,cnt} = '';
                    else
                        labels{k,cnt} = experiments(1).Neuron_head{experiments(1).node_label(k,1),1};
                    end
                end
                cnt = cnt + 1;
            end
        end
        top_labels = calculate_top_labels(labels);
        
        load(data{1,f})
        for n = 1:size(top_labels,2)
            text(mu_r(n,1),mu_r(n,2),top_labels(n).top_labels{1,1},'FontSize',6, 'Color', [1,1,1])
        end
        
        fig = gcf;
        fig.InvertHardcopy = 'off';
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf,[input_dir{1,f},'\rgb_img_RelPos_Col_non_pca.png'])
    end
end

function top_labels = calculate_top_labels(labels)
    top_labels = struct();
    for i = 1:size(labels,1)
        uniq_labels = unique(labels(i,:));
        count = [];
        for j = 1:size(uniq_labels,2)
            count(1,j) = sum(strcmp(uniq_labels(1,j),labels(i,:)));
        end
        [sort_count,sort_index] = sort(count,'descend');
        top_labels(i).top_labels = uniq_labels(1,sort_index(1,1:min(5,size(sort_index,2))));
    end
end