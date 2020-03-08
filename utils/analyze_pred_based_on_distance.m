%%%% function to analyze alternate strain data based on distance of target
%%%% neuron from assigned label.

function analyze_pred_based_on_distance()
    in_direc = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_1r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_3r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_4r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_5r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_6r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_7r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_8r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_9r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_10r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_11r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_14r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_19r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_20r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_22r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_23r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_24r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_28r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_29r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_33r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_34r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_36r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_37r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_38r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_39r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_50r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_52r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_54r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_55r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_57r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_58r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_59r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_60r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_62r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_65r_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_67r_numLand0\'};
    
    data = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_1r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_3r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_4r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_5r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_6r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_7r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_8r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_9r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_10r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_11r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_14r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_19r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_20r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_22r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_23r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_24r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_28r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_29r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_33r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_34r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_36r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_37r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_38r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_39r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_50r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_52r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_54r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_55r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_57r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_58r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_59r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_60r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_62r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_65r.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_67r.mat'};

    nearest_neighbor_dist = [];    
    for i = 1:size(data,2)
        load(data{1,i})
        mu_r_norm = (mu_r - repmat(min(mu_r),size(mu_r,1),1))./repmat(max(mu_r) - min(mu_r),size(mu_r,1),1);
        dist_mat = sqrt(repmat(diag(mu_r_norm*mu_r_norm'),1,size(mu_r,1)) + repmat(diag(mu_r_norm*mu_r_norm')',size(mu_r,1),1) - 2*mu_r_norm*mu_r_norm');
        [sort_dist,sort_index] = sort(dist_mat,2,'ascend');
        nearest_neighbor_dist = [nearest_neighbor_dist;sort_dist(:,2:6)];
    end
    mean_nearest_neighbor_dist = mean(nearest_neighbor_dist);
    median_nearest_neighbor_dist = median(nearest_neighbor_dist);
    
    target_dist = [];
    target_match = [];
    for f = 1:size(in_direc,2)
        fileList = dir(in_direc{1,f});
        load(data{1,f})
        mu_r_norm = (mu_r - repmat(min(mu_r),size(mu_r,1),1))./repmat(max(mu_r) - min(mu_r),size(mu_r,1),1);
        landmark_list = marker_to_neuron_map;
        
        labels = {};
        cnt = 1;
        for i = 1:size(fileList,1)
            name = fileList(i).name;
            if ~strcmp(name,'.') && ~strcmp(name,'..')
                load([in_direc{1,f},name])

                for k = 1:size(mu_r,1)
                    labels{k,cnt} = experiments(1).Neuron_head{experiments(1).node_label(k,1),1};
                end
                cnt = cnt + 1;
            end
        end
        top_labels = calculate_top_labels(labels);
        n = 1;
        top_labels_within_n_cell_dia = calculate_top_labels_within_n_cell_dia(n,top_labels,mu_r_norm,mean_nearest_neighbor_dist);
        top_labels_within_n_nn = calculate_top_labels_within_n_nn(n,top_labels,mu_r_norm,mean_nearest_neighbor_dist);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%% For AX2164, IF only one of URXL,URXR is present in the stack, it was hard
%         %%% to tell whether the neuron was URXL/URXR. So make results
%         %%% indifferent to it
%         ind_URXL = find(strcmp('URXL',marker_name));
%         ind_URXR = find(strcmp('URXR',marker_name));
%         if isempty(ind_URXL) && ~isempty(ind_URXR)
%             marker_name{ind_URXR,1} = 'URX';
%             for i = 1:size(top_labels(ind_URXR).top_labels,2)
%                 if strcmp('URXL',top_labels(ind_URXR).top_labels{1,i})
%                     top_labels(ind_URXR).top_labels{1,i} = 'URX';
%                 elseif strcmp('URXR',top_labels(ind_URXR).top_labels{1,i})
%                     top_labels(ind_URXR).top_labels{1,i} = 'URX';
%                 end
%             end
%         elseif ~isempty(ind_URXL) && isempty(ind_URXR)
%             marker_name{ind_URXL,1} = 'URX';
%             for i = 1:size(top_labels(ind_URXL).top_labels,2)
%                 if strcmp('URXL',top_labels(ind_URXL).top_labels{1,i})
%                     top_labels(ind_URXL).top_labels{1,i} = 'URX';
%                 elseif strcmp('URXR',top_labels(ind_URXL).top_labels{1,i})
%                     top_labels(ind_URXL).top_labels{1,i} = 'URX';
%                 end
%             end
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% for AML 5, since RIGL, AVG and RIGR could not be identified manually don't
        %%% assess prediction for these neurons
        ind_RIGL = find(strcmp('RIGL',marker_name));
        ind_RIGR = find(strcmp('RIGR',marker_name));
        ind_AVG = find(strcmp('AVG',marker_name));
%         top_labels([ind_RIGL;ind_RIGR;ind_AVG]) = [];
        marker_name([ind_RIGL;ind_RIGR;ind_AVG],:) = [];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k = 1:size(marker_name,1)
            match_index = [];
            curr_neuron = marker_name{k,1};
            for n = 1:size(top_labels,2)
                if ~isempty(find(strcmp(curr_neuron,top_labels(n).top_labels(1,1:min(1,size(top_labels(n).top_labels,2))))))
                    match_index = [match_index;n];
                end
            end
            if ~isempty(match_index)
                dist_to_curr_neuron = sqrt(repmat(diag(mu_r_norm(marker_to_neuron_map(k,1),:)*mu_r_norm(marker_to_neuron_map(k,1),:)'),1,size(match_index,1)) + repmat(diag(mu_r_norm(match_index,:)*mu_r_norm(match_index,:)')',1,1) - 2*mu_r_norm(marker_to_neuron_map(k,1),:)*mu_r_norm(match_index,:)');
                [sort_dist,sort_index] = sort(dist_to_curr_neuron,2,'ascend');
                target_dist = [target_dist;sort_dist(1,1)];
                target_match = [target_match;[f,k,1]];
            else
                target_dist = [target_dist;NaN];
                target_match = [target_match;[f,k,0]];
            end
        end
    end
    for i = 1:size(target_dist,1)
        if (target_match(i,3) ~= 0)
            if target_dist(i,1) < mean_nearest_neighbor_dist(1,1)
                target_match(i,4) = 1;
                target_match(i,5) = 0;
                target_match(i,6) = 0;
                target_match(i,7) = 0;
            elseif target_dist(i,1) < mean_nearest_neighbor_dist(1,2)
                target_match(i,4) = 0;
                target_match(i,5) = 1;
                target_match(i,6) = 0;
                target_match(i,7) = 0;
            elseif target_dist(i,1) < mean_nearest_neighbor_dist(1,3)
                target_match(i,4) = 0;
                target_match(i,5) = 0;
                target_match(i,6) = 1;
                target_match(i,7) = 0;
            else
                target_match(i,4) = 0;
                target_match(i,5) = 0;
                target_match(i,6) = 0;
                target_match(i,7) = 1;
            end
        else
            target_match(i,4) = 0;
            target_match(i,5) = 0;
            target_match(i,6) = 0;
            target_match(i,7) = 0;
        end
%         if (target_match(i,3) ~= 0)
%             if target_dist(i,1) < mean_nearest_neighbor_dist(1,1)
%                 target_match(i,4) = 1;
%                 target_match(i,5) = 0;
%                 target_match(i,6) = 0;
%                 target_match(i,7) = 0;
%             elseif target_dist(i,1) < 2*mean_nearest_neighbor_dist(1,1)
%                 target_match(i,4) = 0;
%                 target_match(i,5) = 1;
%                 target_match(i,6) = 0;
%                 target_match(i,7) = 0;
%             elseif target_dist(i,1) < 3*mean_nearest_neighbor_dist(1,1)
%                 target_match(i,4) = 0;
%                 target_match(i,5) = 0;
%                 target_match(i,6) = 1;
%                 target_match(i,7) = 0;
%             else
%                 target_match(i,4) = 0;
%                 target_match(i,5) = 0;
%                 target_match(i,6) = 0;
%                 target_match(i,7) = 1;
%             end
%         else
%             target_match(i,4) = 0;
%             target_match(i,5) = 0;
%             target_match(i,6) = 0;
%             target_match(i,7) = 0;
%         end
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

function top_labels = calculate_top_labels_within_n_cell_dia(k,labels,mu_r_norm,mean_nearest_neighbor_dist)
    top_labels = struct();
    dist_mat = sqrt(repmat(diag(mu_r_norm*mu_r_norm'),1,size(mu_r_norm,1)) + repmat(diag(mu_r_norm*mu_r_norm')',size(mu_r_norm,1),1) - 2*mu_r_norm*mu_r_norm');
    for i = 1:size(mu_r_norm,1)
        within_n_cell_dia = find(dist_mat(i,:)<k*mean_nearest_neighbor_dist(1,1));
        dummy = {};
        for n = 1:size(within_n_cell_dia,2)
            dummy = cat(2,dummy,labels(within_n_cell_dia(1,n)).top_labels(1,1:min(5,size(labels(within_n_cell_dia(1,n)).top_labels,2))));
        end
        uniq_dummy = unique(dummy);
        count = [];
        for j = 1:size(uniq_dummy,2)
            count(1,j) = sum(strcmp(uniq_dummy(1,j),dummy));
        end
        [sort_count,sort_index] = sort(count,2,'descend');
        top_within_n_cell_dia = uniq_dummy(:,sort_index(1:min(5,size(sort_index,2))));
        top_labels(i).top_labels = uniq_dummy;
    end
end

function top_labels = calculate_top_labels_within_n_nn(k,labels,mu_r_norm,mean_nearest_neighbor_dist)
    top_labels = struct();
    dist_mat = sqrt(repmat(diag(mu_r_norm*mu_r_norm'),1,size(mu_r_norm,1)) + repmat(diag(mu_r_norm*mu_r_norm')',size(mu_r_norm,1),1) - 2*mu_r_norm*mu_r_norm');
    for i = 1:size(mu_r_norm,1)
        within_n_nn = find(dist_mat(i,:)<mean_nearest_neighbor_dist(1,k));
        dummy = {};
        for n = 1:size(within_n_nn,2)
            dummy = cat(2,dummy,labels(within_n_nn(1,n)).top_labels(1,1:min(5,size(labels(within_n_nn(1,n)).top_labels,2))));
        end
        uniq_dummy = unique(dummy);
        count = [];
        for j = 1:size(uniq_dummy,2)
            count(1,j) = sum(strcmp(uniq_dummy(1,j),dummy));
        end
        [sort_count,sort_index] = sort(count,2,'descend');
        top_within_n_nn = uniq_dummy(:,sort_index(1:min(5,size(sort_index,2))));
        top_labels(i).top_labels = uniq_dummy;
    end
end
        