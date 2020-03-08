%%%% function to analyze if variability in neuron positions affect
%%%% prediction accuracy

function positionVariability_vs_predictionAccuracy()
% load refernce atlas data
load('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\data_neuron_relationship.mat')
X_rot_norm = (X_rot - min(X_rot))/(max(X_rot) - min(X_rot));
Y_rot_norm = (Y_rot - min(Y_rot))/(max(Y_rot) - min(Y_rot));
Z_rot_norm = (Z_rot - min(Z_rot))/(max(Z_rot) - min(Z_rot));

% create neuropal annotated compiled data
compiled_data = compile_annotated_data();

all_annotated_neurons = {};
for i = 1:size(compiled_data,2)
    compiled_data(i).X_rot_norm = (compiled_data(i).X_rot - min(compiled_data(i).X_rot))/(max(compiled_data(i).X_rot) - min(compiled_data(i).X_rot));
    compiled_data(i).Y_rot_norm = (compiled_data(i).Y_rot - min(compiled_data(i).Y_rot))/(max(compiled_data(i).Y_rot) - min(compiled_data(i).Y_rot));
    compiled_data(i).Z_rot_norm = (compiled_data(i).Z_rot - min(compiled_data(i).Z_rot))/(max(compiled_data(i).Z_rot) - min(compiled_data(i).Z_rot));
    all_annotated_neurons = cat(1,all_annotated_neurons,compiled_data(i).marker_name);
end
uniq_annotated_neurons = unique(all_annotated_neurons);

annotated_neuron_pos = struct();
for i = 1:size(uniq_annotated_neurons,1)
    annotated_neuron_pos(i).name = uniq_annotated_neurons{i,1};
    annotated_neuron_pos(i).X = [];
    annotated_neuron_pos(i).Y = [];
    annotated_neuron_pos(i).Z = [];
    for n = 1:size(compiled_data,2)
        index_curr_neuron = find(strcmp(uniq_annotated_neurons{i,1},compiled_data(n).marker_name));
        if ~isempty(index_curr_neuron)
            annotated_neuron_pos(i).X = cat(1,annotated_neuron_pos(i).X,compiled_data(n).X_rot_norm(compiled_data(n).marker_index(index_curr_neuron,1),:));
            annotated_neuron_pos(i).Y = cat(1,annotated_neuron_pos(i).Y,compiled_data(n).Y_rot_norm(compiled_data(n).marker_index(index_curr_neuron,1),:));
            annotated_neuron_pos(i).Z = cat(1,annotated_neuron_pos(i).Z,compiled_data(n).Z_rot_norm(compiled_data(n).marker_index(index_curr_neuron,1),:));
        end
    end
    
    index_curr_neuron_in_atlas = find(strcmp(uniq_annotated_neurons{i,1},Neuron_head));
    annotated_neuron_pos(i).X_atlas = X_rot_norm(index_curr_neuron_in_atlas,1);
    annotated_neuron_pos(i).Y_atlas = Y_rot_norm(index_curr_neuron_in_atlas,1);
    annotated_neuron_pos(i).Z_atlas = Z_rot_norm(index_curr_neuron_in_atlas,1);
    
    annotated_neuron_pos(i).X_var = sum((annotated_neuron_pos(i).X - repmat(annotated_neuron_pos(i).X_atlas,size(annotated_neuron_pos(i).X,1),1)).^2);
    annotated_neuron_pos(i).Y_var = sum((annotated_neuron_pos(i).Y - repmat(annotated_neuron_pos(i).Y_atlas,size(annotated_neuron_pos(i).Y,1),1)).^2);
    annotated_neuron_pos(i).Z_var = sum((annotated_neuron_pos(i).Z - repmat(annotated_neuron_pos(i).Z_atlas,size(annotated_neuron_pos(i).Z,1),1)).^2);
    
    pos = [annotated_neuron_pos(i).X,annotated_neuron_pos(i).Y,annotated_neuron_pos(i).Z];
    pos_atlas = [annotated_neuron_pos(i).X_atlas,annotated_neuron_pos(i).Y_atlas,annotated_neuron_pos(i).Z_atlas];
    annotated_neuron_pos(i).var = sum(repmat(diag(pos*pos'),1,size(pos_atlas,1)) + repmat(diag(pos_atlas*pos_atlas')',size(pos,1),1) - 2*pos*pos_atlas');
    
    annotated_neuron_pos(i).mean_X_var = (mean(annotated_neuron_pos(i).X) - annotated_neuron_pos(i).X_atlas).^2;
    annotated_neuron_pos(i).mean_Y_var = (mean(annotated_neuron_pos(i).Y) - annotated_neuron_pos(i).Y_atlas).^2;
    annotated_neuron_pos(i).mean_Z_var = (mean(annotated_neuron_pos(i).Z) - annotated_neuron_pos(i).Z_atlas).^2;
    annotated_neuron_pos(i).mean_var = annotated_neuron_pos(i).mean_X_var + annotated_neuron_pos(i).mean_Y_var + annotated_neuron_pos(i).mean_Z_var;
    annotated_neuron_pos(i).mean_var_sqrt = sqrt(annotated_neuron_pos(i).mean_X_var + annotated_neuron_pos(i).mean_Y_var + annotated_neuron_pos(i).mean_Z_var);
end

%%% for figure - deviation from atlas position for each neuron
deviation = [];
deviation_ap = [];
deviation_lr = [];
deviation_dv = [];
for i = 1:size(annotated_neuron_pos,2)
    for n = 1:size(annotated_neuron_pos(i).X,1)
        deviation = [deviation;sqrt(sum(([annotated_neuron_pos(i).X(n,1),annotated_neuron_pos(i).Y(n,1),annotated_neuron_pos(i).Z(n,1)] - [annotated_neuron_pos(i).X_atlas,annotated_neuron_pos(i).Y_atlas,annotated_neuron_pos(i).Z_atlas]).^2)),i];
        deviation_ap = [deviation_ap;abs(annotated_neuron_pos(i).X(n,1) - annotated_neuron_pos(i).X_atlas),i];
        deviation_lr = [deviation_lr;abs(annotated_neuron_pos(i).Y(n,1) - annotated_neuron_pos(i).Y_atlas),i];
        deviation_dv = [deviation_dv;abs(annotated_neuron_pos(i).Z(n,1) - annotated_neuron_pos(i).Z_atlas),i];
    end
end
%%% for figure - inter-nuclear distance between all neurons
top_5_dist = [];
top_5_dist_ap = [];
top_5_dist_lr = [];
top_5_dist_dv = [];
for i = 1:size(compiled_data,2)
    pos = [compiled_data(i).X_rot_norm,compiled_data(i).Y_rot_norm,compiled_data(i).Z_rot_norm];
    inter_nuc_dist = sqrt((repmat(diag(pos*pos'),1,size(pos,1)) + repmat(diag(pos*pos')',size(pos,1),1) - 2*pos*pos'));
    [sort_dist,sort_index] = sort(inter_nuc_dist,2,'ascend');
    top_5_dist = [top_5_dist;cat(2,sort_dist(:,2:6),repmat(i,size(pos,1),1),[1:1:size(pos,1)]')];
    
    pos_ap = [compiled_data(i).X_rot_norm];
    inter_nuc_dist = sqrt((repmat(diag(pos_ap*pos_ap'),1,size(pos_ap,1)) + repmat(diag(pos_ap*pos_ap')',size(pos_ap,1),1) - 2*pos_ap*pos_ap'));
    [sort_dist,sort_index] = sort(inter_nuc_dist,2,'ascend');
    top_5_dist_ap = [top_5_dist_ap;cat(2,sort_dist(:,2:11),repmat(i,size(pos_ap,1),1),[1:1:size(pos_ap,1)]')];
    pos_lr = [compiled_data(i).Y_rot_norm];
    inter_nuc_dist = sqrt((repmat(diag(pos_lr*pos_lr'),1,size(pos_lr,1)) + repmat(diag(pos_lr*pos_lr')',size(pos_lr,1),1) - 2*pos_lr*pos_lr'));
    [sort_dist,sort_index] = sort(inter_nuc_dist,2,'ascend');
    top_5_dist_lr = [top_5_dist_lr;cat(2,sort_dist(:,2:11),repmat(i,size(pos_lr,1),1),[1:1:size(pos_lr,1)]')];
    pos_dv = [compiled_data(i).Z_rot_norm];
    inter_nuc_dist = sqrt((repmat(diag(pos_dv*pos_dv'),1,size(pos_dv,1)) + repmat(diag(pos_dv*pos_dv')',size(pos_dv,1),1) - 2*pos_dv*pos_dv'));
    [sort_dist,sort_index] = sort(inter_nuc_dist,2,'ascend');
    top_5_dist_dv = [top_5_dist_dv;cat(2,sort_dist(:,2:11),repmat(i,size(pos_dv,1),1),[1:1:size(pos_dv,1)]')];
    
    
end

rel_pos_accuracy = get_rel_pos_accuracy(uniq_annotated_neurons);

end

function rel_pos_accuracy = get_rel_pos_accuracy(uniq_annotated_neurons)
%     addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\'))
    in_direc = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_multi_3D_data_20190706_5_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_multi_3D_data_20190706_21_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_multi_3D_data_20190710_2_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_multi_3D_data_20190710_5_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_multi_3D_data_20190710_8_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_multi_3D_data_20190710_10_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_multi_3D_data_20190710_22_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_multi_3D_data_20190710_24_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_multi_3D_data_20190710_27_numLand0\'};
    
    data = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos_Col\data_20190706_5.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos_Col\data_20190706_21.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos_Col\data_20190710_2.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos_Col\data_20190710_5.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos_Col\data_20190710_8.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos_Col\data_20190710_10.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos_Col\data_20190710_22.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos_Col\data_20190710_24.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos_Col\data_20190710_27.mat'};
    
    rel_pos_accuracy = struct();
    for i = 1:size(uniq_annotated_neurons,1)
        rel_pos_accuracy(i).name = uniq_annotated_neurons{i,1};
        rel_pos_accuracy(i).match_matrix = [];
    end
    for f = 1:size(in_direc,2)
        fileList = dir(in_direc{1,f});
        load(data{1,f})
        landmark_list = marker_index;
%         landmark_list = [1:1:size(marker_index,1)]';
        
        labels = {};
        cnt = 1;
        for i = 1:size(fileList,1)
            name = fileList(i).name;
            if ~strcmp(name,'.') && ~strcmp(name,'..')
                load([in_direc{1,f},name])

                landmark_to_be_pred = landmark_list;
                for k = 1:size(landmark_to_be_pred,1)
                    if experiments(1).node_label(landmark_to_be_pred(k,1),1) == 0
                        labels{k,cnt} = '';
                    else
                        labels{k,cnt} = experiments(1).Neuron_head{experiments(1).node_label(landmark_to_be_pred(k,1),1),1};
                    end
                end
                cnt = cnt + 1;
            end
        end
        top_labels = calculate_top_labels(labels);
        true_labels = marker_name;
        for i = 1:size(marker_name,1)
            index = find(strcmp(marker_name{i,1},uniq_annotated_neurons));
            if ~isempty(index)
                if strcmp(top_labels(i).top_labels{1,1},marker_name{i,1})
                    rel_pos_accuracy(index).match_matrix = cat(1,rel_pos_accuracy(index).match_matrix,1);
                else
                    rel_pos_accuracy(index).match_matrix = cat(1,rel_pos_accuracy(index).match_matrix,0);
                end
            else
                rel_pos_accuracy(index).match_matrix = cat(1,rel_pos_accuracy(index).match_matrix,NaN);
            end
        end
    end
    for i = 1:size(rel_pos_accuracy,2)
        rel_pos_accuracy(i).accuracy = nansum(rel_pos_accuracy(i).match_matrix(:,1))/size(rel_pos_accuracy(i).match_matrix(:,1),1);
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
