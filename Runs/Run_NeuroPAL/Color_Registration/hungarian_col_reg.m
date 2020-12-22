%%%% function for Col + Reg based annotation of Neuropal data

function pred_struct = hungarian_col_reg()
% addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Joint_tracking\TPS-RPM'))
% addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GLMD_Demo'))
addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\CPD2'))

data = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\5\data_20190706_5.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\21\data_20190706_21.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\2\data_20190710_2.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\5\data_20190710_5.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\8\data_20190710_8.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\10\data_20190710_10.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\22\data_20190710_22.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\24\data_20190710_24.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\27\data_20190710_27.mat'};

load('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\Col_Reg\all_col_int')
load('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\data_neuron_relationship')
addpath(genpath('/gpfs/pace1/project/pchbe2/schaudhary9/Annotation/NeuroPAL/RelPos_Reg/CPD2'))

numLabelRemove = 50;
pred_struct = struct();
cnt = 1;
for f = 1:size(data,2)
    load(data{1,f})
    
    % generate axis
    mu_r_centered = mu_r - repmat(mean(mu_r),size(mu_r,1),1);
    [coeff,score,latent] = pca(mu_r_centered);
    PA = coeff(:,axes_param(1,1))';
    PA = PA/norm(PA);
    LR = coeff(:,axes_param(1,2))';
    LR = LR/norm(LR);
    DV = coeff(:,axes_param(1,3))';
    DV = DV/norm(DV);

    OLLL_neuron = marker_index(find(strcmp('OLLL',marker_name)),1);
    OLLR_neuron = marker_index(find(strcmp('OLLR',marker_name)),1);
    AVAL_neuron = marker_index(find(strcmp('AVAL',marker_name)),1);
    AVAR_neuron = marker_index(find(strcmp('AVAR',marker_name)),1);
    if ~isempty(OLLL_neuron)
        if ~isempty(AVAL_neuron)
            if (mu_r_centered(OLLL_neuron,:)-mu_r_centered(AVAL_neuron,:))*PA' < 0
                PA = -PA;
            end
        else
            if (mu_r_centered(OLLL_neuron,:)-mu_r_centered(AVAR_neuron,:))*PA' < 0
                PA = -PA;
            end
        end
    else
        if ~isempty(AVAL_neuron)
            if (mu_r_centered(OLLR_neuron,:)-mu_r_centered(AVAL_neuron,:))*PA' < 0
                PA = -PA;
            end
        else
            if (mu_r_centered(OLLR_neuron,:)-mu_r_centered(AVAR_neuron,:))*PA' < 0
                PA = -PA;
            end
        end
    end

    ALA_neuron = marker_index(find(strcmp('ALA',marker_name)),1);
    RMDDL_neuron = marker_index(find(strcmp('RMDDL',marker_name)),1);
    RMDDR_neuron = marker_index(find(strcmp('RMDDR',marker_name)),1);
    if ~isempty(RMDDL_neuron)
        if (mu_r_centered(RMDDL_neuron,:)-mu_r_centered(ALA_neuron,:))*DV' < 0
            DV = -DV;
        end
    else
        if (mu_r_centered(RMDDR_neuron,:)-mu_r_centered(ALA_neuron,:))*DV' < 0
            DV = -DV;
        end
    end
    if cross(PA,LR)*DV' < 0
        LR = -LR;
    end

    % take neurons coordinates to AP, LR, DV axis
    X = mu_r_centered*PA';
    Y = mu_r_centered*LR';
    Z = mu_r_centered*DV';
    
%     lambda_reg_1 = 1:1:10;
%     lambda_reg_2 = 0.1:0.05:0.5;
%     lambda_col_1 = 1:1:10;
%     lambda_col_2 = 0.1:0.1:0.5;
    
    lambda_reg_1 = 10;
    lambda_reg_2 = 0.1;
    lambda_col_1 = 5;
    lambda_col_2 = 0.1;
    
    top_correct_frac = [];
    top_3_correct_frac = [];
    top_5_correct_frac = [];
    labels = {};
    for k = 1:100
        load('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\data_neuron_relationship')
        %%% randomly remove subset of labels
        anterior_index = find(ganglion(:,1) == 1);
        lateral_index = find(ganglion(:,1) == 2);
        ventral_index = find(ganglion(:,1) == 3);
        num_anterior_remove = round(size(anterior_index,1)/size(ventral_index,1)*numLabelRemove/(size(anterior_index,1)/size(ventral_index,1)+size(lateral_index,1)/size(ventral_index,1)+1));
        num_lateral_remove = round(size(lateral_index,1)/size(ventral_index,1)*numLabelRemove/(size(anterior_index,1)/size(ventral_index,1)+size(lateral_index,1)/size(ventral_index,1)+1));
        num_ventral_remove = numLabelRemove - num_anterior_remove - num_lateral_remove;

        remove_anterior = anterior_index(randperm(size(anterior_index,1),num_anterior_remove),:);
        remove_lateral = lateral_index(randperm(size(lateral_index,1),num_lateral_remove),:);
        remove_ventral = ventral_index(randperm(size(ventral_index,1),num_ventral_remove),:);
        remove_index = [remove_anterior;remove_lateral;remove_ventral];
        Neuron_head(remove_index,:) = [];
        DV_matrix(remove_index,:) = [];
        DV_matrix(:,remove_index) = [];
        geo_dist(remove_index,:) = [];
        geo_dist(:,remove_index) = [];
        LR_matrix(remove_index,:) = [];
        LR_matrix(:,remove_index) = [];
        PA_matrix(remove_index,:) = [];
        PA_matrix(:,remove_index) = [];
        X_rot(remove_index,:) = [];
        Y_rot(remove_index,:) = [];
        Z_rot(remove_index,:) = [];
        X_rot_norm(remove_index,:) = [];
        
        reg_node_pot = -reg_based_node_potential(X,Y,Z,X_rot,Y_rot,Z_rot,lambda_reg_1,lambda_reg_2);
        col_node_pot = -col_based_node_potential(all_col_int,data_int,Neuron_head,lambda_col_1,lambda_col_2);
        [assignment,cost] = munkres(reg_node_pot+col_node_pot);
        for n = 1:size(assignment,2)
            if or(assignment(1,n) == 0,isnan(assignment(1,n)))
                labels{n,k} = '';
            else
                labels{n,k} = Neuron_head{assignment(1,n),1};
            end
        end
    end
    
    top_labels = calculate_top_labels(labels);
    true_labels = marker_name;
    top_correct = 0;
    top_3_correct = 0;
    top_5_correct = 0;
    for tr_lb = 1:size(true_labels,1)
        if ~isempty(find(strcmp(true_labels(tr_lb,1),top_labels(marker_index(tr_lb,1)).top_labels(1,1))))
            top_correct = top_correct + 1;
        end
        if ~isempty(find(strcmp(true_labels(tr_lb,1),top_labels(marker_index(tr_lb,1)).top_labels(1,1:min(3,size(top_labels(marker_index(tr_lb,1)).top_labels,2))))))
            top_3_correct = top_3_correct + 1;
        end
        if ~isempty(find(strcmp(true_labels(tr_lb,1),top_labels(marker_index(tr_lb,1)).top_labels(1,1:min(5,size(top_labels(marker_index(tr_lb,1)).top_labels,2))))))
            top_5_correct = top_5_correct + 1;
        end
    end
    top_correct_frac = [top_correct_frac;top_correct/size(true_labels,1)];
    top_3_correct_frac = [top_3_correct_frac;top_3_correct/size(true_labels,1)];
    top_5_correct_frac = [top_5_correct_frac;top_5_correct/size(true_labels,1)];

    pred_struct(cnt).name = data{1,f};
    pred_struct(cnt).lambda_reg_1 = lambda_reg_1;
    pred_struct(cnt).lambda_reg_2 = lambda_reg_2;
    pred_struct(cnt).lambda_col_1 = lambda_col_1;
    pred_struct(cnt).lambda_col_2 = lambda_col_2;
    pred_struct(cnt).top_correct_frac = top_correct_frac;
    pred_struct(cnt).top_3_correct_frac = top_3_correct_frac;
    pred_struct(cnt).top_5_correct_frac = top_5_correct_frac;
    cnt = cnt + 1;
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

% %%%% maximum accuracy for each data across parameters
% name = {};
% score = [];
% for i = 1:size(pred_struct,2)
%     name = cat(1,name,pred_struct(i).name);
%     score = [score;pred_struct(i).correct_frac];
% end
% 
% max_correct_frac = zeros(size(data,2),1);
% for f = 1:size(data,2)
%     index = find(strcmp(data(1,f),name));
%     max_correct_frac(f,1) = max(score(index,:));
% end