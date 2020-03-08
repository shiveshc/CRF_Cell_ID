%%%% function for Col + Reg based annotation of Neuropal data

addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Joint_tracking\TPS-RPM'))
addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GLMD_Demo'))
addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\CPD2'))

compiled_data = compile_annotated_data();
all_annotated_neurons = {};
for i = 1:size(compiled_data,2)
    all_annotated_neurons = cat(1,all_annotated_neurons,compiled_data(i).marker_name);
end
uniq_annotated_neurons = unique(all_annotated_neurons);

correct_fraction = [];
reg_col_accuracy = struct();
for i = 1:size(uniq_annotated_neurons,1)
    reg_col_accuracy(i).name = uniq_annotated_neurons{i,1};
    reg_col_accuracy(i).match_matrix = [];
end

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
    
    lambda_reg_1 = 1;
    lambda_reg_2 = 0.1;
    lambda_col_1 = 5;
    lambda_col_2 = 0.1;
    
    for a = 1:size(lambda_reg_1,2)
        for b = 1:size(lambda_reg_2,2)
            for c = 1:size(lambda_col_1,2)
                for d = 1:size(lambda_col_2,2)
                    reg_node_pot = -reg_based_node_potential(X,Y,Z,X_rot,Y_rot,Z_rot,lambda_reg_1(1,a),lambda_reg_2(1,b));
%                     col_node_pot = -col_based_node_potential(all_col_int,data_int,Neuron_head,lambda_col_1(1,c),lambda_col_2(1,d));
                    [assignment,cost] = munkres(reg_node_pot);

                    %%%% accuracy of hungarian based matching
                    correct = 0;
                    for i = 1:size(marker_name,1)
                        tr_label = marker_name{i,1};
                        curr_label = Neuron_head{assignment(1,marker_index(i,1),1)};
                        index = find(strcmp(marker_name{i,1},uniq_annotated_neurons));
                        if strcmp(tr_label,curr_label)
                            reg_col_accuracy(index).match_matrix = cat(1,reg_col_accuracy(index).match_matrix,1);
                            correct = correct + 1;
                        else
                            reg_col_accuracy(index).match_matrix = cat(1,reg_col_accuracy(index).match_matrix,0);
                        end
                    end
                end
            end
        end
    end
    correct_fraction = [correct_fraction;correct/size(marker_name,1)];
end

for i = 1:size(reg_col_accuracy,2)
    reg_col_accuracy(i).accuracy = nansum(reg_col_accuracy(i).match_matrix(:,1))/size(reg_col_accuracy(i).match_matrix(:,1),1);
end