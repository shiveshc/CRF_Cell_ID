%%%% compare correlation between label consistency score and prediction
%%%% accuracy


function score_matrix = correlate_label_consistency_score_with_prediction()
load('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\data_neuron_relationship_annotation_updated.mat')

in_direc = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_updAtlas\Results_multi_3D_data_20190706_5_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_updAtlas\Results_multi_3D_data_20190706_21_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_updAtlas\Results_multi_3D_data_20190710_2_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_updAtlas\Results_multi_3D_data_20190710_5_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_updAtlas\Results_multi_3D_data_20190710_8_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_updAtlas\Results_multi_3D_data_20190710_10_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_updAtlas\Results_multi_3D_data_20190710_22_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_updAtlas\Results_multi_3D_data_20190710_24_numLand0\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos\Results_updAtlas\Results_multi_3D_data_20190710_27_numLand0\'};
    
data = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos_Col\data_20190706_5.mat',...
'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos_Col\data_20190706_21.mat',...
'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos_Col\data_20190710_2.mat',...
'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos_Col\data_20190710_5.mat',...
'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos_Col\data_20190710_8.mat',...
'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos_Col\data_20190710_10.mat',...
'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos_Col\data_20190710_22.mat',...
'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos_Col\data_20190710_24.mat',...
'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RelPos_Col\data_20190710_27.mat'};

score_matrix = [];
for f = 1:size(in_direc,2)
    fileList = dir(in_direc{1,f});
    load(data{1,f});
        for n = 1:size(fileList,1)
            name = fileList(n).name;
            if ~strcmp('.',name) && ~strcmp('..',name)
                load([in_direc{1,f},fileList(n).name]);
        
                tr_label = marker_name;
                pred_label = experiments(1).node_label(:,1);
                for j = 1:size(tr_label,1)
                    if strcmp(tr_label{j,1},experiments(1).Neuron_head{pred_label(marker_index(j,1),1),1})
                        correct = 1;
                    else
                        correct = 0;
                    end
                
                    %%% create 'geodist_r' as it is not saved in the results
                    K = 6;
                    pos = mu_r;
                    euc_dist = repmat(diag(pos*pos'),1,size(pos,1)) + repmat(diag(pos*pos')',size(pos,1),1) - 2*pos*pos';
                    [sort_euc_dist,sort_index] = sort(euc_dist,2);
                    adj = zeros(size(pos,1),size(pos,1));
                    for i = 1:size(adj,1)
                        adj(i,sort_index(i,2:K+1)) = 1;
                    end
                    adj = max(adj,adj');
                    G = graph(adj);
                    geo_dist_r = distances(G);

                    [bin_pos_score, angle_score, geo_score] = consistency_score(j,pred_label,pos(:,1),pos(:,2),pos(:,3),X_rot,Y_rot,Z_rot,PA_matrix,LR_matrix,DV_matrix,geo_dist,geo_dist_r,experiments(1).Neuron_head,Neuron_head,angle_vec_atlas);
                    score_matrix = [score_matrix;correct,bin_pos_score, angle_score, geo_score];
                end
            end
        end
    end
end

function [bin_pos_score, angle_score, geo_score] = consistency_score(i,node_label,X,Y,Z,X_rot,Y_rot,Z_rot,PA_matrix,LR_matrix,DV_matrix,geo_dist,geo_dist_r,saved_neuron_head,Neuron_head,angle_vec_atlas)
    
    PA_consistency_score = 0;
    LR_consistency_score = 0;
    DV_consistency_score = 0;
    geodist_consistency_score = 0;
    angle_consistency_score = 0;
    cnt = 0;
    i_label = find(strcmp(saved_neuron_head{node_label(i,1),1},Neuron_head));
        for j = 1:size(node_label,1)
            if i ~= j && node_label(j,1) ~= 0
                j_label = find(strcmp(saved_neuron_head{node_label(j,1),1},Neuron_head));
                if X(i,1) < X(j,1)
                    PA_consistency_score = PA_consistency_score + PA_matrix(i_label,j_label);
                else
                    PA_consistency_score = PA_consistency_score + PA_matrix(j_label,i_label);
                end

                if Y(i,1) < Y(j,1)
                    LR_consistency_score = LR_consistency_score + LR_matrix(i_label,j_label);
                else
                    LR_consistency_score = LR_consistency_score + LR_matrix(j_label,i_label);
                end

                if Z(i,1) < Z(j,1)
                    DV_consistency_score = DV_consistency_score + DV_matrix(i_label,j_label);
                else
                    DV_consistency_score = DV_consistency_score + DV_matrix(j_label,i_label);
                end

                geodist_consistency_score = geodist_consistency_score + exp(-(geo_dist(i_label,j_label)-geo_dist_r(i,j)).^2);

                angle_consistency_score = angle_consistency_score + get_angle_consistency_score(X_rot,Y_rot,Z_rot,X,Y,Z,i,j,i_label,j_label,angle_vec_atlas);
                
                cnt = cnt + 1;
            end
        end
        
    bin_pos_score = (PA_consistency_score + LR_consistency_score + DV_consistency_score)/(3*cnt);
    geo_score = geodist_consistency_score/cnt;
    angle_score = angle_consistency_score/cnt;
end

function score = get_angle_consistency_score(X_rot,Y_rot,Z_rot,X,Y,Z,i,j,i_label,j_label,angle_vec_atlas)
    
    p_prime = [X(i,1),Y(i,1),Z(i,1)];
    q_prime = [X(j,1),Y(j,1),Z(j,1)];
    angle_matrix = angle_vec_atlas*((p_prime - q_prime)'/norm(p_prime-q_prime));
    angle_matrix = reshape(angle_matrix,size(X_rot,1),size(X_rot,1));
    score = angle_matrix(i_label,j_label);
end    