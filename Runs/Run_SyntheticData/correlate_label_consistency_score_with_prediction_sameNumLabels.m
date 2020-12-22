%%%% compare correlation between label consistency score and prediction
%%%% accuracy


function score_matrix = correlate_label_consistency_score_with_prediction_sameNumLabels()
load('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\PACE\AtlasTest\data_neuron_relationship.mat')

results_dir = 'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\PACE\AtlasTest\Results\rand130_randLandmarks_sameNumLabels\LB';

score_matrix = [];
fileList = dir(results_dir);
for i = 1:size(fileList,1)
    name = fileList(i).name;
    if ~strcmp('.',name) && ~strcmp('..',name)
        load([results_dir,'\',fileList(i).name]);
        
        tr_label = experiments(1).keep_neurons_from_atlas;
        pred_label = experiments(1).node_label(:,1);
        for j = 1:size(pred_label,1)
            if pred_label(j,1) ~= 0
                if experiments(1).keep_neurons_from_atlas(pred_label(j,1),1) == tr_label(j,1)
                    correct = 1;
                else
                    correct = 0;
                end
                
                %%% create 'geodist_r' as it is not saved in the results
                K = 6;
                pos = [experiments(1).X,experiments(1).Y,experiments(1).Z];
                euc_dist = repmat(diag(pos*pos'),1,size(pos,1)) + repmat(diag(pos*pos')',size(pos,1),1) - 2*pos*pos';
                [sort_euc_dist,sort_index] = sort(euc_dist,2);
                adj = zeros(size(pos,1),size(pos,1));
                for i = 1:size(adj,1)
                    adj(i,sort_index(i,2:K+1)) = 1;
                end
                adj = max(adj,adj');
                G = graph(adj);
                geo_dist_r = distances(G);

                [bin_pos_score, angle_score, geo_score] = consistency_score(j,experiments(1).keep_neurons_from_atlas,pred_label,pos(:,1),pos(:,2),pos(:,3),X_rot,Y_rot,Z_rot,PA_matrix,LR_matrix,DV_matrix,geo_dist,geo_dist_r,Neuron_head);
                score_matrix = [score_matrix;correct,bin_pos_score, angle_score, geo_score];
            end
        end
    end
end
end

function [bin_pos_score, angle_score, geo_score] = consistency_score(i,keep_neurons_from_atlas,node_label,X,Y,Z,X_rot,Y_rot,Z_rot,PA_matrix,LR_matrix,DV_matrix,geo_dist,geo_dist_r,Neuron_head)
    
    PA_consistency_score = 0;
    LR_consistency_score = 0;
    DV_consistency_score = 0;
    geodist_consistency_score = 0;
    angle_consistency_score = 0;
    cnt = 0;
    i_label = keep_neurons_from_atlas(node_label(i,1),1);
        for j = 1:size(node_label,1)
            if i ~= j && node_label(j,1) ~= 0
                j_label = keep_neurons_from_atlas(node_label(j,1),1);
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

                angle_consistency_score = angle_consistency_score + get_angle_consistency_score(X_rot,Y_rot,Z_rot,X,Y,Z,i,j,i_label,j_label,Neuron_head);
                
                cnt = cnt + 1;
            end
        end
        
    bin_pos_score = (PA_consistency_score + LR_consistency_score + DV_consistency_score)/(3*cnt);
    geo_score = geodist_consistency_score/cnt;
    angle_score = angle_consistency_score/cnt;
end

function score = get_angle_consistency_score(X_rot,Y_rot,Z_rot,X,Y,Z,i,j,i_label,j_label,Neuron_head)
    
    p_prime = [X(i,1),Y(i,1),Z(i,1)];
    q_prime = [X(j,1),Y(j,1),Z(j,1)];
    
    p = [X_rot(i_label,1),Y_rot(i_label,1),Z_rot(i_label,1)];
    q = [X_rot(j_label,1),Y_rot(j_label,1),Z_rot(j_label,1)];
    
    score = (p-q)*(p_prime - q_prime)'/(norm(p-q)*norm(p_prime-q_prime));
end
    
