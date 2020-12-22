function [PA_score,LR_score,DV_score,geodist_score,tot_score] = consistency_scores(nNodes,curr_labels,X,Y,Z,PA_matrix,LR_matrix,DV_matrix,geo_dist,geo_dist_r)
    PA_consistency_score = zeros(nNodes,1);
    LR_consistency_score = zeros(nNodes,1);
    DV_consistency_score = zeros(nNodes,1);
    geodist_consistency_score = zeros(nNodes,1);
    for i = 1:nNodes
        i_label = curr_labels(i,1);
        for j = 1:nNodes
            if i ~= j
                j_label = curr_labels(j,1);
                if X(i,1) < X(j,1)
                    PA_consistency_score(i,1) = PA_consistency_score(i,1) + PA_matrix(i_label,j_label);
                else
                    PA_consistency_score(i,1) = PA_consistency_score(i,1) + PA_matrix(j_label,i_label);
                end

                if Y(i,1) < Y(j,1)
                    LR_consistency_score(i,1) = LR_consistency_score(i,1) + LR_matrix(i_label,j_label);
                else
                    LR_consistency_score(i,1) = LR_consistency_score(i,1) + LR_matrix(j_label,i_label);
                end

                if Z(i,1) < Z(j,1)
                    DV_consistency_score(i,1) = DV_consistency_score(i,1) + DV_matrix(i_label,j_label);
                else
                    DV_consistency_score(i,1) = DV_consistency_score(i,1) + DV_matrix(j_label,i_label);
                end

                geodist_consistency_score(i,1) = geodist_consistency_score(i,1) + exp(-0.1*(geo_dist(i_label,j_label)-geo_dist_r(i,j)).^2);
            end
        end
    end
    max_consistency_score = (nNodes - 1); %%% for fully connected CRF model (nNodes - 1 edges)
    PA_score = PA_consistency_score/max_consistency_score;
    LR_score = LR_consistency_score/max_consistency_score;
    DV_score = DV_consistency_score/max_consistency_score;
    geodist_score = geodist_consistency_score/max_consistency_score;
    tot_score = (PA_consistency_score + LR_consistency_score + DV_consistency_score)/(3*max_consistency_score);
 end