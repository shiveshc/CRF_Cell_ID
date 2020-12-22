function node_label = duplicate_labels(curr_labels,X,Y,Z,PA_matrix,LR_matrix,DV_matrix,geo_dist,geo_dist_r,lambda_geo,clamped_neurons)

    % identify duplicate label nodes
    node_label = [curr_labels(:,1),zeros(size(curr_labels,1),2)];
    uniq_label = unique(node_label(:,1));
    for i = 1:size(uniq_label,1)
        curr_label_index = find(node_label(:,1) == uniq_label(i,1));
        if size(curr_label_index,1) > 1
            node_label(curr_label_index,2) = 1;
        end
    end
    
    % calculate score for nodes with duplicate labels
    for i = 1:size(curr_labels,1)
        if node_label(i,2) == 1
            node_label(i,3) = consistency_score(i,node_label(i,1),node_label,X,Y,Z,PA_matrix,LR_matrix,DV_matrix,geo_dist,geo_dist_r,lambda_geo);
        end
    end
    
    % visualize nodes with duplicate labels
%     visualize_duplicates(thisimage_r,mu,node_label,Neuron_head)
    
    % select label for duplicate label nodes
    for i = 1:size(uniq_label,1)
        curr_label_index = find(node_label(:,1) == uniq_label(i,1));
        if size(curr_label_index,1) > 1
            curr_label_score = node_label(curr_label_index,3);
            [sort_score,sort_index] = sort(curr_label_score);
            % fix label of highest score node and reassign remaining nodes,
            % also check landmark neurons
            clamped_neuron_present = clamped_neurons(find(ismember(clamped_neurons,curr_label_index)),1);
            if isempty(clamped_neuron_present)
                node_label(curr_label_index(sort_index(1:end-1,:)),1) = 0;
            else
                node_label(curr_label_index(not(curr_label_index == clamped_neuron_present),1),1) = 0;
            end
            % reassign all duplicate nodes
%             node_label(curr_label_index(sort_index(1:end,:)),1) = 0;
        end
    end
end

function score = consistency_score(i,i_label,node_label,X,Y,Z,PA_matrix,LR_matrix,DV_matrix,geo_dist,geo_dist_r,lambda_geo)
    PA_consistency_score = 0;
    LR_consistency_score = 0;
    DV_consistency_score = 0;
    geodist_consistency_score = 0;
        for j = 1:size(node_label,1)
            if i ~= j && node_label(j,2) == 0
                j_label = node_label(j,1);
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

                geodist_consistency_score = geodist_consistency_score + exp(-lambda_geo*(geo_dist(i_label,j_label)-geo_dist_r(i,j)).^2);
            end
        end
    max_consistency_score = size(find(node_label(:,2) == 0),1); %%% for node i wrt all other non-duplicate nodes
    score = (PA_consistency_score + LR_consistency_score + DV_consistency_score + geodist_consistency_score)/(4*max_consistency_score);
end