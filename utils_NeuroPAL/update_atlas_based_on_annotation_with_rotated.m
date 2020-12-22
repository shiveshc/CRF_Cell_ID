%%%% function to updata relative position relationships based on annotated
%%%% data using both straight worm and rotated worm data (update atlas)


%%%%%%%%%%%%%% load data_neuron_relationship.mat
load('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\data_neuron_relationship.mat')

%%%%%%%%%%%%%% run compiled_data.m to create compiled annotated data
compiled_data = compile_annotated_data_with_rotated;


%%%%%%%%%%%%%% create average positional relationship matrices based on annotated data
all_annotated_neurons = {};
for i = 1:size(compiled_data,2)
    all_annotated_neurons = cat(1,all_annotated_neurons,compiled_data(i).marker_name);
end
uniq_annotated_neurons = unique(all_annotated_neurons);

%%%% remove annotated neurons that are not present in the atlas. This is because
%%%% we do not have data for their true positional relationships
rem_index = [];
for i = 1:size(uniq_annotated_neurons,1)
    if isempty(find(strcmp(uniq_annotated_neurons{i,1},Neuron_head)))
        rem_index = [rem_index;i];
    end
end
uniq_annotated_neurons(rem_index,:) = [];

PA_consistency_matrix = zeros(size(uniq_annotated_neurons,1),size(uniq_annotated_neurons,1));
DV_consistency_matrix = zeros(size(uniq_annotated_neurons,1),size(uniq_annotated_neurons,1));
LR_consistency_matrix = zeros(size(uniq_annotated_neurons,1),size(uniq_annotated_neurons,1));
%%% count matrices to keep track of in how many data sets neuron pair
%%% was observed
PA_consistency_count = zeros(size(uniq_annotated_neurons,1),size(uniq_annotated_neurons,1));
DV_consistency_count = zeros(size(uniq_annotated_neurons,1),size(uniq_annotated_neurons,1));
LR_consistency_count = zeros(size(uniq_annotated_neurons,1),size(uniq_annotated_neurons,1));

% for angle feature calculation based on annotated data, we are going to
% store an array of p-q vectors for all annotated neurons
angle_vec_data = zeros(size(uniq_annotated_neurons,1)*size(uniq_annotated_neurons,1),3);
angle_count_data = zeros(size(uniq_annotated_neurons,1)*size(uniq_annotated_neurons,1),1);
cnt = 1;
for m = 1:size(uniq_annotated_neurons,1)
    neuron_1 = uniq_annotated_neurons{m,1};
    for n = m+1:size(uniq_annotated_neurons,1)
        neuron_2 = uniq_annotated_neurons{n,1};
        for d = 1:size(compiled_data,2)
            neuron_1_index = compiled_data(d).marker_index(find(strcmp(neuron_1,compiled_data(d).marker_name)),1);
            neuron_2_index = compiled_data(d).marker_index(find(strcmp(neuron_2,compiled_data(d).marker_name)),1);
            if ~isempty(neuron_1_index) && ~isempty(neuron_2_index)
                if compiled_data(d).X_rot(neuron_1_index) < compiled_data(d).X_rot(neuron_2_index)
                    PA_consistency_matrix(m,n) = PA_consistency_matrix(m,n) + 1;
                else
                    PA_consistency_matrix(n,m) = PA_consistency_matrix(n,m) + 1;
                end
                PA_consistency_count(m,n) = PA_consistency_count(m,n) + 1;
                PA_consistency_count(n,m) = PA_consistency_count(n,m) + 1;
                
                if compiled_data(d).Y_rot(neuron_1_index) < compiled_data(d).Y_rot(neuron_2_index)
                    LR_consistency_matrix(m,n) = LR_consistency_matrix(m,n) + 1;
                else
                    LR_consistency_matrix(n,m) = LR_consistency_matrix(n,m) + 1;
                end
                LR_consistency_count(m,n) = LR_consistency_count(m,n) + 1;
                LR_consistency_count(n,m) = LR_consistency_count(n,m) + 1;
                
                if compiled_data(d).Z_rot(neuron_1_index) < compiled_data(d).Z_rot(neuron_2_index)
                    DV_consistency_matrix(m,n) = DV_consistency_matrix(m,n) + 1;
                else
                    DV_consistency_matrix(n,m) = DV_consistency_matrix(n,m) + 1;
                end
                DV_consistency_count(m,n) = DV_consistency_count(m,n) + 1;
                DV_consistency_count(n,m) = DV_consistency_count(n,m) + 1;
            end
        end
    end
end
for m = 1:size(uniq_annotated_neurons,1)
    neuron_1 = uniq_annotated_neurons{m,1};
    for n = 1:size(uniq_annotated_neurons,1)
        neuron_2 = uniq_annotated_neurons{n,1};
        for d = 1:size(compiled_data,2)
            neuron_1_index = compiled_data(d).marker_index(find(strcmp(neuron_1,compiled_data(d).marker_name)),1);
            neuron_2_index = compiled_data(d).marker_index(find(strcmp(neuron_2,compiled_data(d).marker_name)),1);
            if ~isempty(neuron_1_index) && ~isempty(neuron_2_index)
                p = [compiled_data(d).X_rot(neuron_1_index),compiled_data(d).Y_rot(neuron_1_index),compiled_data(d).Z_rot(neuron_1_index)];
                q = [compiled_data(d).X_rot(neuron_2_index),compiled_data(d).Y_rot(neuron_2_index),compiled_data(d).Z_rot(neuron_2_index)];
                ange_vec_data_index = sub2ind([size(uniq_annotated_neurons,1),size(uniq_annotated_neurons,1)],m,n);
                angle_vec_data(ange_vec_data_index,:) = angle_vec_data(ange_vec_data_index,:) + [p - q]/norm([p-q]);
                angle_count_data(ange_vec_data_index,1) = angle_count_data(ange_vec_data_index,1) + 1;
            end
        end
    end
end
mean_PA_consistency_matrix = PA_consistency_matrix./PA_consistency_count;
mean_LR_consistency_matrix = LR_consistency_matrix./LR_consistency_count;
mean_DV_consistency_matrix = DV_consistency_matrix./DV_consistency_count;
mean_angle_vec_data = angle_vec_data./repmat(angle_count_data,1,3);

%%% fill average PA, LR, DV information from data into PA, LR, DV matrices
PA_matrix_data = PA_matrix;
LR_matrix_data = LR_matrix;
DV_matrix_data = DV_matrix;
angle_vec_atlas = zeros(size(Neuron_head,1)*size(Neuron_head,1),3);
cnt = 1;
for m = 1:size(Neuron_head,1)
    neuron_1 = find(strcmp(Neuron_head{m,1},uniq_annotated_neurons));
    for n = 1:size(Neuron_head,1)
        neuron_2 = find(strcmp(Neuron_head{n,1},uniq_annotated_neurons));
        angle_vec_atlas_index = sub2ind([size(Neuron_head,1),size(Neuron_head,1)],m,n); 
        if ~isempty(neuron_1) && ~isempty(neuron_2)
            angle_vec_data_index = sub2ind([size(uniq_annotated_neurons,1),size(uniq_annotated_neurons,1)],neuron_1,neuron_2);
            if isnan(mean_PA_consistency_matrix(neuron_1,neuron_2))
            else
                PA_matrix_data(m,n) = mean_PA_consistency_matrix(neuron_1,neuron_2);
            end
            
            if isnan(mean_LR_consistency_matrix(neuron_1,neuron_2))
            else
                LR_matrix_data(m,n) = mean_LR_consistency_matrix(neuron_1,neuron_2);
            end
            
            if isnan(mean_DV_consistency_matrix(neuron_1,neuron_2))
            else
                DV_matrix_data(m,n) = mean_DV_consistency_matrix(neuron_1,neuron_2);
            end
            
            if isnan(mean_angle_vec_data(angle_vec_data_index,1))
                p = [X_rot(m,1),Y_rot(m,1),Z_rot(m,1)];
                q = [X_rot(n,1),Y_rot(n,1),Z_rot(n,1)]; 
                angle_vec_atlas(angle_vec_atlas_index,:) = [p-q]/norm([p-q]);
            else
                angle_vec_atlas(angle_vec_atlas_index,:) = mean_angle_vec_data(angle_vec_data_index,:);
            end
        else
            p = [X_rot(m,1),Y_rot(m,1),Z_rot(m,1)];
            q = [X_rot(n,1),Y_rot(n,1),Z_rot(n,1)]; 
            angle_vec_atlas(angle_vec_atlas_index,:) = [p-q]/norm([p-q]);
        end
    end
end

%%% save new atlas file based on annotated data
DV_matrix = DV_matrix_data;
LR_matrix = LR_matrix_data;
PA_matrix = PA_matrix_data;
save('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\data_neuron_relationship_annotation_updated_latonly','DV_matrix','eigval','eigvec','ganglion','geo_dist','LR_matrix','LR_neurons','LR_neurons_matrix','Neuron_head','PA_matrix','VC_neurons','X_rot','X_rot_norm','Y_rot','Z_rot','angle_vec_atlas')
