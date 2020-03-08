%%%% function to compare AP, LR, DV relationship among neurons

%%%% Inputs - 
%%%% 1. 'compiled_data' variable obtained by running
%%%%    'compile_annotated_data.m'
%%%% 2. PA, LR, DV matrix in the atlas. PA and DV matrix can be used from
%%%%    3d atlas or 2d atlas. 2D atlas can be obtained by running
%%%%    'create_2D_atlas.m'. LR matrix must be obtained from 3D data.
%%%% 3. Neuron_head_2D and Neuron_head variable

all_annotated_neurons = {};
for i = 1:size(compiled_data,2)
    all_annotated_neurons = cat(1,all_annotated_neurons,compiled_data(i).marker_name);
end
uniq_annotated_neurons = unique(all_annotated_neurons);

%%%% remove annotated neurons that are present in 2D atlas. This is because
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
angle_matrix = zeros(size(uniq_annotated_neurons,1),size(uniq_annotated_neurons,1));
%%% count matrices to keep track of in how many data sets neuron pair
%%% was observed
PA_consistency_count = zeros(size(uniq_annotated_neurons,1),size(uniq_annotated_neurons,1));
DV_consistency_count = zeros(size(uniq_annotated_neurons,1),size(uniq_annotated_neurons,1));
LR_consistency_count = zeros(size(uniq_annotated_neurons,1),size(uniq_annotated_neurons,1));
angle_count = zeros(size(uniq_annotated_neurons,1),size(uniq_annotated_neurons,1));

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
                
                p = [compiled_data(d).X_rot(neuron_2_index),compiled_data(d).Z_rot(neuron_2_index)];
                q = [compiled_data(d).X_rot(neuron_1_index),compiled_data(d).Z_rot(neuron_1_index)];
                angle = acosd((p-q)*[1;0]/(norm(p-q)));
                angle_matrix(m,n) = angle_matrix(m,n) + angle;
                angle_matrix(n,m) = angle_matrix(n,m) + 180 - angle;
%                 if angle > 90    
%                     angle_matrix(m,n) = angle_matrix(m,n) + 180 - angle;
%                     angle_matrix(n,m) = angle_matrix(n,m) + 180 - angle;
%                 else
%                     angle_matrix(m,n) = angle_matrix(m,n) + angle;
%                     angle_matrix(n,m) = angle_matrix(n,m) + angle;
%                 end
                angle_count(m,n) = angle_count(m,n) + 1;
                angle_count(n,m) = angle_count(n,m) + 1;
                
            end
        end
    end
end
mean_PA_consistency_matrix = PA_consistency_matrix./PA_consistency_count;
mean_LR_consistency_matrix = LR_consistency_matrix./LR_consistency_count;
mean_DV_consistency_matrix = DV_consistency_matrix./DV_consistency_count;
mean_angle_matrix = angle_matrix./angle_count;

true_PA_matrix = NaN(size(uniq_annotated_neurons,1),size(uniq_annotated_neurons,1));
true_LR_matrix = NaN(size(uniq_annotated_neurons,1),size(uniq_annotated_neurons,1));
true_DV_matrix = NaN(size(uniq_annotated_neurons,1),size(uniq_annotated_neurons,1));
true_angle_matrix = NaN(size(uniq_annotated_neurons,1),size(uniq_annotated_neurons,1));
for m = 1:size(uniq_annotated_neurons,1)
    neuron_1 = uniq_annotated_neurons{m,1};
    neuron_1_index = find(strcmp(neuron_1,Neuron_head));
    for n = m+1:size(uniq_annotated_neurons,1)
        neuron_2 = uniq_annotated_neurons{n,1};
        neuron_2_index = find(strcmp(neuron_2,Neuron_head));
        if PA_consistency_count(m,n) == 0
            true_PA_matrix(m,n) = NaN;
            true_PA_matrix(n,m) = NaN;
        else
            if PA_matrix(neuron_1_index,neuron_2_index) == 1
                true_PA_matrix(m,n) = 1;
                true_PA_matrix(n,m) = 0;
            else
                true_PA_matrix(n,m) = 1;
                true_PA_matrix(m,n) = 0;
            end
        end
        if DV_consistency_count(m,n) == 0
            true_DV_matrix(m,n) = NaN;
            true_DV_matrix(n,m) = NaN;
        else
            if DV_matrix(neuron_1_index,neuron_2_index) == 1
                true_DV_matrix(m,n) = 1;
                true_DV_matrix(n,m) = 0;
            else
                true_DV_matrix(n,m) = 1;
                true_DV_matrix(m,n) = 0;
            end
        end
        if angle_count(m,n) == 0
            true_angle_matrix(m,n) = NaN;
            true_angle_matrix(n,m) = NaN;
        else
            p = [X_rot(neuron_2_index,1),Z_rot(neuron_2_index,1)];
            q = [X_rot(neuron_1_index,1),Z_rot(neuron_1_index,1)];
            angle = acosd((p-q)*[1;0]/(norm(p-q)));
            true_angle_matrix(m,n) = angle;
            true_angle_matrix(n,m) = 180 -angle;
%             if angle > 90
%                 true_angle_matrix(m,n) = 180-angle;
%                 true_angle_matrix(n,m) = 180-angle;
%             else
%                 true_angle_matrix(m,n) = angle;
%                 true_angle_matrix(n,m) = angle;
%             end
        end
    end
end
for m = 1:size(uniq_annotated_neurons,1)
    neuron_1 = uniq_annotated_neurons{m,1};
    neuron_1_index = find(strcmp(neuron_1,Neuron_head));
    for n = m+1:size(uniq_annotated_neurons,1)
        neuron_2 = uniq_annotated_neurons{n,1};
        neuron_2_index = find(strcmp(neuron_2,Neuron_head));
        if LR_consistency_count(m,n) == 0
            true_LR_matrix(m,n) = NaN;
            true_LR_matrix(n,m) = NaN;
        else
            if LR_matrix(neuron_1_index,neuron_2_index) == 1
                true_LR_matrix(m,n) = 1;
                true_LR_matrix(n,m) = 0;
            else
                true_LR_matrix(n,m) = 1;
                true_LR_matrix(m,n) = 0;
            end
        end
    end
end
% %%%%%% sort all matrices by ganglion
% annotated_ganglion = zeros(size(uniq_annotated_neurons,1),1);
% for m = 1:size(uniq_annotated_neurons,1)
%     neuron_1 = uniq_annotated_neurons{m,1};
%     annotated_ganglion(m,1) = ganglion_2D(find(strcmp(neuron_1,Neuron_head_2D)),1);
% end
% [sort_annotated_ganglion,sort_index] = sort(annotated_ganglion);
% true_PA_matrix = true_PA_matrix(sort_index,sort_index);
% true_LR_matrix = true_LR_matrix(sort_index,sort_index);
% true_DV_matrix = true_DV_matrix(sort_index,sort_index);
% true_angle_matrix = true_angle_matrix(sort_index,sort_index);
% mean_PA_consistency_matrix = mean_PA_consistency_matrix(sort_index,sort_index);
% mean_LR_consistency_matrix = mean_LR_consistency_matrix(sort_index,sort_index);
% mean_DV_consistency_matrix = mean_DV_consistency_matrix(sort_index,sort_index);
% mean_angle_matrix = mean_angle_matrix(sort_index,sort_index);

%%%%%% sort all matrices by position
annotated_X = zeros(size(uniq_annotated_neurons,1),1);
annotated_Y = zeros(size(uniq_annotated_neurons,1),1);
annotated_Z = zeros(size(uniq_annotated_neurons,1),1);
for m = 1:size(uniq_annotated_neurons,1)
    neuron_1 = uniq_annotated_neurons{m,1};
    annotated_X(m,1) = X_rot(find(strcmp(neuron_1,Neuron_head)),1);
    annotated_Y(m,1) = Y_rot(find(strcmp(neuron_1,Neuron_head)),1);
    annotated_Z(m,1) = Z_rot(find(strcmp(neuron_1,Neuron_head)),1);
end
[sort_X,sort_index_X] = sort(annotated_X);
true_PA_matrix = true_PA_matrix(sort_index_X,sort_index_X);
mean_PA_consistency_matrix = mean_PA_consistency_matrix(sort_index_X,sort_index_X);
[sort_Y,sort_index_Y] = sort(annotated_Y);
true_LR_matrix = true_LR_matrix(sort_index_Y,sort_index_Y);
mean_LR_consistency_matrix = mean_LR_consistency_matrix(sort_index_Y,sort_index_Y);
[sort_Z,sort_index_Z] = sort(annotated_Z);
true_DV_matrix = true_DV_matrix(sort_index_Z,sort_index_Z);
mean_DV_consistency_matrix = mean_DV_consistency_matrix(sort_index_Z,sort_index_Z);

diff_PA_matrix = true_PA_matrix - mean_PA_consistency_matrix;
diff_LR_matrix = true_LR_matrix - mean_LR_consistency_matrix;
diff_DV_matrix = true_DV_matrix - mean_DV_consistency_matrix;
diff_angle_matrix = true_angle_matrix - mean_angle_matrix;
%%%%%%%%%%%%%%%%%%%%%%%%%% For figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% relationships across all data
data_PA_relationship = [];
data_LR_relationship = [];
data_DV_relationship = [];
data_angle_relationship = [];
true_PA_relationship = [];
true_LR_relationship = [];
true_DV_relationship = [];
true_angle_relationship = [];
for m = 1:size(uniq_annotated_neurons,1)
    for n = 1:size(uniq_annotated_neurons,1)
        data_PA_relationship = [data_PA_relationship;mean_PA_consistency_matrix(m,n)];
        data_LR_relationship = [data_LR_relationship;mean_LR_consistency_matrix(m,n)];
        data_DV_relationship = [data_DV_relationship;mean_DV_consistency_matrix(m,n)];
        data_angle_relationship = [data_angle_relationship;mean_angle_matrix(m,n)];
        
        true_PA_relationship = [true_PA_relationship;true_PA_matrix(m,n)];
        true_LR_relationship = [true_LR_relationship;true_LR_matrix(m,n)];
        true_DV_relationship = [true_DV_relationship;true_DV_matrix(m,n)];
        true_angle_relationship = [true_angle_relationship;true_angle_matrix(m,n)];
    end
end

data_diff_PA_relationship = [];
data_diff_LR_relationship = [];
data_diff_DV_relationship = [];
data_diff_angle_relationship = [];
for m = 1:size(uniq_annotated_neurons,1)
    for n = m+1:size(uniq_annotated_neurons,1)
        data_diff_PA_relationship = [data_diff_PA_relationship;diff_PA_matrix(m,n)];
        data_diff_LR_relationship = [data_diff_LR_relationship;diff_LR_matrix(m,n)];
        data_diff_DV_relationship = [data_diff_DV_relationship;diff_DV_matrix(m,n)];
        data_diff_angle_relationship = [data_diff_angle_relationship;diff_angle_matrix(m,n)];
    end
end

%%%% individual data consistency accuracy
% for d = 1:size(compiled_data,2)
%     curr_