%%%%% function to compare absolute position variability of neurons
%%%%% annotated in NeuroPAL strain

%%%% Inputs - 
%%%% 1. 'compiled_data' variable obtained by running
%%%%    'compile_annotated_data.m'.

all_annotated_neurons = {};
for i = 1:size(compiled_data,2)
    all_annotated_neurons = cat(1,all_annotated_neurons,compiled_data(i).marker_name);
end
uniq_annotated_neurons = unique(all_annotated_neurons);

neuron_positions_ap = NaN(size(compiled_data,2),size(uniq_annotated_neurons,1));
neuron_positions_dv = NaN(size(compiled_data,2),size(uniq_annotated_neurons,1));
neuron_positions_lr = NaN(size(compiled_data,2),size(uniq_annotated_neurons,1));

for i = 1:size(compiled_data,2)
    %%% most anterior neuron
    index_M3L = compiled_data(i).marker_index(find(strcmp('M3L',compiled_data(i).marker_name)),1);
    index_M3R = compiled_data(i).marker_index(find(strcmp('M3R',compiled_data(i).marker_name)),1);
    %%% most posterior neuron
    index_VD2 = compiled_data(i).marker_index(find(strcmp('VD2',compiled_data(i).marker_name)),1);
    
%     %%%% create normalized coordinates
%     if ~isempty(index_M3L) && ~isempty(index_M3R)
%         if compiled_data(i).X_rot(index_M3L,1) > compiled_data(i).X_rot(index_M3R,1)
%             X_norm = (compiled_data(i).X_rot - compiled_data(i).X_rot(index_VD2,1))/(compiled_data(i).X_rot(index_M3L,1) - compiled_data(i).X_rot(index_VD2,1));
%         else
%             X_norm = (compiled_data(i).X_rot - compiled_data(i).X_rot(index_VD2,1))/(compiled_data(i).X_rot(index_M3R,1) - compiled_data(i).X_rot(index_VD2,1));
%         end
%     elseif ~isempty(index_M3L)
%         X_norm = (compiled_data(i).X_rot - compiled_data(i).X_rot(index_VD2,1))/(compiled_data(i).X_rot(index_M3L,1) - compiled_data(i).X_rot(index_VD2,1));
%     else
%         X_norm = (compiled_data(i).X_rot - compiled_data(i).X_rot(index_VD2,1))/(compiled_data(i).X_rot(index_M3R,1) - compiled_data(i).X_rot(index_VD2,1));
%     end
    
    X_norm = (compiled_data(i).X_rot-min(compiled_data(i).X_rot))/(max(compiled_data(i).X_rot)-min(compiled_data(i).X_rot));
    Y_norm = (compiled_data(i).Y_rot-min(compiled_data(i).Y_rot))/(max(compiled_data(i).Y_rot)-min(compiled_data(i).Y_rot));
    Z_norm = (compiled_data(i).Z_rot-min(compiled_data(i).Z_rot))/(max(compiled_data(i).Z_rot)-min(compiled_data(i).Z_rot));
    
    for n = 1:size(neuron_positions_ap,2)
        curr_annotated_neuron_index = compiled_data(i).marker_index(find(strcmp(uniq_annotated_neurons{n,1},compiled_data(i).marker_name)),1);
        if ~isempty(curr_annotated_neuron_index)
            neuron_positions_ap(i,n) = X_norm(curr_annotated_neuron_index,1);
            neuron_positions_dv(i,n) = Z_norm(curr_annotated_neuron_index,1);
            neuron_positions_lr(i,n) = Y_norm(curr_annotated_neuron_index,1);
        end
    end
end

indexed_positions = [];
for i = 1:size(neuron_positions_ap,1)
    indexed_positions = cat(1,indexed_positions,cat(2,[1:1:size(uniq_annotated_neurons,1)]',neuron_positions_ap(i,:)',neuron_positions_lr(i,:)',neuron_positions_dv(i,:)'));
end