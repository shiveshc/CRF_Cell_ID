%%%% function to generate node potentials based on color distributions
%%%% of NeuroPAL annotated

function node_pot = col_based_node_potential(all_col_int,data_int,Neuron_head)

    all_int_names = {};
    for i = 1:size(all_col_int,2)
        all_int_names = cat(1,all_int_names,all_col_int(i).name);
    end
    
    node_pot = zeros(size(data_int,2),size(Neuron_head,1));
    for i = 1:size(data_int,2)
        curr_int_R = mean(double(data_int(i).R_int));
        curr_int_G = mean(double(data_int(i).G_int));
        curr_int_B = mean(double(data_int(i).B_int));
%         curr_hsv = rgb2hsv([curr_int_R,curr_int_G,curr_int_B]/max([curr_int_R,curr_int_G,curr_int_B]));
        curr_hsv = mean(rgb2hsv(cat(3,data_int(i).R_int,data_int(i).G_int,data_int(i).B_int)));
        curr_hsv = cat(2,curr_hsv(1,1,1),curr_hsv(1,1,2),curr_hsv(1,1,3));
        for n = 1:size(Neuron_head,1)
            state_neuron = Neuron_head{n,1};
            if ~isempty(find(strcmp(state_neuron,all_int_names)))
                all_int_index = find(strcmp(state_neuron,all_int_names));
                all_int_R = double(all_col_int(all_int_index).R_int);
                all_int_G = double(all_col_int(all_int_index).G_int);
                all_int_B = double(all_col_int(all_int_index).B_int);
%                 all_hsv = rgb2hsv([all_int_R,all_int_G,all_int_B]./repmat(max([all_int_R,all_int_G,all_int_B],[],2),1,3));
                all_hsv = rgb2hsv(cat(3,all_col_int(all_int_index).R_int,all_col_int(all_int_index).G_int,all_col_int(all_int_index).B_int));
                all_hsv = cat(2,all_hsv(:,1,1),all_hsv(:,1,2),all_hsv(:,1,3));
                %%%% dist in RGB space
                dist = mahal([curr_int_R,curr_int_G,curr_int_B],[all_int_R,all_int_G,all_int_B]);
%                 dist = sum(([curr_int_R,curr_int_G,curr_int_B] - mean([all_int_R,all_int_G,all_int_B])).^2);
%                 dist = [curr_int_R,curr_int_G,curr_int_B]*mean([all_int_R,all_int_G,all_int_B])'/(norm(mean([all_int_R,all_int_G,all_int_B]))*norm([curr_int_R,curr_int_G,curr_int_B]));
                
%                 %%%% dist in HSV space
%                 dist = mahal(curr_hsv(:,1:2),all_hsv(:,1:2));

                node_pot(i,n) = 10*exp(-dist^2);
            else
                node_pot(i,n) = 0.1;
            end
        end
    end
    node_pot(find(node_pot<0.01)) = 0.01;
end