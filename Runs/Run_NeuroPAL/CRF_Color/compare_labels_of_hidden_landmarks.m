function score = compare_labels_of_hidden_landmarks(curr_labels,marker_index,marker_name,landmarks,Neuron_head)

correct = 0;
cnt = 0;
for i = 1:size(marker_name,1)
    if isempty(find(marker_index(i,1) == landmarks))
        landmark_label = marker_name{i,1};
        assigned_label = Neuron_head{curr_labels(marker_index(i,1),1),1};
        if strcmp(landmark_label,assigned_label)
            correct = correct + 1;
        end
        cnt = cnt + 1;
        
%         figure,imshow(max(mat2gray(thisimage_r),[],3))
%         hold on
%         scatter(mu_r(landmark_to_neuron_map(i,1),2),mu_r(landmark_to_neuron_map(i,1),1),'r','LineWidth', 2)
%         title(['true:',landmark_label,', assigned:',assigned_label])
        
%         %%% create confidence map
%         landmark_index = find(strcmp(Neuron_head,landmark_label));
%         [sort_belief,sort_index] = sort(conserved_nodeBel(:,landmark_index),'Descend');
%         max_prob = max(sort_belief(1:10,:));
%         min_prob = min(sort_belief(1:10,:));
%         figure,imshow(max(mat2gray(thisimage_r),[],3))
%         hold on
%         for k = 1:10
%             prob = (sort_belief(k,1)-min_prob)/(max_prob-min_prob);
%             color = prob*[1,0,0] + (1-prob)*[1,1,1];
%             scatter(mu_r(sort_index(k,1),2),mu_r(sort_index(k,1),1),'c',color,'LineWidth', 2)
%         end
%         scatter(mu_r(landmark_to_neuron_map(i,1),2),mu_r(landmark_to_neuron_map(i,1),1),'b','LineWidth', 2)
    end
end
score = correct/cnt;

