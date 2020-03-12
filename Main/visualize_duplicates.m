%%%% function to visualize handling of nodes with duplicate labels.
function visualize_duplicates(thisimage_r,mu,node_label,Neuron_head)

    uniq_label = unique(node_label(:,1));
    for i = 1:size(uniq_label,1)
        curr_label_index = find(node_label(:,1) == uniq_label(i,1));
        if size(curr_label_index,1) > 1
            figure,imshow(max(mat2gray(thisimage_r),[],3))
            hold on
            for k = 1:size(curr_label_index,1)
                scatter(mu(curr_label_index(k),2),mu(curr_label_index(k),1),'r','LineWidth', 2)
                text(mu(curr_label_index(k),2),mu(curr_label_index(k),1)+10,num2str(node_label(curr_label_index(k),3)),'Color','white','FontSize',7)
            end
            [sort_curr_score,sort_index] = sort(node_label(curr_label_index,3),'descend');
            scatter(mu(curr_label_index(sort_index(1,1)),2),mu(curr_label_index(sort_index(1,1)),1),'g','LineWidth', 2)
            caxis([0,0.5])
            title(Neuron_head{node_label(curr_label_index(1,1),1)})
        end
    end % break here to visualize each duplicate
end