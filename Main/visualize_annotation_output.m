% function to visualize labels. Score on top of figure is total consistency
% score
function visualize_annotation_output(thisimage_r,mu,curr_labels,Neuron_head,i)
    if ~isempty(i)
        figure,imshow(max(mat2gray(thisimage_r),[],3), 'border', 'tight')
        hold on
        caxis([0, 0.4])
        scatter(mu(i,1),mu(i,2),'r','LineWidth', 2)
        text(mu(i,1) + 3, mu(i,2) +3, Neuron_head{curr_labels(i,1)}, 'Color', 'w', 'Fontsize', 8)
    else
        figure,imshow(max(mat2gray(thisimage_r),[],3), 'border', 'tight')
        hold on
        caxis([0, 0.4])
        for n = 1:size(mu, 1)
            scatter(mu(n,1),mu(n,2),'.r')
            text(mu(n,1) + 3, mu(n,2) +3, Neuron_head{curr_labels(n,1)}, 'Color', 'w', 'Fontsize', 6)
        end
    end
end