%%%% function to visualize landmark channel to rfp image mapping result

function visualize_landmark_to_neuron_map(thisimage_r,thisimage_c,mu_r,mu_c,landmark_to_neuron_map,landmark_names)
    
    img = cat(2,mat2gray(thisimage_r),max(mat2gray(thisimage_c),mat2gray(thisimage_r)));
    figure,imshow(max(img,[],3),[])
    caxis([0,0.3])
    hold on
    
    for i = 1:size(landmark_to_neuron_map,1)
        c = rand(1,3);
        scatter(mu_r(landmark_to_neuron_map(i,1),2),mu_r(landmark_to_neuron_map(i,1),1),'o','MarkerFaceColor',c)
        text(mu_r(landmark_to_neuron_map(i,1),2),mu_r(landmark_to_neuron_map(i,1),1),landmark_names{i},'Color',c,'FontSize',14)
    end
end

    