%%%% function to visualize landmark channel to rfp image mapping result

function visualize_marker_to_neuron_map(thisimage_r,thisimage_g,mu_r,marker_to_neuron_map,marker_name)
    
    img = cat(2,mat2gray(thisimage_r),mat2gray(thisimage_g));
    figure,imshow(max(img,[],3),[])
    caxis([0,0.3])
    hold on
    
    for i = 1:size(marker_to_neuron_map,1)
        c = rand(1,3);
        scatter(mu_r(marker_to_neuron_map(i,1),1),mu_r(marker_to_neuron_map(i,1),2),'o','MarkerFaceColor',c)
        text(mu_r(marker_to_neuron_map(i,1),1),mu_r(marker_to_neuron_map(i,1),2),marker_name{i},'Color',c,'FontSize',14)
    end
end

    