function plot_all_labels_on_img(thisimage_r,mu_r,curr_labels,Neuron_head)

    scaled_r = imadjustn(mat2gray(max(thisimage_r,[],3)),[0,0.4],[0,1],1);
    temp_img = -max(scaled_r,[],3);
    temp_img = (temp_img - min(temp_img(:)))/(max(temp_img(:))-min(temp_img(:)));
    
    figure, imshow(fliplr(temp_img)), hold on
    for i = 1:size(curr_labels,1)
        text(size(thisimage_r,2) - mu_r(i,1),mu_r(i,2),Neuron_head{curr_labels(i,1)},'FontSize',6, 'Color', [0.5,0,0], 'Rotation', 330)
    end
end
        