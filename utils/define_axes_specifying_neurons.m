%%%%% function to specify axes neurons (A_neuron, P_neuron etc..) which
%%%%% help in defining consistent PA, LR, DV axes based on PCA or
%%%%% optimzation.

function axes_neurons_to_neuron_map = define_axes_specifying_neurons(thisimage_r,thisimage_marker,mu_r,mu_marker,marker_to_neuron_map)
    input_string = {'(A)-P neuron','A-(P) neuron','(L)-R neuron','L-(R) neuron','(D)-V neuron','D-(V) neuron'};
    A_neuron = [];P_neuron = [];L_neuron = [];R_neuron = [];D_neuron = [];V_neuron = [];
    locs = [];
    buttons = [];
    ind_marker_img_used = zeros(6,1);
    for i = 1:6
        if ~isempty(thisimage_r)
%             img = cat(2,max(mat2gray(thisimage_r),[],3),max(mat2gray(thisimage_marker),[],3));
            img = max(mat2gray(thisimage_r),[],3);
            figure,imshow(img,[])
            caxis([0,0.4])
            title(['Click on ',input_string{1,i},'. If not required then right click.'])
            ['select ',input_string{1,i},'. If not present than right click.']
            [x,y,button] = ginput(1);
            close all
            if x > size(thisimage_r,2)
                ind_marker_img_used(i,1) = 1;
                x = x - size(thisimage_r,2);
            end
            locs = [locs;x,y];
            buttons = [buttons;button];
        else
            img = max(mat2gray(thisimage_marker),[],3);
            figure,imshow(img,[])
            caxis([0,0.3])
            ['select ',input_string{1,i},'. If not present than right click.']
            [x,y,button] = ginput(1);
            close all
            ind_marker_img_used(i,1) = 1;
            locs = [locs;x,y];
            buttons = [buttons;button];
        end
        
    end
    %%%%% create axes to neuron map. If red image is used (ind_marker_img_used = 0) and not left click
    %%%%% then map neuron using red image. If marker image used then map to
    %%%%% markers first and then use marker_to_neuron_map variable to map
    %%%%% back to neurons in red image
    if ~isempty(thisimage_r)
        for i = 1:6
            if ind_marker_img_used(i,1) == 0
                dist_mat = repmat(diag(locs(i,:)*locs(i,:)'),1,size(mu_r(:,1:2),1)) + repmat(diag(mu_r(:,1:2)*mu_r(:,1:2)')',size(locs(i,:),1),1) - 2*locs(i,:)*mu_r(:,1:2)';
                [sort_dist,sort_index] = sort(dist_mat,2);
                axes_neurons_to_neuron_map(i,1) = sort_index(1,1);
            else
                dist_mat = repmat(diag(locs(i,:)*locs(i,:)'),1,size(mu_marker(:,1:2),1)) + repmat(diag(mu_marker(:,1:2)*mu_marker(:,1:2)')',size(locs(i,:),1),1) - 2*locs(i,:)*mu_marker(:,1:2)';
                [sort_dist,sort_index] = sort(dist_mat,2);
                axes_neurons_to_neuron_map(i,1) = marker_to_neuron_map(sort_index(1,1),1);
            end
        end
        axes_neurons_to_neuron_map(find(buttons(:,1) == 3),:) = 0;
    else
        dist_mat = repmat(diag(locs*locs'),1,size(mu_marker(:,1:2),1)) + repmat(diag(mu_marker(:,1:2)*mu_marker(:,1:2)')',size(locs,1),1) - 2*locs*mu_marker(:,1:2)';
        [sort_dist,sort_index] = sort(dist_mat,2);
        axes_neurons_to_neuron_map = sort_index(:,1);
        axes_neurons_to_neuron_map(find(buttons(:,1) == 3),:) = 0;
    end
end