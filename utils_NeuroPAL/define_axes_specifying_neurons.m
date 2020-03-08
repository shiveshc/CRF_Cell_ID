%%%%% function to specify axes neurons (A_neuron, P_neuron etc..) which
%%%%% help in defining consistent PA, LR, DV axes based on PCA or
%%%%% optimzation.

function axes_neurons_to_neuron_map = define_axes_specifying_neurons(thisimage_r,thisimage_g,thisimage_b,mu_marker)
    input_string = {'A_neuron','P_neuron','L_neuron','R_neuron','D_neuron','V_neuron'};
    A_neuron = [];P_neuron = [];L_neuron = [];R_neuron = [];D_neuron = [];V_neuron = [];
    locs = [];
    buttons = [];
    ind_marker_img_used = zeros(6,1);
    for i = 1:6
        g_norm = imadjust(max(mat2gray(thisimage_g),[],3),[0;0.7],[0;1],0.8);
        b_norm = imadjust(max(mat2gray(thisimage_b),[],3),[0;0.3],[0;1],0.8);
        r_norm = imadjust(max(mat2gray(thisimage_r),[],3),[0;0.4],[0;1],0.8);
        rgb_img = cat(3,r_norm,g_norm,b_norm);
        % rgb_img(30:31,330:345,:) = 1;
        figure,imshow(rgb_img)
        ['select ',input_string{1,i},'. If not present than right click.']
        [x,y,button] = ginput(1);
        close all
        locs = [locs;x,y];
        buttons = [buttons;button];
        
    end
    %%%%% create axes to neuron map. If not left click, map the selected
    %%%%% points to mu_r
    dist_mat = repmat(diag(locs*locs'),1,size(mu_marker(:,1:2),1)) + repmat(diag(mu_marker(:,1:2)*mu_marker(:,1:2)')',size(locs,1),1) - 2*locs*mu_marker(:,1:2)';
    [sort_dist,sort_index] = sort(dist_mat,2);
    axes_neurons_to_neuron_map = sort_index(:,1);
    axes_neurons_to_neuron_map(find(buttons(:,1) == 3),:) = 0;
end