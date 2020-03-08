function [x_thresh_new,y_thresh_new,z_thresh_new,index_thresh_new] = remove_close_neurons(x_thresh,y_thresh,z_thresh,index_thresh,img,dist_param,dist_param_z)
        
    dist_mat = squareform(pdist([x_thresh,y_thresh]));
    dist_mat = dist_mat + 1000000*eye(size(dist_mat));
    close = dist_mat <= nthroot(2*dist_param,2);
    dist_mat_z = squareform(pdist([z_thresh],'minkowski',1));
    dist_mat_z = dist_mat_z + 1000*eye(size(dist_mat_z));
    close_z = dist_mat_z <= dist_param_z;

    [neu1,neu2] = ind2sub(size(dist_mat),find(close.*close_z));
    remove_list = unique(neu1);
    [img_int,I] = sort(img(index_thresh(remove_list)));
    remove_list = remove_list(I);
    if ~isempty(remove_list)
        x_thresh(remove_list(1,1)) = [];
        y_thresh(remove_list(1,1)) = [];
        z_thresh(remove_list(1,1)) = [];
        index_thresh(remove_list(1,1)) = [];
    end
        
    
%     figure,imshow(max(img,[],3),[])
%     hold on
%     scatter(y_thresh,x_thresh,'.r')
    while ~isempty(remove_list)
        dist_mat = squareform(pdist([x_thresh,y_thresh]));
        dist_mat = dist_mat + 1000000*eye(size(dist_mat));
        close = dist_mat <= nthroot(2*dist_param,2);
        dist_mat_z = squareform(pdist([z_thresh],'minkowski',1));
        dist_mat_z = dist_mat_z + 1000*eye(size(dist_mat_z));
        close_z = dist_mat_z <= dist_param_z;

        [neu1,neu2] = ind2sub(size(dist_mat),find(close.*close_z));
        remove_list = unique(neu1);
        [img_int,I] = sort(img(index_thresh(remove_list)));
        remove_list = remove_list(I);
        
        if ~isempty(remove_list)
%             scatter(y_thresh(remove_list(1,1)),x_thresh(remove_list(1,1)),'.y')
            x_thresh(remove_list(1,1)) = [];
            y_thresh(remove_list(1,1)) = [];
            z_thresh(remove_list(1,1)) = [];
            index_thresh(remove_list(1,1)) = [];
        end
    end
    x_thresh_new = x_thresh;
    y_thresh_new = y_thresh;
    z_thresh_new = z_thresh;
    index_thresh_new = index_thresh;
end