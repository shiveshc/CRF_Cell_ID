%%% function for pre-processing images before identity labelling. This
%%% function is very simialr to another version of
%%% 'preprocess_landmark_data.m'. This script is to be used for segmenting
%%% individual stacks of worms whereas the other script can be used for
%%% preparing files for annotation whole brain videos after tracking.
%%% Performs - 1.segmentation of RFP and CyOFP images
%%%            2.preparing labeled_img
addpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\WholeBrainCode_PACE')
out_folder = 'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\20190414_AML70_unc47_gcy32_CyOFP\9\';

%%%%% segment red-channel
    cmap = rand(130+1,3);
    cmap(1,:) = [0,0,0];
    full_img = thisimage_r(:,:,:,1);

    for i = 1:1;
        temp = double(full_img(:,:,:,i));
        %%% remove background in images
        temp_back_corr = temp;
        t_index = find(ones(size(temp,1),size(temp,2)));
        [t_x,t_y] = ind2sub([size(temp,1),size(temp,2)],t_index);
        t_data = [t_x,t_y,ones(size(t_x,1),1)];
        for z = 1:size(temp,3)
            curr = temp(:,:,z);
            coeff = (t_data'*t_data)\t_data'*curr(t_index);
            curr_back = reshape(t_data*coeff,size(curr,1),size(curr,2));
            temp_back_corr(:,:,z) = temp(:,:,z) - curr_back;
        end
        
%         temp = mat2gray(temp);
%         temp_gray = pcafilt_imgvol(temp,coeff);
%         temp_gray = mat2gray(temp);
        
        % Bandpass filtering
%         MIJ.createImage('img',temp,true);
%         MIJ.run('Bandpass Filter...', 'filter_large=50 filter_small=3 suppress=None tolerance=5 saturate process');
%         temp_bandpass = mat2gray(MIJ.getImage('img'));
%         MIJ.run('Close All')
%         for k = 1:size(temp,3)
%             temp_bandpass(:,:,k) = gaussianbpf(mat2gray(temp(:,:,k)),1,5);
%         end
%         init_psf = ones(5,5,3);
%         [temp_gray_deblur,P] = deconvblind(temp_gray,init_psf);
        
        
        % 3d image
%         temp_med = ordfilt3(temp_gray(:,:,1:15),'med',3);
%         temp_med_gauss = imgaussfilt3(temp_med,1);        
        for z = 1:size(temp,3)
            temp_med(:,:,z) = medfilt2(temp(:,:,z));            
            temp_med_gauss(:,:,z) = imgaussfilt(temp_med(:,:,z),1);
%             temp_log(:,:,z) = imfilter(temp_med_gauss(:,:,z),lapGaufilt,'replicate');
        end
        
        temp_med_gauss = mat2gray(temp_med_gauss);
        curr_img = temp_med_gauss;
        curr_img = imhmax(curr_img,0.01);
        loc_max = imregionalmax(curr_img);
        index = find(loc_max);
        [x,y,z] = ind2sub(size(loc_max),index);
        peaks = curr_img(index);
        thresh_val = prctile(peaks,50);
        peak_thresh =  peaks(peaks>thresh_val);
        x_thresh = x(peaks>thresh_val);
        y_thresh = y(peaks>thresh_val);
        z_thresh = z(peaks>thresh_val); 
        index_thresh = index(peaks>thresh_val);
        [x_thresh_new,y_thresh_new,z_thresh_new,index_thresh_new] = remove_close_neurons(x_thresh,y_thresh,z_thresh,index_thresh,curr_img);
        figure,imshow(max(temp_med_gauss,[],3),[])
        caxis([0,0.6])
        hold on
        scatter(y_thresh,x_thresh,'.g')
        scatter(y_thresh_new,x_thresh_new,'.r')
        
        [comp_weights,mu,sigma] = em_segmentation_3d_gpu(curr_img,index_thresh_new,size(index_thresh_new,1),[],[],5,1,cmap);
        
        dummy = struct();
        dummy.comp_weights = comp_weights;
        dummy.mu = mu;
        dummy.sigma = sigma;
        dummy.index_thresh = index_thresh_new;
        seg_struct.(['f',num2str(i)]) = dummy;
    end
    
    seg_type = 3; % segmentation done on 2d image or 3d image
    segNorm1 = false(size(full_img));
    allCenters = cell(1,size(segNorm1,4));
    numObjFound = zeros(size(segNorm1,4),1);
    for i = 1:size(segNorm1,4)
        temp_img = zeros(size(segNorm1,1),size(segNorm1,2),size(segNorm1,3));
        temp_img_filt = temp_img;
        temp_img(seg_struct.(['f',num2str(i)]).index_thresh) = 1;
        for m = 1:size(seg_struct.(['f',num2str(i)]).index_thresh,1)
            [p,q,r] = ind2sub(size(temp_img),seg_struct.(['f',num2str(i)]).index_thresh(m,1));
            temp_img_filt(max(p-10,1):min(p+10,size(temp_img,1)),max(q-10,1):min(q+10,size(temp_img,2)),max(r-2,1):min(r+2,size(temp_img,3))) = 1;
        end
        index_filt = find(temp_img_filt);
        [x,y,z] = ind2sub(size(temp_img),index_filt);
        img_data = [x,y,z];
        
        mu = seg_struct.(['f',num2str(i)]).mu;
        sigma = seg_struct.(['f',num2str(i)]).sigma;
        comp_weights = seg_struct.(['f',num2str(i)]).comp_weights;
        
        likelihood = zeros(size(img_data,1),size(comp_weights,2));
        posterior = zeros(size(img_data,1),size(comp_weights,2));
        prob_cutoff = zeros(size(img_data,1),size(comp_weights,2));
        for k = 1:size(comp_weights,2)
            likelihood(:,k) = mvnpdf(img_data,mu(k,:),sigma(:,:,k));
            posterior(:,k) = comp_weights(k)*mvnpdf(img_data,mu(k,:),sigma(:,:,k));
            prob_cutoff(:,k) = 1/sqrt((2*pi)^3*det(sigma(:,:,k)))*exp(-1/2*1);
        end
        posterior = posterior./repmat(sum(posterior,2),1,size(comp_weights,2));
    
        likelihood = likelihood.*(likelihood>= prob_cutoff);
        posterior = posterior.*(likelihood>= prob_cutoff);
        [max_posterior,max_comp] = max(posterior,[],2);
        img_pos = zeros(size(segNorm1(:,:,:,i)));
        img_pos(index_filt) = max_comp.*(max_posterior > 0.6);
        segNorm1(:,:,:,i) = img_pos;
        
        dummy = struct();
        for k = 1:size(comp_weights,2)
            Area = size(find(img_pos == k),1);
            Centroid = mu(k,:);
            pixelIndex = find(img_pos == k);
            [x,y,z] = ind2sub(size(img_pos),pixelIndex);
            PixelList = [x,y,z];
            dummy(k).Area = Area;
            dummy(k).Centroid = Centroid;
            dummy(k).PixelList = PixelList;
        end
        allCenters{1,i} = dummy;
        numObjFound(i,1) = size(comp_weights,2);
    end
    
    labeled_img = ones(size(segNorm1(:,:,:,:)));
    for i = 1:size(allCenters,2)
        for n = 1:size(allCenters{i},2)
            curr_pixelList = allCenters{i}(n).PixelList;
            if ~isempty(curr_pixelList)
                for k = 1:size(curr_pixelList,1)
                    labeled_img(curr_pixelList(k,1),curr_pixelList(k,2),curr_pixelList(k,3),i) = n+1;
                end
            end
        end
    end

    cmap = rand(150,3);
    cmap(1,:) = [0,0,0];
    for i = 1:size(labeled_img,4)
        name = strcat([out_folder,'labeled_img_r.tif']);
        indexed_img_to_tiff(labeled_img(:,:,:,i),cmap,name)
    end
    
mu_r = mu;
labeled_img_r = labeled_img;

%%%% segment CyOFP-channel (landmarks)
    cmap = rand(130+1,3);
    cmap(1,:) = [0,0,0];
    full_img = thisimage_c;

    for i = 1:1;
        temp = double(full_img(:,:,:,i));
        %%% remove background in images
        temp_back_corr = temp;
        t_index = find(ones(size(temp,1),size(temp,2)));
        [t_x,t_y] = ind2sub([size(temp,1),size(temp,2)],t_index);
        t_data = [t_x,t_y,ones(size(t_x,1),1)];
        for z = 1:size(temp,3)
            curr = temp(:,:,z);
            coeff = (t_data'*t_data)\t_data'*curr(t_index);
            curr_back = reshape(t_data*coeff,size(curr,1),size(curr,2));
            temp_back_corr(:,:,z) = temp(:,:,z) - curr_back;
        end
               
        for z = 1:size(temp,3)
            temp_med(:,:,z) = medfilt2(temp(:,:,z));            
            temp_med_gauss(:,:,z) = imgaussfilt(temp_med(:,:,z),1);
%             temp_log(:,:,z) = imfilter(temp_med_gauss(:,:,z),lapGaufilt,'replicate');
        end
        
        temp_med_gauss = mat2gray(temp_med_gauss);
        curr_img = temp_med_gauss;
%         curr_img = imhmax(curr_img,0.01);
        loc_max = imregionalmax(curr_img);
        index = find(loc_max);
        [x,y,z] = ind2sub(size(loc_max),index);
        peaks = curr_img(index);
        thresh_val = prctile(peaks,99.987);
        peak_thresh =  peaks(peaks>thresh_val);
        x_thresh = x(peaks>thresh_val);
        y_thresh = y(peaks>thresh_val);
        z_thresh = z(peaks>thresh_val); 
        index_thresh = index(peaks>thresh_val);
        [x_thresh_new,y_thresh_new,z_thresh_new,index_thresh_new] = remove_close_neurons(x_thresh,y_thresh,z_thresh,index_thresh,curr_img);
        figure,imshow(max(temp_med_gauss,[],3),[])
        caxis([0,0.6])
        hold on
        scatter(y_thresh,x_thresh,'.g')
        scatter(y_thresh_new,x_thresh_new,'.r')
        
        [comp_weights,mu,sigma] = em_segmentation_3d_gpu(curr_img,index_thresh_new,size(index_thresh_new,1),[],[],5,1,cmap);
        
        dummy = struct();
        dummy.comp_weights = comp_weights;
        dummy.mu = mu;
        dummy.sigma = sigma;
        dummy.index_thresh = index_thresh_new;
        seg_struct.(['f',num2str(i)]) = dummy;
    end
    
    seg_type = 3; % segmentation done on 2d image or 3d image
    segNorm1 = false(size(full_img));
    allCenters = cell(1,size(segNorm1,4));
    numObjFound = zeros(size(segNorm1,4),1);
    for i = 1:size(segNorm1,4)
        temp_img = zeros(size(segNorm1,1),size(segNorm1,2),size(segNorm1,3));
        temp_img_filt = temp_img;
        temp_img(seg_struct.(['f',num2str(i)]).index_thresh) = 1;
        for m = 1:size(seg_struct.(['f',num2str(i)]).index_thresh,1)
            [p,q,r] = ind2sub(size(temp_img),seg_struct.(['f',num2str(i)]).index_thresh(m,1));
            temp_img_filt(max(p-10,1):min(p+10,size(temp_img,1)),max(q-10,1):min(q+10,size(temp_img,2)),max(r-2,1):min(r+2,size(temp_img,3))) = 1;
        end
        index_filt = find(temp_img_filt);
        [x,y,z] = ind2sub(size(temp_img),index_filt);
        img_data = [x,y,z];
        
        mu = seg_struct.(['f',num2str(i)]).mu;
        sigma = seg_struct.(['f',num2str(i)]).sigma;
        comp_weights = seg_struct.(['f',num2str(i)]).comp_weights;
        
        likelihood = zeros(size(img_data,1),size(comp_weights,2));
        posterior = zeros(size(img_data,1),size(comp_weights,2));
        prob_cutoff = zeros(size(img_data,1),size(comp_weights,2));
        for k = 1:size(comp_weights,2)
            likelihood(:,k) = mvnpdf(img_data,mu(k,:),sigma(:,:,k));
            posterior(:,k) = comp_weights(k)*mvnpdf(img_data,mu(k,:),sigma(:,:,k));
            prob_cutoff(:,k) = 1/sqrt((2*pi)^3*det(sigma(:,:,k)))*exp(-1/2*1);
        end
        posterior = posterior./repmat(sum(posterior,2),1,size(comp_weights,2));
    
        likelihood = likelihood.*(likelihood>= prob_cutoff);
        posterior = posterior.*(likelihood>= prob_cutoff);
        [max_posterior,max_comp] = max(posterior,[],2);
        img_pos = zeros(size(segNorm1(:,:,:,i)));
        img_pos(index_filt) = max_comp.*(max_posterior > 0.6);
        segNorm1(:,:,:,i) = img_pos;
        
        dummy = struct();
        for k = 1:size(comp_weights,2)
            Area = size(find(img_pos == k),1);
            Centroid = mu(k,:);
            pixelIndex = find(img_pos == k);
            [x,y,z] = ind2sub(size(img_pos),pixelIndex);
            PixelList = [x,y,z];
            dummy(k).Area = Area;
            dummy(k).Centroid = Centroid;
            dummy(k).PixelList = PixelList;
        end
        allCenters{1,i} = dummy;
        numObjFound(i,1) = size(comp_weights,2);
    end
    
    labeled_img = ones(size(segNorm1(:,:,:,:)));
    for i = 1:size(allCenters,2)
        for n = 1:size(allCenters{i},2)
            curr_pixelList = allCenters{i}(n).PixelList;
            if ~isempty(curr_pixelList)
                for k = 1:size(curr_pixelList,1)
                    labeled_img(curr_pixelList(k,1),curr_pixelList(k,2),curr_pixelList(k,3),i) = n+1;
                end
            end
        end
    end

    cmap = rand(150,3);
    cmap(1,:) = [0,0,0];
    for i = 1:size(labeled_img,4)
        name = strcat([out_folder,'labeled_img_c.tif']);
        indexed_img_to_tiff(labeled_img(:,:,:,i),cmap,name)
    end
    
mu_c = mu;
labeled_img_c = labeled_img;

%%%%%%% manually label landmark neurons by taking user input %%%%%%%
landmark_names_orig = cell(size(mu_c,1),1);
for i = 1:size(mu_c,1)
    imshow(max(labeled_img_c,[],3),cmap)
    hold on
    scatter(mu_c(i,2),mu_c(i,1),'*y')
    curr_name = input('Landmark name?');
    landmark_names_orig{i,1} = curr_name;
    close all
end

%%%%%%% map CyOFP channel neurons to RFP neurons %%%%%%%
not_landmarks = find(strcmp(landmark_names_orig,''));
mu_c_final = mu_c;
mu_c_final(not_landmarks,:) = [];
landmark_names_orig(not_landmarks,:) = [];
dist_mat = repmat(diag(mu_c_final*mu_c_final'),1,size(mu_r,1)) + repmat(diag(mu_r*mu_r')',size(mu_c_final,1),1) - 2*mu_c_final*mu_r';
[sort_dist,sort_index] = sort(dist_mat,2);
landmark_to_neuron_map_orig = sort_index(:,1);
visualize_landmark_to_neuron_map(thisimage_r,thisimage_c,mu_r,mu_c,landmark_to_neuron_map_orig,landmark_names_orig)


axes_neurons_to_neuron_map = [];
ind_PCA = 1; %%% set to 1 if want to generate axes using PCA else set 0
specify_PA = 0; %%% set to 1 if for non-pca axes generation - want to define PA axis using PCA
landmark_names = cat(1,landmark_names_orig,add_landmarks);
landmark_to_neuron_map = cat(1,landmark_to_neuron_map_orig,add_landmark_to_neuron_map);
save([inp_folder,data_name],'labeled_img_c','labeled_img_r','landmark_names','landmark_to_neuron_map','mu_c','mu_r','thisimage_c','thisimage_r','cmap','ind_PCA','specify_PA','axes_neurons_to_neuron_map')