%%% function for pre-processing images before identity labelling. This
%%% function is very simialr to another version of
%%% 'preprocess_landmark_data.m'. This script is to be used for segmenting
%%% individual stacks of worms whereas the other script can be used for
%%% preparing files for annotation whole brain videos after tracking.
%%% Performs - 1.segmentation of RFP and CyOFP images
%%%            2.preparing labeled_img

function preprocess_data(out_folder,data_name,img_1,varargin)
    if ~exist('img_1')
        disp("Please provide the image in which cells are to be identified. Exiting.")
        return
    end
    
    num_imgs = length(varargin);
    
    %%%%% segment channel in which cells are to identified
    disp("Segmenting main image channel ... ")
    cmap = rand(130+1,3);
    cmap(1,:) = [0,0,0];
    full_img = img_1;
    
    % remove noise
    temp = double(full_img);
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
    
    % determine thresh parameter
    accept_thresh = 'n';
    while ~strcmp(accept_thresh,'y')
        [thresh_param,accept_thresh,index_thresh_new] = determine_thresh_param(curr_img);
    end
    
    % segment stack
    [comp_weights,mu,sigma] = em_segmentation_3d_gpu(curr_img,index_thresh_new,size(index_thresh_new,1),[],[],5,1,cmap);
    
    % store segmentation information
    dummy = struct();
    dummy.comp_weights = comp_weights;
    dummy.mu = mu;
    dummy.sigma = sigma;
    dummy.index_thresh = index_thresh_new;
    seg_struct.(['f',num2str(1)]) = dummy;
    
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
    close all
    cmap = rand(150,3);
    cmap(1,:) = [0,0,0];
    for i = 1:size(labeled_img,4)
        name = strcat([out_folder,'\segmented_img_r.tif']);
        indexed_img_to_tiff(labeled_img(:,:,:,i),cmap,name)
    end
    mu_r = mu;
    labeled_img_r = labeled_img;

    %%%% segment other channels that may have landmarks
    mu_other_channels = {};
    labeled_img_other_channels = {};
    for chan = 1:num_imgs
        disp(['Segmenting image channel ',num2str(chan+1)])
        cmap = rand(130+1,3);
        cmap(1,:) = [0,0,0];
        full_img = varargin{1,chan};

        % remove background in images
        temp = double(full_img);
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
        
        % determine thresh parameter
        accept_thresh = 'n';
        while ~strcmp(accept_thresh,'y')
            [thresh_param,accept_thresh,index_thresh_new] = determine_thresh_param(curr_img);
        end
        
        % segment image stack
        [comp_weights,mu,sigma] = em_segmentation_3d_gpu(curr_img,index_thresh_new,size(index_thresh_new,1),[],[],5,1,cmap);
        
        % store segmentation information
        dummy = struct();
        dummy.comp_weights = comp_weights;
        dummy.mu = mu;
        dummy.sigma = sigma;
        dummy.index_thresh = index_thresh_new;
        seg_struct.(['f',num2str(i)]) = dummy;
    
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
        close all
        cmap = rand(150,3);
        cmap(1,:) = [0,0,0];
        for i = 1:size(labeled_img,4)
            name = strcat([out_folder,'\segmented_img_',num2str(chan+1),'.tif']);
            indexed_img_to_tiff(labeled_img(:,:,:,i),cmap,name)
        end
        
        mu_other_channels{1,chan} = mu;
        labeled_img_other_channels{1,chan} = labeled_img;
    end
    disp('Segmentation finished')

    %%%%%%% manually label landmark neurons by taking user input %%%%%%%
    mu_to_map = [mu_r(:,2),mu_r(:,1)];
    landmark_channels = input("Enter which channels to use for specifying landmarks e.g [2,4] else enter blank (single quotes) -");
    if isempty(landmark_channels)
        mu_c = [];
        landmark_names_orig = {};
        landmark_to_neuron_map_orig = [];
    else
        mu_c = [];
        landmark_names_orig = {};
        landmark_to_neuron_map_orig = [];
        for i = 1:size(landmark_channels,2)
            curr_channel = landmark_channels(1,i);
            mu_curr_channel = [];
            if curr_channel == 1
                curr_mu = mu_r;
                channel_done = 'n';
                while strcmp(channel_done,'n')
                    imshow(max(mat2gray(img_1),[],3))
                    caxis([0,0.3])
                    hold on
                    scatter(curr_mu(:,2),curr_mu(:,1),'.r')
                    title(["Click on a landmark cell, then enter its name on terminal e.g. 'RMEL'"])
                    ["Click on a landmark cell, then enter its name on terminal e.g. 'RMEL'"]
                    [x,y,button] = ginput(1);
                    close all
                    mu_curr_channel = [mu_curr_channel;x,y]; %%% check the coordinate system
                    curr_name = input("Enter name of the selected landmark e.g. 'RMEL' -");
                    landmark_names_orig = cat(1,landmark_names_orig,curr_name);
                    channel_done = input("If done with this channel, enter 'y' -");
                end
                dist_mat = repmat(diag(mu_curr_channel*mu_curr_channel'),1,size(mu_to_map(:,1:2),1)) + repmat(diag(mu_to_map(:,1:2)*mu_to_map(:,1:2)')',size(mu_curr_channel,1),1) - 2*mu_curr_channel*mu_to_map(:,1:2)';
                [sort_dist,sort_index] = sort(dist_mat,2);
                landmark_to_neuron_map_orig = [landmark_to_neuron_map_orig;sort_index(:,1)];
                mu_c = [mu_c;mu_r(sort_index(:,1),:)];
            else
                curr_mu = mu_other_channels{1,curr_channel-1};
                channel_done = 'n';
                while strcmp(channel_done,'n')
                    imshow(max(mat2gray(varargin{1,curr_channel-1}),[],3))
                    caxis([0,0.3])
                    hold on
                    scatter(curr_mu(:,2),curr_mu(:,1),'.r')
                    title(["Click on a landmark cell, then enter its name on terminal e.g. 'RMEL'"])
                    ["Click on a landmark cell, then enter its name on terminal e.g.'RMEL'"]
                    [x,y,button] = ginput(1);
                    close all
                    mu_curr_channel = [mu_curr_channel;x,y]; %%% check the coordinate system
                    curr_name = input("Enter name of the selected landmark e.g. 'RMEL' -");
                    landmark_names_orig = cat(1,landmark_names_orig,curr_name);
                    channel_done = input("If done with this channel, enter 'y' -");
                end
                dist_mat = repmat(diag(mu_curr_channel*mu_curr_channel'),1,size(mu_to_map(:,1:2),1)) + repmat(diag(mu_to_map(:,1:2)*mu_to_map(:,1:2)')',size(mu_curr_channel,1),1) - 2*mu_curr_channel*mu_to_map(:,1:2)';
                [sort_dist,sort_index] = sort(dist_mat,2);
                landmark_to_neuron_map_orig = [landmark_to_neuron_map_orig;sort_index(:,1)];
                mu_c = [mu_c;mu_r(sort_index(:,1),:)];
            end
        end
    end
    
%     visualize_landmark_to_neuron_map(thisimage_r,thisimage_c,mu_r,mu_c,landmark_to_neuron_map_orig,landmark_names_orig)

    %%%%%% define axes specifying neurons
    axes_neurons_to_neuron_map = define_axes_specifying_neurons(img_1,[mu_r(:,2),mu_r(:,1),mu_r(:,3)]);
    
    %%%%%%%%%%%%%% define axes param based on PCA
    mu_r_centered = mu_r - repmat(mean(mu_r),size(mu_r,1),1);
    [coeff,score,latent] = pca(mu_r_centered);
    axes_param = input("Enter PCA coefficients for specifying AP, LR and DV axes e.g [1,2,3] or [1,3,2] -");
    PA = coeff(:,axes_param(1,1))';
    PA = PA/norm(PA);
    LR = coeff(:,axes_param(1,2))';
    LR = LR/norm(LR);
    DV = coeff(:,axes_param(1,3))';
    DV = DV/norm(DV);
    A_neuron = axes_neurons_to_neuron_map(1,1);
    P_neuron = axes_neurons_to_neuron_map(2,1);
    L_neuron = axes_neurons_to_neuron_map(3,1);
    R_neuron = axes_neurons_to_neuron_map(4,1);
    D_neuron = axes_neurons_to_neuron_map(5,1);
    V_neuron = axes_neurons_to_neuron_map(6,1);
    if (mu_r_centered(A_neuron,:)-mu_r_centered(P_neuron,:))*PA' < 0
        PA = -PA;
    end
    if D_neuron ~= 0 && V_neuron ~= 0
        if (mu_r_centered(V_neuron,:)-mu_r_centered(D_neuron,:))*DV' < 0
            DV = -DV;
        end
        if cross(DV,PA)*LR' < 0
            LR = -LR;
        end
    else
        if (mu_r_centered(R_neuron,:)-mu_r_centered(L_neuron,:))*LR' < 0
            LR = -LR;
        end
        if cross(PA,LR)*DV' < 0
            DV = -DV;
        end
    end
    figure,scatter3(mu_r_centered(:,1),mu_r_centered(:,2),mu_r_centered(:,3),'.r')
    hold on
    plot3([0,50*PA(1,1)],[0,50*PA(1,2)],[0,50*PA(1,3)],'b','LineWidth',2.5)
    plot3([0,20*LR(1,1)],[0,20*LR(1,2)],[0,20*LR(1,3)],'g','LineWidth',2.5)
    plot3([0,20*DV(1,1)],[0,20*DV(1,2)],[0,20*DV(1,3)],'k','LineWidth',2.5)

    %%%%%% set axes defining parameters and save data
    ind_PCA = 1; %%% set to 1 if want to generate axes using PCA else set 0
    specify_PA = 0; %%% set to 1 if for non-pca axes generation - want to define PA axis using PCA
    landmark_names = landmark_names_orig;
    landmark_to_neuron_map = landmark_to_neuron_map_orig;
    save([out_folder,'\',data_name],'mu_r','mu_other_channels','labeled_img_r','labeled_img_other_channels','landmark_names','landmark_to_neuron_map','mu_c','img_1','varargin','cmap','ind_PCA','specify_PA','axes_neurons_to_neuron_map')
end