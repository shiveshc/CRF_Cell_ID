%%%%% first half code of GMM segmentation. Run to decide thresh_param and
%%%%% remove_close_neurons_distance param.

 cmap = rand(130+1,3);
    cmap(1,:) = [0,0,0];
    full_img = thisimage_r(:,:,:,1);

    for i = 1:1
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
%         curr_img = imhmax(curr_img,0.001);
        loc_max = imregionalmax(curr_img);
        index = find(loc_max);
        [x,y,z] = ind2sub(size(loc_max),index);
        peaks = curr_img(index);
        thresh_val = prctile(peaks,thresh_param);
        peak_thresh =  peaks(peaks>thresh_val);
        x_thresh = x(peaks>thresh_val);
        y_thresh = y(peaks>thresh_val);
        z_thresh = z(peaks>thresh_val); 
        index_thresh = index(peaks>thresh_val);
        [x_thresh_new,y_thresh_new,z_thresh_new,index_thresh_new] = remove_close_neurons(x_thresh,y_thresh,z_thresh,index_thresh,curr_img,dist_param,dist_param_z);
        figure,imshow(max(temp_med_gauss,[],3),[])
        caxis([0,0.6])
        hold on
        scatter(y_thresh,x_thresh,'.g')
        scatter(y_thresh_new,x_thresh_new,'.r')
    end