%%%% expectation maximization for segmentation
% based on Pan K. Gaussian mixtures for intensity modeling of spots in microscopy. 2010 7th IEEE International Symposium on Biomedical Imaging: From Nano to Macro, ISBI 2010 - Proceedings. 2010. pp. 121–124. doi:10.1109/ISBI.2010.5490398
% input index_thresh is index of local peaks in image

function [comp_weights,mu,sigma] = em_segmentation_3d_gpu(img,index_thresh,num_comp,mu,sigma,max_iter,change_mu,cmap)
    %%% specify x,y,z resolution, if not available set to 1
    x_res = 1;
    y_res = 1;
    z_res = 1;
    
    index = 1:1:size(img,1)*size(img,2)*size(img,3);
    [x,y,z] = ind2sub(size(img),index');
    img_data = [x,y,z];
    I_norm = img(index')/sum(sum(sum(img)));
    I_em = zeros(size(I_norm));
    
    % get data for local regions around peaks
    if ~isempty(index_thresh)
        temp_img_filt = zeros(size(img));
        for m = 1:size(index_thresh,1)
            [p,q,r] = ind2sub(size(temp_img_filt),index_thresh(m,1));
            temp_img_filt(max(p-10,1):min(p+10,size(img,1)),max(q-10,1):min(q+10,size(img,2)),max(r-2,1):min(r+2,size(img,3))) = 1;
        end
        index_filt = find(temp_img_filt);
        img_data = img_data(index_filt,:);
        if ~isempty(x_res)
            img_data = img_data.*repmat([y_res,x_res,z_res],size(img_data,1),1);
        end
        I_norm = img(index_filt)/sum(sum(sum(img(index_filt))));
    end
    size(img_data)
    if nargin < 3
        num_comp = size(index_thresh,1);
    elseif isempty(num_comp)
        num_comp = size(index_thresh,1);
    end
    if nargin < 4
        [x,y,z] = ind2sub(size(img),index_thresh);
        mu = [x,y,z].*repmat([y_res,x_res,z_res],size(x,1),1);
    elseif isempty(mu)
        [x,y,z] = ind2sub(size(img),index_thresh);
        mu = [x,y,z].*repmat([y_res,x_res,z_res],size(x,1),1);
    end
    if nargin < 5
        sigma = repmat(x_res*5*eye(size(img_data,2),size(img_data,2)),1,1,num_comp);
        sigma(3,3,:) = z_res*1;
    elseif  isempty(sigma)
        sigma = repmat(x_res*5*eye(size(img_data,2),size(img_data,2)),1,1,num_comp);
        sigma(3,3,:) = z_res*1;
    end
    if nargin < 5
        max_iter = 100;
    elseif isempty(max_iter)
        max_iter = 100;
    end
    % principle curvature values
%     [K,H,k1,k2] = curvature(img,find(ones(size(img))));
    
    % Start EM
    iter = 1;
    comp_weights = 1/num_comp*ones(1,num_comp);
    sig_img_data = gpuArray(cat(2,reshape(img_data(:,1),1,1,size(img_data,1)),reshape(img_data(:,2),1,1,size(img_data,1)),reshape(img_data(:,3),1,1,size(img_data,1))));
    [mu,sigma,comp_weights] = em_step(img,index_filt,img_data,I_norm,sig_img_data,iter,max_iter,mu,sigma,comp_weights,change_mu,cmap);
    [likelihood,old_error] = calc_likelihood(img_data,I_norm,mu,sigma,comp_weights,num_comp);
    
%     % First split and merge step
%     [sm_mu,sm_sigma,sm_comp_weights,sm_num_comp] = split_merge_em_3d(img,index_filt,img_data,I_norm,sig_img_data,k2,iter,max_iter,mu,sigma,comp_weights,1,old_error,cmap);
%     [new_likelihood,new_error] = calc_likelihood(img_data,I_norm,sm_mu,sm_sigma,sm_comp_weights,sm_num_comp);
%     if new_error < old_error
%         mu = sm_mu;
%         sigma = sm_sigma;
%         comp_weights = sm_comp_weights;
%     end
%             
%     % Iterate split and merge step
%     while old_error > new_error
%         old_error = new_error;
%         [sm_mu,sm_sigma,sm_comp_weights,sm_num_comp] = split_merge_em_3d(img,index_filt,img_data,I_norm,sig_img_data,iter,max_iter,mu,sigma,comp_weights,1,old_error,cmap);
%         [new_likelihood,new_error] = calc_likelihood(img_data,I_norm,sm_mu,sm_sigma,sm_comp_weights,sm_num_comp);
% 
%         if new_error < old_error
%             mu = sm_mu;
%             sigma = sm_sigma;
%             comp_weights = sm_comp_weights;
%         else
%             break
%         end
%     end
        
end

function [mu,sigma,comp_weights] = em_step(img,index_filt,img_data,I_norm,sig_img_data,iter,max_iter,mu,sigma,comp_weights,change_mu,cmap)
    num_comp = size(comp_weights,2);
    spatial_contribution = zeros(size(img_data,1),num_comp);
    old_error = 1000;
    new_error = 999;
    old_likelihood = -1000;
    new_likelihood = -999;
    while  abs(new_error - old_error) > 10^-3 && new_error < old_error && iter <= max_iter
        old_error = new_error;
        old_likelihood = new_likelihood;
        for k = 1:num_comp
            spatial_contribution(:,k) = comp_weights(1,k)*mvnpdf(img_data,mu(k,:),sigma(:,:,k));
        end
        P = spatial_contribution./repmat(sum(spatial_contribution,2),1,num_comp);
        contribution = P.*repmat(I_norm,1,num_comp);
        % Treat NaN value in contribution (pixels which have zero spatial
        % contribution to all components)
        contribution(isnan(contribution)) = 0;
        comp_weights = sum(contribution);
        for k = 1:num_comp
            if change_mu == 1
                mu(k,:) = sum(img_data.*repmat(contribution(:,k),1,size(mu,2)))/sum(contribution(:,k));
            end
            sig_mu_data = gpuArray(repmat(mu(k,:),1,1,size(img_data,1)));
            sig = pagefun(@mtimes, pagefun(@transpose,sig_img_data - sig_mu_data),(sig_img_data - sig_mu_data));
            sigma(:,:,k) = gather((sum(pagefun(@times,sig,repmat(reshape(contribution(:,k),1,1,size(img_data,1)),size(sigma,1),size(sigma,2),1)),3))/sum(contribution(:,k)));
        end
        
        [new_likelihood,new_error] = calc_likelihood(img_data,I_norm,mu,sigma,comp_weights,num_comp);
        
        ['iter:',num2str(iter),', likelihood:',num2str(new_likelihood),', error:',num2str(new_error)]
        iter = iter + 1;
    end
end

function [likelihood,error] = calc_likelihood(img_data,I_norm,mu,sigma,comp_weights,num_comp)
    spatial_contribution = zeros(size(img_data,1),num_comp);
    for k = 1:num_comp
        spatial_contribution(:,k) = comp_weights(1,k)*mvnpdf(img_data,mu(k,:),sigma(:,:,k));
    end
    posterior = spatial_contribution./repmat(sum(spatial_contribution,2),1,num_comp);
    posterior(find(isinf(posterior))) = 0;
    posterior(find(isnan(posterior))) = 0;
    temp = log(spatial_contribution);
    temp(find(isinf(temp))) = 0;
    temp(find(isnan(temp))) = 0;
    likelihood = sum(sum(temp.*posterior,2).*I_norm);
    for k = 1:num_comp
        spatial_contribution(:,k) = comp_weights(1,k)*mvnpdf(img_data(:,1:3),mu(k,1:3),sigma(1:3,1:3,k));
    end
    error = sum(abs(sum(spatial_contribution,2) - I_norm));
end