%%%%%% function second half of segmentation
[comp_weights,mu,sigma] = em_segmentation_3d_gpu(curr_img,index_thresh_new,size(index_thresh_new,1),[],[],5,1,cmap);      
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

cmap = rand(150,3);
cmap(1,:) = [0,0,0];
for i = 1:size(labeled_img,4)
    name = strcat([input_2{1,1},'\labeled_img_r.tif']);
    indexed_img_to_tiff(labeled_img(:,:,:,i),cmap,name)
end

mu_r = mu;
labeled_img_r = labeled_img;