%%%% function to visualize multi-label prediction results for alternate
%%%% strain and make video. Example of inputs is hown below
%%%% Inputs - 
%%%% 1. results_direc
%%%% 2. input_img_dir
%%%% 3. input_data_dir
%%%% 4. data

%%%%% e.g of inputs
% results_direc = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_36r_numLand0',...
%     'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\Results_data_20190809_37r_numLand0'};
% input_img_dir = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190809_AML5_array\36\36.tif',...
% 	'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190809_AML5_array\37\37.tif'};
% input_data_dir = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190809_AML5_array\36',...
% 	'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190809_AML5_array\37'};
% data{1,1} = 'data_20190809_1r';
% data{1,2} = 'data_20190809_3r';

function visualize_multi_label_results(results_direc,input_img_dir,input_data_dir,data,d)

load('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\data_neuron_relationship_annotation_updated.mat')

%% read images and annotation data
input_1 = {input_img_dir{1,d}};
input_2 = input_data_dir{1,d};
curr_data = data{1,d};
readCyOFP_SD
% 2 channels therefore splitting full_img into 2
thisimage_g = full_img(:,:,1:size(full_img,3)/2);
thisimage_r = full_img(:,:,size(full_img,3)/2+1:end);
load([input_2,'\',curr_data])
% %% visualize labels
% neuron = 17;
% figure,imshow(max(mat2gray(thisimage_r),[],3))
% hold on
% scatter(mu_r(neuron,1),mu_r(neuron,2),'r','LineWidth', 2)
% caxis([0,0.2])
% for i = 1:size(label_struct(neuron).top_5,2)
%     text(30 ,10+20*i,[label_struct(neuron).top_5{1,i},' - ',num2str(label_struct(neuron).top_5_count(1,i))],'Color',[1,0,0],'FontSize',7)
% end

%% create annotation variable
results_dir1 = results_direc{1,d};
fileList1 = dir(results_dir1);
label_struct = struct();

% create 'geodist_r' as it is not saved in the results
K = 6;
pos = mu_r;
euc_dist = repmat(diag(pos*pos'),1,size(pos,1)) + repmat(diag(pos*pos')',size(pos,1),1) - 2*pos*pos';
[sort_euc_dist,sort_index] = sort(euc_dist,2);
adj = zeros(size(pos,1),size(pos,1));
for i = 1:size(adj,1)
    adj(i,sort_index(i,2:K+1)) = 1;
end
adj = max(adj,adj');
G = graph(adj);
geo_dist_r = distances(G);

% generate X, Y, Z from mu_r
A_neuron = axes_neurons_to_neuron_map(1,1);
P_neuron = axes_neurons_to_neuron_map(2,1);
L_neuron = axes_neurons_to_neuron_map(3,1);
R_neuron = axes_neurons_to_neuron_map(4,1);
D_neuron = axes_neurons_to_neuron_map(5,1);
V_neuron = axes_neurons_to_neuron_map(6,1);
mu_r_centered = mu_r - repmat(mean(mu_r),size(mu_r,1),1);
[coeff,score,latent] = pca(mu_r_centered);
PA = coeff(:,axes_param(1,1))';
PA = PA/norm(PA);
if ind_PCA == 1
    PA = coeff(:,axes_param(1,1))';
    PA = PA/norm(PA);
    LR = coeff(:,axes_param(1,2))';
    LR = LR/norm(LR);
    DV = coeff(:,axes_param(1,3))';
    DV = DV/norm(DV);
    
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
else
    if L_neuron ~= 0 && R_neuron ~= 0
        LR = mu_r_centered(R_neuron,:) - mu_r_centered(L_neuron,:); % LR axis based on L,R
        LR = LR/norm(LR);
        % fun = @(x)-(x(1)*coeff(1,1) + x(2)*coeff(2,1) + x(3)*coeff(3,1))^2/(x(1)^2 + x(2)^2 + x(3)^2);
        fun = @(x)-(x(1)*LR(1,1) + x(2)*LR(1,2) + LR(1,3))^2/(x(1)^2 + x(2)^2 + 1^2);
        Aeq = PA(1,1:2);
        beq = -PA(1,3);
        % x0 = coeff(:,1);
        x0 = LR(1,1:2);
        PA = fmincon(fun,x0,[],[],Aeq,beq);   % PA axis (perperndicular to LR and in direction of PC1)
        PA = [PA,1];
        PA = PA/norm(PA);

        % DV axis (perperndicular to LR and PA axis)
        A = [PA(1,1:2);LR(1,1:2)];
        b = [-PA(1,3);-LR(1,3)];
        DV = inv(A'*A)*A'*b;           
        DV = [DV',1];
        DV = DV/norm(DV);
        
        if (mu_r_centered(A_neuron,:)-mu_r_centered(P_neuron,:))*PA' < 0
            PA = -PA;
        end
        if cross(PA,LR)*DV' < 0
            DV = -DV;
        end
    end
end
% take neurons coordinates to AP, LR, DV axis
X = mu_r_centered*PA';
Y = mu_r_centered*LR';
Z = mu_r_centered*DV';

cnt = 0;
for i = 1:size(fileList1,1)
    name = fileList1(i).name;
    if ~strcmp(name,'.') && ~strcmp(name,'..')
        load([results_dir1,'\',name])
        cnt = cnt + 1;
        for n = 1:size(experiments(1).node_label(:,1),1)
            label_struct(n).label{1,i-2} = experiments(1).Neuron_head{experiments(1).node_label(n,1),1};
            [bin_pos_score, angle_score, geo_score] = consistency_score(n,experiments(1).node_label,X,Y,Z,X_rot,Y_rot,Z_rot,PA_matrix,LR_matrix,DV_matrix,geo_dist,geo_dist_r,experiments(1).Neuron_head,Neuron_head,angle_vec_atlas);
            label_struct(n).binPosScore(1,i-2) = bin_pos_score;
            label_struct(n).angleScore(1,i-2) = angle_score;
            label_struct(n).geoScore(1,i-2) = geo_score;
        end
    end
end
%%%% calculate frequency of top 5 labels assigned
for n = 1:size(label_struct,2)
    uniq_labels = unique(label_struct(n).label);
    label_struct(n).uniq_label = uniq_labels;
    count = zeros(1,size(uniq_labels,2));
    for i = 1:size(uniq_labels,2)
        count(1,i) = sum(strcmp(uniq_labels{1,i},label_struct(n).label));
    end
    label_struct(n).count = count;
    [sort_count,sort_index] = sort(count,'descend');
    label_struct(n).top_5 = {uniq_labels{1,sort_index(1,1:min(5,size(sort_index,2)))}};
    label_struct(n).top_5_count = count(1,sort_index(1,1:min(5,size(sort_index,2))));
    
    top_binPosScore = [];
    top_angleScore = [];
    top_geoScore = [];
    for t = 1:size(label_struct(n).top_5,2)
        curr_label = label_struct(n).top_5{1,t};
        curr_index = find(strcmp(curr_label,label_struct(n).label));
        top_binPosScore = [top_binPosScore,mean(label_struct(n).binPosScore(:,curr_index))];
        top_angleScore = [top_angleScore,mean(label_struct(n).angleScore(:,curr_index))];
        top_geoScore = [top_geoScore,mean(label_struct(n).geoScore(:,curr_index))];
    end
    label_struct(n).top_binPosScore = top_binPosScore;
    label_struct(n).top_angleScore = top_angleScore;
    label_struct(n).top_geoScore = top_geoScore;
    
    %%% optional - sort according to one score
    [~,sort_index] = sort(label_struct(n).top_angleScore,'descend');
    label_struct(n).top_5 = label_struct(n).top_5(:,sort_index);
    label_struct(n).top_5_count = label_struct(n).top_5_count(:,sort_index);
    label_struct(n).top_binPosScore = label_struct(n).top_binPosScore(:,sort_index);
    label_struct(n).top_angleScore = label_struct(n).top_angleScore(:,sort_index);
    label_struct(n).top_geoScore = label_struct(n).top_geoScore(:,sort_index);
end

% %%% calculate label-consistency scores for top labels
% % create 'geodist_r' as it is not saved in the results
% K = 6;
% pos = mu_r;
% euc_dist = repmat(diag(pos*pos'),1,size(pos,1)) + repmat(diag(pos*pos')',size(pos,1),1) - 2*pos*pos';
% [sort_euc_dist,sort_index] = sort(euc_dist,2);
% adj = zeros(size(pos,1),size(pos,1));
% for i = 1:size(adj,1)
%     adj(i,sort_index(i,2:K+1)) = 1;
% end
% adj = max(adj,adj');
% G = graph(adj);
% geo_dist_r = distances(G);
% 
% for n = 1:size(label_struct,2)
%     binPosScore = zeros(1,size(label_struct(i).top_5,2));
%     angleScore = zeros(1,size(label_struct(i).top_5,2));
%     geoScore = zeros(1,size(label_struct(i).top_5,2));
%     for k = 1:size(label_struct(n).top_5,2)
%         [bin_pos_score, angle_score, geo_score] = consistency_score(n,label_struct(n).top_5,mu_r(:,1),mu_r(:,3),mu_r(:,3),X_rot,Y_rot,Z_rot,PA_matrix,LR_matrix,DV_matrix,geo_dist,geo_dist_r,experiments(1).Neuron_head,Neuron_head,angle_vec_atlas);
%         binPosScore(1,k) = bin_pos_score; 
%         angleScore(1,k) = angle_score;
%         geoScore(1,k) = geo_score;
%     end
%     label_struct(k).binPosScore = binPosScore;
%     label_struct(k).angleScore = angleScore;
%     label_struct(k).geoScore = geoScore;
% end

%% create video
out_direc = input_2;
vid_name = [curr_data,'_pred_updAtlas_v2.avi'];
underscore_index = strfind(curr_data,'_');
txt_on_vid = curr_data(1,underscore_index(1,end)+1:end);
[sort_mu,sort_index] = sort(mu_r(:,1),'ascend');
writerObj = VideoWriter([out_direc,'\',vid_name]); % Name it.
writerObj.FrameRate = 2;
open(writerObj);
for n = 1:size(mu_r,1)
    curr_neuron = sort_index(n,1);
    img_r = imadjust(max(mat2gray(thisimage_r),[],3),[0,0.4],[0,1]);
    img_g = imadjust(max(mat2gray(thisimage_g),[],3),[0,0.3],[0,1]);
    img = cat(3,img_r,max(cat(3,img_r,img_g),[],3),img_r);
    img(30:31,410:424,:) = 1;               % scale bar for 5 um (0.33 um/pixel for opterra)
    imshow(img,'Border','tight')
    hold on
    scatter(mu_r(curr_neuron,1),mu_r(curr_neuron,2),'r','LineWidth', 2)
    for k = 1:size(label_struct(curr_neuron).top_5,2)
        text(30 ,200+20*k,[label_struct(curr_neuron).top_5{1,k}],'Color',[1,1,1],'FontSize',7)
%         text(80 ,200+20*k,[' - '],'Color',[1,1,1],'FontSize',7)
%         text(100 ,200+20*k,[num2str(label_struct(curr_neuron).top_5_count(1,k))],'Color',[1,1,1],'FontSize',7)
%         text(140 ,200+20*k,[num2str(round(label_struct(curr_neuron).top_binPosScore(1,k),2))],'Color',[1,1,1],'FontSize',7)
        text(100 ,200+20*k,[num2str(round(label_struct(curr_neuron).top_angleScore(1,k),2))],'Color',[1,1,1],'FontSize',7)
%         text(220 ,200+20*k,[num2str(round(label_struct(curr_neuron).top_geoScore(1,k),2))],'Color',[1,1,1],'FontSize',7)
    end
    text(200 ,10,['data 0809 AML5',' ',txt_on_vid],'Color',[0,1,1],'FontSize',7)
    set(gcf,'Color','black','InvertHardCopy','off')
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
end
close all
close(writerObj);
end

function [bin_pos_score, angle_score, geo_score] = consistency_score(i,node_label,X,Y,Z,X_rot,Y_rot,Z_rot,PA_matrix,LR_matrix,DV_matrix,geo_dist,geo_dist_r,saved_neuron_head,Neuron_head,angle_vec_atlas)
    
    PA_consistency_score = 0;
    LR_consistency_score = 0;
    DV_consistency_score = 0;
    geodist_consistency_score = 0;
    angle_consistency_score = 0;
    cnt = 0;
    i_label = find(strcmp(saved_neuron_head{node_label(i,1),1},Neuron_head));
        for j = 1:size(node_label,1)
            if i ~= j && node_label(j,1) ~= 0
                j_label = find(strcmp(saved_neuron_head{node_label(j,1),1},Neuron_head));
                if X(i,1) < X(j,1)
                    PA_consistency_score = PA_consistency_score + PA_matrix(i_label,j_label);
                else
                    PA_consistency_score = PA_consistency_score + PA_matrix(j_label,i_label);
                end

                if Y(i,1) < Y(j,1)
                    LR_consistency_score = LR_consistency_score + LR_matrix(i_label,j_label);
                else
                    LR_consistency_score = LR_consistency_score + LR_matrix(j_label,i_label);
                end

                if Z(i,1) < Z(j,1)
                    DV_consistency_score = DV_consistency_score + DV_matrix(i_label,j_label);
                else
                    DV_consistency_score = DV_consistency_score + DV_matrix(j_label,i_label);
                end

                geodist_consistency_score = geodist_consistency_score + exp(-(geo_dist(i_label,j_label)-geo_dist_r(i,j)).^2);

                angle_consistency_score = angle_consistency_score + get_angle_consistency_score(X_rot,Y_rot,Z_rot,X,Y,Z,i,j,i_label,j_label,angle_vec_atlas);
                
                cnt = cnt + 1;
            end
        end
        
    bin_pos_score = (PA_consistency_score + LR_consistency_score + DV_consistency_score)/(3*cnt);
    geo_score = geodist_consistency_score/cnt;
    angle_score = angle_consistency_score/cnt;
end

function score = get_angle_consistency_score(X_rot,Y_rot,Z_rot,X,Y,Z,i,j,i_label,j_label,angle_vec_atlas)
    
    p_prime = [X(i,1),Y(i,1),Z(i,1)];
    q_prime = [X(j,1),Y(j,1),Z(j,1)];
    angle_matrix = angle_vec_atlas*((p_prime - q_prime)'/norm(p_prime-q_prime));
    angle_matrix = reshape(angle_matrix,size(X_rot,1),size(X_rot,1));
    score = angle_matrix(i_label,j_label);
end