%%%% function to visualize multi-label prediction results for alternate
%%%% strain and make video. Example of inputs is hown below
%%%% Inputs - 
%%%% 1. results_direc
%%%% 2. input_img_dir
%%%% 3. input_data_dir
%%%% 4. data

%%%%% e.g of inputs
% input_img_dir = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190809_AML5_array\36\36.tif',...
% 	'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190809_AML5_array\37\37.tif'};
% input_data_dir = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190809_AML5_array\36',...
% 	'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190809_AML5_array\37'};
% data{1,1} = 'data_20190809_1r';
% data{1,2} = 'data_20190809_3r';

function visualize_label_results_marker_only(input_img_dir,input_data_dir,data,d,landmark_list,match_matrix,match_name,match_score)

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

%% define prediction results
%%%% for registration based matching, prediction results for any data are
%%%% defined using landmark list and match_matrix variable. These two
%%%% variables are defined by running 'registration_based_matching.m'
landmark_list = landmark_list;
match_matrix = match_matrix(31,:);
match_name = match_name(31,:);

%%%% for Rel.Pos - prediction results are defined using node_label and keep_neurons
%%%% these variables are defined by running CRF function
Neuron_head = match_score(31).Neuron_head;
node_label = match_score(31).node_label;

%% create video
out_direc = input_2;
vid_name = [curr_data,'_pred_marker_only.avi'];
underscore_index = strfind(curr_data,'_');
txt_on_vid = curr_data(1,underscore_index(1,end)+1:end);
[sort_mu,sort_index] = sort(mu_marker(:,1),'ascend');
writerObj = VideoWriter([out_direc,'\',vid_name]); % Name it.
writerObj.FrameRate = 2;
open(writerObj);
for n = 1:size(mu_marker,1)
    curr_neuron = sort_index(n,1);
    curr_neuron_reg = find(strcmp(marker_name{curr_neuron,1},landmark_list));
%     img_r = imadjust(max(mat2gray(thisimage_r),[],3),[0,0.4],[0,1]);
    img_r = max(zeros(size(thisimage_r)),[],3);
    img_g = imadjust(max(mat2gray(thisimage_g),[],3),[0,0.3],[0,1]);
    img = cat(3,img_r,max(cat(3,img_r,img_g),[],3),img_r);
    img(30:31,410:424,:) = 1;               % scale bar for 5 um (0.33 um/pixel for opterra)
    imshow(img,'Border','tight')
    hold on
    
%     if or(strcmp(marker_name{curr_neuron,1},'RIGR'),strcmp(marker_name{curr_neuron,1},'AVG'))
%     else
        scatter(mu_marker(curr_neuron,1),mu_marker(curr_neuron,2),'r','LineWidth', 2)
        text(30 ,50,'Predicted Labels','Color',[1,1,1],'FontSize',7)
        text(30 ,50+20,'CRF','Color',[1,1,1],'FontSize',7)
        text(90 ,50+20,'Registration','Color',[1,1,1],'FontSize',7)
        text(75 ,50+20,'|','Color',[1,1,1],'FontSize',7)
        
        if node_label(curr_neuron,1) == 0
            text(30 ,50+40,'','Color',[1,1,1],'FontSize',7)
        else
            if strcmp(Neuron_head{node_label(curr_neuron,1),1},marker_name{curr_neuron,1})
                text(30 ,50+40,Neuron_head{node_label(curr_neuron,1),1},'Color',[1,1,1],'FontSize',7)
            else
                text(30 ,50+40,Neuron_head{node_label(curr_neuron,1),1},'Color',[1,0,0],'FontSize',7)
            end
        end
        
        text(75 ,50+40,'|','Color',[1,1,1],'FontSize',7)
        
        if match_matrix(1,curr_neuron_reg) == 1 % for data 18
            text(90 ,50+40,match_name{1,curr_neuron_reg},'Color',[1,1,1],'FontSize',7)
        else
            text(90 ,50+40,match_name{1,curr_neuron_reg},'Color',[1,0,0],'FontSize',7)
        end
        set(gcf,'Color','black','InvertHardCopy','off')
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
%     end
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