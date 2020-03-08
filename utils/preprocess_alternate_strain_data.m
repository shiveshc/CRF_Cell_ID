%%% function for pre-processing images of alternate strains before identity
%%% labelling.
%%% 
%%% Performs - 1.segmentation of RFP (pan-neuronal channel is available)
%%%            2. Reads marker_names and marker_files for the current image
%%%            3. create data that goes as input to
%%%            'annotation_CRF_alternate_strain_multi_3d.m'


%%% Inputs - 1. input_1 folder specifies location of tif images
%%%          2. input_2 folder specifies the location where data file will
%%%             be saved and locations of marker files

input_1 = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190809_AML5_array\29\29.tif'};
input_2 = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190809_AML5_array\29'};

%%%%%%%%%%%%%%%% read images
% addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\functions_alternateStrains'))
% readCyOFP_SD
% % 2 channels therefore splitting full_img into 2
% thisimage_g = full_img(:,:,1:size(full_img,3)/2);
% thisimage_r = full_img(:,:,size(full_img,3)/2+1:end);
% 
% %%%%%%%%%%%%%%% segment red channel
% %%%%% segment red-channel
% thresh_param = 99.83;
% dist_param = 16;
% dist_param_z = 4;
% segment_red_first_half
% segment_red_second_half % break point here
% mu_r = [mu_r(:,2),mu_r(:,1),mu_r(:,3)];
% 
% %%%%%%%%%%%%%%%%%%%% read marker_names and markers file
[X,Y,Z,marker_name,marker_index] = read_marker_files(input_2{1,1});
mu_marker = [X,Y,Z];
figure,scatter3(mu_r(:,1),mu_r(:,2),mu_r(:,3),'.r')
% %%%%%%%%%%%%%%%%% correct for image artifacats
% ind_remove = find(mu_r(:,3) < 8);
% mu_r(ind_remove,:) = [];

%%%%%%%%%%%%%%%%%%%%%%% create marker to neuron map
dist_mat = repmat(diag(mu_marker*mu_marker'),1,size(mu_r,1)) + repmat(diag(mu_r*mu_r')',size(mu_marker,1),1) - 2*mu_marker*mu_r';
[sort_dist,sort_index] = sort(dist_mat,2);
ind_change = find(sort_dist(:,1) > 64);
mu_r = [mu_r;mu_marker(ind_change,:)];
dist_mat = repmat(diag(mu_marker*mu_marker'),1,size(mu_r,1)) + repmat(diag(mu_r*mu_r')',size(mu_marker,1),1) - 2*mu_marker*mu_r';
[sort_dist,sort_index] = sort(dist_mat,2);
marker_to_neuron_map = sort_index(:,1);
visualize_marker_to_neuron_map(thisimage_r,thisimage_g,mu_r,marker_to_neuron_map,marker_name)

%%%%%%%%%%%%%%%%%%% define A-P, L-R, D-V neuron in red channel
axes_neurons_to_neuron_map = define_axes_specifying_neurons(thisimage_r,thisimage_g,mu_r,mu_marker,marker_to_neuron_map);

%%%%%%%%%%%%%% define axes param based on PCA
mu_r_centered = mu_marker - repmat(mean(mu_marker),size(mu_marker,1),1);
mu_r_centered = mu_r - repmat(mean(mu_r),size(mu_r,1),1);
[coeff,score,latent] = pca(mu_r_centered);
PA = coeff(:,1)';
PA = PA/norm(PA);
LR = coeff(:,3)';
LR = LR/norm(LR);
DV = coeff(:,2)';
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
figure,imshow(mat2gray(max(thisimage_r,[],3)))
caxis([0,0.4])
hold on
line([mean(mu_r(:,1)),mean(mu_r(:,1))+100],[mean(mu_r(:,2)),mean(mu_r(:,2))+100*PA(1,2)/PA(1,1)],'Color','c','LineWidth',3.5)
line([mean(mu_r(:,1)),mean(mu_r(:,1))+50],[mean(mu_r(:,2)),mean(mu_r(:,2))+50*DV(1,2)/DV(1,1)],'Color','r','LineWidth',3.5)
figure,scatter3(mu_r_centered(:,1),mu_r_centered(:,2),mu_r_centered(:,3),'.r')
hold on
plot3([0,50*PA(1,1)],[0,50*PA(1,2)],[0,50*PA(1,3)],'b','LineWidth',2.5)
plot3([0,20*LR(1,1)],[0,20*LR(1,2)],[0,20*LR(1,3)],'g','LineWidth',2.5)
plot3([0,20*DV(1,1)],[0,20*DV(1,2)],[0,20*DV(1,3)],'k','LineWidth',2.5)
axes_param = [1,3,2]; % define axes_param here
ind_PCA = 1;

%%%%%%%%%%%%%%%%% save data for full image annotation
save([input_2{1,1},'\data_20190809_29r_temp'],'mu_r','marker_to_neuron_map','marker_name','marker_index','axes_neurons_to_neuron_map','axes_param','ind_PCA')

%%%%%%%%%%%%%%%%%%% define A-P, L-R, D-V neuron in marker channel
axes_neurons_to_neuron_map = define_axes_specifying_neurons([],thisimage_g,[],mu_marker,[]);
ind_PCA = 1;
save([input_2{1,1},'\data_20190809_67g'],'mu_marker','marker_name','marker_index','axes_neurons_to_neuron_map','axes_param','ind_PCA')