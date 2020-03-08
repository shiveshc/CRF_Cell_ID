%%% function for pre-processing images of alternate strains before identity
%%% labelling. If data is collected using Opeterra
%%% 
%%% Performs - 1.segmentation of RFP (pan-neuronal channel is available)
%%%            2. Reads marker_names and marker_files for the current image
%%%            3. create data that goes as input to
%%%            'annotation_CRF_alternate_strain_multi_3d.m'


%%% Inputs - 1. input_1 folder specifies location of tif images
%%%          2. input_2 folder specifies the location where data file will
%%%             be saved and locations of marker files. KEEP IT as input_2

input_1 = {'D:\PhD\Imaging\20190724_AX2164_OH10690\ZSeries-07242019-1708-326\',...
           'D:\PhD\Imaging\20190724_AX2164_OH10690\ZSeries-07242019-1708-327\'};
input_4 = {'ZSeries-07242019-1708-326_Cycle00001_Ch3_000001.ome.tif',...
           'ZSeries-07242019-1708-327_Cycle00001_Ch2_000001.ome.tif'};
input_3 = 53;

input_2 = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190724_AX2164_OH10690\Z326_327'};

%%%%%%%%%%%%%%%% read images
% addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\functions_alternateStrains'))
in_direc = input_1{1,1};
out_direc = [input_2{1,1},'\'];
img_name = input_4{1,1};
num_zplane = input_3 - 4;
readtif_Brucker_for_preprocess_data
thisimage_r = thisimage{1,1};

in_direc = input_1{1,2};
out_direc = [input_2{1,1},'\'];
img_name = input_4{1,2};
num_zplane = input_3 - 4;
readtif_Brucker_for_preprocess_data
thisimage_g = thisimage{1,1};

%%%%%%%%%%%%%%% segment red channel
%%%%% segment red-channel
thresh_param = 99.89;
dist_param = 12;
dist_param_z = 5;
segment_red_first_half
segment_red_second_half % break point here
mu_r = [mu_r(:,2),mu_r(:,1),mu_r(:,3)];

%%%%%%%%%%%%%%%%%%%% read marker_names and markers file
[X,Y,Z,marker_name,marker_index] = read_marker_files(input_2{1,1});
mu_marker = [X,Y,Z];
figure,scatter3(mu_r(:,1),mu_r(:,2),mu_r(:,3),'.r')
%%%%%%%%%%%%%%%%% correct for image artifacats
% ind_remove = find(mu_r(:,3) < 8);
% mu_r(ind_remove,:) = [];

%%%%%%%%%%%%%%%%%%%%%% create marker to neuron map
mu_r = [mu_r;mu_marker];
dist_mat = repmat(diag(mu_marker*mu_marker'),1,size(mu_r,1)) + repmat(diag(mu_r*mu_r')',size(mu_marker,1),1) - 2*mu_marker*mu_r';
[sort_dist,sort_index] = sort(dist_mat,2);
marker_to_neuron_map = sort_index(:,1);
visualize_marker_to_neuron_map(thisimage_r,thisimage_g,mu_r,marker_to_neuron_map,marker_name)

%%%%%%%%%%%%%%%%%%% define A-P, L-R, D-V neuron in red channel
axes_neurons_to_neuron_map = define_axes_specifying_neurons(thisimage_r,thisimage_g,mu_r,mu_marker,marker_to_neuron_map);

%%%%%%%%%%%%%% define axes param based on PCA
mu_r_centered = mu_r - repmat(mean(mu_r),size(mu_r,1),1);
[coeff,score,latent] = pca(mu_r_centered);
PA = coeff(:,1)';
PA = PA/norm(PA);
LR = coeff(:,3)';
LR = LR/norm(LR);
DV = coeff(:,2)';
DV = DV/norm(DV);
figure,scatter3(mu_r_centered(:,1),mu_r_centered(:,2),mu_r_centered(:,3),'.r')
hold on
plot3([0,50*PA(1,1)],[0,50*PA(1,2)],[0,50*PA(1,3)],'b','LineWidth',2.5)
plot3([0,20*LR(1,1)],[0,20*LR(1,2)],[0,20*LR(1,3)],'g','LineWidth',2.5)
plot3([0,20*DV(1,1)],[0,20*DV(1,2)],[0,20*DV(1,3)],'k','LineWidth',2.5)
axes_param = [1,3,2]; % define axes_param here
ind_PCA = 1;

%%%%%%%%%%%%%%%%% save data for full image annotation
save([input_2{1,1},'\data_20190724_326_327r'],'mu_r','marker_to_neuron_map','marker_name','marker_index','axes_neurons_to_neuron_map','axes_param','ind_PCA')

%%%%%%%%%%%%%%%%%%% define A-P, L-R, D-V neuron in marker channel
axes_neurons_to_neuron_map = define_axes_specifying_neurons([],thisimage_g,[],mu_marker,[]);
ind_PCA = 0;
save([input_2{1,1},'\data_20190724_326_327g'],'mu_marker','marker_to_neuron_map','marker_name','marker_index','axes_neurons_to_neuron_map','axes_param','ind_PCA')