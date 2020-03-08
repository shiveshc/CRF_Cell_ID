%%%%% function to compare annotation accuracy of registration vs CRF method
%%%%% for alternate strains.
%%%%% registration here is 3D done with 3D atlas
%%%%% Methods used for registration
%%%%% 1. CPD
%%%%% 2. GLMD (Global and Local Mixture Distance based (GLMD) Non-rigid
%%%%%    Point Set Matching)
%%%%% 3. TPS-RPM


addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Joint_tracking\TPS-RPM'))
addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GLMD_Demo'))
addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\CPD2'))

data = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_1g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_3g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_4g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_5g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_6g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_7g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_8g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_9g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_10g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_11g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_14g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_19g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_20g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_22g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_23g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_24g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_28g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_29g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_33g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_34g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_36g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_37g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_38g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_39g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_50g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_52g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_54g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_55g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_57g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_58g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_59g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_60g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_62g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_65g.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos\data_20190809_67g.mat'};

%%%%%%%%%%%% creat annotated neurons list
all_marker_names = {};
for i = 1:size(data,2)
    load(data{1,i})
    all_marker_names = cat(1,all_marker_names,marker_name);
end
landmark_list = unique(all_marker_names);


%%%%%%%%%%%%%%%%%% define GLTP_EM parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
param.omega = 0.1; param.K = 5; param.alpha = 3; param.beta = 2; param.lambda = 3; param.all_feat = 1; param.plot_iter = 1; param.norm = 1; param.denorm = 0;

%%%%%%%%%%%%%%%%%%%%%%% define CPD parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.method='nonrigid'; % use nonrigid registration
opt.beta=1;            % the width of Gaussian kernel (smoothness)
opt.lambda=3;          % regularization weight

opt.viz=0;              % DON't show every iteration
opt.outliers=0.3;       % Noise weight
opt.fgt=0;              % do not use FGT (default)
opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
opt.corresp=1;          % compute correspondence vector at the end of registration (not being estimated by default)

opt.max_it=100;         % max number of iterations
opt.tol=1e-10;          % tolerance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

match_matrix = NaN(size(data,2),size(landmark_list,1));
for i = 1:size(data,2)
%     load('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\data_neuron_relationship_2D_3D.mat')
    load('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\data_neuron_relationship.mat')
    source = [X_rot,Y_rot,Z_rot];
    
    load(data{1,i})
    
    % generate axis
%     mu_r_centered = mu_r - repmat(mean(mu_r),size(mu_r,1),1);
    mu_r_centered = mu_marker - repmat(mean(mu_marker),size(mu_marker,1),1);
    [coeff,score,latent] = pca(mu_r_centered);
    PA = coeff(:,axes_param(1,1))';
    PA = PA/norm(PA);
    LR = coeff(:,axes_param(1,2))';
    LR = LR/norm(LR);
    DV = coeff(:,axes_param(1,3))';
    DV = DV/norm(DV);
    
    A_neuron = axes_neurons_to_neuron_map(1,1);P_neuron = axes_neurons_to_neuron_map(2,1);
    L_neuron = axes_neurons_to_neuron_map(3,1);R_neuron = axes_neurons_to_neuron_map(4,1);
    D_neuron = axes_neurons_to_neuron_map(5,1);V_neuron = axes_neurons_to_neuron_map(6,1);
    
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

    % take neurons coordinates to AP, LR, DV axis
    X = mu_r_centered*PA';
    Y = mu_r_centered*LR';
    Z = mu_r_centered*DV';
    
    %%% remove all labels that are not present in the annotated neurons. This
    %%% is the same case matching to same number of labels in atlas simulations
%     keep_neurons = {'BAGL';'BAGR';'URXL';'URXR';'ASEL';'ASER';'AFDL';'AFDR';'AQR'}; %AX2164
    keep_neurons = {'IL2R';'IL1R';'IL1VR';'IL1VL';'IL1DR';'IL1DL';'IL2L';'IL1L';'AIBR';'AIBL';'AIZR';'AIZL';'ASGR';'ASGL';'RIVR';'RIVL';'RIGR';'AVG'}; % AML5
    ind_remove = 1:1:size(Neuron_head,1);
    ind_keep = [];
    for k = 1:size(marker_name,1)
        ind_keep = [ind_keep;find(strcmp(marker_name{k,1},Neuron_head))];
    end
    ind_remove(:,ind_keep) = [];
    Neuron_head(ind_remove,:) = [];
    DV_matrix(ind_remove,:) = [];
    DV_matrix(:,ind_remove) = [];
    geo_dist(ind_remove,:) = [];
    geo_dist(:,ind_remove) = [];
    LR_matrix(ind_remove,:) = [];
    LR_matrix(:,ind_remove) = [];
    PA_matrix(ind_remove,:) = [];
    PA_matrix(:,ind_remove) = [];
    ganglion(ind_remove,:) = [];
    X_rot(ind_remove,:) = [];
    Y_rot(ind_remove,:) = [];
    Z_rot(ind_remove,:) = [];
    X_rot_norm(ind_remove,:) = [];
    
    % simulate missing neurons (extrachromosomal arrays)
    rng(rng_missing_neurons(1,f))
    index_neurons_missing = randperm(size(mu_marker,1),num_missing_neurons(1,d));
    X(index_neurons_missing,:) = [];
    Y(index_neurons_missing,:) = [];
    Z(index_neurons_missing,:) = [];
    marker_index(index_neurons_missing,:) = [];
    marker_name(index_neurons_missing,:) = [];
    target = [X,Y,Z];
    
    source = [X_rot,Y_rot,Z_rot];
    [Transform, C] = cpd_register(source, target, opt);
%     Transform = GLTP_EM(source',target',param);
    for n = 1:size(landmark_list,1)
        if ~isempty(find(strcmp(landmark_list{n,1},marker_name)))
            curr_neuron = find(strcmp(landmark_list{n,1},marker_name));
            CPD_match = Neuron_head{C(curr_neuron,1),1};
            if strcmp(landmark_list{n,1},CPD_match)
                match_matrix(i,n) = 1;
            else
                match_matrix(i,n) = 0;
            end
        end
    end
end

%%%% calculate correct prediction fraction
% NaN in match matrix are landmarks not present in the data. 0 are not
% matched correctly and 1 are matched correctly.
% match_matrix(:,[7,17]) = [];
correct_fraction = zeros(size(match_matrix,1),1);
for i = 1:size(match_matrix,1)
    num_landmarks = size(match_matrix,2) - sum(isnan(match_matrix(i,:)));
    correct = nansum(match_matrix(i,:));
    correct_fraction(i,1) = correct/num_landmarks;
end