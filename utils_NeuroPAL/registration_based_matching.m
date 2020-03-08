%%%%% function to compare annotation accuracy of registration vs CRF method
%%%%% for NeuroPAL strain.
%%%%% registration here is 3D done with 3D atlas
%%%%% Methods used for registration
%%%%% 1. CPD
%%%%% 2. GLMD (Global and Local Mixture Distance based (GLMD) Non-rigid
%%%%%    Point Set Matching)
%%%%% 3. TPS-RPM


addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Joint_tracking\TPS-RPM'))
addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GLMD_Demo'))
addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\CPD2'))

data = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\5\data_20190706_5.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\21\data_20190706_21.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\2\data_20190710_2.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\5\data_20190710_5.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\8\data_20190710_8.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\10\data_20190710_10.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\22\data_20190710_22.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\24\data_20190710_24.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\27\data_20190710_27.mat'};

%%%%%%%%%%%% creat annotated neurons list
all_marker_names = {};
for i = 1:size(data,2)
    load(data{1,i})
    all_marker_names = cat(1,all_marker_names,marker_name);
end
landmark_list = unique(all_marker_names);

load('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\data_neuron_relationship.mat')
source = [X_rot,Y_rot,Z_rot];

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

match_top_matrix = NaN(size(data,2),size(landmark_list,1));
match_top_3_matrix = NaN(size(data,2),size(landmark_list,1));
match_top_5_matrix = NaN(size(data,2),size(landmark_list,1));
for i = 1:size(data,2)
    load(data{1,i})
    
    % generate axis
    mu_r_centered = mu_r - repmat(mean(mu_r),size(mu_r,1),1);
    [coeff,score,latent] = pca(mu_r_centered);
    PA = coeff(:,axes_param(1,1))';
    PA = PA/norm(PA);
    LR = coeff(:,axes_param(1,2))';
    LR = LR/norm(LR);
    DV = coeff(:,axes_param(1,3))';
    DV = DV/norm(DV);

    OLLL_neuron = marker_index(find(strcmp('OLLL',marker_name)),1);
    OLLR_neuron = marker_index(find(strcmp('OLLR',marker_name)),1);
    AVAL_neuron = marker_index(find(strcmp('AVAL',marker_name)),1);
    AVAR_neuron = marker_index(find(strcmp('AVAR',marker_name)),1);
    if ~isempty(OLLL_neuron)
        if ~isempty(AVAL_neuron)
            if (mu_r_centered(OLLL_neuron,:)-mu_r_centered(AVAL_neuron,:))*PA' < 0
                PA = -PA;
            end
        else
            if (mu_r_centered(OLLL_neuron,:)-mu_r_centered(AVAR_neuron,:))*PA' < 0
                PA = -PA;
            end
        end
    else
        if ~isempty(AVAL_neuron)
            if (mu_r_centered(OLLR_neuron,:)-mu_r_centered(AVAL_neuron,:))*PA' < 0
                PA = -PA;
            end
        else
            if (mu_r_centered(OLLR_neuron,:)-mu_r_centered(AVAR_neuron,:))*PA' < 0
                PA = -PA;
            end
        end
    end

    ALA_neuron = marker_index(find(strcmp('ALA',marker_name)),1);
    RMDDL_neuron = marker_index(find(strcmp('RMDDL',marker_name)),1);
    RMDDR_neuron = marker_index(find(strcmp('RMDDR',marker_name)),1);
    if ~isempty(RMDDL_neuron)
        if (mu_r_centered(RMDDL_neuron,:)-mu_r_centered(ALA_neuron,:))*DV' < 0
            DV = -DV;
        end
    else
        if (mu_r_centered(RMDDR_neuron,:)-mu_r_centered(ALA_neuron,:))*DV' < 0
            DV = -DV;
        end
    end
    if cross(PA,LR)*DV' < 0
        LR = -LR;
    end

    % take neurons coordinates to AP, LR, DV axis
    X = mu_r_centered*PA';
    Y = mu_r_centered*LR';
    Z = mu_r_centered*DV';
    
    target = [X,Y,Z];
    [Transform, C] = cpd_register(source, target, opt);
%     Transform = GLTP_EM(source',target',param);
    dist = repmat(diag(Transform.T*Transform.T'),1,size(Transform.X,1)) + repmat(diag(Transform.X*Transform.X')',size(Transform.T,1),1) - 2*Transform.T*Transform.X';
    [sort_dist,sort_index] = sort(dist,2,'ascend');
    for n = 1:size(landmark_list,1)
        if ~isempty(find(strcmp(landmark_list{n,1},marker_name)))
            curr_neuron = marker_index(find(strcmp(landmark_list{n,1},marker_name)),1);
            CPD_match = Neuron_head{C(curr_neuron,1),1};
            if strcmp(landmark_list{n,1},CPD_match)
                match_top_matrix(i,n) = 1;
            else
                match_top_matrix(i,n) = 0;
            end
            top_3_labels = Neuron_head(sort_index(curr_neuron,1:3),:);
            if ~isempty(find(strcmp(landmark_list{n,1},top_3_labels)))
                match_top_3_matrix(i,n) = 1;
            else
                match_top_3_matrix(i,n) = 0;
            end
            top_5_labels = Neuron_head(sort_index(curr_neuron,1:5),:);
            if ~isempty(find(strcmp(landmark_list{n,1},top_5_labels)))
                match_top_5_matrix(i,n) = 1;
            else
                match_top_5_matrix(i,n) = 0;
            end
        end
    end
end

%%%% calculate correct prediction fraction
% NaN in match matrix are landmarks not present in the data. 0 are not
% matched correctly and 1 are matched correctly.
correct_fraction = zeros(size(match_top_matrix,1),1);
correct_fraction_top_3 = zeros(size(match_top_matrix,1),1);
correct_fraction_top_5 = zeros(size(match_top_matrix,1),1);
for i = 1:size(match_top_matrix,1)
    num_landmarks = size(match_top_matrix,2) - sum(isnan(match_top_matrix(i,:)));
    correct = nansum(match_top_matrix(i,:));
    correct_fraction(i,1) = correct/num_landmarks;
    correct_top_3 = nansum(match_top_3_matrix(i,:));
    correct_fraction_top_3(i,1) = correct_top_3/num_landmarks;
    correct_top_5 = nansum(match_top_5_matrix(i,:));
    correct_fraction_top_5(i,1) = correct_top_5/num_landmarks;
end