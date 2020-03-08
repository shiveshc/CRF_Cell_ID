%%%%% function to compare annotation accuracy of registration vs CRF method
%%%%% for neuropal strain
%%%%% registration here is 2D done with 2D atlas
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

%%%% load 2D atlas
%%% run create_2D_atlas.m
source = [X_rot_2D,Y_rot_2D];

%%%%%%%%%%%%%%%%%%%%%%% define CPD parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.method='nonrigid'; % use nonrigid registration
opt.beta=1;            % the width of Gaussian kernel (smoothness)
opt.lambda=3;          % regularization weight

opt.viz=1;              % DON't show every iteration
opt.outliers=0.3;       % Noise weight
opt.fgt=0;              % do not use FGT (default)
opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
opt.corresp=1;          % compute correspondence vector at the end of registration (not being estimated by default)

opt.max_it=100;         % max number of iterations
opt.tol=1e-10;          % tolerance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

match_matrix = NaN(size(data,2),size(landmark_list,1));
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
    
    target = [X,Z];
    [Transform, C] = cpd_register(source, target, opt);
    
    %%%%% match one side only atlas
    for n = 1:size(landmark_list,1)
        if ~isempty(find(strcmp(landmark_list{n,1},marker_name)))
            curr_neuron = marker_index(find(strcmp(landmark_list{n,1},marker_name)),1);
            CPD_match = Neuron_head_2D{C(curr_neuron,1),1};
            CPD_match_char = convertStringsToChars(CPD_match);
            curr_neuron_name_char = landmark_list{n,1};
            curr_neuron_name = convertCharsToStrings(landmark_list{n,1});
            if strcmp(curr_neuron_name_char(1,end),'R') || strcmp(curr_neuron_name_char(1,end),'L')
                if strcmp(curr_neuron_name_char(1,1:end-1),CPD_match_char(1,1:end-1))
                    match_matrix(i,n) = 1;
                else
                    match_matrix(i,n) = 0;
                end
            elseif strcmp(curr_neuron_name,CPD_match)
                match_matrix(i,n) = 1;
            else
                match_matrix(i,n) = 0;
            end
        end
    end
    
%     %%%%% match both side atlas
%     for n = 1:size(landmark_list,1)
%         if ~isempty(find(strcmp(landmark_list{n,1},marker_name)))
%             curr_neuron = marker_index(find(strcmp(landmark_list{n,1},marker_name)),1);
%             CPD_match = Neuron_head_2D{C(curr_neuron,1),1};
%             if strcmp(landmark_list{n,1},CPD_match)
%                 match_matrix(i,n) = 1;
%             else
%                 match_matrix(i,n) = 0;
%             end
%         end
%     end
end

%%%% calculate correct prediction fraction
% NaN in match matrix are landmarks not present in the data. 0 are not
% matched correctly and 1 are matched correctly.
correct_fraction = zeros(size(match_matrix,1),1);
for i = 1:size(match_matrix,1)
    num_landmarks = size(match_matrix,2) - sum(isnan(match_matrix(i,:)));
    correct = nansum(match_matrix(i,:));
    correct_fraction(i,1) = correct/num_landmarks;
end