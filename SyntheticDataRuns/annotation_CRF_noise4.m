%%%% function for automatic annotation of one stack
%%%% Changes - 
%%%% 1. Geodesic distance based edge potentials
%%%% 2. log linear hard edge potentials
%%%% 3. Handle duplicates
%%%%    3a. Reassign all duplicate nodes
%%%%    3b. Form graph structure of all nodes, change potential so that
%%%%    unassigned nodes can be assigned unassigned labels, clamp potential
%%%%    for assigned nodes
%%%% 4. Hide landmarks and check their predicted identities
%%%% 5. Node potential based on normalized distance along PA
%%%% 6. Relative angle based edge-potentials

function annotation_CRF_noise3(numLandmarks,num)
rng shuffle
%%% load neuron relationship data (Atlas)
load('/gpfs/pace1/project/pchbe2/schaudhary9/Annotation/Atlas/data_neuron_relationship.mat')
% load('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\PACE\AtlasTest\data_neuron_relationship.mat')

%%% generate random data from atlas by adding random perturbations to
%%% neuron locations
mu_r = [X_rot,Y_rot,Z_rot];
distances_AP = sqrt(repmat(diag(X_rot*X_rot'),1,size(X_rot,1)) + repmat(diag(X_rot*X_rot')',size(X_rot,1),1) - 2*X_rot*X_rot');
distances_LR = sqrt(repmat(diag(Y_rot*Y_rot'),1,size(Y_rot,1)) + repmat(diag(Y_rot*Y_rot')',size(Y_rot,1),1) - 2*Y_rot*Y_rot');
distances_DV = sqrt(repmat(diag(Z_rot*Z_rot'),1,size(Z_rot,1)) + repmat(diag(Z_rot*Z_rot')',size(Z_rot,1),1) - 2*Z_rot*Z_rot');
% figure,hist(distances_AP) % check distances to decide perturbation sigma
% figure,hist(distances_LR)
% figure,hist(distances_DV)
mu_r_perturbed = mu_r + [0.05*randn(size(mu_r,1),1), 0.021*randn(size(mu_r,1),1), 0.020*randn(size(mu_r,1),1)];

%%% randomly select neurons (to mimic real image scenario)
keep_neurons_from_atlas = randperm(size(mu_r,1),130)';

% take neurons coordinates to AP, LR, DV axis
X = mu_r_perturbed(keep_neurons_from_atlas,1);
Y = mu_r_perturbed(keep_neurons_from_atlas,2);
Z = mu_r_perturbed(keep_neurons_from_atlas,3);
X_norm = (X-min(X))/(max(X)-min(X));

addpath(genpath('/gpfs/pace1/project/pchbe2/schaudhary9/Annotation/Atlas/UGM'))
% addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Joint_tracking\UGM'))
% create spatial neighborhood
K = 6;
pos = [mu_r_perturbed(keep_neurons_from_atlas,1),mu_r_perturbed(keep_neurons_from_atlas,2),mu_r_perturbed(keep_neurons_from_atlas,3)];
euc_dist = repmat(diag(pos*pos'),1,size(pos,1)) + repmat(diag(pos*pos')',size(pos,1),1) - 2*pos*pos';
[sort_euc_dist,sort_index] = sort(euc_dist,2);
adj = zeros(size(X,1),size(X,1));
for i = 1:size(adj,1)
    adj(i,sort_index(i,2:K+1)) = 1;
end
adj = max(adj,adj');
G = graph(adj);
geo_dist_r = distances(G);
[sOut,tOut] = findedge(G);
% figure,scatter3(X(:,1),Y(:,1),Z(:,1),30,'.r')
% hold on
% for i = 1:size(sOut,1)
%     plot3([X(sOut(i,1));X(tOut(i,1))],[Y(sOut(i,1));Y(tOut(i,1))],[Z(sOut(i,1));Z(tOut(i,1))],'k')
% end
% adj_lfs = exp(-geo_dist_r.^2/(2*max(max(geo_dist_r))));
% D = diag(sum(adj_lfs,2));
% L = D - adj_lfs;
% [eigvec_r,eigval_r] = eig(L);

% create node and edge potential
adj = ones(size(X,1),size(X,1)); % fully connected graph structure of CRF
adj = adj - diag(diag(adj));
nStates = size(Neuron_head,1);
nNodes = size(X,1);
edgeStruct = UGM_makeEdgeStruct(adj,nStates);

% node_pot =  ones(nNodes,nStates);
loc_sigma = 0.2;
node_pot = zeros(nNodes,nStates);
for i = 1:nNodes
    node_pot(i,:) = diag(exp(-((ones(size(X_rot_norm))*X_norm(i,1) - X_rot_norm)*(ones(size(X_rot_norm))*X_norm(i,1) - X_rot_norm)')/(2*loc_sigma^2)))';
end

% VC_neurons_index = find(VC_neurons(:,1) == 1);
% VC_ind_matrix = zeros(nStates,nStates);
% VC_ind_matrix(VC_neurons_index,VC_neurons_index) = 1;
% edge_pot = zeros(nStates,nStates,edgeStruct.nEdges);
% 
% dist_PA = repmat(diag(pos(:,1)*pos(:,1)'),1,size(pos(:,1),1)) + repmat(diag(pos(:,1)*pos(:,1)')',size(pos(:,1),1),1) - 2*pos(:,1)*pos(:,1)';
% dist_PA(find(dist_PA<0)) = 0;
% dist_PA = sqrt(dist_PA);
% dist_PA = dist_PA/(0.1); % dist_PA/4 = 5 => dist_PA = 20 => 20 pixels radius will reach 0.99 value
% 
% dist_LR = repmat(diag(pos(:,2)*pos(:,2)'),1,size(pos(:,2),1)) + repmat(diag(pos(:,2)*pos(:,2)')',size(pos(:,2),1),1) - 2*pos(:,2)*pos(:,2)';
% dist_LR(find(dist_LR<0)) = 0;
% dist_LR = sqrt(dist_LR);
% dist_LR = dist_LR/(0.1); % dist_LR/2 = 5 => dist_LR = 10 => 10 pixels radius will reach 0.99 value
% 
% 
% dist_DV = repmat(diag(pos(:,3)*pos(:,3)'),1,size(pos(:,3),1)) + repmat(diag(pos(:,3)*pos(:,3)')',size(pos(:,3),1),1) - 2*pos(:,3)*pos(:,3)';
% dist_DV(find(dist_DV<0)) = 0;
% dist_DV = sqrt(dist_DV);
% dist_DV = dist_DV/(0.1); %(0.00001*std(dist_DV(find(triu(dist_DV,1)>0))));
% 
% PA_matrix2 = PA_matrix;
% LR_matrix2 = LR_matrix;
% DV_matrix2 = DV_matrix;
% PA_matrix2(find(PA_matrix == 0)) = -1;
% LR_matrix2(find(LR_matrix == 0)) = -1;
% DV_matrix2(find(DV_matrix == 0)) = -1;

lambda_PA = 0;
lambda_LR = 0;
lambda_DV = 0;
lambda_geo = 0;
lambda_angle = 1;
for i = 1:edgeStruct.nEdges
    node1 = edgeStruct.edgeEnds(i,1);
    node2 = edgeStruct.edgeEnds(i,2);
    angle_matrix = get_relative_angles(X_rot,Y_rot,Z_rot,X,Y,Z,node1,node2);
    if X(node1,1) < X(node2,1)
        if Y(node1,1) < Y(node2,1)
            if Z(node1,1) < Z(node2,1)
                pot = exp(lambda_PA*PA_matrix).*exp(lambda_DV*DV_matrix).*exp(lambda_LR*LR_matrix).*exp(exp(-lambda_geo*(geo_dist-geo_dist_r(node1,node2)).^2)).*exp(lambda_angle*angle_matrix);
            else
                pot = exp(lambda_PA*PA_matrix).*exp(lambda_DV*DV_matrix').*exp(lambda_LR*LR_matrix).*exp(exp(-lambda_geo*(geo_dist-geo_dist_r(node1,node2)).^2)).*exp(lambda_angle*angle_matrix);
            end
        else
            if Z(node1,1) < Z(node2,1)
                pot = exp(lambda_PA*PA_matrix).*exp(lambda_DV*DV_matrix).*exp(lambda_LR*LR_matrix').*exp(exp(-lambda_geo*(geo_dist-geo_dist_r(node1,node2)).^2)).*exp(lambda_angle*angle_matrix);
            else
                pot = exp(lambda_PA*PA_matrix).*exp(lambda_DV*DV_matrix').*exp(lambda_LR*LR_matrix').*exp(exp(-lambda_geo*(geo_dist-geo_dist_r(node1,node2)).^2)).*exp(lambda_angle*angle_matrix);
            end
        end
    else
        if Y(node1,1) < Y(node2,1)
            if Z(node1,1) < Z(node2,1)
                pot = exp(lambda_PA*PA_matrix').*exp(lambda_DV*DV_matrix).*exp(lambda_LR*LR_matrix).*exp(exp(-lambda_geo*(geo_dist-geo_dist_r(node1,node2)).^2)).*exp(lambda_angle*angle_matrix);
            else
                pot = exp(lambda_PA*PA_matrix').*exp(lambda_DV*DV_matrix').*exp(lambda_LR*LR_matrix).*exp(exp(-lambda_geo*(geo_dist-geo_dist_r(node1,node2)).^2)).*exp(lambda_angle*angle_matrix);
            end
        else
            if Z(node1,1) < Z(node2,1)
                pot = exp(lambda_PA*PA_matrix').*exp(lambda_DV*DV_matrix).*exp(lambda_LR*LR_matrix').*exp(exp(-lambda_geo*(geo_dist-geo_dist_r(node1,node2)).^2)).*exp(lambda_angle*angle_matrix);
            else
                pot = exp(lambda_PA*PA_matrix').*exp(lambda_DV*DV_matrix').*exp(lambda_LR*LR_matrix').*exp(exp(-lambda_geo*(geo_dist-geo_dist_r(node1,node2)).^2)).*exp(lambda_angle*angle_matrix);
            end
        end
    end
    pot(find(pot<0.01)) = 0.001; %  small potential of incompatible matches
    pot = pot - diag(diag(pot)) + 0.001*eye(size(pot,1)); 
    edge_pot(:,:,i) = pot;
end

% randomly select landmarks
num_landmarks = numLandmarks;
rand_selection = randperm(size(X,1),num_landmarks);
clamped = zeros(nNodes,1);
for i = 1:size(rand_selection,2)
    clamped(rand_selection(1,i),1) = keep_neurons_from_atlas(rand_selection(1,i),1);
end

[nodeBel,edgeBel,logZ] = UGM_Infer_Conditional(node_pot,edge_pot,edgeStruct,clamped,@UGM_Infer_LBP);
conserved_nodeBel = nodeBel; %node belief matrix to maintain marginal probabilities after clamping in subsequent steps
% optimal_decode = UGM_Decode_Conditional(node_pot,edge_pot,edgeStruct,clamped,@UGM_Decode_LBP);
[sort_nodeBel,nodeBel_sort_index] = sort(nodeBel,2,'descend');
curr_labels = nodeBel_sort_index(:,1);
[PA_score,LR_score,DV_score,geodist_score,tot_score] = consistency_scores(nNodes,curr_labels,X,Y,Z,PA_matrix,LR_matrix,DV_matrix,geo_dist,geo_dist_r);
landmark_match_score = compare_labels_of_hidden_landmarks(curr_labels,rand_selection,X,keep_neurons_from_atlas,conserved_nodeBel);

%%% handle duplicate assignments 
orig_state_array = [1:1:size(Neuron_head,1)]';
clamped_neurons = rand_selection';
node_label = duplicate_labels(curr_labels,X,Y,Z,PA_matrix,LR_matrix,DV_matrix,geo_dist,geo_dist_r,mu_r,Neuron_head,lambda_geo,clamped_neurons);
cnt = 2;
while find(node_label(:,1) == 0)
    assigned_nodes = find(node_label(:,1) ~= 0);
    assigned_labels = node_label(node_label(:,1) ~= 0,1);
    unassigned_nodes = find(node_label(:,1) == 0);
    
%     node_pot =  ones(nNodes,nStates);
    loc_sigma = 0.2;
    node_pot = zeros(nNodes,nStates);
    for i = 1:nNodes
        node_pot(i,:) = diag(exp(-((ones(size(X_rot_norm))*X_norm(i,1) - X_rot_norm)*(ones(size(X_rot_norm))*X_norm(i,1) - X_rot_norm)')/(2*loc_sigma^2)))';
    end
    node_pot(unassigned_nodes,assigned_labels) = 0;
    node_pot(find(node_pot<0.01)) = 0.001;
    edge_pot = zeros(nStates,nStates,edgeStruct.nEdges);
    for i = 1:size(edgeStruct.edgeEnds,1)
        node1 = edgeStruct.edgeEnds(i,1);
        node2 = edgeStruct.edgeEnds(i,2);
        angle_matrix = get_relative_angles(X_rot,Y_rot,Z_rot,X,Y,Z,node1,node2);
        if X(node1,1) < X(node2,1)
            if Y(node1,1) < Y(node2,1)
                if Z(node1,1) < Z(node2,1)
                    pot = exp(lambda_PA*PA_matrix).*exp(lambda_DV*DV_matrix).*exp(lambda_LR*LR_matrix).*exp(exp(-lambda_geo*(geo_dist-geo_dist_r(node1,node2)).^2)).*exp(lambda_angle*angle_matrix);
                else
                    pot = exp(lambda_PA*PA_matrix).*exp(lambda_DV*DV_matrix').*exp(lambda_LR*LR_matrix).*exp(exp(-lambda_geo*(geo_dist-geo_dist_r(node1,node2)).^2)).*exp(lambda_angle*angle_matrix);
                end
            else
                if Z(node1,1) < Z(node2,1)
                    pot = exp(lambda_PA*PA_matrix).*exp(lambda_DV*DV_matrix).*exp(lambda_LR*LR_matrix').*exp(exp(-lambda_geo*(geo_dist-geo_dist_r(node1,node2)).^2)).*exp(lambda_angle*angle_matrix);
                else
                    pot = exp(lambda_PA*PA_matrix).*exp(lambda_DV*DV_matrix').*exp(lambda_LR*LR_matrix').*exp(exp(-lambda_geo*(geo_dist-geo_dist_r(node1,node2)).^2)).*exp(lambda_angle*angle_matrix);
                end
            end
        else
            if Y(node1,1) < Y(node2,1)
                if Z(node1,1) < Z(node2,1)
                    pot = exp(lambda_PA*PA_matrix').*exp(lambda_DV*DV_matrix).*exp(lambda_LR*LR_matrix).*exp(exp(-lambda_geo*(geo_dist-geo_dist_r(node1,node2)).^2)).*exp(lambda_angle*angle_matrix);
                else
                    pot = exp(lambda_PA*PA_matrix').*exp(lambda_DV*DV_matrix').*exp(lambda_LR*LR_matrix).*exp(exp(-lambda_geo*(geo_dist-geo_dist_r(node1,node2)).^2)).*exp(lambda_angle*angle_matrix);
                end
            else
                if Z(node1,1) < Z(node2,1)
                    pot = exp(lambda_PA*PA_matrix').*exp(lambda_DV*DV_matrix).*exp(lambda_LR*LR_matrix').*exp(exp(-lambda_geo*(geo_dist-geo_dist_r(node1,node2)).^2)).*exp(lambda_angle*angle_matrix);
                else
                    pot = exp(lambda_PA*PA_matrix').*exp(lambda_DV*DV_matrix').*exp(lambda_LR*LR_matrix').*exp(exp(-lambda_geo*(geo_dist-geo_dist_r(node1,node2)).^2)).*exp(lambda_angle*angle_matrix);
                end
            end
        end
        
        if node_label(node1,1) == 0 && node_label(node2,1) == 0 % unassigned-unassigned nodes
            pot(assigned_labels,assigned_labels) = 0;
        elseif node_label(node1,1) == 0 && node_label(node2,1) ~= 0 % unassigned-assigned nodes
            pot(assigned_labels,:) = 0;
        elseif node_label(node1,1) ~= 0 && node_label(node2,1) == 0 % assigned-unassigned nodes
            pot(:,assigned_labels) = 0;
        else
        end 
        pot(find(pot<0.01)) = 0.001; %  small potential of incompatible matches
        pot = pot - diag(diag(pot)) + 0.001*eye(size(pot,1));
        edge_pot(:,:,i) = pot;
    end
    
    clamped = zeros(nNodes,1);
    clamped(assigned_nodes) = assigned_labels;
    
    [nodeBel,edgeBel,logZ] = UGM_Infer_Conditional(node_pot,edge_pot,edgeStruct,clamped,@UGM_Infer_LBP);
    conserved_nodeBel(unassigned_nodes,:) = nodeBel(unassigned_nodes,:);
    [sort_nodeBel,nodeBel_sort_index] = sort(nodeBel,2,'descend');
    
    curr_labels = nodeBel_sort_index(:,1);
    [PAscore,LRscore,DVscore,geodistscore,totscore] = consistency_scores(nNodes,curr_labels,X,Y,Z,PA_matrix,LR_matrix,DV_matrix,geo_dist,geo_dist_r);
    PA_score(:,cnt) = PAscore;
    LR_score(:,cnt) = LRscore;
    DV_score(:,cnt) = DVscore;
    geodist_score(:,cnt) = geodistscore;
    tot_score(:,cnt) = totscore;
    landmarkMatchScore = compare_labels_of_hidden_landmarks(curr_labels,rand_selection,X,keep_neurons_from_atlas,conserved_nodeBel);
    landmark_match_score(:,cnt) = landmarkMatchScore;
    
    node_label = duplicate_labels(curr_labels,X,Y,Z,PA_matrix,LR_matrix,DV_matrix,geo_dist,geo_dist_r,mu_r,Neuron_head,lambda_geo,clamped_neurons);
    cnt = cnt + 1;
end
%%% save experiments results
experiments_file = ['/gpfs/pace1/project/pchbe2/schaudhary9/Annotation/Atlas/Results/Results_130rand_206label_randLandmarks/RP_noise3/experiments_atlas_',num2str(numLandmarks),'_',num2str(num),'.mat'];
% experiments_file = ['C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\PACE\AtlasTest\experiments_atlas_',num2str(numLandmarks),'.mat'];
if exist(experiments_file)
    load(experiments_file)
    num_exp = size(experiments,2);
    
    experiments(num_exp+1).K = K;
    experiments(num_exp+1).lambda_PA = lambda_PA;
    experiments(num_exp+1).lambda_LR = lambda_LR;
    experiments(num_exp+1).lambda_DV = lambda_DV;
    experiments(num_exp+1).lambda_geo = lambda_geo;
    experiments(num_exp+1).PA_score = PA_score;
    experiments(num_exp+1).LR_score = LR_score;
    experiments(num_exp+1).DV_score = DV_score;
    experiments(num_exp+1).geodist_score = geodist_score;
    experiments(num_exp+1).tot_score = tot_score;
    experiments(num_exp+1).landmark_match_score = landmark_match_score;
    experiments(num_exp+1).num_landmarks = num_landmarks;
    experiments(num_exp+1).loc_sigma = loc_sigma;
    experiments(num_exp+1).lambda_angle = lambda_angle;
    experiments(num_exp+1).node_label = node_label;
%     experiments(num_exp+1).mu_r = mu_r;
%     experiments(num_exp+1).Neuron_head = Neuron_head;
    experiments(num_exp+1).landmarks_used = rand_selection;
    experiments(num_exp+1).keep_neurons_from_atlas = keep_neurons_from_atlas;
%     experiments(num_exp+1).landmark_names = landmark_names;
%     experiments(num_exp+1).landmark_to_neuron_map = landmark_to_neuron_map;
    save(experiments_file,'experiments')
else
    experiments = struct();
    experiments(1).K = K;
    experiments(1).lambda_PA = lambda_PA;
    experiments(1).lambda_LR = lambda_LR;
    experiments(1).lambda_DV = lambda_DV;
    experiments(1).lambda_geo = lambda_geo;
    experiments(1).PA_score = PA_score;
    experiments(1).LR_score = LR_score;
    experiments(1).DV_score = DV_score;
    experiments(1).geodist_score = geodist_score;
    experiments(1).tot_score = tot_score;
    experiments(1).landmark_match_score = landmark_match_score;
    experiments(1).num_landmarks = num_landmarks;
    experiments(1).loc_sigma = loc_sigma;
    experiments(1).lambda_angle = lambda_angle;
    experiments(1).node_label = node_label;
%     experiments(1).mu_r = mu_r;
%     experiments(1).Neuron_head = Neuron_head;
    experiments(1).landmarks_used = rand_selection;
    experiments(1).keep_neurons_from_atlas = keep_neurons_from_atlas;
%     experiments(1).landmark_names = landmark_names;
%     experiments(1).landmark_to_neuron_map = landmark_to_neuron_map;
    save(experiments_file,'experiments')
end