%%%% function to generate node potentials based on registration of data to
%%%% atlas

function node_pot = reg_based_node_potential(X,Y,Z,X_rot,Y_rot,Z_rot)

% addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\CPD2'))
% addpath(genpath('/gpfs/pace1/project/pchbe2/schaudhary9/Annotation/NeuroPAL/RelPos_Reg/CPD2'))
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

source = [X_rot,Y_rot,Z_rot];
target = [X,Y,Z];
[Transform, C] = cpd_register(source, target, opt);
reg_source = [Transform.X(:,1),Transform.X(:,2),Transform.X(:,3)];
reg_target = [Transform.T(:,1),Transform.T(:,2),Transform.T(:,3)];

reg_dist = repmat(diag(reg_target*reg_target'),1,size(reg_source,1)) + repmat(diag(reg_source*reg_source')',size(reg_target,1),1) - 2*reg_target*reg_source';
node_pot = 1*exp(-reg_dist/(2*0.1));
node_pot(find(node_pot<0.001)) = 0.001; %  small potential of incompatible matches

end