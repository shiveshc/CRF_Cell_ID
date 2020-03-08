%%% function for pre-processing images of NeuroPAL strains before identity
%%% labelling.
%%% 
%%% Performs - 
%%%            1. Reads marker_names and marker_files for the current image
%%%            2. create data that goes as input to
%%%            'annotation_CRF_neuroPAL_multi_3d.m'


%%% Inputs - 1. input_1 folder specifies location of tif images
%%%          2. input_2 folder specifies the location where data file will
%%%             be saved and locations of marker files

input_1 = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\21'};
input_2 = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\21'};

%%%%%%%%%%%%%%%% read images
curr_files = dir([input_2{1,1}]);
for n = 1:size(curr_files,1)
    curr_name = curr_files(n).name;
    if ~(strcmp('.',curr_name)) && ~(strcmp('..',curr_name)) && isfolder([input_2{1,1},'\',curr_name]) && isempty(regexp(curr_name,'index','match'))
        z_plane_list = dir([input_2{1,1},'\',curr_name]);
        z_plane_names = {z_plane_list(:).name};
        rem_index = [find(strcmp(z_plane_names(:),'.')),find(strcmp(z_plane_names(:),'..'))];
        z_plane_names(rem_index) = [];
        num_z_planes = numel(z_plane_names);

        img_size = size(imread([input_2{1,1},'\',curr_name,'\',z_plane_names{1}]));

        full_img = zeros(img_size(1),img_size(2),num_z_planes,1,'uint16');

        for j = 1:num_z_planes
%                 full_img(:,:,num_z_planes - j + 1,1) = imread([in_direc{1,i},'\',curr_name,'\',z_plane_names{j}]);
            full_img(:,:,j,1) = imread([input_2{1,1},'\',curr_name,'\',z_plane_names{j}]);
        end

        if ~isempty(regexp(curr_name,'[rR]','match'))
        elseif ~isempty(regexp(curr_name,'[bB]','match')) 
            thisimage_b = full_img;
        elseif ~isempty(regexp(curr_name,'[cC]','match')) 
            thisimage_g = full_img;
        elseif ~isempty(regexp(curr_name,'[mM]','match'))
            thisimage_r = full_img;
        end
    end
end

%%%%%%%%%%%%%%%%%%%% read marker_names and markers file
% [X,Y,Z,marker_name,marker_index] = read_marker_files(input_2{1,1});
[X,Y,Z] = read_marker_files_wo_marker_name(input_2{1,1});
mu_r = [X,Y,Z];
ind_PCA = 0;
axes_param = [1,3,2];
specify_PA = 0;
%%%%%%%%%%%%%%%%%%%%%% create marker to neuron map


%%%%%%%%%%%%%%%%%%% define A-P, L-R, D-V neuron in red channel
%%%%% not needed if ind_PCA = 1.
axes_neurons_to_neuron_map = define_axes_specifying_neurons(thisimage_r,thisimage_g,thisimage_b,mu_r);
A_neuron = axes_neurons_to_neuron_map(1,1);P_neuron = axes_neurons_to_neuron_map(2,1);L_neuron = axes_neurons_to_neuron_map(3,1);
R_neuron = axes_neurons_to_neuron_map(4,1);D_neuron = axes_neurons_to_neuron_map(5,1);V_neuron = axes_neurons_to_neuron_map(6,1);

%%%%%%%%%%%%%% define axes param based on PCA
mu_r_centered = mu_r - repmat(mean(mu_r),size(mu_r,1),1);
[coeff,score,latent] = pca(mu_r_centered);
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
    if L_neuron~=0 && R_neuron~=0
        if (mu_r_centered(R_neuron,:)-mu_r_centered(L_neuron,:))*LR' < 0
            LR = -LR;
        end
        if cross(PA,LR)*DV' < 0
            DV = -DV;
        end
    else
        if (mu_r_centered(V_neuron,:)-mu_r_centered(D_neuron,:))*DV' < 0
            DV = -DV;
        end
        if cross(DV,PA)*LR' < 0
            LR = -LR;
        end
    end
    figure,scatter3(mu_r_centered(:,1),mu_r_centered(:,2),mu_r_centered(:,3),'.r')
    hold on
    plot3([0,50*PA(1,1)],[0,50*PA(1,2)],[0,50*PA(1,3)],'b','LineWidth',2.5)
    plot3([0,20*LR(1,1)],[0,20*LR(1,2)],[0,20*LR(1,3)],'g','LineWidth',2.5)
    plot3([0,20*DV(1,1)],[0,20*DV(1,2)],[0,20*DV(1,3)],'k','LineWidth',2.5)
%     axes_param = [1,3,2]; % define axes_param here
else
    if specify_PA
        PA = coeff(:,1)';
        PA = PA/norm(PA);
        if L_neuron ~= 0 && R_neuron ~= 0
            LR = mu_r_centered(R_neuron,:) - mu_r_centered(L_neuron,:); % LR axis based on L,R
            LR = LR/norm(LR);
            % fun = @(x)-(x(1)*coeff(1,1) + x(2)*coeff(2,1) + x(3)*coeff(3,1))^2/(x(1)^2 + x(2)^2 + x(3)^2);
            fun = @(x)-(x(1)*LR(1,1) + x(2)*LR(1,2) + LR(1,3))^2/(x(1)^2 + x(2)^2 + 1^2);
            Aeq = PA(1,1:2);
            beq = -PA(1,3);
            % x0 = coeff(:,1);
            x0 = LR(1,1:2);
            LR = fmincon(fun,x0,[],[],Aeq,beq);   % PA axis (perperndicular to LR and in direction of PC1)
            LR = [LR,1];
            LR = LR/norm(LR);

            % DV axis (perperndicular to LR and PA axis)
            A = [PA(1,1:2);LR(1,1:2)];
            b = [-PA(1,3);-LR(1,3)];
            DV = inv(A'*A)*A'*b;           
            DV = [DV',1];
            DV = DV/norm(DV);

            if (mu_r_centered(A_neuron,:)-mu_r_centered(P_neuron,:))*PA' < 0
                PA = -PA;
            end
            if (mu_r_centered(R_neuron,:)-mu_r_centered(L_neuron,:))*LR' < 0
                LR = -LR;
            end
            if cross(PA,LR)*DV' < 0
                DV = -DV;
            end
        end
    else
        A_neuron = axes_neurons_to_neuron_map(1,1);P_neuron = axes_neurons_to_neuron_map(2,1);L_neuron = axes_neurons_to_neuron_map(3,1);
        R_neuron = axes_neurons_to_neuron_map(4,1);D_neuron = axes_neurons_to_neuron_map(5,1);V_neuron = axes_neurons_to_neuron_map(6,1);
        LR = mu_r_centered(R_neuron,:) - mu_r_centered(L_neuron,:); % LR axis based on L,R
        LR = LR/norm(LR);
        if A_neuron ~= 0 && P_neuron ~= 0
%             PA = mu_r_centered(A_neuron,:) - mu_r_centered(P_neuron,:); % LR axis based on L,R
%             PA = PA/norm(PA);
            PA = coeff(:,1)';
            PA = PA/norm(PA);
            % fun = @(x)-(x(1)*coeff(1,1) + x(2)*coeff(2,1) + x(3)*coeff(3,1))^2/(x(1)^2 + x(2)^2 + x(3)^2);
            fun = @(x)-(x(1)*PA(1,1) + x(2)*PA(1,2) + PA(1,3))^2/(x(1)^2 + x(2)^2 + 1^2);
            Aeq = LR(1,1:2);
            beq = -LR(1,3);
            % x0 = coeff(:,1);
            x0 = PA(1,1:2);
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
    figure,scatter3(mu_r_centered(:,1),mu_r_centered(:,2),mu_r_centered(:,3),'.r')
    hold on
    plot3([0,50*PA(1,1)],[0,50*PA(1,2)],[0,50*PA(1,3)],'b','LineWidth',2.5)
    plot3([0,20*LR(1,1)],[0,20*LR(1,2)],[0,20*LR(1,3)],'g','LineWidth',2.5)
    plot3([0,20*DV(1,1)],[0,20*DV(1,2)],[0,20*DV(1,3)],'k','LineWidth',2.5)
%     axes_param = [1,3,2]; % define axes_param here
end

ind_rotate = 0;
if ind_rotate == 1
    spec_LR_vec = mu_r_centered(axes_neurons_to_neuron_map(4,1),:) - mu_r_centered(axes_neurons_to_neuron_map(3,1),:);
    %%% proj this vector in LR-DV plane by substracting PA component
    spec_LR_proj = spec_LR_vec - (spec_LR_vec*PA')*PA;
    spec_LR_proj = spec_LR_proj/(norm(spec_LR_proj));
    %%% calc angle between proj vector and LR axis
    angle = acos((spec_LR_proj*LR')/(norm(spec_LR_proj)));
    if cross(LR,spec_LR_proj)*PA' > 0
        sign_angle = 1;
    else
        sign_angle = 0;
    end
    %%% rotate LR and DV axes by this angle
    if sign_angle
        LR_new = LR*[1, cos(angle), sin(angle); 0, -sin(angle), cos(angle); 0,0,0];
        DV_new = DV*[1, cos(angle), sin(angle); 0, -sin(angle), cos(angle); 0,0,0];
    else
        LR_new = LR*[1, cos(-angle), sin(-angle); 0, -sin(-angle), cos(-angle); 0,0,0];
        DV_new = DV*[1, cos(-angle), sin(-angle); 0, -sin(-angle), cos(-angle); 0,0,0];
    end
    
    LR_new = LR_new/norm(LR_new);
    DV_new = DV_new/norm(DV_new);
    hold on
    plot3([0,20*LR_new(1,1)],[0,20*LR_new(1,2)],[0,20*LR_new(1,3)],'g--','LineWidth',2.5)
    plot3([0,20*DV_new(1,1)],[0,20*DV_new(1,2)],[0,20*DV_new(1,3)],'k--','LineWidth',2.5)
end

%%%%%%%%%%%%%%%%% get intensity distributions of neurons %%%%%%%%%%%%
data_int = [];
for an = 1:size(X,1)
    X_range = max(round(X(an,1)) - 2,1):1:min(round(X(an,1))+2,size(thisimage_r,2));
    Y_range = max(round(Y(an,1)) - 2,1):1:min(round(Y(an,1))+2,size(thisimage_r,1));
    Z_range = max(round(Z(an,1)) - 1,1):1:min(round(Z(an,1))+1,size(thisimage_r,3));
    [x,y,z] = meshgrid(Y_range,X_range,Z_range);
    pixels = sub2ind(size(thisimage_r),x(:),y(:),z(:));
    R_int = thisimage_r(pixels);
    G_int = thisimage_g(pixels);
    B_int = thisimage_b(pixels);

    data_int(an).R_int = R_int;
    data_int(an).G_int = G_int;
    data_int(an).B_int = B_int;
end

%%%%%%%%%%%%%%%%% save data for full image annotation
% save([input_2{1,1},'\data_20190710_21_v2'],'mu_r','marker_index','marker_name','axes_param','ind_PCA','data_int')
save([input_2{1,1},'\data_20190720_13_v2'],'mu_r','axes_neurons_to_neuron_map','ind_PCA','specify_PA','data_int','axes_param')
