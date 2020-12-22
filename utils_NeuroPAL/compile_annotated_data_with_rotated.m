%%%%%% NeuroPAL validation
%%%%% function to compile all NeuroPAL annotated data - *inculuding rotated
%%%%% worms*
%%%%% 'compiled_data' variable can be used to compare variation in position
%%%%% of neurons marked in NeuroPAL strains and relative positions of those neurons

function compiled_data = compile_annotated_data_with_rotated()

    %%%% first compile straight worm data
    in_direc = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\5',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\21',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\2',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\5',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\8',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\10',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\22',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\24',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\27'};

%     in_direc = {};

    compiled_data = struct();

    for i = 1:size(in_direc,2)
        %%% read marker files
        [X,Y,Z,marker_name,marker_index] = read_marker_files(in_direc{1,i});

        %%%% generate axes and rotated coordinates
        mu = [X,Y,Z];
        mu_centered = mu - repmat(mean(mu),size(mu,1),1);
        [coeff,score,latent] = pca(mu_centered);
        PA = coeff(:,1)';
        PA = PA/norm(PA);
        DV = coeff(:,2)';
        DV = DV/norm(DV);
        LR = coeff(:,3)';
        LR = LR/norm(LR);

        OLLL_neuron = marker_index(find(strcmp('OLLL',marker_name)),1);
        OLLR_neuron = marker_index(find(strcmp('OLLR',marker_name)),1);
        AVAL_neuron = marker_index(find(strcmp('AVAL',marker_name)),1);
        AVAR_neuron = marker_index(find(strcmp('AVAR',marker_name)),1);
        if ~isempty(OLLL_neuron)
            if ~isempty(AVAL_neuron)
                if (mu_centered(OLLL_neuron,:)-mu_centered(AVAL_neuron,:))*PA' < 0
                    PA = -PA;
                end
            else
                if (mu_centered(OLLL_neuron,:)-mu_centered(AVAR_neuron,:))*PA' < 0
                    PA = -PA;
                end
            end
        else
            if ~isempty(AVAL_neuron)
                if (mu_centered(OLLR_neuron,:)-mu_centered(AVAL_neuron,:))*PA' < 0
                    PA = -PA;
                end
            else
                if (mu_centered(OLLR_neuron,:)-mu_centered(AVAR_neuron,:))*PA' < 0
                    PA = -PA;
                end
            end
        end

        ALA_neuron = marker_index(find(strcmp('ALA',marker_name)),1);
        RMDDL_neuron = marker_index(find(strcmp('RMDDL',marker_name)),1);
        RMDDR_neuron = marker_index(find(strcmp('RMDDR',marker_name)),1);
        if ~isempty(RMDDL_neuron)
            if (mu_centered(RMDDL_neuron,:)-mu_centered(ALA_neuron,:))*DV' < 0
                DV = -DV;
            end
        else
            if (mu_centered(RMDDR_neuron,:)-mu_centered(ALA_neuron,:))*DV' < 0
                DV = -DV;
            end
        end
        if cross(PA,LR)*DV' < 0
            LR = -LR;
        end

        X_rot = mu_centered*PA';
        Y_rot = mu_centered*LR';
        Z_rot = mu_centered*DV';
        figure,scatter3(mu_centered(:,1),mu_centered(:,2),mu_centered(:,3),20,'or','filled')
        hold on
        plot3(100*[0,PA(1,1)],100*[0,PA(1,2)],100*[0,PA(1,3)],'b','LineWidth',5)
        plot3(50*[0,LR(1,1)],50*[0,LR(1,2)],30*[0,LR(1,3)],'g','LineWidth',5)
        plot3(50*[0,DV(1,1)],50*[0,DV(1,2)],50*[0,DV(1,3)],'k','LineWidth',5)
        
        %%%%keep only one RIGR neuron
        ind_RIGR = find(strcmp('RIGR',marker_name));
        if size(ind_RIGR,1) > 1
%             X_rot(ind_RIGR(1,1),:) = [];
%             Y_rot(ind_RIGR(1,1),:) = [];
%             Z_rot(ind_RIGR(1,1),:) = [];
            marker_index(ind_RIGR(1,1),:) = [];
            marker_name(ind_RIGR(1,1),:) = [];
        end
        ind_RIGL = find(strcmp('RIGL',marker_name));
        if size(ind_RIGL,1) > 1
%             X_rot(ind_RIGL(1,1),:) = [];
%             Y_rot(ind_RIGL(1,1),:) = [];
%             Z_rot(ind_RIGL(1,1),:) = [];
            marker_index(ind_RIGL(1,1),:) = [];
            marker_name(ind_RIGL(1,1),:) = [];
        end
        
        compiled_data(i).X_rot = X_rot;
        compiled_data(i).Y_rot = Y_rot;
        compiled_data(i).Z_rot = Z_rot;
        compiled_data(i).marker_index = marker_index;
        compiled_data(i).marker_name = marker_name;

%         rgb_img_annotated_data(in_direc{1,i},X,Y,marker_name,marker_index)
    end
    
    % next compile rotated worm data
    % Few rotated datasets were not used to make atlas. used are - 
    % 20190705_27
    % 20190706_16
    % 20190710_21
    % 20190720_4
    % 20190720_8
    % 20190720_9
    % 20190720_10
    % 20190720_11
    % 20190720_13
    if isempty(fieldnames(compiled_data))
        num_non_rotated_data = 0;
    else
        num_non_rotated_data = size(compiled_data,2);
    end
%     in_direc = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190705_OH15495_array\27',...
%         'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\16',...
%         'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\21',...
%         'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\4',...
%         'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\8',...
%         'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\9',...
%         'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\10',...
%         'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\11'};
    
    in_direc = {};
    
    data_direc = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190705_27_pca_ensem_int_colconstnorm.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190706_16_pca_ensem_int_colconstnorm.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190710_21_pca_ensem_int_colconstnorm.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190720_4_pca_ensem_int_colconstnorm.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190720_8_pca_ensem_int_colconstnorm.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190720_9_pca_ensem_int_colconstnorm.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190720_10_pca_ensem_int_colconstnorm.mat',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190720_11_pca_ensem_int_colconstnorm.mat'};

    for i = 1:size(in_direc,2)
        %%% read marker files
        [~,~,~,marker_name,marker_index] = read_marker_files(in_direc{1,i});
        
        %%% load data (needed for ind_PCA and axes_neurons_to_neuron_map)
        load(data_direc{1,i})

        %%%% generate axes and rotated coordinates
        mu = mu_r;
        mu_centered = mu - repmat(mean(mu),size(mu,1),1);
        [coeff,score,latent] = pca(mu_centered);
        if ind_PCA == 1
            PA = coeff(:,axes_param(1,1))';
            PA = PA/norm(PA);
            LR = coeff(:,axes_param(1,2))';
            LR = LR/norm(LR);
            DV = coeff(:,axes_param(1,3))';
            DV = DV/norm(DV);
            A_neuron = axes_neurons_to_neuron_map(1,1);P_neuron = axes_neurons_to_neuron_map(2,1);L_neuron = axes_neurons_to_neuron_map(3,1);
            R_neuron = axes_neurons_to_neuron_map(4,1);D_neuron = axes_neurons_to_neuron_map(5,1);V_neuron = axes_neurons_to_neuron_map(6,1);

            if (mu_centered(A_neuron,:)-mu_centered(P_neuron,:))*PA' < 0
                PA = -PA;
            end
            if L_neuron~=0 && R_neuron~=0
                if (mu_centered(R_neuron,:)-mu_centered(L_neuron,:))*LR' < 0
                    LR = -LR;
                end
                if cross(PA,LR)*DV' < 0
                    DV = -DV;
                end
            else
                if (mu_centered(V_neuron,:)-mu_centered(D_neuron,:))*DV' < 0
                    DV = -DV;
                end
                if cross(DV,PA)*LR' < 0
                    LR = -LR;
                end
            end
        else
            if specify_PA
                A_neuron = axes_neurons_to_neuron_map(1,1);P_neuron = axes_neurons_to_neuron_map(2,1);L_neuron = axes_neurons_to_neuron_map(3,1);
                R_neuron = axes_neurons_to_neuron_map(4,1);D_neuron = axes_neurons_to_neuron_map(5,1);V_neuron = axes_neurons_to_neuron_map(6,1);
                PA = coeff(:,1)';
                PA = PA/norm(PA);
                if L_neuron ~= 0 && R_neuron ~= 0
                    LR = mu_centered(R_neuron,:) - mu_centered(L_neuron,:); % LR axis based on L,R
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

                    if (mu_centered(A_neuron,:)-mu_centered(P_neuron,:))*PA' < 0
                        PA = -PA;
                    end
                    if (mu_centered(R_neuron,:)-mu_centered(L_neuron,:))*LR' < 0
                        LR = -LR;
                    end
                    if cross(PA,LR)*DV' < 0
                        DV = -DV;
                    end
                end
            else
                A_neuron = axes_neurons_to_neuron_map(1,1);P_neuron = axes_neurons_to_neuron_map(2,1);L_neuron = axes_neurons_to_neuron_map(3,1);
                R_neuron = axes_neurons_to_neuron_map(4,1);D_neuron = axes_neurons_to_neuron_map(5,1);V_neuron = axes_neurons_to_neuron_map(6,1);
                LR = mu_centered(R_neuron,:) - mu_centered(L_neuron,:); % LR axis based on L,R
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

                    if (mu_centered(A_neuron,:)-mu_centered(P_neuron,:))*PA' < 0
                        PA = -PA;
                    end
                    if cross(PA,LR)*DV' < 0
                        DV = -DV;
                    end
                end
            end
        end

        X_rot = mu_centered*PA';
        Y_rot = mu_centered*LR';
        Z_rot = mu_centered*DV';
        figure,scatter3(mu_centered(:,1),mu_centered(:,2),mu_centered(:,3),20,'or','filled')
        hold on
        plot3(100*[0,PA(1,1)],100*[0,PA(1,2)],100*[0,PA(1,3)],'b','LineWidth',5)
        plot3(50*[0,LR(1,1)],50*[0,LR(1,2)],30*[0,LR(1,3)],'g','LineWidth',5)
        plot3(50*[0,DV(1,1)],50*[0,DV(1,2)],50*[0,DV(1,3)],'k','LineWidth',5)
        
        %%%%keep only one RIGR neuron
        ind_RIGR = find(strcmp('RIGR',marker_name));
        if size(ind_RIGR,1) > 1
%             X_rot(ind_RIGR(1,1),:) = [];
%             Y_rot(ind_RIGR(1,1),:) = [];
%             Z_rot(ind_RIGR(1,1),:) = [];
            marker_index(ind_RIGR(1,1),:) = [];
            marker_name(ind_RIGR(1,1),:) = [];
        end
        ind_RIGL = find(strcmp('RIGL',marker_name));
        if size(ind_RIGL,1) > 1
%             X_rot(ind_RIGL(1,1),:) = [];
%             Y_rot(ind_RIGL(1,1),:) = [];
%             Z_rot(ind_RIGL(1,1),:) = [];
            marker_index(ind_RIGL(1,1),:) = [];
            marker_name(ind_RIGL(1,1),:) = [];
        end
               
        compiled_data(num_non_rotated_data+i).X_rot = X_rot;
        compiled_data(num_non_rotated_data+i).Y_rot = Y_rot;
        compiled_data(num_non_rotated_data+i).Z_rot = Z_rot;
        compiled_data(num_non_rotated_data+i).marker_index = marker_index;
        compiled_data(num_non_rotated_data+i).marker_name = marker_name;
    end
end