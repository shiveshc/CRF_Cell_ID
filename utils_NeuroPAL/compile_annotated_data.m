%%%%%% NeuroPAL validation
%%%%% function to compile all NeuroPAL annotated data
%%%%% 'compiled_data' variable can be used to compare variation in position
%%%%% of neurons marked in NeuroPAL strains and relative positions of those neurons

function compiled_data = compile_annotated_data()
    in_direc = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\5',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\21',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\2',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\5',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\8',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\10',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\22',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\24',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\27'};

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
            X_rot(ind_RIGR(1,1),:) = [];
            Y_rot(ind_RIGR(1,1),:) = [];
            Z_rot(ind_RIGR(1,1),:) = [];
            marker_index(ind_RIGR(1,1),:) = [];
            marker_name(ind_RIGR(1,1),:) = [];
        end
        ind_RIGL = find(strcmp('RIGL',marker_name));
        if size(ind_RIGL,1) > 1
            X_rot(ind_RIGL(1,1),:) = [];
            Y_rot(ind_RIGL(1,1),:) = [];
            Z_rot(ind_RIGL(1,1),:) = [];
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
end