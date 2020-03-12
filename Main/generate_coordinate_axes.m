%%% function to generate coordinate axes in images.
%%% NOTE - landmarks cells are used to disambiguate axes. This function
%%%        works for the datasets used in the current study. For newer
%%%        strains thus function nneds to be modified.

function [PA,LR,DV] = generate_coordinate_axes(application,mu_r,landmark_to_neuron_map,axes_param,axes_neurons_to_neuron_map,landmark_names,ind_PCA,specify_PA)

switch application
    case or('GeneExpressionAnalysis','MultiCellCalciumImaging')
        A_neuron = axes_neurons_to_neuron_map(1,1);P_neuron = axes_neurons_to_neuron_map(2,1);
        L_neuron = axes_neurons_to_neuron_map(3,1);R_neuron = axes_neurons_to_neuron_map(4,1);
        D_neuron = axes_neurons_to_neuron_map(5,1);V_neuron = axes_neurons_to_neuron_map(6,1);
        mu_r_centered = mu_r - repmat(mean(mu_r),size(mu_r,1),1);
        [coeff,score,latent] = pca(mu_r_centered);
        PA = coeff(:,axes_param(1,1))';
        PA = PA/norm(PA);
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
        else
            if L_neuron ~= 0 && R_neuron ~= 0
                LR = mu_r_centered(R_neuron,:) - mu_r_centered(L_neuron,:); % LR axis based on L,R
                LR = LR/norm(LR);
                % fun = @(x)-(x(1)*coeff(1,1) + x(2)*coeff(2,1) + x(3)*coeff(3,1))^2/(x(1)^2 + x(2)^2 + x(3)^2);
                fun = @(x)-(x(1)*LR(1,1) + x(2)*LR(1,2) + LR(1,3))^2/(x(1)^2 + x(2)^2 + 1^2);
                Aeq = PA(1,1:2);
                beq = -PA(1,3);
                % x0 = coeff(:,1);
                x0 = LR(1,1:2);
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
    case 'NeuroPAL'
        mu_r_centered = mu_r - repmat(mean(mu_r),size(mu_r,1),1);
        [coeff,score,latent] = pca(mu_r_centered);
        if ind_PCA == 1
            PA = coeff(:,axes_param(1,1))';
            PA = PA/norm(PA);
            LR = coeff(:,axes_param(1,2))';
            LR = LR/norm(LR);
            DV = coeff(:,axes_param(1,3))';
            DV = DV/norm(DV);
            A_neuron = axes_neurons_to_neuron_map(1,1);P_neuron = axes_neurons_to_neuron_map(2,1);L_neuron = axes_neurons_to_neuron_map(3,1);
            R_neuron = axes_neurons_to_neuron_map(4,1);D_neuron = axes_neurons_to_neuron_map(5,1);V_neuron = axes_neurons_to_neuron_map(6,1);

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
        else
            if specify_PA
                A_neuron = axes_neurons_to_neuron_map(1,1);P_neuron = axes_neurons_to_neuron_map(2,1);L_neuron = axes_neurons_to_neuron_map(3,1);
                R_neuron = axes_neurons_to_neuron_map(4,1);D_neuron = axes_neurons_to_neuron_map(5,1);V_neuron = axes_neurons_to_neuron_map(6,1);
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
        end
    case 'LandmarkStrain'
        mu_r_centered = mu_r - repmat(mean(mu_r),size(mu_r,1),1);
        [coeff,score,latent] = pca(mu_r_centered);
        if ind_PCA == 1
            PA = coeff(:,1)';
            PA = PA/norm(PA);
            LR = coeff(:,3)';
            LR = LR/norm(LR);
            DV = coeff(:,2)';
            DV = DV/norm(DV);

            RMEL_neuron = landmark_to_neuron_map(find(strcmp('RMEL',landmark_names)),1);
            RMER_neuron = landmark_to_neuron_map(find(strcmp('RMER',landmark_names)),1);
            RMED_neuron = landmark_to_neuron_map(find(strcmp('RMED',landmark_names)),1);
            RMEV_neuron = landmark_to_neuron_map(find(strcmp('RMEV',landmark_names)),1);
            RIS_neuron = landmark_to_neuron_map(find(strcmp('RIS',landmark_names)),1);
            if ~isempty(RIS_neuron)
                if (mu_r_centered(RMER_neuron,:)-mu_r_centered(RIS_neuron,:))*PA' < 0
                    PA = -PA;
                end
            else
                DD_neuron = landmark_to_neuron_map(find(strcmp('DD1',landmark_names)),1);
                if (mu_r_centered(RMER_neuron,:)-mu_r_centered(DD_neuron,:))*PA' < 0
                    PA = -PA;
                end
            end
            if (mu_r_centered(RMER_neuron,:)-mu_r_centered(RMEL_neuron,:))*LR' < 0
                LR = -LR;
            end
            if cross(PA,LR)*DV' < 0
                DV = -DV;
            end
        else
            if specify_PA
                A_neuron = axes_neurons_to_neuron_map(1,1);P_neuron = axes_neurons_to_neuron_map(2,1);L_neuron = axes_neurons_to_neuron_map(3,1);
                R_neuron = axes_neurons_to_neuron_map(4,1);D_neuron = axes_neurons_to_neuron_map(5,1);V_neuron = axes_neurons_to_neuron_map(6,1);
                PA = coeff(:,1)';
                PA = PA/norm(PA);
                if L_neuron ~= 0 && R_neuron ~= 0
                    LR = mu_r_centered(R_neuron,:) - mu_r_centered(L_neuron,:); % LR axis based on L,R neurons
                    LR = LR/norm(LR);

                    fun = @(x)-(x(1)*LR(1,1) + x(2)*LR(1,2) + LR(1,3))^2/(x(1)^2 + x(2)^2 + 1^2);
                    Aeq = PA(1,1:2);
                    beq = -PA(1,3);
                    x0 = LR(1,1:2);
                    LR = fmincon(fun,x0,[],[],Aeq,beq);   % LR axis (perperndicular to PA axis and in direction of LR neurons)
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
                    PA = coeff(:,1)';
                    PA = PA/norm(PA);
                    fun = @(x)-(x(1)*PA(1,1) + x(2)*PA(1,2) + PA(1,3))^2/(x(1)^2 + x(2)^2 + 1^2);
                    Aeq = LR(1,1:2);
                    beq = -LR(1,3);
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
        end
end