%%% function to generate coordinate axes in images.
%%% NOTE - landmarks cells are used to disambiguate axes. This function
%%%        works for the datasets used in the current study. For newer
%%%        strains thus function nneds to be modified.

function [PA,LR,DV] = generate_coordinate_axes(strain,mu_r,landmark_to_neuron_map,axes_param,axes_neurons_to_neuron_map,landmark_names,ind_PCA,specify_PA)

A_neuron = axes_neurons_to_neuron_map(1,1);P_neuron = axes_neurons_to_neuron_map(2,1);
L_neuron = axes_neurons_to_neuron_map(3,1);R_neuron = axes_neurons_to_neuron_map(4,1);
D_neuron = axes_neurons_to_neuron_map(5,1);V_neuron = axes_neurons_to_neuron_map(6,1);
mu_r_centered = mu_r - repmat(mean(mu_r),size(mu_r,1),1);
[coeff,score,latent] = pca(mu_r_centered);

switch strain
    case 'other'
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
        elseif specify_PA == 1
            PA = coeff(:,axes_param(1,1))';
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
                LR = fmincon(fun,x0,[],[],Aeq,beq);   % LR axis (perperndicular to PA and in direction of LR neurons)
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
                if cross(PA,LR)*DV' < 0
                    DV = -DV;
                end
            elseif D_neuron ~= 0 && V_neuron ~= 0
                DV = mu_r_centered(V_neuron,:) - mu_r_centered(D_neuron,:); % DV axis based on D,V neurons
                DV = DV/norm(DV);
                % fun = @(x)-(x(1)*coeff(1,1) + x(2)*coeff(2,1) + x(3)*coeff(3,1))^2/(x(1)^2 + x(2)^2 + x(3)^2);
                fun = @(x)-(x(1)*DV(1,1) + x(2)*DV(1,2) + DV(1,3))^2/(x(1)^2 + x(2)^2 + 1^2);
                Aeq = PA(1,1:2);
                beq = -PA(1,3);
                % x0 = coeff(:,1);
                x0 = DV(1,1:2);
                DV = fmincon(fun,x0,[],[],Aeq,beq);   % DV axis (perperndicular to PA and in direction of DV neurons)
                DV = [DV,1];
                DV = DV/norm(DV);

                % LR axis (perperndicular to PA and DV axis)
                A = [PA(1,1:2);DV(1,1:2)];
                b = [-PA(1,3);-DV(1,3)];
                LR = inv(A'*A)*A'*b;           
                LR = [LR',1];
                LR = LR/norm(LR);

                if (mu_r_centered(A_neuron,:)-mu_r_centered(P_neuron,:))*PA' < 0
                    PA = -PA;
                end
                if cross(DV,PA)*LR' < 0
                    LR = -LR;
                end
            end
        else
            if L_neuron ~= 0 && R_neuron ~= 0
                LR = mu_r_centered(R_neuron,:) - mu_r_centered(L_neuron,:); % LR axis based on L,R
                LR = LR/norm(LR);
                PA = coeff(:,axes_param(1,1))';
                PA = PA/norm(PA);
                % fun = @(x)-(x(1)*coeff(1,1) + x(2)*coeff(2,1) + x(3)*coeff(3,1))^2/(x(1)^2 + x(2)^2 + x(3)^2);
                fun = @(x)-(x(1)*PA(1,1) + x(2)*PA(1,2) + PA(1,3))^2/(x(1)^2 + x(2)^2 + 1^2);
                Aeq = LR(1,1:2);
                beq = -LR(1,3);
                % x0 = coeff(:,1);
                x0 = PA(1,1:2);
                PA = fmincon(fun,x0,[],[],Aeq,beq);   % PA axis (perperndicular to LR and in direction of coeff(:,1))
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
            elseif D_neuron ~= 0 && V_neuron ~= 0
                DV = mu_r_centered(V_neuron,:) - mu_r_centered(D_neuron,:); % DV axis based on D,V neurons
                DV = DV/norm(DV);
                PA = coeff(:,axes_param(1,1))';
                PA = PA/norm(PA);
                % fun = @(x)-(x(1)*coeff(1,1) + x(2)*coeff(2,1) + x(3)*coeff(3,1))^2/(x(1)^2 + x(2)^2 + x(3)^2);
                fun = @(x)-(x(1)*PA(1,1) + x(2)*PA(1,2) + PA(1,3))^2/(x(1)^2 + x(2)^2 + 1^2);
                Aeq = DV(1,1:2);
                beq = -DV(1,3);
                % x0 = coeff(:,1);
                x0 = DV(1,1:2);
                PA = fmincon(fun,x0,[],[],Aeq,beq);   % PA axis (perperndicular to DV and in direction of coeff(:,1))
                PA = [PA,1];
                PA = PA/norm(PA);

                % LR axis (perperndicular to PA and DV axis)
                A = [PA(1,1:2);DV(1,1:2)];
                b = [-PA(1,3);-DV(1,3)];
                LR = inv(A'*A)*A'*b;           
                LR = [LR',1];
                LR = LR/norm(LR);

                if (mu_r_centered(A_neuron,:)-mu_r_centered(P_neuron,:))*PA' < 0
                    PA = -PA;
                end
                if cross(DV,PA)*LR' < 0
                    LR = -LR;
                end
            end
        end
    case 'LandmarkStrain'
        RMEL_neuron = landmark_to_neuron_map(find(strcmp('RMEL',landmark_names)),1);
        RMER_neuron = landmark_to_neuron_map(find(strcmp('RMER',landmark_names)),1);
        RMED_neuron = landmark_to_neuron_map(find(strcmp('RMED',landmark_names)),1);
        RMEV_neuron = landmark_to_neuron_map(find(strcmp('RMEV',landmark_names)),1);
        RIS_neuron = landmark_to_neuron_map(find(strcmp('RIS',landmark_names)),1);
        DD_neuron = landmark_to_neuron_map(find(strcmp('DD',landmark_names)),1);
        VD_neuron = landmark_to_neuron_map(find(strcmp('VD',landmark_names)),1);
        if ind_PCA == 1
            PA = coeff(:,axes_param(1,1))';
            PA = PA/norm(PA);
            LR = coeff(:,axes_param(1,2))';
            LR = LR/norm(LR);
            DV = coeff(:,axes_param(1,3))';
            DV = DV/norm(DV);
            
            if ~isempty(RIS_neuron)
                if (mu_r_centered(RMED_neuron,:)-mu_r_centered(RIS_neuron,:))*PA' < 0
                    PA = -PA;
                end
            elseif ~isempty(DD_neuron)
                if (mu_r_centered(DD_neuron,:)-mu_r_centered(RIS_neuron,:))*PA' < 0
                    PA = -PA;
                end
            elseif ~isempty(VD_neuron)    
                if (mu_r_centered(VD_neuron,:)-mu_r_centered(RIS_neuron,:))*PA' < 0
                    PA = -PA;
                end
            else
                if (mu_r_centered(A_neuron,:)-mu_r_centered(P_neuron,:))*PA' < 0
                    PA = -PA;
                end
            end
            
            if and(~isempty(RMER_neuron),~isempty(RMEL_neuron))
                if (mu_r_centered(RMER_neuron,:)-mu_r_centered(RMER_neuron,:))*LR' < 0
                    LR = -LR;
                end
                if cross(PA,LR)*DV' < 0
                    DV = -DV;
                end
            elseif and(~isempty(RMEV_neuron),~isempty(RMED_neuron))
                if (mu_r_centered(RMEV_neuron,:)-mu_r_centered(RMED_neuron,:))*DV' < 0
                    DV = -DV;
                end
                if cross(DV,PA)*LR' < 0
                    LR = -LR;
                end
            elseif and(R_neuron ~= 0, L_neuron ~= 0)
                if (mu_r_centered(R_neuron,:)-mu_r_centered(L_neuron,:))*LR' < 0
                    LR = -LR;
                end
                if cross(PA,LR)*DV' < 0
                    DV = -DV;
                end
            elseif and(V_neuron ~= 0, D_neuron ~= 0)
                if (mu_r_centered(V_neuron,:)-mu_r_centered(D_neuron,:))*DV' < 0
                    DV = -DV;
                end
                if cross(DV,PA)*LR' < 0
                    LR = -LR;
                end
            end
        elseif specify_PA == 1
            PA = coeff(:,axes_param(1,1))';
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
                LR = fmincon(fun,x0,[],[],Aeq,beq);   % LR axis (perperndicular to PA and in direction of LR neurons)
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
                if cross(PA,LR)*DV' < 0
                    DV = -DV;
                end
            elseif D_neuron ~= 0 && V_neuron ~= 0
                DV = mu_r_centered(V_neuron,:) - mu_r_centered(D_neuron,:); % DV axis based on D,V neurons
                DV = DV/norm(DV);
                % fun = @(x)-(x(1)*coeff(1,1) + x(2)*coeff(2,1) + x(3)*coeff(3,1))^2/(x(1)^2 + x(2)^2 + x(3)^2);
                fun = @(x)-(x(1)*DV(1,1) + x(2)*DV(1,2) + DV(1,3))^2/(x(1)^2 + x(2)^2 + 1^2);
                Aeq = PA(1,1:2);
                beq = -PA(1,3);
                % x0 = coeff(:,1);
                x0 = DV(1,1:2);
                DV = fmincon(fun,x0,[],[],Aeq,beq);   % DV axis (perperndicular to PA and in direction of DV neurons)
                DV = [DV,1];
                DV = DV/norm(DV);

                % LR axis (perperndicular to PA and DV axis)
                A = [PA(1,1:2);DV(1,1:2)];
                b = [-PA(1,3);-DV(1,3)];
                LR = inv(A'*A)*A'*b;           
                LR = [LR',1];
                LR = LR/norm(LR);

                if (mu_r_centered(A_neuron,:)-mu_r_centered(P_neuron,:))*PA' < 0
                    PA = -PA;
                end
                if cross(DV,PA)*LR' < 0
                    LR = -LR;
                end
            end
        else
            if L_neuron ~= 0 && R_neuron ~= 0
                LR = mu_r_centered(R_neuron,:) - mu_r_centered(L_neuron,:); % LR axis based on L,R
                LR = LR/norm(LR);
                PA = coeff(:,axes_param(1,1))';
                PA = PA/norm(PA);
                % fun = @(x)-(x(1)*coeff(1,1) + x(2)*coeff(2,1) + x(3)*coeff(3,1))^2/(x(1)^2 + x(2)^2 + x(3)^2);
                fun = @(x)-(x(1)*PA(1,1) + x(2)*PA(1,2) + PA(1,3))^2/(x(1)^2 + x(2)^2 + 1^2);
                Aeq = LR(1,1:2);
                beq = -LR(1,3);
                % x0 = coeff(:,1);
                x0 = PA(1,1:2);
                PA = fmincon(fun,x0,[],[],Aeq,beq);   % PA axis (perperndicular to LR and in direction of coeff(:,1))
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
            elseif D_neuron ~= 0 && V_neuron ~= 0
                DV = mu_r_centered(V_neuron,:) - mu_r_centered(D_neuron,:); % DV axis based on D,V neurons
                DV = DV/norm(DV);
                PA = coeff(:,axes_param(1,1))';
                PA = PA/norm(PA);
                % fun = @(x)-(x(1)*coeff(1,1) + x(2)*coeff(2,1) + x(3)*coeff(3,1))^2/(x(1)^2 + x(2)^2 + x(3)^2);
                fun = @(x)-(x(1)*PA(1,1) + x(2)*PA(1,2) + PA(1,3))^2/(x(1)^2 + x(2)^2 + 1^2);
                Aeq = DV(1,1:2);
                beq = -DV(1,3);
                % x0 = coeff(:,1);
                x0 = DV(1,1:2);
                PA = fmincon(fun,x0,[],[],Aeq,beq);   % PA axis (perperndicular to DV and in direction of coeff(:,1))
                PA = [PA,1];
                PA = PA/norm(PA);

                % LR axis (perperndicular to PA and DV axis)
                A = [PA(1,1:2);DV(1,1:2)];
                b = [-PA(1,3);-DV(1,3)];
                LR = inv(A'*A)*A'*b;           
                LR = [LR',1];
                LR = LR/norm(LR);

                if (mu_r_centered(A_neuron,:)-mu_r_centered(P_neuron,:))*PA' < 0
                    PA = -PA;
                end
                if cross(DV,PA)*LR' < 0
                    LR = -LR;
                end
            end
        end
end