function [obj_val] = objFun(test,objectiveStruct)
    %% Globals
    global variable_tracking_struct 
    Fopt = [];
    
    %% get intenter's center of mass (not used currently)
    indenter_center_of_mass_sim = [NaN NaN NaN];
    if isfield(objectiveStruct,'febio_spec')
        indenter_center_of_mass_sim = objectiveStruct.febio_spec.Material.material{2}.center_of_mass;
    end
    if isfield(test,'MeshGeometry')
        indenter_center_of_mass_sim = test.MeshGeometry.Indenter.center_of_mass;
    end
    
    
    %% baseline ("experimental") data
    disp_exp = objectiveStruct.disp_exp; % (Nnodes)x3x(Nt)
    pos_exp = objectiveStruct.pos_exp; % (Nnodes)x3x(Nt)
    force_exp = objectiveStruct.force_exp; % 1x(Nt)
    %% trial ("simulated") data
    [~,pos_sim,disp_sim] = getNPosMat(test);        
    force_sim = test.indenter_RB_out.Fz.data;
    indentation_depth_sim = test.indenter_RB_out.z.data-indenter_center_of_mass_sim(3);
    indentation_depth_sim(1) = 0; 
    %% Nodal Weights
    distances = abs(pos_sim(:,1,1));
    [Wn] =getWn(distances,objectiveStruct.indenterRadius);
    Wn = Wn/norm(Wn,2);
    %% Residuals in force
    force_res=[force_exp-force_sim];
    force_res_mag = abs(force_res);
    force_exp_mag = abs(force_exp); % Normalization factor for force_res_mag
    squared_normalized_force_res = (force_res_mag./force_exp_mag).^2;
    % Replace NaN values with 0 (originating from normalizing by 0).
    squared_normalized_force_res(isnan(squared_normalized_force_res))=0; % Ff(delta)@(p_sim,p_exp) eq. (3.1) in paper
    %% Residuals in displacement
    disp_res_mag = vecnorm(disp_exp-disp_sim,2,2); % Not interpolated (synchronised+same mesh)
    disp_exp_mag = vecnorm(disp_exp,2,2);
    squared_normalized_disp_res = squeeze((disp_res_mag./disp_exp_mag).^2); % notice the squeeze(): dim1 - nodes, dim2 - indentation steps
    % Replace NaN values with 0 (normalizing by disp_exp_mag(:,:,1)==0)
    squared_normalized_disp_res(isnan(squared_normalized_disp_res))=0;
    node_weighted_squared_normalized_disp_res = (Wn')*squared_normalized_disp_res; % Fu(delta)@(p_sim,p_exp) eq. (3.2) in paper
    
    %% Set output structure
    obj_val.Ff = squared_normalized_force_res;
    obj_val.Fu = node_weighted_squared_normalized_disp_res;
    
%% old 
%%%%%%%%FORCE %%%%%%%%%%%%%%%
%         
%        
%          %%%%%%%%%%%%%%%%%%%%%%%%%%% "Time interpolation". "Space" interpolation
%         %%%%%%%%%%%%%%%%%%%%%%%%%%% not need currently as both meshes are
%         %%%%%%%%%%%%%%%%%%%%%%%%%%% identical
%         Nn = size(pos_sim,1);
%         Nt = length(disp_exp);
% %         pos_sim_exp = nan(Nn,size(pos_sim,2),Nt);
% %         pos_sim_exp(:,1,:) = interp1(disp_sim,squeeze(pos_sim(:,1,:))',disp_exp,'pchip')';
% %         pos_sim_exp(:,2,:) = interp1(disp_sim,squeeze(pos_sim(:,2,:))',disp_exp,'pchip')';
% %         pos_sim_exp(:,3,:) = interp1(disp_sim,squeeze(pos_sim(:,3,:))',disp_exp,'pchip')';
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %Derive Fopt
%         force_res=[force_exp-force_sim]';
%         
%         positionDev = pos_exp-pos_sim; %not interpolated since mesh is identicle
%         positionDevMag = (squeeze((positionDev(:,1,:).^2+positionDev(:,2,:).^2+positionDev(:,3,:).^2))).^.5; % magnitude of position difference
% %         positionDevMag = (squeeze((positionDev(:,3,:).^2))); %square of magnitude of position difference
%         %% switch objectiveStruct.formulation
%         
%         switch objectiveStruct.formulation
%             case 1
%                 %%
%                 Fopt=sum((force_res).^2);
%                 Fopt_plot=Fopt; %Sum of squared differences
%                 obj_val = Fopt_plot;
%             case 2
%                 %%
%                 Fopt=force_res(:);%(stressDev).^2; %Squared differences
%                 Fopt_plot=sum((force_res).^2);%(stressDev).^2; %Squared differences
%                 obj_val = Fopt_plot;
%             case 3
%                 %%
%     %             eta =0.5;
%                 %%%%%%%%%%
%                 Wt = ones(Nt,1);
% %                 Wt(1) = 0;
% %                 Wt(1:end-1) = 0;
%                 %%%%%%%%%%
%                 %%%%%%%%%% Spacial weights
% %                 Wn = ones(Nn,1).*(1./sqrt((sum(pos_sim(:,:,1).^2,2))));
%     %             Wn = ones(Nn,1).*(1./(sum(pos_sim(:,:,1).^2,2)));
%                     distances = sqrt((sum(pos_sim(:,:,1).^2,2)));
% 
%                 [Wn] =getWn(distances,objectiveStruct.indenterRadius);
%                 %%%%%%%%%%
%                 nodeWeightedPositionDevMag = (Wn.*positionDevMag)/sum(Wn); %A - (Nn)x(Nt) matrix, A_nt = (1/sum(Wn))*[pos_dev_nt*Wn_n]
%                 timeWeightedPositionDevMag = (nodeWeightedPositionDevMag*Wt)/sum(Wt);%B - (Nn)x(1) vector, B_n= A_nt*Wt_t
% 
%                 F2 = sum(timeWeightedPositionDevMag); %contribution from difference in nodal position
%                 F1 = sum((force_res).^2); %contribution from difference in reaction force magnitude
% 
%                 switch variable_tracking_struct.isfirst_flg
%                     case 1 % always 1 on first step
%                         % Set factors according to values of first step
% %                         F1_norm_Factor = 1/F1; 
%                         F1_norm_Factor = 1/F1
%                         F2_norm_Factor = 1/F2;
%                         % Save factors in tracking struct
%                         variable_tracking_struct.F1_norm_Factor = F1_norm_Factor; 
%                         variable_tracking_struct.F2_norm_Factor = F2_norm_Factor;
%                         % Update flag
%                         variable_tracking_struct.isfirst_flg = 0;
%                     case 0
%                         % Retreive factors from tracking struct
%                         F1_norm_Factor = variable_tracking_struct.F1_norm_Factor;
%                         F2_norm_Factor = variable_tracking_struct.F2_norm_Factor;
%                 end
%                 % Apply factors to normalize objective functions
%                 F1_normalized = F1_norm_Factor*F1;
%                 F2_normalized = F2_norm_Factor*F2;
% 
%                 % Merge normalized objective functions and apply penalty
%                 % factor eta
%                 Fopt=(eta)*F1_normalized+(1-eta)*F2_normalized;
%                 Fopt_plot=Fopt;
%                 
%                 F1_opt_plot=(eta)*F1_normalized;
%                 F2_opt_plot=(1-eta)*F2_normalized;    
%                 
%                 obj_val = Fopt_plot;
%             case 4
%                 %%
%                 %%%%%%%%%%
%                 Wt = ones(Nt,1);
%                 Wt(1) = 0;
%     %             Wt(1:end-1) = 0; %wrong convergence
%     %             Wt = [1:Nt]'; % very wrong P_opt
% 
%                 %%%%%%%%%%
%                 %%%%%%%%%% Spacial weights
%                 distances = sqrt((sum(pos_sim(:,:,1).^2,2)));
% 
%     %             Wn = ones(Nn,1).*(1./distances);
%     %              Wn = ones(Nn,1).*(1./sqrt((sum(pos_sim(:,:,1).^2,2))));
%     %             Wn = ones(Nn,1).*(1./(distances.^2));             Wn(Wn==Inf)=0;
%     %             objectiveStruct.indenterRadius
%     %             getWn(distances,objectiveStruct.indenterRadius)
%                 [Wn] =getWn(distances,objectiveStruct.indenterRadius);
%                 %%%%%%%%%%
%                 nodeWeightedPositionDevMag = (Wn.*positionDevMag); %A - (Nn)x(Nt) matrix, A_nt = (1/sum(Wn))*[pos_dev_nt*Wn_n]
%     %             nodeWeightedPositionDevMag = (Wn.*positionDevMag); %A - (Nn)x(Nt) matrix, A_nt = (1/sum(Wn))*[pos_dev_nt*Wn_n]
%                 timeWeightedPositionDevMag = (nodeWeightedPositionDevMag*Wt)/sum(Wt);%B - (Nn)x(1) vector, B_n= A_nt*Wt_t
% 
%     %             timeWeightedForceDev = (force_res.^2).*Wt/sum(Wt);
%                 timeWeightedForceDev = abs(force_res).*Wt/sum(Wt);
% 
%                 F2 = sum(timeWeightedPositionDevMag); %contribution from difference in nodal position
%     %             F1 = sum((force_res).^2); %contribution from difference in reaction force magnitude
%                 F1 = sum(timeWeightedForceDev); %contribution from difference in reaction force magnitude
%                 switch variable_tracking_struct.isfirst_flg
%                     case 1 % always 1 on first step
%                         % Set factors according to values of first step
%                         F1_norm_Factor = 1/F1; 
%                         F2_norm_Factor = 1/F2;
%                         % Save factors in tracking struct
%                         variable_tracking_struct.F1_norm_Factor = F1_norm_Factor; 
%                         variable_tracking_struct.F2_norm_Factor = F2_norm_Factor;
%                         % Update flag
%                         variable_tracking_struct.isfirst_flg = 0;
%                     case 0
%                         % Retreive factors from tracking struct
%                         F1_norm_Factor = variable_tracking_struct.F1_norm_Factor;
%                         F2_norm_Factor = variable_tracking_struct.F2_norm_Factor;
%                 end
%                 % Apply factors to normalize objective functions
%                 F1_normalized = F1_norm_Factor*F1;
%                 F2_normalized = F2_norm_Factor*F2;
% 
%                 % Merge normalized objective functions and apply penalty
%                 % factor eta
%     %             Fopt=(eta)*F1_normalized+(1-eta)*F2_normalized;
%                 Fopt=sqrt([(eta)*F1_normalized;(1-eta)*F2_normalized]);
%                 Fopt_plot=(eta)*F1_normalized+(1-eta)*F2_normalized;
%                 F1_opt_plot=(eta)*F1_normalized;
%                 F2_opt_plot=(1-eta)*F2_normalized;    
%                 
%                 obj_val = Fopt_plot;
%              case 5
%                 %%
%                 %%%%%%%%%%
%                 Wt = ones(Nt,1);
%                 Wt(1) = 0;
% %                 Wt(1:end-1) = 0; %wrong convergence
%     %             Wt = [1:Nt]'; % very wrong P_opt
% 
%                 %%%%%%%%%%
%                 %%%%%%%%%% Spacial weights
%                 distances = sqrt((sum(pos_sim(:,:,1).^2,2)));
% 
%     %             Wn = ones(Nn,1).*(1./distances);
%     %              Wn = ones(Nn,1).*(1./sqrt((sum(pos_sim(:,:,1).^2,2))));
%     %             Wn = ones(Nn,1).*(1./(distances.^2));             Wn(Wn==Inf)=0;
%     %             objectiveStruct.indenterRadius
%     %             getWn(distances,objectiveStruct.indenterRadius)
%                 [Wn] =getWn(distances,objectiveStruct.indenterRadius);
%                 %%%%%%%%%%
%                 nodeWeightedPositionDevMag = (Wn.*positionDevMag); %A - (Nn)x(Nt) matrix, A_nt = (1/sum(Wn))*[pos_dev_nt*Wn_n]
%     %             nodeWeightedPositionDevMag = (Wn.*positionDevMag); %A - (Nn)x(Nt) matrix, A_nt = (1/sum(Wn))*[pos_dev_nt*Wn_n]
% %                 timeWeightedPositionDevMag = (nodeWeightedPositionDevMag*Wt)/sum(Wt);%B - (Nn)x(1) vector, B_n= A_nt*Wt_t
%                 timeWeightedPositionDevMag = (nodeWeightedPositionDevMag*Wt)/norm(Wt);%B - (Nn)x(1) vector, B_n= A_nt*Wt_t
% 
%                 %timeWeightedForceDev = (force_res.^2).*Wt/sum(Wt);
% %                 timeWeightedForceDev = (force_res).*Wt/sum(Wt);
%                                 timeWeightedForceDev = (force_res).*Wt/norm(Wt);
% 
% 
% %                 timeWeightedForceDev = abs(force_res).*Wt/sum(Wt);
% 
%                 F2 = timeWeightedPositionDevMag; %contribution from difference in nodal position
% %                 F2 = sum(timeWeightedPositionDevMag); %contribution from difference in nodal position
%     %             F1 = sum((force_res).^2); %contribution from difference in reaction force magnitude
%                 F1 = timeWeightedForceDev; %contribution from difference in reaction force magnitude
%                 switch variable_tracking_struct.isfirst_flg
%                     case 1 % always 1 on first step
%                         % Set factors according to values of first step
%                         F1_norm_Factor = 1/norm(F1); 
%                         F2_norm_Factor = 1/norm(F2);
%                         
% %                         F1_norm_Factor(isinf(F1_norm_Factor))= 0;
% %                         F2_norm_Factor(isinf(F2_norm_Factor))= 0;
% 
%                         % Save factors in tracking struct
%                         variable_tracking_struct.F1_norm_Factor = F1_norm_Factor; 
%                         variable_tracking_struct.F2_norm_Factor = F2_norm_Factor;
%                         % Update flag
%                         variable_tracking_struct.isfirst_flg = 0;
%                     case 0
%                         % Retreive factors from tracking struct
%                         F1_norm_Factor = variable_tracking_struct.F1_norm_Factor;
%                         F2_norm_Factor = variable_tracking_struct.F2_norm_Factor;
%                 end
%                 % Apply factors to normalize objective functions
%                 F1_normalized = F1_norm_Factor*F1;
%                 F2_normalized = F2_norm_Factor*F2;
% 
%                 % Merge normalized objective functions and apply penalty
%                 % factor eta
%     %             Fopt=(eta)*F1_normalized+(1-eta)*F2_normalized;
%                 Fopt = [eta*F1_normalized;(1-eta)*F2_normalized];
% %                 Fopt=sqrt([(eta)*F1_normalized;(1-eta)*F2_normalized]);
% %                 Fopt_plot=(eta)*F1_normalized+(1-eta)*F2_normalized;
%                 F1_opt_plot=(eta)*norm(F1_normalized);
%                 F2_opt_plot=(1-eta)*norm(F2_normalized);
%                 Fopt_plot = F1_opt_plot+F2_opt_plot;
%                 
%                 obj_val = Fopt_plot;
%             case 6 % Suggested formulation (24/10)
%                 %%
%                 %%%%%%%%%%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % note to self: this combination increased sensitivity in both low
%                 % and high values of eta
%                                 Wt_force = ones(Nt,1); Wt_force(1) = 0;
%                                 Wt_pos = ones(Nt,1); Wt_pos(1) = 0;
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                 Wt_force = (0:(Nt-1))';
% %                 Wt_pos = (0:(Nt-1))';
% %                 Wt_force = ones(Nt,1); Wt_force(1) = 0;
% %                 Wt_pos = ones(Nt,1); Wt_pos(1:end) = 0;
% %                 Wt_disp(1:end
% %                 Wt = ones(Nt,1);
% %                 Wt(1) = 0;
% %                 Wt(1:end-1) = 0; %wrong convergence
%     %             Wt = [1:Nt]'; % very wrong P_opt
% 
%                 %%%%%%%%%%
%                 %%%%%%%%%% Spacial weights
%                 distances = sqrt((sum(pos_sim(:,:,1).^2,2)));
% 
%     %             Wn = ones(Nn,1).*(1./distances);
%     %              Wn = ones(Nn,1).*(1./sqrt((sum(pos_sim(:,:,1).^2,2))));
%     %             Wn = ones(Nn,1).*(1./(distances.^2));             Wn(Wn==Inf)=0;
%     %             objectiveStruct.indenterRadius
%     %             getWn(distances,objectiveStruct.indenterRadius)
%                 [Wn] =getWn(distances,objectiveStruct.indenterRadius);
%                 %%%%%%%%%%
% %                         positionDevMag = (squeeze((positionDev(:,1,:).^2+positionDev(:,2,:).^2+positionDev(:,3,:).^2))); %square of magnitude of position difference
%                 positionDev_fromReference = pos_exp(:,:,:)-pos_exp(:,:,1); % total displacment of reference node from initial state
%                 positionDev_fromReferenceMag = sqrt(squeeze(sum(positionDev_fromReference.^2,2)));
% %                 (squeeze((positionDev_fromReference(:,1,:).^2+positionDev_fromReference(:,2,:).^2+positionDev_fromReference(:,3,:).^2))); %square of magnitude of total displacment of reference node from initial state (%D hat)
%                 
% 
%                 positionDevMag_normalized = positionDevMag./positionDev_fromReferenceMag;
%                 positionDevMag_normalized = positionDevMag_normalized.^2;
%                 
%                 %%%%%%  norm1
%                     nodeWeightedPositionDevMag1 = (Wn'*positionDevMag_normalized)/sum(Wn); 
%                     
%                     nodeWeightedPositionDevMag1(isnan(nodeWeightedPositionDevMag1))=0; % add zero to first element (which is NaN because normalizing factor is 0)
%                 %%%%%%  norm2
%                     nodeWeightedPositionDevMag2 = (Wn'*positionDevMag_normalized)/sqrt(sum(Wn.^2)); 
%                     nodeWeightedPositionDevMag2(isnan(nodeWeightedPositionDevMag2))=0; % add zero to first element (which is NaN because normalizing factor is 0)
%     %             nodeWeightedPositionDevMag = (Wn.*positionDevMag); %A - (Nn)x(Nt) matrix, A_nt = (1/sum(Wn))*[pos_dev_nt*Wn_n]
%                 
%     
%                 timeWeightedPositionDevMag1 = (nodeWeightedPositionDevMag1*Wt_pos)/sum(Wt_pos);%B - (Nn)x(1) vector, B_n= A_nt*Wt_t
%                 timeWeightedPositionDevMag2 = (nodeWeightedPositionDevMag2*Wt_pos)/sum(Wt_pos);%B - (Nn)x(1) vector, B_n= A_nt*Wt_t
% 
% %                 timeWeightedPositionDevMag = (nodeWeightedPositionDevMag*Wt_pos)/sqrt(sum(Wt_pos.^2));%B - (Nn)x(1) vector, B_n= A_nt*Wt_t
%                 %%%% norm2_sq
% %                 nodeWeightedPositionDevMag2_sq = (Wn'*positionDevMag_normalized.^2)/sum(Wn); 
% %                 nodeWeightedPositionDevMag2_sq(isnan(nodeWeightedPositionDevMag1_sq))=0; % add zero to first element (which is NaN because normalizing factor is 0)
%                 
% 
% %                 Wt = ones(Nt,1); Wt(1) = 0;
%                 ForceDevMagNormalized = ((force_res'.^2)./(force_exp.^2)); %Squared force
% %                 ForceDevMagNormalized = (abs(force_res')./abs(force_exp)); % force not squared
%                 ForceDevMagNormalized(isnan(ForceDevMagNormalized))=0;
%                 timeWeightedForceDev = ForceDevMagNormalized*Wt_force/sum(Wt_force);
%                 
%                 
% 
%                 F2 = timeWeightedPositionDevMag1; %contribution from difference in nodal position
%                 F2_2 = timeWeightedPositionDevMag2; %contribution from difference in nodal position
%                 
%     %             F1 = sum((force_res).^2); %contribution from difference in reaction force magnitude
%                 F1 = timeWeightedForceDev; %contribution from difference in reaction force magnitude
% %                 switch variable_tracking_struct.isfirst_flg
% %                     case 1 % always 1 on first step
% %                         % Set factors according to values of first step
% %                         F1_norm_Factor = 1/F1; 
% %                         F2_norm_Factor = 1/F2;
% %                         % Save factors in tracking struct
% %                         variable_tracking_struct.F1_norm_Factor = F1_norm_Factor; 
% %                         variable_tracking_struct.F2_norm_Factor = F2_norm_Factor;
% %                         % Update flag
% %                         variable_tracking_struct.isfirst_flg = 0;
% %                     case 0
% %                         % Retreive factors from tracking struct
% %                         F1_norm_Factor = variable_tracking_struct.F1_norm_Factor;
% %                         F2_norm_Factor = variable_tracking_struct.F2_norm_Factor;
% %                 end
%                 % Apply factors to normalize objective functions
%                 F1_normalized = F1; % force
%                 F2_normalized = F2;
% 
%                 % Merge normalized objective functions and apply penalty
%                 % factor eta
%     %             Fopt=(eta)*F1_normalized+(1-eta)*F2_normalized;
% %                 Fopt=sqrt([(eta)*F1_normalized;(1-eta)*F2_normalized]);
% %                 Fopt_plot=(eta)*F1_normalized+(1-eta)*F2_normalized;
% %                 F1_opt_plot=(eta)*F1_normalized;
% %                 F2_opt_plot=(1-eta)*F2_normalized;    
%                 
%                 
%                 
%                 obj_val.F1 = F1;
%                 obj_val.F2 = F2;
%                 obj_val.F2_2 = F2_2;
%                 
% %                 obj_val.Fopt_plot = Fopt_plot;
%                 
%                 obj_val.nodeWeightedPositionDevMag1 = nodeWeightedPositionDevMag1;
%                 obj_val.nodeWeightedPositionDevMag2 = nodeWeightedPositionDevMag2;                 obj_val.ForceDevMagNormalized = ForceDevMagNormalized;
%         end
end
