function runOptimization(febio_spec,febioAnalysis, runFlag,savePath,MeshGeometry)
        %% Import FEBio results 
        % Importing nodal displacements from a log file
        [time_mat, N_disp_mat,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp)); %Nodal displacements
        time_mat=[0; time_mat(:)]; %Time

        N_disp_mat=N_disp_mat(:,2:end,:);
        sizImport=size(N_disp_mat);
        sizImport(3)=sizImport(3)+1;
        N_disp_mat_n=zeros(sizImport);
        N_disp_mat_n(:,:,2:end)=N_disp_mat;
        N_disp_mat=N_disp_mat_n;
        N_disp_magnitude_mat=squeeze(sqrt(sum(N_disp_mat.^2,2)));
        DN=N_disp_mat(:,:,end);
        DN_magnitude=N_disp_magnitude_mat(:,end);
        V_def=V+DN;
        V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
        X_DEF=V_DEF(:,1,:);
        Y_DEF=V_DEF(:,2,:);
        Z_DEF=V_DEF(:,3,:);
        Fb1 = Fb;
        E2_joined = E3;
        [CF]=vertexToFaceMeasure(Fb1,DN_magnitude);

                    % Importing rigid body data from a log file
        [~, RigidBody_mat,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_rigidBody)); %Nodal displacements
        %     RigidBody_mat=RigidBody_mat(:,2:end);
        %     RigidBody_mat=[zeros(1,size(RigidBody_mat,2)); RigidBody_mat];
        %     RigidBody_force_magnitude=sqrt(sum(RigidBody_mat(:,1:3).^2,2));
        %     RigidBody_torque_magnitude=sqrt(sum(RigidBody_mat(:,4:6).^2,2));
        RigidBody_mat = squeeze(RigidBody_mat)';
        RigidBody_mat(:,1) = []; %remove Indices
        RigidBody_mat = [zeros(1,size(RigidBody_mat,2)); RigidBody_mat]; % Add zero state row
        RigidBody_displacement=abs(RigidBody_mat(:,3)-center_of_mass(3));
        RigidBody_force_magnitude=sqrt(sum(RigidBody_mat(:,1:3).^2,2));
        RigidBody_torque_magnitude=sqrt(sum(RigidBody_mat(:,4:6).^2,2));






         % Importing surface nodes nodal position from a log file
        [time_mat, N_pos_mat,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_pos),1,1); %Nodal position
        N_pos_mat(:,:,1) = febio_spec.Mesh.Nodes{1}.node.VAL(febio_spec.Mesh.NodeSet{2}.node.ATTR.id,:); % Get zero state from febio_spec


        %     
        %     % Importing reaction forces from a log file
        %     [~, R_force_mat,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_force)); %Nodal displacements
        %     
        %     R_force_mat=R_force_mat(:,2:end,:);
        %     sizImport=size(R_force_mat);
        %     sizImport(3)=sizImport(3)+1;
        %     R_force_mat_n=zeros(sizImport);
        %     R_force_mat_n(:,:,2:end)=R_force_mat;
        %     R_force_mat=R_force_mat_n;
        %     indenter_force=squeeze(R_force_mat(E2_joined(1,1),:,:))';
        %     indenter_force_magnitude=sqrt(sum(indenter_force.^2,2));
        %% Plotting the simulated results using |anim8| to visualize and animate
        % deformations

        % Create basic view and store graphics handle to initiate animation
        hf=cFigure; %Open figure
        hp1=gpatch(Fb1,V_def,CF,'none',1); %Add graphics object to animate

        % fake legend
        for ii=1:nSteps
                gpatch([1 2 3],zeros(3),colors(ii,:),'none',faceAlpha2); %Add graphics object to animate
        end

        % title with force result
        ht1=gtitle([num2str(RigidBody_force_magnitude(end),'%.2f') ' N'],fontSize2);
        % plot indenter
        hp2=gpatch(E3,V_def,colors(ii,:),'none',faceAlpha2); %Add graphics object to animate

        axisGeom(gca,fontSize1);
        colormap cool;
        hc=colorbar;
        title(hc,'dispMgn [mm]');
        caxis([0 max(N_disp_magnitude_mat(:))]);
        axis([min(X_DEF(:)) max(X_DEF(:)) min(Y_DEF(:)) max(Y_DEF(:)) min(Z_DEF(:)) max(Z_DEF(:))]);
        camlight headlight;

        % Set up animation features
        animStruct.Time=time_mat; %The time vector
        for qt=1:1:size(N_disp_mat,3) %Loop over time increments
            DN=N_disp_mat(:,:,qt); %Current displacement
            DN_magnitude=N_disp_magnitude_mat(:,qt); %Current displacement magnitude
            V_def=V+DN; %Current nodal coordinates
            [CF]=vertexToFaceMeasure(Fb1,DN_magnitude); %Current color data to use
            % indenter color according to time step
            colorInd=ceil(time_mat(qt));
            colorInd(colorInd==0)=1;

            %Set entries in animation structure
            animStruct.Handles{qt}=[hp1 hp1 hp2 hp2 ht1]; %Handles of objects to animate
            animStruct.Props{qt}={'Vertices','CData','Vertices','FaceColor','String'}; %Properties of objects to animate
            animStruct.Set{qt}={V_def,CF,V_def,colors(colorInd,:),[num2str(RigidBody_force_magnitude(qt),'%.2f') ' N']}; %Property values for to set in order to animate
        end
        anim8(hf,animStruct); %Initiate animation feature
        set(gca,'FontSize',fontSize2);

        legendStrings={'cylinder'};
        for ii=1:nSteps
            legendStrings=[legendStrings, {['indentation ', num2str(ii)]}];
        end
        legend(legendStrings);

         %     addFigureButtons
        drawnow;
        %% "Experimental" results

        must_time_inds=find(sum(time_mat==[0:nSteps],2));% time indices of must points
        force_exp=RigidBody_force_magnitude(must_time_inds); % "experimental results" (FEA results with accurate parameters in the data must points)
        displacement_exp=squeeze(abs(N_disp_mat(E3(1,1),3,must_time_inds)));
        position_exp = N_pos_mat(:,:,must_time_inds);
        %     
        %     %Interpolate to higher sampling
        %     n_sampling=100;
        %     displacement_exp_n = linspace(0,max(displacement_exp),n_sampling);
        %     force_exp_n = interp1(displacement_exp,force_exp,displacement_exp_n,'pchip');
        %     force_exp=force_exp_n;
        %     displacement_exp=displacement_exp_n;
        %     %Add noise
        %     force_exp = awgn(force_exp,snr);
        %     force_exp(force_exp<0)=0;
        %% PLOT force-displacement curve - sinulation vs. experimental
        hf2=cFigure; hold on;
        force_sim=RigidBody_force_magnitude;
        displacement_sim=squeeze(abs(N_disp_mat(E3(1,1),3,:)));
        hl(1)=plot(displacement_exp,force_exp,'bo','lineWidth',lineWidth,'markerSize',2*markerSize2);
        hl(2)=plot(displacement_sim,force_sim,'ro--','lineWidth',lineWidth,'markerSize',markerSize2);
        ht(1)=text(0.1,0.8,['c_t_r_u=' num2str(c1_tru,'%3.3e') ', m_t_r_u=' num2str(m1_tru,'%.3f')  ', k_t_r_u=' num2str(k_tru,'%.3f')],'Units','Normalized','FontSize',fontSize2);
        ht(2)=text(0.1,0.7,['c_i_n_i=' num2str(c1_ini,'%3.3e') ', m_i_n_i=' num2str(m1_ini,'%.3f')  ', k_i_n_i=' num2str(k_ini,'%.3f')],'Units','Normalized','FontSize',fontSize2);
        ht(3)=text(0.1,0.6,['c_s_i_m=' num2str(c1_ini,'%3.3e') ', m_s_i_m=' num2str(m1_ini,'%.3f') ', k_s_i_m=' num2str(k_ini,'%.3f')],'Units','Normalized','FontSize',fontSize2);
    %     ht(4)=text(0.1,0.5,['F_o_p_t=' num2str(NaN,'%3.3e')],'Units','Normalized','FontSize',fontSize2);
        ht(4)=text(0.1,0.5,['F_o_p_t=' num2str(NaN,'%3.3e') ',{F_{opt}}^1=',num2str(NaN,'%3.3e'),',{F_{opt}}^2=',num2str(NaN,'%3.3e') ],'Units','Normalized','FontSize',fontSize2);
        ht(5)=title(['Iteration number ' num2str(0)]);
        xlabel('Displacement [mm]');
        ylabel('Force [N]');
        set(gca,'FontSize',fontSize2);
        legend(hl,{'Experiment','Simulation'},'Location','northwest');
        drawnow;
        %% Create structures for optimization
        % Material structure
        mat_struct.id=1; %Material id
        mat_struct.par_names={'c1','m1','c2','m2','k'}; %Parameter names
        mat_struct.par_values={c1_ini m1_ini c1_ini -m1_ini k_ini}; %Parameter values

        % docNode=set_mat_par_FEBIO(FEB_struct.run_filename,FEB_struct.run_filename,{mat_struct});

        febioAnalysis.disp_on=0;
        febioAnalysis.disp_log_on=0;

        % What should be known to the objective function:
        objectiveStruct.hl=hl(2); % handle for simulated force-displacement curve
        objectiveStruct.ht2=ht(3); % handle for parameter text
        objectiveStruct.ht3=ht(4); % handle for Fopt
        objectiveStruct.ht4=ht(5); % handle for title
        objectiveStruct.force_exp=force_exp;
        objectiveStruct.displacement_exp=displacement_exp;

        objectiveStruct.position_exp=position_exp; %<<<<<<< ADDED THIS

        objectiveStruct.febioAnalysis=febioAnalysis;
        objectiveStruct.febio_spec=febio_spec;
        objectiveStruct.febioFebFileName=febioFebFileName;
        objectiveStruct.mat_struct=mat_struct;
        objectiveStruct.parNormFactors=P; %This will normalize the parameters to ones(size(P))
%         objectiveStruct.Pb_struct.xx_c=P; %Parameter constraining centre
%         objectiveStruct.Pb_struct.xxlim=[[P(1)/100 P(2)/100 P(3)/100]' [P(1)*100 P(2)*100 P(3)*100]']; %Parameter bounds
        objectiveStruct.indenterVertex=E2_joined(1,1);
        objectiveStruct.iterationNumber=0;

        %File names of output files
        output_names.displacement=fullfile(savePath,febioLogFileName_disp);
        output_names.force=fullfile(savePath,febioLogFileName_rigidBody);
        objectiveStruct.run_output_names=output_names;
        objectiveStruct.indenterRadius = indenterRadius;

        objectiveStruct.k_factor_tru = k_factor_tru
        global variable_tracking_struct;%<-create a global for tracking data
        variable_tracking_struct=struct;
        variable_tracking_struct.isfirst_flg = 1;
        variable_tracking_struct.Fopt_previous = 1;
        % variable_tracking_struct.displacement_sim_all=nan(size(displacement_sim,1),maxNumberIterations);
        % variable_tracking_struct.displacement_sim_all=nan(size(displacement_sim,1),maxNumberIterations);
        % variable_tracking_struct.force_sim_all=nan(size(force_sim,1),maxNumberIterations);
        % variable_tracking_struct.P_all=nan(maxNumberIterations,3);
        % variable_tracking_struct.Fopt_all=nan(maxNumberIterations,1);
        % variable_tracking_struct.job_count = 1;
        %% start optimization
        maxNumberIterations = 10;
        parameterTolerance =1e-2;
        rng default
        Pn=P./objectiveStruct.parNormFactors;
        objectiveStruct.method=5; % 1=> fminsearch and Nelder-Mead, 2=> lsqnonlin and Levenberg-Marquardt, 3=> fminsearch and Nelder-Mead (force+displacment),4=> lsqnonlin and Levenberg-Marquardt (force+displacment)
        switch objectiveStruct.method
            case 1 %fminsearch and Nelder-Mead
                OPT_options=optimset('fminsearch'); % 'Nelder-Mead simplex direct search'
                OPT_options = optimset(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
                    'MaxIter',maxNumberIterations,...
                    'TolFun',functionTolerance,...
                    'TolX',parameterTolerance,...
                    'Display',displayTypeIterations,...
                    'FinDiffRelStep',1e-2,...
                    'DiffMaxChange',0.5);
                [Pn_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= fminsearch(@(Pn) objectiveFunctionIFEA(Pn,objectiveStruct),Pn,OPT_options);
            case 2 %lsqnonlin and Levenberg-Marquardt
                OPT_options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
                OPT_options = optimoptions(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
                    'MaxIter',maxNumberIterations,...
                    'TolFun',functionTolerance,...
                    'TolX',parameterTolerance,...
                    'Display',displayTypeIterations,...
                    'FinDiffRelStep',1e-2,...
                    'DiffMaxChange',0.5);
                [Pn_opt,OPT_out.resnorm,OPT_out.residual]= lsqnonlin(@(Pn) objectiveFunctionIFEA(Pn,objectiveStruct),Pn,[],[],OPT_options);
             case 3 %fminsearch and Nelder-Mead
                OPT_options=optimoptions('patternsearch','Cache','on');
                OPT_options = optimoptions(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
                    'MaxIter',maxNumberIterations,...
                    'TolFun',functionTolerance,...
                    'TolX',parameterTolerance,...
                    'Display',displayTypeIterations);
                % patternsearch on C1
                [C1n_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= patternsearch(@(t) objectiveFunctionIFEA([t,1],objectiveStruct,1),Pn(1),[],[],[],[],[Pn(1)*1e-2],[Pn(1)*1e2],[],OPT_options);
                    C1n_temp_meshsize = OPT_out.output.meshsize;
                    variable_tracking_struct.Fopt_previous = OPT_out.fval*variable_tracking_struct.Fopt_previous;
                    variable_tracking_struct.isfirst_flg = 1;
                %                 C1n_bound = max([0,0],C1n_opt+[-OPT_out.output.meshsize, OPT_out.output.meshsize]);
%                 fprintf('C1 bounds(unormalized): [%d,%d]\n',C1n_bound(1)*objectiveStruct.parNormFactors(1),C1n_bound(2)*objectiveStruct.parNormFactors(1));
%                 fprintf('C1 bounds(normalized): [%d,%d]\n',C1n_bound(1),C1n_bound(2));
                % patternsearch on M1
                [M1n_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= patternsearch(@(t) objectiveFunctionIFEA([C1n_opt,t],objectiveStruct,0),Pn(2),[],[],[],[],[Pn(2)*1e-2],[Pn(2)*1e2],[],OPT_options);
%                 M1n_bound = max([0,0],M1n_opt+[-OPT_out.output.meshsize, OPT_out.output.meshsize]);
%                 fprintf('M1 bounds: [%d,%d]',M1n_bound(1)*objectiveStruct.parNormFactors(2),M1n_bound(2)*objectiveStruct.parNormFactors(2));
                M1n_temp_meshsize = OPT_out.output.meshsize;
                variable_tracking_struct.Fopt_previous = OPT_out.fval*variable_tracking_struct.Fopt_previous;
                variable_tracking_struct.isfirst_flg = 1;
                % patternsearch on K
%                 [Kn_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= patternsearch(@(t) objectiveFunctionIFEA([C1n_opt,M1n_opt,t],objectiveStruct,0),Pn(3),[],[],[],[],[Pn(3)*1e-2],[Pn(3)*1e2],[],OPT_options);
%                 Kn_temp_meshsize = OPT_out.output.meshsize;
%                 variable_tracking_struct.Fopt_previous = OPT_out.fval*variable_tracking_struct.Fopt_previous;
%                 variable_tracking_struct.isfirst_flg = 1;
%                 Kn_bound = Kn_opt+[-OPT_out.output.meshsize, OPT_out.output.meshsize];
             
%                 Kn_opt = 1;
                 
%                 [C1n_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= patternsearch(@(t) objectiveFunctionIFEA([t,M1n_opt,Kn_opt],objectiveStruct,1),C1n_opt,[],[],[],[],[C1n_bound(1)],[C1n_bound(2)],[],OPT_options);
% 
%                 [M1n_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= patternsearch(@(t) objectiveFunctionIFEA([C1n_opt,t,Kn_opt],objectiveStruct,0),M1n_opt,[],[],[],[],[M1n_bound(1)],[M1n_bound(2)],[],OPT_options);
                % patternsearch on C1 (#2)
                OPT_options = optimoptions(OPT_options,'InitialMeshSize', 0.5*C1n_temp_meshsize, 'MeshContractionFactor',0.25,'MeshExpansionFactor',1.5);  
                [C1n_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= patternsearch(@(t) objectiveFunctionIFEA([t,M1n_opt],objectiveStruct,1),C1n_opt,[],[],[],[],[Pn(1)*1e-2],[Pn(1)*1e2],[],OPT_options);
                variable_tracking_struct.Fopt_previous = OPT_out.fval*variable_tracking_struct.Fopt_previous;
                variable_tracking_struct.isfirst_flg = 1;

                OPT_options = optimoptions(OPT_options,'InitialMeshSize', 0.5*M1n_temp_meshsize, 'MeshContractionFactor',0.25,'MeshExpansionFactor',1.5);   
                % patternsearch on m1 (#2)
                [M1n_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= patternsearch(@(t) objectiveFunctionIFEA([C1n_opt,t],objectiveStruct,0),M1n_opt,[],[],[],[],[Pn(2)*1e-2],[Pn(2)*1e2],[],OPT_options);
                variable_tracking_struct.Fopt_previous = OPT_out.fval*variable_tracking_struct.Fopt_previous;
                variable_tracking_struct.isfirst_flg = 1;
                % patternsearch on K (#2)
%                 OPT_options = optimoptions(OPT_options,'InitialMeshSize', Kn_temp_meshsize, 'MeshContractionFactor',0.5,'MeshExpansionFactor',2);   
%                 [Kn_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= patternsearch(@(t) objectiveFunctionIFEA([C1n_opt,M1n_opt,t],objectiveStruct,0),Kn_opt,[],[],[],[],[Pn(3)*1e-2],[Pn(3)*1e2],[],OPT_options);
%                 variable_tracking_struct.Fopt_previous = OPT_out.fval*variable_tracking_struct.Fopt_previous;
%                 variable_tracking_struct.isfirst_flg = 1;
                % lsqnonlin to find local minimum
                objectiveStruct.method=4;
                maxNumberIterations = 5;
                OPT_options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
                maxNumberFunctionEvaluations = 50;
                FinDiffRelSteps = ones(size(Pn))*1e-4;
                OPT_options = optimoptions(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
                    'MaxIter',maxNumberIterations,...
                    'TolFun',functionTolerance,...       %'TolX',parameterTolerance,...
                    'Display',displayTypeIterations,...
                    'FinDiffRelStep',FinDiffRelSteps,... %"FiniteDifferenceStepSize" in non legacy
                    'ScaleProblem', 'jacobian',... %<<<<
                    'StepTolerance', 1e-8,... %<<<<
                    'DiffMaxChange',1e-1);
%                 Pn_opt = [C1n_opt,M1n_opt,Kn_opt];
                Pn_opt = [C1n_opt,M1n_opt];
                [Pn_opt,OPT_out.resnorm,OPT_out.residual,OPT_out.exitflag,OPT_out.output]= lsqnonlin(@(Pn) objectiveFunctionIFEA([Pn(1),Pn(2),Pn(1)*k_factor_tru],objectiveStruct,0.15),Pn_opt,[Pn_opt*1e-2],[Pn_opt*1e2],OPT_options);
                variable_tracking_struct.Fopt_previous = OPT_out.resnorm*variable_tracking_struct.Fopt_previous;
%                 variable_tracking_struct.isfirst_flg = 1;
                objectiveStruct.method = 3;

%                 C1n_bound
%                 M1n_bound
%                 Kn_bound
                %%%%%%%% convert to hybrid
%                 objectiveStruct.method = 4;
% 
%                 OPT_options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
%         %         maxNumberFunctionEvaluations = 100;
%                 FinDiffRelSteps = 1e-6;%[1e-4, 1e-4,1e-2];
%                 OPT_options = optimoptions(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
%                     'MaxIter',4,...
%                     'TolFun',functionTolerance,...       %'TolX',parameterTolerance,...
%                     'Display',displayTypeIterations,...
%                     'FinDiffRelStep',FinDiffRelSteps,... %"FiniteDifferenceStepSize" in non legacy
%                     'ScaleProblem', 'none',... %<<<<
%                     'StepTolerance', 1e-3);
%                 %             'StepTolerance', 1e-6)
%         %             'FiniteDifferenceStepSize', [0.01, 0.05, 1e-3],... %<<<
%         %         lambda_0=0.5;
%         %         Pn(end+1)=0.5;
%         %         Pn(end) = [];
%         %         [Pn_opt,OPT_out.resnorm,OPT_out.residual,OPT_out.exitflag,OPT_out.output]= lsqnonlin(@(Pn) objectiveFunctionIFEA([Pn(1),Pn(2),1],objectiveStruct),Pn,[Pn*1e-2],[Pn*1e2],OPT_options); %impose k_ini = k_tru
%                 [Pn_opt,OPT_out.resnorm,OPT_out.residual,OPT_out.exitflag,OPT_out.output]= lsqnonlin(@(Pn) objectiveFunctionIFEA([Pn(1),Pn(2) Kn_opt],objectiveStruct,0.15),[mean(C1n_bound),mean(M1n_bound),mean(Kn_bound)],[C1n_bound(1),M1n_bound(1)],[C1n_bound(2),M1n_bound(2)],OPT_options);
% %                 
%                 objectiveStruct.method = 3;

%                 OPT_options = optimoptions(OPT_options,'MaxIter', 10); 
%                 [C1n_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= patternsearch(@(t) objectiveFunctionIFEA([t,M1n_opt,1],objectiveStruct,1),C1n_opt,[],[],[],[],[C1n_bound(1)],[C1n_bound(2)],[],OPT_options);
%                 [M1n_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= patternsearch(@(t) objectiveFunctionIFEA([C1n_opt,t,1],objectiveStruct,0),M1n_opt,[],[],[],[],[M1n_bound(1)],[M1n_bound(2)],[],OPT_options);
                    %                     Pn_opt = [C1n_opt,M1n_opt,1];

%                 OPT_options=optimset('fminsearch'); % 'Nelder-Mead simplex direct search'
%                 OPT_options = optimset(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
%                     'MaxIter',maxNumberIterations,...
%                     'TolFun',functionTolerance,...
%                     'TolX',parameterTolerance,...
%                     'Display',displayTypeIterations,...
%                     'FinDiffRelStep',1e-3,...
%                     'DiffMaxChange',0.05);
%                 [Pn_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= fminsearch(@(Pn) objectiveFunctionIFEA(Pn,objectiveStruct),Pn,OPT_options);
%                 [Pn_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= fminbnd(@(Pn) objectiveFunctionIFEA(Pn,objectiveStruct),[Pn*1e-2],[Pn*1e2],OPT_options);

%                 [C1n_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= fmincon(@(t) objectiveFunctionIFEA([t,1,1],objectiveStruct,1),Pn(1),[],[],[],[],[Pn(1)*1e-2],[Pn(1)*1e2],[],OPT_options);
%                 [M1n_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= fmincon(@(t) objectiveFunctionIFEA([C1n_opt,t,1],objectiveStruct,0),Pn(2),[],[],[],[],[Pn(2)*1e-2],[Pn(2)*1e2],[],OPT_options);
%                 [C1n_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= fmincon(@(t) objectiveFunctionIFEA([t,M1n_opt,1],objectiveStruct,1),C1n_opt,[],[],[],[],[C1n_opt*0.5],[C1n_opt*1e2],[],OPT_options);
%                 [M1n_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= fmincon(@(t) objectiveFunctionIFEA([C1n_opt,t,1],objectiveStruct,0),M1n_opt,[],[],[],[],[M1n_opt*0.5],[M1n_opt*1e2],[],OPT_options);
                

            case 5 %lsqnonlin and Levenberg-Marquardt
                objectiveStruct.method=5;
                        maxNumberIterations = 10;
                OPT_options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
                maxNumberFunctionEvaluations = 122;
                FinDiffRelSteps = [1e-1, 1e-1,1e-1];
                OPT_options = optimoptions(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
                    'MaxIter',maxNumberIterations,...
                    'TolFun',functionTolerance,...       %'TolX',parameterTolerance,...
                    'Display',displayTypeIterations,...
                    'FinDiffRelStep',FinDiffRelSteps,... %"FiniteDifferenceStepSize" in non legacy
                    'ScaleProblem', 'none',... %<<<<
                    'StepTolerance', 1e-8,... %<<<<
                    'DiffMaxChange',1e-1,...
                    'FiniteDifferenceType', 'central');
                %             'StepTolerance', 1e-6)
        %             'FiniteDifferenceStepSize', [0.01, 0.05, 1e-3],... %<<<
        %         lambda_0=0.5;
        %         Pn(end+1)=0.5;
        %         Pn(end) = [];
        %         [Pn_opt,OPT_out.resnorm,OPT_out.residual,OPT_out.exitflag,OPT_out.output]= lsqnonlin(@(Pn) objectiveFunctionIFEA([Pn(1),Pn(2),1],objectiveStruct),Pn,[Pn*1e-2],[Pn*1e2],OPT_options); %impose k_ini = k_tru
                [Pn_opt,OPT_out.resnorm,OPT_out.residual,OPT_out.exitflag,OPT_out.output]= lsqnonlin(@(Pn) objectiveFunctionIFEA(Pn,objectiveStruct,0.2),Pn,[Pn*1e-2],[Pn*1e2],OPT_options);
%                 variable_tracking_struct.Fopt_previous = OPT_out.resnorm*variable_tracking_struct.Fopt_previous;
%                 variable_tracking_struct.isfirst_flg = 1;

        %         [Pn_opt,OPT_out.resnorm,OPT_out.residual,OPT_out.exitflag,OPT_out.output]= lsqnonlin(@(Pn) objectiveFunctionIFEA(Pn,objectiveStruct),Pn,[-Inf,-Inf,-Inf,0],[Inf,Inf,Inf,1],OPT_options);
        end
        %% Repeat analysis with optimal set
        [Fopt,OPT_stats_out]=objectiveFunctionIFEA(Pn_opt,objectiveStruct,0.2);
        %         Pn_opt(end) = [];
        %% Unnormalize and constrain parameters

        P_opt=Pn_opt.*objectiveStruct.parNormFactors; %Scale back, undo normalization

        %Constraining parameters
        % for q=1:1:numel(P_opt)
        %     [P_opt(q)]=boxconstrain(P_opt(q),objectiveStruct.Pb_struct.xxlim(q,1),objectiveStruct.Pb_struct.xxlim(q,2),objectiveStruct.Pb_struct.xx_c(q));
        % end

        disp_text=sprintf('%6.16e,',P_opt); disp_text=disp_text(1:end-1);
        disp(['P_opt=',disp_text]);
        %% animate the optimization process
        x_max = 0;
        y_max = 0;
        x_len_max = 0;
        y_len_max = 0;
        if ~ishandle(hf2)
            hf2 = figure;
        end
        for i=1:numel(variable_tracking_struct.displacement_sim_all)
            x_max = max(x_max, max(variable_tracking_struct.displacement_sim_all{i}));
            y_max = max(y_max, max(variable_tracking_struct.force_sim_all{i}));
            x_len_max = max(x_len_max,length(variable_tracking_struct.displacement_sim_all{i}));
            y_len_max = max(y_len_max,numel(variable_tracking_struct.force_sim_all{i}));
            num_iterations = i;
        end
        variable_tracking_struct.displacement_sim_all_new = nan(x_len_max,num_iterations);
        variable_tracking_struct.force_sim_all_new = nan(y_len_max,num_iterations);

        % for i=1:numel(variable_tracking_struct.displacement_sim_all)
        %     variable_tracking_struct.displacement_sim_all_new(:,i) = variable_tracking_struct.displacement_sim_all{i};
        %     variable_tracking_struct.force_sim_all_new(:,i) = variable_tracking_struct.force_sim_all{i};
        %     
        % %     variable_tracking_struct.force_sim_all_new.P_all_new = cell2mat(variable_tracking_struct.force_sim_all_new
        % end



        % xlim([0 max(max(variable_tracking_struct.displacement_sim_all_new))]);
        % ylim([0 max(max(variable_tracking_struct.force_sim_all_new))]);
        xlim([0 x_max]);
        ylim([0 y_max]);

        % num_iterations=sum(~isnan(variable_tracking_struct.P_all(:,1)));
        animStruct.Time=1:1:num_iterations; %The time vector
        for qt=1:1:num_iterations %Loop over time increments
            %Set entries in animation structure
            animStruct.Handles{qt}=[hl(2) hl(2) ht(3) ht(4) ht(5)]; %Handles of objects to animate
            animStruct.Props{qt}={'XData','YData','String','String','String'}; %Properties of objects to animate
            animStruct.Set{qt}={variable_tracking_struct.displacement_sim_all{qt},variable_tracking_struct.force_sim_all{qt},...
                ['c_s_i_m=' num2str(variable_tracking_struct.P_all{qt}(1),'%3.3e') ', m_s_i_m=' num2str(variable_tracking_struct.P_all{qt}(2),'%.2f') ', k_s_i_m=' num2str(variable_tracking_struct.P_all{qt}(3),'%.2f')],...
                ['F_o_p_t=' num2str(variable_tracking_struct.Fopt_all{qt},'%3.3e') ',{F_{opt}}^1=',num2str(variable_tracking_struct.F1opt_all{qt},'%3.3e'),',{F_{opt}}^2=',num2str(variable_tracking_struct.F2opt_all{qt},'%3.3e')],...
                ['Iteration number ' num2str(qt)]}; %Property values for to set in order to animate
            %         ['F_o_p_t=' num2str(variable_tracking_struct.Fopt_all{qt},'%3.3e')],...

        end
        anim8(hf2,animStruct); %Initiate animation feature
        drawnow;


        c1_data = [];
        m1_data = [];
        F_opt_data = [];
        figure;
        hold on
        for qt=1:1:num_iterations %Loop over time increments
            c1_data(qt) = variable_tracking_struct.P_all{qt}(1);
            m1_data(qt) = variable_tracking_struct.P_all{qt}(2);
            F_opt_data(qt) = variable_tracking_struct.Fopt_all{qt}(1);
        end

        [~, opt_ind] = min(F_opt_data);
        C_map = colormap(jet(num_iterations));
        for i=1:num_iterations-1
            plot3(c1_data(i:i+1),m1_data(i:i+1), F_opt_data(i:i+1), '-', 'Color', C_map(i,:), 'HandleVisibility', 'Off');
        end
        xlabel('c1_i');
        ylabel('m1_i');
        zlabel('F_{opt}_i');
        % P_opt=Pn_opt.*objectiveStruct.parNormFactors;
        % colormap(C_map);
        scatter_size = 36;
        scatter3(c1_data,m1_data, F_opt_data,scatter_size*ones(1,num_iterations),C_map, 'HandleVisibility' , 'Off'); % All points
        scatter3(c1_data(opt_ind),m1_data(opt_ind), F_opt_data(opt_ind),1.5*scatter_size,C_map(end,:),'filled', 'DisplayName', 'Optimal'); %Optimal point
        scatter3(c1_data(1),m1_data(1), F_opt_data(1),1.5*scatter_size,C_map(1,:),'filled', 'DisplayName', 'Innitial'); %Initial point
        scatter3(c1_tru,m1_tru,0, 'x', 'DisplayName', 'Tru Values'); % tru point
        grid on;
        legend show;
        caxis([0 num_iterations])
        cbar = colorbar;
        cbar.Title.String = 'Itteration #';
        title(OPT_out.output.message);
        disp(['Optimum of objective Fcn: ', num2str(F_opt_data(opt_ind))]);
        disp(['Optimal C1: ', num2str(c1_data(opt_ind))]);
        disp(['Optimal m1: ', num2str(m1_data(opt_ind))]);
        disp(['Optimal k: ', num2str(variable_tracking_struct.P_all{opt_ind}(3))]);
        ax1 = gca;
        ax1.ZScale = 'log';
        zlim([0, 1.01]);
        t =xline(c1_tru, 'k'); yline(m1_tru, 'k'); 

        % disp(['C1_error', num2str(abs((c1_data(opt_ind)-c1_tru)/c1_tru))]);
        % disp(['m1_error', num2str(abs((c1_data(opt_ind)-c1_tru)/c1_tru))]);
        % disp(['m1_error', num2str(abs((c1_data(opt_ind)-c1_tru)/c1_tru))]);
        figure;
        semilogy(F_opt_data);
        xlabel('k (itteration)');
        ylabel('F_{opt}^k')

end

