function figH = MakePlots(test,mat_type,eta_arr,dir_name,normalize_param,colormap_data_field,F_pos_version,save_figures,fig_save_path)  
%% 
desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
path_parts = strsplit(dir_name,filesep);
% group_name = sprintf('(%s)_%s_%s_wieght_norm=%d_param_norm_%d',mat_type,path_parts{end-1},path_parts{end},F_pos_version,normalize_param);
group_name = sprintf('%s_wieght_norm=%d_param_norm_%d',mat_type,F_pos_version,normalize_param);
myGroup = desktop.addGroup(group_name);
desktop.setGroupDocked(group_name, 0);
myDim   = java.awt.Dimension(1, 2);  
% 1: Maximized, 2: Tiled, 3: Floating
desktop.setDocumentArrangement(group_name, 1, myDim)
figH    = gobjects(1, 2);
bakWarn = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

warning(bakWarn);

% finite differences options and data
    global finite_diff_data
    finite_diff_data.n = 8; % truncation error order O(h^n)
    finite_diff_data.leap = 1; % distances between adjacent data points for calculation differences
    [finite_diff_data.N finite_diff_data.D finite_diff_data.B finite_diff_data.K,finite_diff_data.W_xy] = getXYCoeffMat(finite_diff_data.n);
    finite_diff_data.Wtt = secondDerWeights(finite_diff_data.n);
%% Original code
%     normalize_param = 0;
    ax = [];
    figH(1) = figure('WindowStyle', 'normal','NumberTitle', 'off','SizeChangedFcn',@fig1SizeChangeCBack,'units','centimeters');
    figH(1).Position(2) = 0;
    figH(1).Position(3) = [18.1]; figH(1).Position(4) = 17; %[15.6]
    
    axis_font_size = 10;
    label_multiplyer = 1.1;
    title_multiplyer = 1.1;
    axis_label_font_size = axis_font_size*label_multiplyer;
    title_font_size =axis_font_size*title_multiplyer;
    cb_title_font_size = title_font_size;
    tile_layout_font_size = title_font_size;
%     p = uipanel('Position',[.1 .1 .8 .8]);
%     figH(1).Scrollable = 'on'
    N_times = length(test{1}.indenter_RB_out.time)-1; % number of indentation steps not including delta=0
%     N_times = 1;
%     t = tiledlayout(N_times,length(eta_arr), 'Padding', 'compact',...
%         'TileSpacing','none','Units', 'centimeters','OuterPosition',[0.1 0.1 17.9 15.4]);
    t = tiledlayout(N_times,length(eta_arr), 'Padding', 'compact',...
        'TileSpacing','none','Units', 'normalized','OuterPosition',[0 0 1 1]);
    t.InnerPosition(1) = 0.1; t.InnerPosition(3) = 0.8;
    for kk=2:(N_times+1)
        for eta_ind = 1:length(eta_arr) 
            temp_eta = eta_arr(eta_ind);
            P1_arr = [];
            P2_arr = [];
            D1_arr = [];
            D2_arr = [];
            min_OBJ_VAL = Inf;
            
            
            
            
            %% Calculate datapoints coordinates and values for graphic objects 
            N_p1 = test{end}.P1_ind;
            N_p2 = test{end}.P2_ind;
            X = zeros(N_p2,N_p1); Y = zeros(N_p2,N_p1); % Coordinate grids for plotting graphics
            Z = zeros(N_p2,N_p1); % Grid of datapoint values
            for test_ind=1:numel(test)
                Ff =  test{test_ind}.obj_fun_val.Ff;
                Fu = test{test_ind}.obj_fun_val.Fu;
                Ffu = temp_eta*Ff+(1-temp_eta)*Fu;
                if N_times==1 % show only one row of contour plots (average over delta)
                   OBJ_VAL = mean(Ffu);
                else
                   OBJ_VAL = Ffu(kk); % show contour plots at multiple indentation depths delta
                end
                Z(test{test_ind}.P2_ind,test{test_ind}.P1_ind) = (OBJ_VAL);
                
                if OBJ_VAL<min_OBJ_VAL % find the index of the global minimum in each plot
                    optimal_ind(kk-1,eta_ind) = test_ind; % index of global minimum
                    min_OBJ_VAL = OBJ_VAL; % global minimum
                end
                
                P1_arr = sort(unique([P1_arr test{test_ind}.p1]));
                P2_arr = sort(unique([P2_arr test{test_ind}.p2]));
                
                X(test{test_ind}.P2_ind,test{test_ind}.P1_ind) = test{test_ind}.p1;
                Y(test{test_ind}.P2_ind,test{test_ind}.P1_ind) = test{test_ind}.p2;
            end
            optimal_OBJ_VAL(kk-1,eta_ind) = min_OBJ_VAL; % global minimum
            optimal_params{kk-1,eta_ind} = [test{optimal_ind(kk-1,eta_ind)}.p1, test{optimal_ind(kk-1,eta_ind)}.p2];
            optimal_params_temp = optimal_params{kk-1,eta_ind};
           
            if normalize_param
%                 [X,Y] = meshgrid(P1_arr./optimal_params_temp(1),P2_arr./optimal_params_temp(2));
                X = X/optimal_params_temp(1);
                Y = Y/optimal_params_temp(2);
            else % convert to easier untis
                % Scale axis for showing KPa instead of GPA
                if any(strcmp(mat_type, {'OG','OM','MR','NH-Lamme','NH'}))
                    X = X*1000;
                    optimal_params_temp(1) = optimal_params_temp(1)*1000;
                end
                if any(strcmp(mat_type, {'MR','NH-Lamme'}))%scale also P2
                    Y = Y*1000;
                    optimal_params_temp(2) = optimal_params_temp(2)*1000;
                end
            end
            
            
%% old
                %                 if N_times==1
% %                     F_force = test{test_ind}.obj_fun_val.mean(test{test_ind}.obj_fun_val.ForceDevMagNormalized(2:end));
%                     F_force = mean(F_force);
% %                 else 
% %                     F_force = test{test_ind}.obj_fun_val.ForceDevMagNormalized(kk);
%                 end
%                 
%                 
%                 switch F_pos_version
%                     case 1
%                         if N_times==1
%                             F_pos = mean(test{test_ind}.obj_fun_val.nodeWeightedPositionDevMag1(2:end));
%                         else
%                             F_pos = test{test_ind}.obj_fun_val.nodeWeightedPositionDevMag1(kk);
%                         end
%                     case 2
%                         if N_times==1
%                             F_pos = mean(test{test_ind}.obj_fun_val.nodeWeightedPositionDevMag2(2:end));
%                         else
%                             F_pos = test{test_ind}.obj_fun_val.nodeWeightedPositionDevMag2(kk);
%                         end
%                 end
%                 OBJ_VAL = temp_eta*sqrt(test{test_ind}.obj_fun_val.F1)+(1-temp_eta)*sqrt(test{test_ind}.obj_fun_val.F2); %sqrt
%                 OBJ_VAL = temp_eta*test{test_ind}.obj_fun_val.F1+(1-temp_eta)*test{test_ind}.obj_fun_val.F2; %no sqrt

                
%                 OBJ_VAL = temp_eta*F_force+(1-temp_eta)*F_pos; %no sqrt
%                     OBJ_VAL = temp_eta*sqrt(F_force)+(1-temp_eta)*sqrt(F_pos); % sqrt

%                                         OBJ_VAL = sqrt(temp_eta*F_force+(1-temp_eta)*F_pos); %no sqrt
%                     if any(strcmp(mat_type, {'NH-Lamme'}))
%                         Z(test{test_ind}.P2_ind,test{test_ind}.P1_ind*test{test_ind}.P2_ind) = OBJ_VAL;
%                     else
                
%                     end
%%
          
%% old
%             for i=1:numel(test)
%                 switch mat_type
%                     case 'NH-Lamme' % irrelevant
%                         mu_arr = [2:2:38]*1e-3;
%                         lambda_arr = linspace(1.8462,1862,19)*1e-3;
%                         test{i}.mu = mu_arr(test{i}.P1_ind);
%                         test{i}.lambda = lambda_arr(test{i}.P2_ind);
%                         P1_arr = sort(unique([P1_arr test{i}.mu]));
%                         P2_arr = sort(unique([P2_arr test{i}.lambda]));
%                         optimal_params_temp(1) = mu_arr(test{optimal_ind(kk-1,eta_ind)}.P1_ind);
%                         optimal_params_temp(2) = lambda_arr(test{optimal_ind(kk-1,eta_ind)}.P2_ind);
%                          
%                     otherwise
%                         P1_arr = sort(unique([P1_arr test{i}.p1]));
%                         P2_arr = sort(unique([P2_arr test{i}.p2]));
%                 end
%             end
            %% draw settings
            alpha = 0.8;
            draw_grad_arrows = 0;
            %% draw
            [dfdx,dfdy] = gradient(Z); % gradient 
            grad_mag = sqrt(dfdx.^2+dfdy.^2); % gradient magnitude
            ax{kk-1,eta_ind} = nexttile();
            hold on;
            switch colormap_data_field
                case 'grad_mag_surface'
                % norm(grad(fval))
                    c_h = surf(X,Y,Z,grad_mag, 'FaceAlpha', alpha, 'EdgeColor', 'none');
                    fig_title_str = sprintf('Colormap data field: %s','Objective function gradient magnitude');
                case 'Fval_surface'
                % Fval surface
                    s_h = surf(X,Y,Z, 'FaceColor', 'interp','FaceAlpha',alpha, 'EdgeColor', 'k');
                    fig_title_str = sprintf('Colormap data field: %s','Objective function value');
                case 'Fval_contour'
                % Fval contour
                    [~,c_h{kk-1,eta_ind}] = contourf(ax{kk-1,eta_ind},X,Y,Z,'Tag','c_h');
                    fig_title_str = sprintf('Colormap data field: %s','Objective function value');
%                      [c_h2{kk-1,eta_ind}] = scatter3(reshape(X,[],1),reshape(Y,[],1),'MarkerEdgeColor','k','SizeData',2);
%                      [c_h2{kk-1,eta_ind}] = scatter3(reshape(X,[],1),reshape(Y,[],1),reshape(Z,[],1), 2,reshape(Z,[],1));
                case 'Fval_contour3'
                % Fval contour
                    [~,c_h{kk-1,eta_ind}] = contour3(X,Y,Z,10);
                    fig_title_str = sprintf('Colormap data field: %s','Objective function value');
                case 'grad_mag_contour'
                % norm(gradient) contour
                    [~,c_h{kk-1,eta_ind}] = contourf(X,Y,grad_mag,50);
                    fig_title_str = sprintf('Colormap data field: %s','Objective function gradient magnitude');
                case 'scatter3Fval'
                     [c_h{kk-1,eta_ind}] = scatter3(reshape(X,[],1),reshape(Y,[],1),reshape(Z,[],1), 6,reshape(Z,[],1));
                    fig_title_str = sprintf('Colormap data field: %s','Objective function value');
            end
            if draw_grad_arrows
                quiver(X,Y,-dfdx,-dfdy, 'r');
            end
            axis tight
            if kk==2
                title(sprintf('\\eta=%g',eta_arr(eta_ind)));
            end



            axis_metadata.eta = temp_eta;
            axis_metadata.T = kk;
            axis_metadata.i = kk-1;
            axis_metadata.j = eta_ind;
            ax{kk-1,eta_ind}.UserData = axis_metadata;
            if normalize_param
%                     xlabel('{\boldmath$\frac{c}{c^\ast}$}', 'interpreter', 'latex');
%                     ylabel('{\boldmath$\frac{m}{m^\ast}$}', 'interpreter', 'latex');
                true_scatter_h{kk-1,eta_ind} = scatter(test{optimal_ind(kk-1,eta_ind)}.p1/optimal_params_temp(1),test{optimal_ind(kk-1,eta_ind)}.p2/optimal_params_temp(2),'MarkerFaceColor','w', 'MarkerEdgeColor', 'k', 'DisplayName', '"True" parameters','Tag','True_param_scatter');
%                     axis equal
            else
%                     xlabel('{\boldmath$c$}', 'interpreter', 'latex');
%                     ylabel('{\boldmath$m$}', 'interpreter', 'latex');
                true_scatter_h{kk-1,eta_ind} = scatter(optimal_params_temp(1),optimal_params_temp(2),'MarkerFaceColor','w', 'MarkerEdgeColor', 'k', 'DisplayName', '"True" parameters','Tag','True_param_scatter');
               
            end
            xlim(ax{kk-1,eta_ind},[min(X,[],'all'),max(X,[],'all')]);
            ylim(ax{kk-1,eta_ind},[min(Y,[],'all'),max(Y,[],'all')]);
        end
    end
    % determine and assign proper title for figure
%     switch colormap_data_field
%         case 'grad_mag_surface'
%             fig_title_str = sprintf('Colormap data field: %s','Objective function gradient magnitude');
%         case 'Fval_surface'
%             fig_title_str = sprintf('Colormap data field: %s','Objective function value');
%         case 'Fval_contour'
%             fig_title_str = sprintf('Colormap data field: %s','Objective function value');
%         case 'Fval_contour3'
%         % Fval contour
%             fig_title_str = sprintf('Colormap data field: %s','Objective function value');
%         case 'grad_mag_contour'
%         % norm(gradient) contour
%             fig_title_str = sprintf('Colormap data field: %s','Objective function gradient magnitude');
%         case 'scatter3Fval'
%             fig_title_str = sprintf('Colormap data field: %s','Objective function value');
%     end
    set(gcf, 'Name', sprintf('%s',fig_title_str));
    % add colorbar
    cb = colorbar;
    cb.Layout.Tile = 'east';
    cb.Title.Interpreter = 'latex';
    cb.Title.String ='$\bar{F}_{fu}$';
    cb.Title.FontSize = cb_title_font_size;
    cb.ButtonDownFcn = @cb_callback;
   %% Adjust contourf plots, run after opening a figure
    % Controls
    min_val = 0;
    switch F_pos_version
        case 1
            max_val = 0.5; % max value in auto mode is 2.25, can lower this to truncate higher values
        case 2 
            max_val = 0.5; % max value in auto mode is 2.25, can lower this to truncate higher values
    end
    df = (max_val-min_val)/10; % increment between isolines

    
    %%%%%%%% prevrious
%     levels_arr = 0:df:max_val; 
%     levels_arr_cbar = [0:df:(max_val+df)];
%     cb_ticks = [levels_arr_cbar+df/2].^2;
%     levels_arr = levels_arr.^2;
%     levels_arr_cbar = levels_arr_cbar.^2;
%     
%     r = diff(levels_arr_cbar)./levels_arr_cbar(end);
%     cmap = parula(length(levels_arr));
%     new_cmap = [];
%     
%     for ii=1:length(r)
%         temp_map = ones(round(100*length(levels_arr)*r(ii)),1)*cmap(ii,:);
%         new_cmap = [new_cmap; temp_map];
%         y_tick_str{ii} = sprintf('%d%%-%d%%', df*100*(ii-1),df*100*ii);
%     end
%     y_tick_str{end} = sprintf('>%d%%',max_val*100);
%     colormap(new_cmap);
%     cb.Ticks = cb_ticks;
%     cb.TickLabels =y_tick_str;
%     
    
    
    %% Custom colormap and colorbar
    %%%%%%%% new
    levels_arr = 0:df:max_val; 
    levels_arr_cbar = [0:df:(max_val+df)];
    cb_ticks = [levels_arr_cbar+df/2].^2;
    levels_arr = levels_arr.^2;
    levels_arr_cbar = levels_arr_cbar.^2;
    
    r = diff(levels_arr_cbar)./levels_arr_cbar(end);
    cmap = parula(length(levels_arr));
    new_cmap = [];
    
    for ii=1:length(r)
        temp_map = ones(round(100*length(levels_arr)*r(ii)),1)*cmap(ii,:);
        new_cmap = [new_cmap; temp_map];
        y_tick_str{ii} = sprintf('%d', df*100*(ii-1));
    end
    y_tick_str{1} = [];
    y_tick_str{end+1} = sprintf('>%d',max_val*100);
    colormap(new_cmap);
    cb.Ticks = levels_arr_cbar;
    cb.TickLabels =y_tick_str;
    cb.Title.String ='$\bar{F}_{fu}$ (\%)';
    %% Arrange axis
    for i=1:size(ax,1)
        for j=1:size(ax,2) 
            optimal_params_temp = optimal_params{i,j};
            % Scale axis for showing KPa instead of GPA
            if any(strcmp(mat_type, {'OG','OM','MR','NH-Lamme','NH'}))
                P1_arr = P1_arr*1000;
                optimal_params_temp(1) = optimal_params_temp(1)*1000;
            end
            if any(strcmp(mat_type, {'MR','NH-Lamme'}))%scale also P2
                P2_arr = P2_arr*1000;
                optimal_params_temp(2) = optimal_params_temp(2)*1000;
            end
            ax{i,j}.FontSize = axis_font_size;
        %     xlabel(ax(i),'$c/c^\ast$', 'interpreter', 'latex');
            %% Set yabel
            if j==1
                switch mat_type
                    case {'OG','OM'}
                        if normalize_param
                            ylabel(ax{i,j},'$m/m^\ast$', 'interpreter', 'latex');
                        else
                            ylabel(ax{i,j},'$m$', 'interpreter', 'latex');
                        end
                    case 'MR' 
                         if normalize_param
                            ylabel(ax{i,j},'$C_{01}/C_{01}^\ast$', 'interpreter', 'latex');
                        else
                            ylabel(ax{i,j},'$C_{10}$ (KPa)', 'interpreter', 'latex');
                         end
                     case 'NH-Lamme' 
                         if normalize_param
                            ylabel(ax{i,j},'$\lambda/\lambda^\ast$', 'interpreter', 'latex');
                        else
                            ylabel(ax{i,j},'$\lambda$', 'interpreter', 'latex');
                         end
                     case 'NH' 
                         if normalize_param
                            ylabel(ax{i,j},'$\nu/\nu^\ast$', 'interpreter', 'latex');
                        else
                            ylabel(ax{i,j},'$\nu$', 'interpreter', 'latex');
                        end
                end
            else 
                ax{i,j}.YLabel = [];
                ax{i,j}.YTickLabel = [];
            end
            %% Set xlabel
            if i<size(ax,1)
                xlabel(ax{i,j},[]);
                ax{i,j}.XTickLabel = [];
            else
                switch mat_type
                    case {'OG','OM'}
                        if normalize_param
                            xlabel(ax{i,j},'$c/c^\ast$', 'interpreter', 'latex');
                        else
                            xlabel(ax{i,j},'$c$ (KPa)', 'interpreter', 'latex');
                        end
                    case 'MR' 
                         if normalize_param
                            xlabel(ax{i,j},'$C_{10}/C_{10}^\ast$', 'interpreter', 'latex');
                        else
                            xlabel(ax{i,j},'$C_{10}$ (KPa)', 'interpreter', 'latex');
                         end
                     case 'NH-Lamme' 
                         if normalize_param
                            xlabel(ax{i,j},'$\mu/\mu^\ast$', 'interpreter', 'latex');
                        else
                            xlabel(ax{i,j},'$\mu$ (KPa)', 'interpreter', 'latex');
                         end
                     case 'NH' 
                         if normalize_param
                            xlabel(ax{i,j},'$E/E^\ast$', 'interpreter', 'latex');
                        else
                            xlabel(ax{i,j},'$E$ (KPa)', 'interpreter', 'latex');
                         end
                end
            end
        %     colorbar(ax(i),'off');
%             c_h{i,j} = findobj(ax{i,j}.Children,'Type','contour'); % handle for countourf object.
            %% update plots
            c_h{i,j}.LevelList = levels_arr; % assign new LevelList
            ax{i,j}.CLim = [levels_arr_cbar(1) levels_arr_cbar(end)]; % update Color limmit range
%             c_h{i,j}.LineColor = 'k'; %'k' % hide isolines
%             c_h{i,j}.Fill = 'off';
            ax{i,j}.TickLength(1) = 0.025;
            ax{i,j}.LineWidth = 1;
            ax{i,j}.Layer = 'top';
        %     ax(i).XTickLabelMode = 'manual'
            ax{i,j}.XLabel.FontSize = axis_label_font_size;
        %     ax(i).XTick = ax(i).XLim(1):0.2:ax(i).XLim(2);
        %     ax(i).XTick = 0:0.2:2;
            ax{i,j}.Title.FontSize = title_font_size;
            ax{i,j}.YLabel.FontSize = axis_label_font_size;
        %     ax(i).YTick = 0.2:0.2:1.8
            ax{i,j}.XTickLabelRotation = 30;
            ax{i,j}.YTickLabelRotation = 30;
            ax{i,j}.Box = 'on';
            axis(ax{i,j},'square');
    %         ax{i,j}.XLim = [1 15]/8
    
    
            %% Hessian calculation
            Hessian_temp = getHessian(c_h{i,j}.XData,c_h{i,j}.YData,c_h{i,j}.ZData,optimal_ind(i,j),test);
%             Hessian_temp2 = getHessian(c_h{i,j}.XData,c_h{i,j}.YData,c_h{i,j}.ZData,optimal_ind(i,j),test,4);
%             Hessian_temp3 = getHessian(c_h{i,j}.XData,c_h{i,j}.YData,c_h{i,j}.ZData,optimal_ind(i,j),test,8);
        
            c_h{i,j}.UserData.Hessian = Hessian_temp;
            c_h{i,j}.UserData.Hessian_tilde = [1 Hessian_temp(1,2)/sqrt(Hessian_temp(1,1)*Hessian_temp(2,2));
                                              Hessian_temp(2,1)/sqrt(Hessian_temp(2,2)*Hessian_temp(1,1)) 1];
%             uncertainty_limits = getUncertaintyLimits(c_h{i,j}.XData,c_h{i,j}.YData,c_h{i,j}.ZData,optimal_ind(i,j),test)
            %% uncertainty limits
            if normalize_param
                uncertainty_limits = getUncertaintyLimits(c_h{i,j},[1,1]);
            else
                uncertainty_limits = getUncertaintyLimits(c_h{i,j},optimal_params_temp);
            end
            c_h{i,j}.UserData.uncertainty_limits = uncertainty_limits;
            %% Plotting principle curvatures directions and marks for fin-diff
            [V,D] = eigs(c_h{i,j}.UserData.Hessian);
            dx = (ax{i,j}.XLim(2)-ax{i,j}.XLim(1));
            dy = (ax{i,j}.YLim(2)-ax{i,j}.YLim(1));
            quiver_r = 0.2;
            alpha_1 = quiver_r/sqrt((V(1,1)/dx)^2+(V(2,1)/dy)^2);
            alpha_2 = quiver_r/sqrt((V(1,2)/dx)^2+(V(2,2)/dy)^2);
            p1_scale = 1;%(ax{i,j}.XLim(2)-ax{i,j}.XLim(1))/10;
            p2_scale = 1;%;(ax{i,j}.YLim(2)-ax{i,j}.YLim(1))/10;
            
            
            q_scale_factor = sqrt(dx^2+dy^2)/10;
            p1_scale = q_scale_factor; p2_scale = q_scale_factor;
            if normalize_param
                xc_quiver = 1;
                yc_quiver = 1;
            else
                xc_quiver = optimal_params_temp(1);
                yc_quiver = optimal_params_temp(2);
            end
            q1 = quiver(ax{i,j},xc_quiver,yc_quiver,V(1,1)*alpha_1,V(2,1)*alpha_1,1,...
                'ShowArrowHead', 'off');
            q1_r = quiver(ax{i,j},xc_quiver,yc_quiver,-V(1,1)*alpha_1,-V(2,1)*alpha_1,1,...
                'ShowArrowHead', 'off','Color',q1.Color);
            q2 = quiver(ax{i,j},xc_quiver,yc_quiver,V(1,2)*alpha_2,V(2,2)*alpha_2,1,...
                'ShowArrowHead', 'off');
            q2_r = quiver(ax{i,j},xc_quiver,yc_quiver,-V(1,2)*alpha_2,-V(2,2)*alpha_2,1,...
                'ShowArrowHead', 'off','Color', q2.Color);
            
            ax{i,j}.UserData.fin_diff_marks = markPointsForDerivatives(ax{i,j},c_h{i,j}.XData,c_h{i,j}.YData,optimal_ind(i,j),test);
            
            ax{i,j}.UserData.quiver_h = [q1 q1_r q2 q2_r];
%             ax{i,j}.Subtitle. Interpreter = 'latex';
            %% Hessians data textbox
            show_H_data = 0;
            if show_H_data
                text_box{i,j} = makeTextbox(ax{i,j},c_h{i,j}); set(text_box{i,j},'Visible', 'Off');
            else
                text_box{i,j} = [];
            end
            %% Quadratic approximation fcontours
            if normalize_param
                theta_fun = @(x,y) ([x;y]-[1;1]);
            else
                theta_fun = @(x,y) ([x;y]-optimal_params_temp(1:2)');
            end
            ellipsoid_contour{i,j} = fcontour(ax{i,j}, @(x,y) optimal_OBJ_VAL(i,j)+0.5*dot((theta_fun(x,y)'*Hessian_temp)',theta_fun(x,y)),'--', 'LevelList', optimal_OBJ_VAL(i,j)+levels_arr, 'Visible', 'off','tag','Contours');
            ellipsoid_contour_white{i,j} = fcontour(ax{i,j}, @(x,y) optimal_OBJ_VAL(i,j)+0.5*dot((theta_fun(x,y)'*Hessian_temp)',theta_fun(x,y)),'--w', 'LevelList', optimal_OBJ_VAL(i,j)+0.05^2, 'Visible', 'on','Tag','WhiteLevels');
        end
    end
%     t.TileSpacing = 'tight';
    set(t.YLabel,'Interpreter', 'latex','String',{'$\leftarrow$ Indentation depth, $\delta$'},'FontSize', tile_layout_font_size);
%     set(t.XLabel,'Interpreter', 'latex','String',{'$\eta \rightarrow$'},'FontSize', tile_layout_font_size);

    set(t.Title, 'Interpreter', 'latex','String',{['\textbf{(a)} Objective function for the ',mat_type,' model']},...
        'FontSize', tile_layout_font_size,'Color', 'k');
    if ~isequal(text_box,cell(size(text_box)))
        adjustTextbox(ax,text_box);
    end
%     rect = annotation(gcf,'rectangle', t.Position);
%     rect_out = annotation(gcf,'rectangle', t.OuterPosition);

%     addTransperacy(c_h,0.5); %add transperacy to contour plots (optional) [https://www.mathworks.com/matlabcentral/answers/514830-adjusting-the-transparency-of-a-contour-plot-using-a-gradient-of-alpha-values]
%     figH(1).UserData.rect = rect;
%     figH(1).UserData.rect_out = rect_out;
    figH(1).UserData.ax = ax;
    figH(1).UserData.t = t;
    figH(1).UserData.c_h = c_h;
    figH(1).UserData.ellipsoid_contour = ellipsoid_contour;
    figH(1).UserData.ellipsoid_contour_white = ellipsoid_contour_white;
    figH(1).UserData.text_box = text_box;
    figH(1).UserData.cb = cb;
    figH(1).UserData.mat_type = mat_type;
    figH(1).UserData.optimal_params = optimal_params;
    figH(1).UserData.optimal_params_temp = optimal_params_temp;
    figH(1).UserData.eta_arr = eta_arr;
    figH(1).UserData.optimal_ind = optimal_ind;
    figH(1).UserData.dir_name = dir_name;
    figH(1).UserData.normalize_param = normalize_param;
    figH(1).UserData.colormap_data_field = colormap_data_field;
    figH(1).UserData.F_pos_version = F_pos_version;
    figH(1).UserData.indentation_depths = test{1}.pos_out.z.data(intersect(find(test{1}.pos_out.y.data(:,2)==0),find(test{1}.pos_out.x.data(:,2)==0)),:);
    figH(1).UserData.fig_type = 'ContourPlots';
%     fig_name = sprintf('(%s)_%s_%s_wieght_norm=%d_param_norm_%d',mat_type,path_parts{end-1},path_parts{end},F_pos_version,normalize_param);
    fig_name = strcat('Contour plots [',group_name,']');
    figH(1).UserData.fig_name = fig_name;
    figH(1).UserData.group_name = group_name;
    figH(1).Name = fig_name;
    adjustPlotsToolbar(figH(1));

    %% second figure (heatmaps)
    figH = open_stat_fig(figH)
    
    drawnow;
    pause(0.02);
%     set(get(handle(figH(2)), 'javaframe'), 'GroupName', group_name);
   figH(2).UserData.figH1 = figH(1);
   figH(2).UserData.fig_type = 'HessianData';
   figH(2).UserData.fig_name = figH(2).Name;
    %% save figures
    figH(1).UserData.t.XLabel.String = [];
    if strcmp(save_figures, 'yes')
        for ii=1:length(figH)
            set(figH(ii), 'WindowStyle', 'normal','Renderer', 'painters')
            pause(1);
            brighten(figH(ii), 0.2);
            %%%%%%%% This must be done after(!) brighten() is called:
            condH_ax = findobj(figH(ii),'Tag','condH','Type','Axes');
            if ~isempty(condH_ax)
                condH_ax.Colormap = flip(condH_ax.Colormap,1); 
            end
            %%%%%%%%
            %%axes
            if ii==1
                if normalize_param
                    x_tick_arr = [0.5 1,1.5];
                    switch mat_type
                        case 'NH' 
                            y_tick_arr =[0.8 1,1.2];
                        otherwise % {'OG,'OM','MR','NH-Lamme','NH'}
                             y_tick_arr = x_tick_arr;
                    end
                    for kk = 1:numel(figH(ii).UserData.ax)
                        set(figH(ii).UserData.ax{kk},'YTickLabelRotation',0,'YTick', y_tick_arr,'XTickLabelRotation',0,'XTick',x_tick_arr);
                    end
                end
            end
            %% ADJUST PLOTS
            lineWidth = 1;
            level = 0.1^2;
            set(findobj(figH(ii),'Tag','WhiteLevels'),'LineColor', 'm','LineWidth',lineWidth*0.8,'LevelList',level);
            set(findobj('Tag','HighlightLevels'),'LineColor', 'w','LineWidth',lineWidth*0.8,'LevelList',level);
            %% QUIVER
            set(findobj('type', 'quiver'),'Visible', 'On','LineWidth',lineWidth*1.5,'linestyle',':','Color','g');
            %% Contours
            true_param_scatter = findobj('Tag','True_param_scatter');
            for mm=1:length(true_param_scatter) 
                uistack(true_param_scatter(mm),'top'); 
            end
            t = figH(1).UserData.t;
            if ii==1
                if isequal(t.GridSize,[1 1])
                    t.YLabel.String = [];
                    t.Title.String = [];
                    ax{1}.Title.String = [];
                    figH(1).UserData.cb.Title.String = '$\bar{F}_{f}^{tot}$ (\%)';
                    figH(1).Position(3) = 8.1;
                    figH(1).Position(4) = 7.3;
                    drawnow();
                    string_lbl = inputdlg('caption index','get image caption label',1,{'(a)'});
                    subfig_lbl = annotation(figH(1),'textbox','String', string_lbl{1}, 'FontSize',8,'FontName','Times New Roman',...
                            'VerticalAlignment','bottom','HorizontalAlignment','center','Margin',0,'EdgeColor', 'none');
                    subfig_lbl.Position(1:2) = [0.5, 0]; subfig_lbl.Position(1) = subfig_lbl.Position(1)-0.5*subfig_lbl.Position(3);
                else
                   figH(1).UserData.t.XLabel.String = [];
                end
            end
            drawnow(); 
            %% adjust tile innerposition (remove spaces and center the plot)
            if ii==1
                removeSpaces(figH(ii).UserData.t,figH(ii))
            end
            %% save
    %       print('-painters',figs(ii).fig_h,fullfile(selpath,figs(ii).f_name(1:end-4)),'-dsvg');
            print('-painters',figH(ii),fullfile(fig_save_path,figH(ii).UserData.fig_name),'-dsvg');
            savefig(figH(ii),fullfile(fig_save_path,group_name));
    %         hgsave(figH(1).fullfile(fig_save_path,fig_name));    
        end
    end
end


%% Auxiliary functions 
function [mu,lambda] = Young2Lamme(E,nu)
    mu = E./(2*(1+nu));
    lambda = nu*E./((1+nu).*(1-2*nu));
end

% caluclates hessian from objective function evaluations (Z matrix of contourf)
function [H] = getHessian(X,Y,Z,optimal_ind,test)
    global finite_diff_data
    ii = test{optimal_ind}.P2_ind;
    jj = test{optimal_ind}.P1_ind;
    
    r = finite_diff_data.n/2;
    leap = finite_diff_data.leap;

    
    try
        h_p1 = X(ii,jj+leap)-X(ii,jj); % delta_p1
        h_p2 = Y(ii+leap,jj)-Y(ii,jj); % delta_p2
        
        H11 = Z(ii,(jj-r*leap):leap:(jj+r*leap))*finite_diff_data.Wtt'/(h_p1^2);
        H22 = Z((ii-r*leap):leap:(ii+r*leap),jj)'*finite_diff_data.Wtt'/(h_p2^2);
        H12 = sum((finite_diff_data.W_xy.*Z((ii-r*leap):leap:(ii+r*leap),(jj-r*leap):leap:(jj+r*leap))), 'all')/(h_p1*h_p2);
        H = [H11, H12; H12, H22];
        
        
        
        
    catch ME % optimal parameter set is on the boundary of the parameters range
        H = zeros(2,2);
    end
    %% old
%     delta_n = 2;
    % https://v8doc.sas.com/sashtml/ormp/chap5/sect28.htm
    % horizontal + vertical
%     F_0_0 = Z(ii,jj); % F_obj at p*
%     F_0_P = Z(ii+delta_n,jj); % F_obj at (p1*+delta_p1,p2*)
%     F_0_PP = Z(ii+2*delta_n,jj); % F_obj at (p1*+2delta_p1,p2*)
%     F_0_M = Z(ii-delta_n,jj); % F_obj at (p1*-delta_p1,p2*)
%     F_0_MM = Z(ii-2*delta_n,jj); %...
%     F_P_0 = Z(ii, jj+delta_n);
%     F_PP_0 = Z(ii, jj+2*delta_n);
%     F_M_0 = Z(ii, jj-delta_n);
%     F_MM_0 = Z(ii, jj-2*delta_n);
%     % diagonals
%     F_P_P = Z(ii+delta_n, jj+delta_n);
%     F_P_M = Z(ii-delta_n, jj+delta_n);
%     F_M_P = Z(ii+delta_n, jj-delta_n);
%     F_M_M = Z(ii-delta_n, jj-delta_n);
%     F_PP_PP = Z(ii+2*delta_n, jj+2*delta_n);
%     F_MM_MM = Z(ii-2*delta_n, jj-2*delta_n);
%     F_PP_MM = Z(ii-2*delta_n, jj+2*delta_n);
%     F_MM_PP = Z(ii+2*delta_n, jj-2*delta_n);
%     
%     % second diagonals
%         % right side 
%     F_P_PP = Z(ii+2*delta_n, jj+delta_n);
%     F_P_MM = Z(ii-2*delta_n, jj+delta_n);
%     F_PP_P = Z(ii+delta_n, jj+2*delta_n);
%     F_PP_M = Z(ii-delta_n, jj+2*delta_n);
%         % left side 
%     F_M_PP = Z(ii+2*delta_n, jj-delta_n);
%     F_M_MM = Z(ii-2*delta_n, jj-delta_n);
%     F_MM_P = Z(ii+delta_n, jj-2*delta_n);
%     F_MM_M = Z(ii-delta_n, jj-2*delta_n);
%     
%     
%     
%     h_p1 = X(ii,jj+delta_n)-X(ii,jj); % delta_p1
%     h_p2 = Y(ii+delta_n,jj)-Y(ii,jj); % delta_p2
%     
%     H11 = (-(F_PP_0+F_MM_0)+16*(F_P_0+F_M_0)-30*F_0_0)/(12*h_p1*h_p1);
%     H22 = (-(F_0_PP+F_0_MM)+16*(F_0_P+F_0_M)-30*F_0_0)/(12*h_p2*h_p2);
%     
%     
%     
%     
% %     H12 = (F_P_P-F_P_M-F_M_P+F_M_M)/(4*h_p1*h_p2); % previous +O(hx*hy)
%     
%     n1 = F_P_MM+F_PP_M+F_MM_P+F_M_PP;
%     n2 = F_M_MM+F_MM_M+F_P_PP+F_PP_P;
%     n3 = F_PP_MM+F_MM_PP-F_MM_MM-F_PP_PP;
%     n4 = F_M_M+F_P_P-F_P_M-F_M_P;
%     
%     
%     H12 = (8*n1-8*n2-n3+64*n4)/(144*h_p1*h_p2); % new +O(hx^2*hy^2)   http://www.holoborodko.com/pavel/2014/11/04/computing-mixed-derivatives-by-finite-differences/
%     
% %     H11 = (-Z(ii,jj+2*delta_n)+16*Z(ii,jj+delta_n)-30*Z(ii,jj)+16*Z(ii,jj-delta_n)-Z(ii,jj-2*delta_n))/(12*(X(ii,jj+1*delta_n)-X(ii,jj))^2);
% %     H22 = (-Z(ii+2*delta_n,jj)+16*Z(ii+delta_n,jj)-30*Z(ii,jj)+16*Z(ii-delta_n,jj)-Z(ii-2*delta_n,jj))/(12*(Y(ii+1*delta_n,jj)-Y(ii,jj))^2);
% % %     H22 = (Z(ii+2*delta_n,jj)-Z(ii,jj)-Z(ii,jj)+Z(ii-2*delta_n,jj))/(4*(Y(ii+1*delta_n,jj)-Y(ii,jj))*(Y(ii+1*delta_n,jj)-Y(ii,jj)));
% % %     H11 = (Z(ii,jj+2*delta_n)-Z(ii,jj)-Z(ii,jj)+Z(ii,jj-2*delta_n))/(4*(X(ii,jj+1*delta_n)-X(ii,jj))*(X(ii,jj+1*delta_n)-X(ii,jj)));
% %     H12 = (Z(ii+1*delta_n,jj+1*delta_n)-Z(ii-1*delta_n,jj+1*delta_n)-Z(ii+1*delta_n,jj-1*delta_n)+Z(ii-1*delta_n,jj-1*delta_n))/(4*(X(ii,jj+1*delta_n)-X(ii,jj))*(Y(ii+1*delta_n,jj)-Y(ii,jj)));
%     
%     
% %     H12 = -(Z(ii,jj+delta_n)+Z(ii,jj-delta_n)+Z(ii+delta_n,jj)+Z(ii-delta_n,jj)-2*Z(ii,jj)-Z(ii+delta_n,jj+delta_n)-Z(ii-delta_n,jj-delta_n))/(2*(X(ii,jj+1*delta_n)-X(ii,jj))*(Y(ii+1*delta_n,jj)-Y(ii,jj)));
%     H21 = H12;
end 



function [uncert_lim] = getUncertaintyLimits(c_h,optimal_params)
%%%%%%%%%%%%   ->  uncert_lim = [height, x_min, x_max, y_min, y_max, bounding_flag];
% returns the bounding rectangular dimensions for different objective function values defined by level_list

%     level_list = linspace(0,2,51);
    level_list = c_h.LevelList;
    optimal_params = [1,1];
%     M = c_h.ContourMatrix;
%     ii=1;
% uncert_lim = [];
    ax = c_h.Parent;
    [M_temp,c_h_temp] = contour(ax,c_h.XData,c_h.YData,c_h.ZData,'w', 'fill', 'Off', 'LevelList', level_list,'Tag','HighlightLevels');
%     delete(c_h_temp);
    ii=1;
    uncert_lim = [];
    while ii<size(M_temp,2)
        current_level = M_temp(1,ii);
        N_pts = M_temp(2,ii);
        x_pts = M_temp(1,(ii+1):(ii+N_pts));
        y_pts = M_temp(2,(ii+1):(ii+N_pts));
        if ~isempty(uncert_lim)
            if current_level == uncert_lim(end,1)
                % duplicated level ->merge limits
                uncert_lim(end,1:5) = [uncert_lim(end,1), min([uncert_lim(end,2),min(x_pts)]),max([uncert_lim(end,3),max(x_pts)]),min([uncert_lim(end,4),min(y_pts)]),max([uncert_lim(end,5),max(y_pts)])];
            else
                % new level -> create new entry
                uncert_lim(end+1,1:5) = [current_level, min(x_pts),max(x_pts),min(y_pts),max(y_pts)];
            end
        else
            % first level -> create new entry
                uncert_lim(end+1,1:5) = [current_level, min(x_pts),max(x_pts),min(y_pts),max(y_pts)];
        end
        uncert_lim(end,6) = ((uncert_lim(end,2)<optimal_params(1))&&(optimal_params(1)<uncert_lim(end,3)))&&((uncert_lim(end,4)<optimal_params(2))&&(optimal_params(2)<uncert_lim(end,5)));
        ii = ii+N_pts+1;
%         if ii>=size(M_temp,2)
%             break;
%         end
    end
%     while ii<=size(M,2)
%         current_level = c_h(ii);
%         if current_level
% 
%     end
end
function [diff_circ_h] = markPointsForDerivatives(ax,X,Y,optimal_ind,test)
    global finite_diff_data;
    try
        ii = test{optimal_ind}.P2_ind;
        jj = test{optimal_ind}.P1_ind;
        r = finite_diff_data.n/2;
        leap = finite_diff_data.leap;
        XData = X((ii-r*leap):leap:(ii+r*leap),(jj-r*leap):leap:(jj+r*leap));
        YData = Y((ii-r*leap):leap:(ii+r*leap),(jj-r*leap):leap:(jj+r*leap));
        C = 'k';
        diff_circ_h = {};
        for ii=1:size(XData,1)
            for jj=1:size(XData,2)
                diff_circ_h{end+1} = scatter(ax,XData(ii,jj),YData(ii,jj)',5,...
                    'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'w', 'Visible', 'Off','Tag','finDifMarker');
            end
        end
    catch ME
        diff_circ_h = [];
    end
end


function [text_box_temp] = makeTextbox(ax_temp,c_h_temp)
% Hessians data textbox
    temp_str = {['det($H_{\hat{\theta}}$)=',sprintf('%.2e',det(c_h_temp.UserData.Hessian))];...
        ['cond($H_{\hat{\theta}}$)=',sprintf('%.2f',cond(c_h_temp.UserData.Hessian))];...
        ['det($\tilde{H_{\hat{\theta}}}$)=',sprintf('%.2f',det(c_h_temp.UserData.Hessian_tilde))]};
    text_box_temp = annotation('textbox', getTextboxDim(ax_temp), 'String',temp_str,...
        'Interpreter', 'Latex', 'Margin',5,'Color', 'w', 'BackgroundColor', 'k','FaceAlpha', 0.2,'EdgeColor', 'none', 'HorizontalAlignment', 'left');
    
%      set(text_box{max_l_ii,max_l_jj},...
%                     'Color', 'w', 'BackgroundColor', 'k','FaceAlpha', 0.2, 'FontSize', font_size,...
%                     'EdgeColor', 'none');
%     
end



% sets position of textboxes in axes 
function adjustTextbox(ax,text_box) 
    width_ratio = .9;
    max_length = 0;
    max_l_ii = [];
    max_l_jj = [];
    for ii=1:size(ax,1)
        for jj=1:size(ax,2)
            if max(strlength(text_box{ii,jj}.String)-[17;16;24])>=max_length
                max_l_ii = ii;
                max_l_jj = jj;
                max_length = max(strlength(text_box{ii,jj}.String)-[17;16;24]);
            end
        end
    end
    
    font_size = 14;
    drawnow
    set(text_box{max_l_ii,max_l_jj},'FontSize', font_size, 'FitBoxToText', 'On');
% drawnow
    
%         set(text_box{max_l_ii,max_l_jj}, 'FitBoxToText', 'On');
%         drawnow;
      pause(1e-4)
%       if text_box{max_l_ii,max_l_jj}.Position(3)>ax{max_l_ii,max_l_jj}.Position(3)*width_ratio
        font_size = max([text_box{max_l_ii,max_l_jj}.FontSize*((ax{max_l_ii,max_l_jj}.Position(3)*width_ratio)/text_box{max_l_ii,max_l_jj}.Position(3)),5]);
        text_box{max_l_ii,max_l_jj}.FontSize = font_size;
%         box_height = 0;
%         text_box{max_l_ii,max_l_jj}.Position(3) = target_pos(3);\
      

%       end
      drawnow
          drawnow limitrate nocallbacks;

%       text_box{max_l_ii,max_l_jj}.FitBoxToText
%        drawnow
%        drawnow limitrate nocallbacks;
       box_width = text_box{max_l_ii,max_l_jj}.Position(3);
       box_height = text_box{max_l_ii,max_l_jj}.Position(4);
       distance_from_bottom = text_box{max_l_ii,max_l_jj}.Position(2)-ax{max_l_ii,max_l_jj}.Position(2);
       total_gap = ax{max_l_ii,max_l_jj}.Position(4)-text_box{max_l_ii,max_l_jj}.Position(4);
      for ii=1:size(ax,1)
        for jj=1:size(ax,2)
                
%             target_pos = getTextboxDim(ax{ii,jj});
            text_box{ii,jj}.FontSize = font_size;
            text_box{ii,jj}.Position =[ax{ii,jj}.Position(1) ax{ii,jj}.Position(2)+total_gap box_width, box_height] ;
%             set(text_box{ii,jj}, 'FitBoxToText', 'On');
%             text_box{ii,jj}.Position(3) = box_width;
        end
      end
      drawnow();
%     set(text_box, 'FitBoxToText', 'On')
%     
%     for ii=1:size(ax,1)
%         for jj=1:size(ax,2)
% %             str_length = strlength(text_box{ii,jj}.String)-[17;16;24];
% %             old_units = text_box{ii,jj}.FontUnits;
% %             text_box{ii,jj}.FontUnits = 'normalized';
% %                 font_size = ax{ii,jj}.Position(3)/max(str_length)*2.5
% % %             text_box{ii,jj}.FontUnits = old_units
%             text_box{ii,jj}.FontUnits = 'points';
%             font_size = 14;
%             target_pos = getTextboxDim(ax{ii,jj});
% %             set(text_box{ii,jj},'Position',target_pos);
%               set(text_box{ii,jj},'Position',target_pos,...
%                     'Color', 'w', 'BackgroundColor', 'k','FaceAlpha', 0.2, 'FontSize', font_size,...
%                     'EdgeColor', 'none');
%               set(text_box{ii,jj}, 'FitBoxToText', 'On');
%               drawnow();
%               if text_box{ii,jj}.Position(3)>target_pos(3)
%                   text_box{ii,jj}.FontSize = text_box{ii,jj}.FontSize*(target_pos(3)/text_box{ii,jj}.Position(3))
%               end
%               
% %             if text_box{ii,jj}.Position(3)>target_pos(3)
% %                 set(text_box{ii,jj}, 'FitBoxToText', 'Off');
% %                 set(text_box{ii,jj},'Position', target_pos);
% % %                 pause(0.001);
% %             end
% %             while (text_box{ii,jj}.Position(3)>target_pos(3))&&(text_box{ii,jj}.FontSize>=5)
% %                 text_box{ii,jj}.FontSize = text_box{ii,jj}.FontSize-1;
% % 
% %                 target_pos = getTextboxDim(ax{ii,jj});
% %             pause(1e-8);
% %             end
%             set(text_box{ii,jj},'Position',target_pos);
%         end
%     end
end
%     
function [dim] = getTextboxDim(ax)
%     switch ax.Parent.PositionConstraint
%         case 'innerposition'
%             p = ax.OuterPosition;
%         case 'outerposition'
%             p = ax.Position;
%     end
    p = ax.Position;
    x0 = 0;
    y0 = 0;
    h = 0.5;
    w = 1;
    dim = [p(1)+x0*p(3), p(2)+y0*p(4) w*p(3) h*p(4)];
end

function [test_ind] = paramInd2testInd(P1_ind_target, P2_ind_target,test)
    for test_ind = 1:numel(test)
        if isequal([test(test_ind).P1_ind,test(test_ind).P2_ind], [P1_ind_target, P2_ind_target])
            return
        end
    end
    error('Error calculating finite difference, test with [P1_ind,P2_ind]==[%d,%d] not found', P1_ind_target,P2_ind_target);
end

function adjustPlotsToolbar(fig_h)
tb = uitoolbar(fig_h);
% axis_equal_tt = uitoggletool(tb);
% axis_square_tt = uitoggletool(tb);
arrange_textbox_pt = uipushtool(tb, 'Interruptible', 'off');

m = uimenu('Text', 'Options');
mitem_equal = uimenu(m, 'Text', 'Axis equal', 'Checked','off','MenuSelectedFcn', @mitem_equal_cback);
mitem_det_h = uimenu(m, 'Text', 'det(H) plot', 'Checked','off','MenuSelectedFcn', @mitem_det_h_cback);
mitem_show_marker = uimenu(m, 'Text', 'Show markers (finite-difference)', 'Checked','off','MenuSelectedFcn', @mitem_show_marker_cback);

mitem_remove_spaces = uimenu(m, 'Text', 'Remove padding', 'Checked','off','MenuSelectedFcn', {@removeSpaces,fig_h.UserData.t,fig_h});


% axis_equal_tt.ClickedCallback = {@axis_equal_tt_cback,axis_square_tt,fig_h};
% axis_square_tt.ClickedCallback = {@axis_square_tt_cback,axis_equal_tt,fig_h};
% arrange_textbox_pt.ClickedCallback = {@arrange_textbox_pt_cback,fig_h};
% arrange_textbox_pt.ClickedCallback = {@removeSpaces,fig_h.UserData.t,fig_h};
end

function mitem_equal_cback(src,event)
    fig_h = src.Parent.Parent;
    ax = fig_h.UserData.ax;
    text_box = fig_h.UserData.text_box;  
    c_h = fig_h.UserData.c_h;
    switch src.Checked
        case 'on' % revert to axis square
            axis_mode = 'square';
            src.Checked = 'off';
        case 'off' % turn on axis equal
            axis_mode = 'equal';
            src.Checked = 'on';
    end        
    for ii=1:size(ax,1)
        for jj=1:size(ax,2)
            axis(ax{ii,jj},axis_mode);
            xlim(ax{ii,jj}, [min(c_h{ii,jj}.XData,[],'all') max(c_h{ii,jj}.XData,[],'all')]);
            ylim(ax{ii,jj}, [min(c_h{ii,jj}.YData,[],'all') max(c_h{ii,jj}.YData,[],'all')]);
        end
    end  
    if ~isequal(text_box,cell(size(text_box)))
        adjustTextbox(ax,text_box);
    end
end


function mitem_show_marker_cback(src,event)
    fig_h = src.Parent.Parent;
%     ax = fig_h.UserData.ax;
%     text_box = fig_h.UserData.text_box;  
%     c_h = fig_h.UserData.c_h;
    markers_h = findobj('Tag', 'finDifMarker');
    switch src.Checked
        case 'on' % revert to visible 'off'
            set(markers_h, 'Visible', 'Off');
            src.Checked = 'off';
        case 'off' % set visible 'on'
            set(markers_h, 'Visible', 'On');
            src.Checked = 'on';
    end        
end


% not in use
function mitem_det_h_cback(src,event)
    fig_h = src.Parent.Parent;
    ax = fig_h.UserData.ax;
    text_box = fig_h.UserData.text_box;  
    c_h = fig_h.UserData.c_h;
    temp_fig = figure('Name',strcat(fig_h.UserData.fig_name,'_(Hessian data)'));
    t2 = tiledlayout(2,3, 'Padding', 'tight','TileSpacing','tight','Units', 'centimeters');
    t2.InnerPosition(3) = 12;
    
    det_hessian_tol = det(c_h{end,end}.UserData.Hessian)*1e-3;
    
    t_ind = 1;
    im_ax(t_ind) = nexttile(); % H11
    for ii=1:size(ax,1)
        for jj=1:size(ax,2)   
            C(ii,jj) = c_h{ii,jj}.UserData.Hessian(1,1);
        end
    end
    im_h(t_ind) = image(im_ax(t_ind),C,'CDataMapping','scaled');
    colorbar(im_ax(t_ind));
    im_ax(t_ind).Subtitle.Interpreter = 'latex';
    im_ax(t_ind).Subtitle.String = '$H_{\hat{\theta},11}$';
    t_ind = t_ind+1;
    
    im_ax(t_ind) = nexttile(); % H22
    for ii=1:size(ax,1)
        for jj=1:size(ax,2)   
            C(ii,jj) = c_h{ii,jj}.UserData.Hessian(2,2);
        end
    end
    im_h(t_ind) = image(im_ax(t_ind),C,'CDataMapping','scaled');
    colorbar(im_ax(t_ind));
    im_ax(t_ind).Subtitle.Interpreter = 'latex';
    im_ax(t_ind).Subtitle.String = '$H_{\hat{\theta},22}$';
    t_ind = t_ind+1;
    
    im_ax(t_ind) = nexttile(); % max(lambda)
    for ii=1:size(ax,1)
        for jj=1:size(ax,2)   
            C(ii,jj) = eigs(c_h{ii,jj}.UserData.Hessian,1);
        end
    end
    im_h(t_ind) = image(im_ax(t_ind),C,'CDataMapping','scaled');
    colorbar(im_ax(t_ind));
    im_ax(t_ind).Subtitle.Interpreter = 'latex';
    im_ax(t_ind).Subtitle.String = '$\rho(H_{\hat{\theta}}$)';
    t_ind = t_ind+1;
    
    im_ax(t_ind) = nexttile(); % det(H)
    for ii=1:size(ax,1)
        for jj=1:size(ax,2)   
            
            C(ii,jj) = det(c_h{ii,jj}.UserData.Hessian);
            A(ii,jj) = 1;
            if det(c_h{ii,jj}.UserData.Hessian)<det_hessian_tol
                C(ii,jj) = NaN;
                A(ii,jj) = 0;
            end
        end
    end
    im_h(t_ind) = image(im_ax(t_ind),C,'CDataMapping','scaled', 'AlphaData',A);
    colorbar(im_ax(t_ind));
    im_ax(t_ind).Subtitle.Interpreter = 'latex';
    im_ax(t_ind).Subtitle.String = 'det($H_{\hat{\theta}}$)';
    t_ind = t_ind+1;
    
    im_ax(t_ind) = nexttile(); % cond(H)
    for ii=1:size(ax,1)
        for jj=1:size(ax,2)   
            C(ii,jj) = cond(c_h{ii,jj}.UserData.Hessian);
            A(ii,jj) = 1;
            if det(c_h{ii,jj}.UserData.Hessian)<det_hessian_tol
                C(ii,jj) = Inf;
                A(ii,jj) = 0;
            end
        end
    end
    im_h(t_ind) = image(im_ax(t_ind),C,'CDataMapping','scaled', 'AlphaData',A);
    colorbar(im_ax(t_ind));
    im_ax(t_ind).Subtitle.Interpreter = 'latex';
    im_ax(t_ind).Subtitle.String = 'cond($H_{\hat{\theta}}$)';
    im_ax(t_ind).Colormap = flip(im_ax(t_ind).Colormap,1);

    t_ind = t_ind+1;
    
    im_ax(t_ind) = nexttile(); % det(Hessian_tilde)
    for ii=1:size(ax,1)
        for jj=1:size(ax,2)   
            C(ii,jj) = det(c_h{ii,jj}.UserData.Hessian_tilde);
            A(ii,jj) = 1;
            if det(c_h{ii,jj}.UserData.Hessian)<det_hessian_tol
                C(ii,jj) = NaN;
                A(ii,jj) = 0;
            end
        end
    end
    im_h(t_ind) = image(im_ax(t_ind),C,'CDataMapping','scaled', 'AlphaData',A);
    colorbar(im_ax(t_ind));
    im_ax(t_ind).Subtitle.Interpreter = 'latex';
    im_ax(t_ind).Subtitle.String = 'det($\tilde{H_{\hat{\theta}}}$)';
%     im_ax(t_ind).Colormap = flip(im_ax(t_ind).Colormap,1);
    t_ind = t_ind+1;
    
    for kk = 1:length(im_ax)
        set(im_ax(kk).XLabel,'Interpreter', 'latex','String',{'$\eta \rightarrow$'},'FontSize', 10)
        set(im_ax(kk).YLabel,'Interpreter', 'latex','String',{'$\leftarrow$ Indentation depth (mm)'},'FontSize', 10)
        set(im_ax(kk),'XTick', 1:length(fig_h.UserData.eta_arr), 'XTickLabel', strsplit(num2str(fig_h.UserData.eta_arr,3))')
        set(im_ax(kk),'YTick', 1:length(fig_h.UserData.indentation_depths(2:end)), 'YTickLabel', strsplit(num2str(fig_h.UserData.indentation_depths(2:end),3))')
%         im_ax(kk).XTickLabel = num2str(fig_h.UserData.eta_arr,3)';
        im_ax(kk).Subtitle.Interpreter = 'latex';
    end
    axis(im_ax, 'square');
    figH.UserData.t2 = t2;
end




function arrange_textbox_pt_cback(src,event,fig_h)
    ax = fig_h.UserData.ax;
    text_box = fig_h.UserData.text_box;
    if ~isequal(text_box,cell(size(text_box)))
        adjustTextbox(ax,text_box);
    end
    disp('done');
end

function fig1SizeChangeCBack(src,callbackdata)

    if isfield(src.UserData, 'ax')&&isfield(src.UserData,'text_box')
        text_box = src.UserData.text_box;
        if ~isequal(text_box,cell(size(text_box)))
            adjustTextbox(src.UserData.ax,src.UserData.text_box);
        end
    end
    
%     if isfield(src.UserData, 'rect')&&isfield(src.UserData,'t')
%         rect = src.UserData.rect;
%         rect_out = src.UserData.rect_out;
%         t = src.UserData.t;
        
        
%         t.PositionConstraint = 'OuterPosition';
% %                 t.OuterPosition = [0.1 0.1 0.3 0.8];
% 
%         a = src.Position(3);
%         b = 100;
%         c = 50;
%         switch t.PositionConstraint
%             case 'innerposition'
%                 t.Position(1) = b/a;
%                 t.Position(3) = 1-(b+c)/a;
%                 rect.Position = t.InnerPosition;
%                 rect_out.Position = t.OuterPosition;
%             case 'outerposition'
% %                 rect.Position = t.OuterPosition;
%             
%         end
%         t.Position(1)-t.OuterPosition(1);
%         if (t.Position(1)-t.OuterPosition(1))<0.1
% %             t.Position(3) = t.Position(3)-(t.Position(1)-t.OuterPosition(1)-0.1);
%             t.PositionConstraint = 'innerposition';
% %             t.Position(3) = 1-(b+c)/a-abs(t.Position(1)-t.OuterPosition(1)-0.1)
% 1
%             t.Position(3) = t.Position(4)/src.Position(4)*a
%         end    
%     end
% %     t.OuterPosition
% %     t.Units
% t.Position(3)
end

function cb_callback(src,callbackdata)
    prompt = {'Colormap levels:'};
    dlgtitle = 'Colormap levels:';
    dims = [1 35];
    definput = {'0:0.1:0.5','hsv'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    levels_arr = eval(answer{1});
    cmap = colormap(parula(length(levels_arr)));
    fig_h = src.Parent.Parent;
    c_h = fig_h.UserData.c_h;
    ellipsoid_contour = fig_h.UserData.ellipsoid_contour;
    ax = fig_h.UserData.ax;
    for ii=1:size(c_h,1)
        for jj=1:size(c_h,2)
            c_h{ii,jj}.LevelList = levels_arr;
            ellipsoid_contour{ii,jj}.LevelList = levels_arr;
            ax{ii,jj}.CLim = [levels_arr(1), levels_arr(end)];
            ax{ii,jj}.Colormap = cmap;
        end
    end
%     src.
end

function addTransperacy(c_h,alpha)
    for ii=1:size(c_h,1)
        for jj=1:size(c_h,2)
            % This is the secret that 'keeps' the transparency.
            eventFcn = @(srcObj, e) updateTransparency(srcObj);
            addlistener(c_h{ii,jj}, 'MarkedClean', eventFcn);
        end
    end
    % Elsewhere in script, a separate file, or another method of your class.
    function updateTransparency(contourObj)
        contourFillObjs = contourObj.FacePrims;
        for i = 1:length(contourFillObjs)
            % Have to set this. The default is 'truecolor' which ignores alpha.
            contourFillObjs(i).ColorType = 'truecoloralpha';
            % The 4th element is the 'alpha' value. First 3 are RGB. Note, the
            % values expected are in range 0-255.
            contourFillObjs(i).ColorData(4) = alpha*255;
        end
    end
end