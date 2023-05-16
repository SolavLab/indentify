function [figH] = open_stat_fig(figH)
    fig_h = figH(1);
    c_h = fig_h.UserData.c_h;
    cb_h_location = 'EastOutside';
    CScale  = {'log','log','linear'};
    figH(2) = figure('WindowStyle', 'normal','Name',strcat('Hessian Statistics [',figH(1).UserData.group_name,']'),'NumberTitle', 'off','units','centimeters');
    figH(2).Position(3) = figH(1).Position(3); figH(2).Position(4) = 5.1;
    t2 = tiledlayout(1,3, 'Padding', 'none','TileSpacing','loose','Units',...
        'normalized','OuterPosition',[0 0 1 1]); % for presentation
    t.InnerPosition(1) = 0.1; t.InnerPosition(3) = 0.8;
    t2.YLabel.FontSize = figH(1).UserData.t.YLabel.FontSize;
    t2.XLabel.FontSize = figH(1).UserData.t.XLabel.FontSize;
    t2.Title.FontSize = figH(1).UserData.t.Title.FontSize;
    det_hessian_tol = 0;
    t_ind = 1;
    ax = fig_h.UserData.ax;
    N_rows = size(ax,1);
    N_cols = size(ax,2);
    [X,Y] = meshgrid(0.5:(0.5+N_rows),0.5:(0.5+N_cols));
    %% det(H) heatmap
    im_ax(t_ind) = nexttile(); % det(H)
    C = nan(N_rows,N_cols);
    for ii=1:N_rows
        for jj=1:N_cols       
            if det(c_h{ii,jj}.UserData.Hessian)>det_hessian_tol
                C(ii,jj) = det(c_h{ii,jj}.UserData.Hessian);
            end
        end
    end
    % Plotting
    invalid_mask = isnan(C); % mask of invalid hessians (blank squares)
    [invalid_rows, invalid_cols] = find(invalid_mask); % and their locations
    im_h(t_ind) = imagesc(im_ax(t_ind),C,'CDataMapping','scaled', 'AlphaData',~invalid_mask); % draw heatmap
    % cross-out blank squares
    im_ax(t_ind).NextPlot = 'add'; %needs to be after plotting image
    plotCross(im_ax(t_ind),X,Y,invalid_rows,invalid_cols)
    
    % Arrange plot
    im_ax(t_ind).Tag ='detH';
    im_ax(t_ind).XAxisLocation = 'top';
    im_ax(t_ind).ColorScale = CScale{t_ind};
    cb_h(t_ind) = colorbar(im_ax(t_ind), 'Location',cb_h_location);
    im_ax(t_ind).Title.Interpreter = 'latex';
    im_ax(t_ind).Title.String = '\textbf{(b)} $\textit{M1}=det($\boldmath${\bar{H}^\ast}$)';
    im_ax(t_ind).Title.FontName = 'Times New Roman';
    temp_CLim = [min(C,[],'all','omitnan'), max(C,[],'all','omitnan')];
    im_ax(t_ind).CLim = temp_CLim;

    % proceed to next plot 
    t_ind = t_ind+1;
    %% cond(H) heatmap
    im_ax(t_ind) = nexttile(); 
    C = nan(N_rows,N_cols);
    for ii=1:N_rows
        for jj=1:N_cols   
            if det(c_h{ii,jj}.UserData.Hessian)>det_hessian_tol
                C(ii,jj) = cond(c_h{ii,jj}.UserData.Hessian);
            end
        end
    end
    % Plotting
    invalid_mask = isnan(C); % mask of invalid hessians (blank squares)
    [invalid_rows, invalid_cols] = find(invalid_mask); % and their locations
    im_h(t_ind) = imagesc(im_ax(t_ind),C,'CDataMapping','scaled', 'AlphaData',~invalid_mask); % draw heatmap
    % cross-out blank squares
    im_ax(t_ind).NextPlot = 'add'; %needs to be after plotting image
    plotCross(im_ax(t_ind),X,Y,invalid_rows,invalid_cols)
    
    % Arrange plot
    im_ax(t_ind).Tag ='condH';
    im_ax(t_ind).XAxisLocation = 'top';
    im_ax(t_ind).ColorScale = CScale{t_ind};
    im_ax(t_ind).Colormap = flip(im_ax(t_ind).Colormap,1); 
    cb_h(t_ind)=colorbar(im_ax(t_ind), 'Location',cb_h_location);
    im_ax(t_ind).Title.Interpreter = 'latex';
    im_ax(t_ind).Title.String = '\textbf{(c)} $\textit{M2}=cond($\boldmath${\bar{H}^\ast}$)';
    im_ax(t_ind).Title.FontName = 'Times New Roman';
    temp_CLim = [min(C,[],'all','omitnan'), max(C,[],'all','omitnan')];
    im_ax(t_ind).CLim = temp_CLim;
    
    % Proceed to next plot 
    t_ind = t_ind+1;
    %% det(Hessian_tilde) heatmap
    im_ax(t_ind) = nexttile(); 
    C = nan(N_rows,N_cols);
    for ii=1:N_rows
        for jj=1:N_cols   
            if det(c_h{ii,jj}.UserData.Hessian)>det_hessian_tol
                C(ii,jj) = det(c_h{ii,jj}.UserData.Hessian_tilde);
            end
        end
    end
    % Plotting
    invalid_mask = isnan(C); % mask of invalid hessians (blank squares)
    [invalid_rows, invalid_cols] = find(invalid_mask); % and their locations
    im_h(t_ind) = imagesc(im_ax(t_ind),C,'CDataMapping','scaled', 'AlphaData',~invalid_mask); % draw heatmap
    % cross-out blank squares
    im_ax(t_ind).NextPlot = 'add'; %needs to be after plotting image
    plotCross(im_ax(t_ind),X,Y,invalid_rows,invalid_cols)

    % Arrange plot
    im_ax(t_ind).Tag ='detH_tilde';
    im_ax(t_ind).XAxisLocation = 'top';
    im_ax(t_ind).ColorScale = CScale{t_ind};
    cb_h(t_ind)=colorbar(im_ax(t_ind), 'Location',cb_h_location);
    im_ax(t_ind).Title.Interpreter = 'latex';
    im_ax(t_ind).Title.String = '\textbf{(d)} $\textit{M3}=det($\boldmath${\tilde{H}^\ast}$)';
    im_ax(t_ind).Title.FontName = 'Times New Roman';
    temp_CLim = [min(C,[],'all','omitnan'), max(C,[],'all','omitnan')];
    im_ax(t_ind).CLim = temp_CLim;

    % Proceed to next plot 
    t_ind = t_ind+1;
    %% final adjusting of plots
    for kk = 1:length(im_ax)
        im_ax(kk).YAxis.TickLength = [0 0];
        im_ax(kk).XAxis.TickLength = [0 0];
        set(im_ax(kk).Subtitle,'Interpreter', 'latex','String',{'$\eta$'},'FontSize', 11,...
                'HorizontalAlignment','right','Position', [0.5 0.3 0],'FontWeight', 'bold');
        if kk==1
            set(im_ax(kk).YLabel,'Interpreter', 'latex','String',{'$\delta$ (mm)'},'FontSize', 11);
            set(im_ax(kk),'YTick', 1:length(fig_h.UserData.indentation_depths(2:end)), 'YTickLabel', strsplit(num2str(abs(fig_h.UserData.indentation_depths(2:end)),3))')
        else 
            set(im_ax(kk),'YTick',[]);
        end
        set(im_ax(kk),'XTick', 1:length(fig_h.UserData.eta_arr), 'XTickLabel', strsplit(num2str(fig_h.UserData.eta_arr,3))','XTickLabelRotation',0)
        im_ax(kk).FontSize = figH(1).UserData.ax{1,1}.FontSize;
        im_ax(kk).Subtitle.Interpreter = 'latex';
    
        % Plot gridlines
        hold(im_ax(kk),'On')
        mesh(im_ax(kk),X,Y,zeros(size(X)),'FaceColor', 'None','EdgeColor','k','LineWidth',1);
        drawnow();

        %% Arrange colorbar ticks
        tick_label_opt = 1;
        cb_h(kk).TicksMode = 'manual';
        if length(cb_h(kk).Ticks)<4
            minor_ticks = cb_h(kk).Ruler.MinorTickValues;
            a = minor_ticks(1); 
            b = minor_ticks(end);
            [~, c_ind] = min(abs(10^(0.75*log10(a)+0.25*log10(b))-minor_ticks)); c = minor_ticks(c_ind);
            [~, d_ind] = min(abs(10^(0.5*log10(a)+0.5*log10(b))-minor_ticks)); d = minor_ticks(d_ind);
            [~, e_ind] = min(abs(10^(0.25*log10(a)+0.75*log10(b))-minor_ticks)); e = minor_ticks(e_ind);
            new_ticks = [a c d e b];
            tol = 1/10;
            mask = abs(log10(new_ticks)-log10(cb_h(kk).Ruler.TickValues'))<=((log10(b)-log10(a))*tol);
            [~,mask_col] = find(mask);
            new_ticks(unique(mask_col)) = [];
            cb_h(kk).Ticks = unique(sort([cb_h(kk).Ruler.TickValues new_ticks]));
            switch tick_label_opt
                case 1
                    if (cb_h(kk).Ticks(1)<=1e-1)||(cb_h(kk).Ticks(1)>=1e1)||(cb_h(kk).Ticks(end)/cb_h(kk).Ticks(1)>=100)
                        cb_h(kk).Ruler.TickLabelFormat = '%.e';
                        if (cb_h(kk).Ticks(end)/cb_h(kk).Ticks(1)<=10000)
                            exponent_factor = floor(log10(cb_h(kk).Ticks(1)));
                            for ii=1:length(cb_h(kk).Ticks)
                                new_lbls{ii} = num2str(cb_h(kk).Ticks(ii)/(10^(exponent_factor)));
                           end
                           cb_h(kk).TickLabels = new_lbls;
                           cb_h(kk).Title.String = ['\times10^{',num2str(exponent_factor),'}'];
                           cb_h(kk).Ruler.TickLabelFormat = '%g';
                        end
                    end
                case 2
                    exponent_factor = floor(log10(cb_h(kk).Ticks(1)));
                    if exponent_factor~=0
                        for ii=1:length(cb_h(kk).Ticks)
                            new_lbls{ii} = num2str(cb_h(kk).Ticks(ii)/(10^(exponent_factor)));
                        end
                        cb_h(kk).TickLabels = new_lbls;
                        cb_h(kk).Title.String = ['\times10^{',num2str(exponent_factor),'}'];
                        cb_h(kk).Ruler.TickLabelFormat = '%g';
                        cb_h(kk).Title.HorizontalAlignment = 'left';
                    end
            end
            
           
        end
    end
    %% Arrange figure and store graphic object's handles to UserData
    set(cb_h,'Location','EastOutside');
    axis(im_ax, 'tight');
    axis(im_ax, 'square');
    drawnow;
    t2.InnerPosition(1) = 0.1; t2.InnerPosition(3) = 0.8; 
    drawnow;
    figH(2).UserData.im_ax = im_ax;
    figH(2).UserData.cb_h = cb_h;
    figH(2).UserData.im_h = im_h;
    figH(2).UserData.t2 = t2;
end

function plotCross(ax_h,X,Y,invalid_rows,invalid_cols)
        for kk=1:length(invalid_rows)
            ii=invalid_rows(kk);
            jj=invalid_cols(kk);
            plot(ax_h,[X(ii,jj) X(ii,jj)+1],[Y(ii,jj) Y(ii,jj)+1],'Color', 'k', 'LineWidth',0.5);
            plot(ax_h,[X(ii,jj) X(ii,jj)+1],[Y(ii,jj)+1 Y(ii,jj)],'Color', 'k', 'LineWidth',0.5);
        end
end
%% 
% _*indentify footer text*_ 
% 
% License: <https://github.com/SolavLab/indentify/blob/main/LICENSE>
% 
% indentify: An open-source project for exploring the identifiability of 
% soft-tissue material parameters from noninvasive indentation test and
% inverse finite-element analysis.
% 
% Copyright (C) 2022 Zohar Oddes, Dana Solav, and the indentify contributors
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.
% 
% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.