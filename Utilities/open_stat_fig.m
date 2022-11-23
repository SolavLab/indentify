function [figH] = open_stat_fig(figH)
    fig_h = figH(1);
    ax = fig_h.UserData.ax;
    c_h = fig_h.UserData.c_h;
    cb_h_location = 'EastOutside';
    CScale  = {'log','log','linear'};
%     CScale  = {'linear','linear','linear'};
    
    figH(2) = figure('WindowStyle', 'normal','Name',strcat('Hessian Statistics [',figH(1).UserData.group_name,']'),'NumberTitle', 'off','units','centimeters');
    figH(2).Position(3) = figH(1).Position(3); figH(2).Position(4) = 5.1;
    
%     t = tiledlayout(2,3, 'Padding', 'tight','TileSpacing','tight');
    t2 = tiledlayout(1,3, 'Padding', 'none','TileSpacing','loose','Units',...
        'normalized','OuterPosition',[0 0 1 1]); % for presentation
    t.InnerPosition(1) = 0.1; t.InnerPosition(3) = 0.8;
%     t2.InnerPosition = [0.1 0 0.8 0.75];
%     t2.PositionConstraint = 'OuterPosition'
%     t2.InnerPosition(3) = 12;
    t2.YLabel.FontSize = figH(1).UserData.t.YLabel.FontSize;
    t2.XLabel.FontSize = figH(1).UserData.t.XLabel.FontSize;
    t2.Title.FontSize = figH(1).UserData.t.Title.FontSize;
%     det_hessian_tol = det(c_h{end,end-1}.UserData.Hessian)*1e-1;
    det_hessian_tol = 0;
    t_ind = 1;
    %% det(H) heatmap
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
%     mean_val = mean(C,'all','omitnan');
%     temp_mask = zeros(size(ax))
%     for ii=1:size(ax,1)
%         for jj=1:size(ax,2) 
%             if C(ii,jj)<mean_val*1e-1
%                 A(ii,jj) = 0;
%                 temp_mask(ii,jj) = 1;
%                 temp_mask_txt(ii,jj) = C(ii,jj);
%                 C(ii,jj) = NaN;
%                 A(ii,jj) = 0;
%             end
%         end
%     end
    im_h(t_ind) = imagesc(im_ax(t_ind),C,'CDataMapping','scaled', 'AlphaData',A);
%         im_h(t_ind) = heatmap(C);
        
    im_ax(t_ind).NextPlot = 'add'; %needs to be after plotting image
    im_ax(t_ind).Tag ='detH';
    im_ax(t_ind).XAxisLocation = 'top';
    im_ax(t_ind).ColorScale = CScale{t_ind};
    cb_h(t_ind) = colorbar(im_ax(t_ind), 'Location',cb_h_location);
    im_ax(t_ind).Title.Interpreter = 'latex';
    im_ax(t_ind).Title.String = '\textbf{(b)} $\textit{M1}=det($\boldmath${\bar{H}^\ast}$)';
    im_ax(t_ind).Title.FontName = 'Times New Roman';
    
    temp_CLim = [min(C,[],'all','omitnan'), max(C,[],'all','omitnan')];
    
    
    
    im_ax(t_ind).CLim = temp_CLim;
    t_ind = t_ind+1;
    %% cond(H) heatmap
    im_ax(t_ind) = nexttile(); 
    for ii=1:size(ax,1)
        for jj=1:size(ax,2)   
            C(ii,jj) = cond(c_h{ii,jj}.UserData.Hessian);
            A(ii,jj) = 1;
            if det(c_h{ii,jj}.UserData.Hessian)<det_hessian_tol
                C(ii,jj) = NaN;
                A(ii,jj) = 0;
            end
        end
    end
    im_h(t_ind) = imagesc(im_ax(t_ind),C,'CDataMapping','scaled', 'AlphaData',A);
    im_ax(t_ind).NextPlot = 'add'; %needs to be after plotting image
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
    t_ind = t_ind+1;
    %% det(Hessian_tilde) heatmap
    im_ax(t_ind) = nexttile(); 
    
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
    im_h(t_ind) = imagesc(im_ax(t_ind),C,'CDataMapping','scaled', 'AlphaData',A);
    im_ax(t_ind).NextPlot = 'add'; %needs to be after plotting image
    im_ax(t_ind).Tag ='detH_tilde';
    im_ax(t_ind).XAxisLocation = 'top';

    im_ax(t_ind).ColorScale = CScale{t_ind};

    cb_h(t_ind)=colorbar(im_ax(t_ind), 'Location',cb_h_location);
    im_ax(t_ind).Title.Interpreter = 'latex';
    im_ax(t_ind).Title.String = '\textbf{(d)} $\textit{M3}=det($\boldmath${\tilde{H}^\ast}$)';
    im_ax(t_ind).Title.FontName = 'Times New Roman';
%     im_ax(t_ind).Colormap = flip(im_ax(t_ind).Colormap,1);
%     im_ax(t_ind).CLim(2) =min([max(C,[],'all'),1]);

    temp_CLim = [min(C,[],'all','omitnan'), max(C,[],'all','omitnan')];
    
    
    
    im_ax(t_ind).CLim = temp_CLim;
    t_ind = t_ind+1;
    %% final adjusting of plots
    for kk = 1:length(im_ax)
        
        im_ax(kk).YAxis.TickLength = [0 0];
        im_ax(kk).XAxis.TickLength = [0 0];
        set(im_ax(kk).Subtitle,'Interpreter', 'latex','String',{'$\eta$'},'FontSize', 11,...
                'HorizontalAlignment','right','Position', [0.5 0.3 0],'FontWeight', 'bold');
        if kk==1
%             set(im_ax(kk).Subtitle,'Interpreter', 'latex','String',{'$\eta \rightarrow$'},'FontSize', 11,...
%                 'HorizontalAlignment','right','Position', [0 0.4101 0]);
%             set(im_ax(kk).Subtitle,'Interpreter', 'latex','String',{'$\eta$'},'FontSize', 11,...
%                 'HorizontalAlignment','right','Position', [0 0.4101 0]);
            
%             set(im_ax(kk).YLabel,'Interpreter', 'latex','String',{'$\leftarrow$Indentation depth (mm)'},'FontSize', 10)
            set(im_ax(kk).YLabel,'Interpreter', 'latex','String',{'$\delta$ (mm)'},'FontSize', 11);
%             set(t2.YLabel,'String',{'$\leftarrow$ Indentation depth'}, 'Interpreter', 'latex','FontSize', 11);
            set(im_ax(kk),'YTick', 1:length(fig_h.UserData.indentation_depths(2:end)), 'YTickLabel', strsplit(num2str(abs(fig_h.UserData.indentation_depths(2:end)),3))')
        else 
            set(im_ax(kk),'YTick',[]);
        end
        set(im_ax(kk),'XTick', 1:length(fig_h.UserData.eta_arr), 'XTickLabel', strsplit(num2str(fig_h.UserData.eta_arr,3))','XTickLabelRotation',0)
%         im_ax(kk).XTick = [];
%         im_ax(kk).XTickLabel = num2str(fig_h.UserData.eta_arr,3)';
        im_ax(kk).FontSize = figH(1).UserData.ax{1,1}.FontSize;
        im_ax(kk).Subtitle.Interpreter = 'latex';
    
        hold(im_ax(kk),'On')
        [X,Y] = meshgrid(0.5:5.5,0.5:5.5);
%         im_ax(kk).CLimMode = 'manual';
        mesh(im_ax(kk),X,Y,zeros(size(X)),'FaceColor', 'None','EdgeColor','k','LineWidth',1);
%         im_ax(kk).CLimMode = 'auto'
                drawnow();
%         for ii=1:size(ax,1)
%             for jj=1:size(ax,2) 
%                 if temp_mask(ii,jj)
%                     txt_h(ii,jj) = text(im_ax(1),X(ii,jj)+0.5,Y(ii,jj)+0.5,sprintf('%.e',temp_mask_txt(ii,jj)),'HorizontalAlignment','center','FontSize', 8);
%                 end
%             end
%         end
                
%         cb_h(kk).Ticks = 10.^[-4:6];
       
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
            opt = 2;
            switch opt
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
                        else
                    %                     lbls = cb_h(kk).TickLabels;
                    %                     for ii=1:numel(lbls)
                    %                         if strcmp(lbls{ii}(end-2:end),'{0}')
                    %                             lbls{ii} = lbls{ii}(1:end-12);
                    %                         end
                    %                     end
                    %                    cb_h(kk).TickLabels = lbls
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
    set(cb_h,'Location','EastOutside');
    
    
%     set(im_ax,'XAxisLocation', 'bottom');
    axis(im_ax, 'tight');
    axis(im_ax, 'square');
    drawnow;
%     figH(2).Position(4) = figH(2).Position(4)*im_ax(1).OuterPosition(4)+1;
    t2.InnerPosition(1) = 0.1; t2.InnerPosition(3) = 0.8; 

    drawnow;
    figH(2).UserData.im_ax = im_ax;
    figH(2).UserData.cb_h = cb_h;
    figH(2).UserData.im_h = im_h;
    figH(2).UserData.t2 = t2;
end