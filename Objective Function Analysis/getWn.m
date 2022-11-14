function [Wn] =getWn(distances,indenterRadius)

% 1
% 2
% 3
    OPT = 2;
    Wn = zeros(size(distances));
    switch OPT
        case 1
            ROI_bounds_normalized = [1 4]; %boundaries in units of indenterRadius 
            ROI_bounds = ROI_bounds_normalized*indenterRadius;% absolute boundaries
            nbis = 5; 
            edges = linspace(ROI_bounds(1),ROI_bounds(2),nbis);

        %     disp(['ROI_bounds=', num2str(ROI_bounds)]);
            [N,edges] = histcounts(distances,edges);
            %                 norm_factor = zeros(size(distances));
            Wn = zeros(size(distances));
            for i=1:length(distances)
                for j=1:length(edges)-1
                    if (edges(j)<=distances(i))&&(distances(i)<edges(j+1))
        %                 Wn(i) = (1/distances(i)^2)*(1/N(j));
        %                     Wn(i) = (1/N(j))*1/(distances(i)^2);
                        Wn(i) = (1/N(j));
                        break;
                    end
                end
            end
        case 2
            Wn(distances>=indenterRadius) = 1./distances(distances>=indenterRadius);
%             Wn(distances>=indenterRadius) = distances(distances>=indenterRadius);
% Wn(distances>=indenterRadius) = distances(distances>=indenterRadius).^(2/3);
%             for i=1:length(distances)
% %                 if distances(i)>=indenterRadius*(1+1/3)
% %                 if distances(i)>=indenterRadius*(1)
%                 if distances(i)>=indenterRadius*(1)
% 
%                     Wn(i) = 1/distances(i);
%                 end
%             end
        case 3
            for i=1:length(distances)
                if distances(i)>=indenterRadius
                    Wn(i) = 1/(distances(i)^2);
                end
            end
        case 4
            for i=1:length(distances)
                if distances(i)>=1.5*indenterRadius
                    Wn(i) = 1;
                end
            end
    end
end