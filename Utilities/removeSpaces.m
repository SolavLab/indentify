% removes spaces between tiles of tiledlatout t and aligns its
% OuterPosition to the figure window
function removeSpaces(varargin)
   switch numel(varargin)
       case 2
           t = varargin{1};
           fig = varargin{2};
       case 4
           src = varargin{1};
           event = varargin{2};
           t = varargin{3};
           fig = varargin{4};
   end
   % Retrieve original units of t and fig
   fig_units = fig.Units;
   t_units = t.Units; 
   % adjust units
   fig.Units = 'centimeters';
   t.Units = 'normalized';
   %
   t.TileSpacing = 'none';
   drawnow;
   t.OuterPosition = [0 0 1 1];
    if t.InnerPosition(3)*fig.Position(3)>t.InnerPosition(4)*fig.Position(4)
        t.InnerPosition(3) = t.InnerPosition(4)*fig.Position(4)/fig.Position(3);
    else
        t.InnerPosition(4) = t.InnerPosition(3)*fig.Position(3)/fig.Position(4);
    end
    drawnow;
    t.PositionConstraint = 'outerposition';
%     drawnow;
    if (t.OuterPosition(3)>1)&&(t.OuterPosition(4)>1)
        t.OuterPosition(3:4) = [1 1];
        removeSpaces(t,fig);
    end
    if t.OuterPosition(3)>1
        t.OuterPosition(3) = 1;
        removeSpaces(t,fig);
    end
    if t.OuterPosition(4)>1
        t.OuterPosition(4) = 1;
        removeSpaces(t,fig);
    end
%     t.OuterPosition(1:2) = [0,0];
    t.OuterPosition(1) = 0.5*(1-t.OuterPosition(3));
    t.OuterPosition(2) = 0.5*(1-t.OuterPosition(4));
    drawnow;
%     t.InnerPosition(1) = 0.1; t.InnerPosition(3) = 0.8;
    % restore original units
    fig.Units = fig_units;
    t.Units = t_units;
end