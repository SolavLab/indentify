function Wn = getWeight(pos_data,varargin)
%% Initialize variables and parameters
%Get nodal weights with optional indenter obstruction consideration
%getWeight(pos_data) gets the appropriate weights for all nodes
%getWeight(pos_data, type,nodeList) will use data to block out obscured nodes

p = inputParser; % Create an input parser object

n = size(pos_data, 1); % Determine the number of nodes
defaultNodeList = true(1, n); % Default nodeList as a logical array of ones
defaultType = 1; %Distance from indenter is factor

% Define the input scheme
addOptional(p,'type', defaultType);
addOptional(p, 'nodeList', defaultNodeList); % Add parameter for node list
parse(p,varargin{:}); % Parse the input arguments

% Extract results
nodeList = p.Results.nodeList; % Get the list of nodes to consider
type = p.Results.type;


rel_pos = pos_data(:,:,1)';

% Find average distance between nodes to asses density 
dMat = pdist(rel_pos);
dMat(dMat==0) = NaN;
node_density = mean(min(dMat)); % Calculate the average minimum distance
% Add a small fraction of the node density to the distances to prevent
% infinite values when two points coincide. This fraction (node_density/10)
% is small enough not to significantly affect the overall distance values,
% but large enough to prevent computational errors due to division by zero.

distances = vecnorm(rel_pos)+node_density/10;

Wn = zeros(1,n);

%% Calculate nodal weights
switch type
    case 1 %Distance from indenter is factor
        Wn(nodeList) = 1./distances(nodeList);
        
    case 2 %Distance from x symmetry line
        Wn(nodeList) = 1./((0.1*rel_pos(1,nodeList).^2+rel_pos(2,nodeList).^2).^(0.5)+node_density/10);

    case 3 %Distance from y symmetry line
        Wn(nodeList) = 1./((0.1*rel_pos(2,nodeList).^2+rel_pos(1,nodeList).^2).^(0.5)+node_density/10);

end
Wn = Wn(nodeList);
end