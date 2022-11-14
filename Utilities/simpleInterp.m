%% A weighted average between elements of v1 and v2 (0<=eta<=1)
% Ffu = simpleInterp(Ff,Fu,eta)
function [val] = simpleInterp(v1,v2,eta)
    val = eta*v1+(1-eta)*v2;
end