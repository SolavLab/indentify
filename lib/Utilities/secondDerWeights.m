% n - order of error (+O(dt^n))
% W - function weights from left (minus) to right (plus)

% df2dx2 = [f(t-n*dt), f(t-(n+1)*dt,.... f(t), ...., f(t+(n-1)*dt),
% f(t+n*dt)]*W'/(dt^2)
function [W] = secondDerWeights(n)
    switch n
        case  2
            W = [1, -2, 1];
        case 4
            W = [-1/12, 4/3, -5/2, 4/3, -1/12];
        case 6
            W = [1/90, -3/20, 3/2, -49/18, 3/2, -3/20, 1/90];
        case 8
            W = [-1/560, 8/315, -1/5, 8/5, -205/72, 8/5, -1/5, 8/315, -1/560];
    end
end
