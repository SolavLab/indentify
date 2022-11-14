function [N,D,B,K,W] = getXYCoeffMat(n_accuracy)
    syms x y hx hy f(x,y)
    n=n_accuracy;
    df_dx_dy(x,y) = diff_t(diff_t(f,y,hy,n),x,hx,n);
    pos_arr = -(n/2):1:(n/2);
    for ii=1:length(pos_arr)
        for jj=1:length(pos_arr)
            W(ii,jj) = double(diff(df_dx_dy(x,y)*hx*hy,...
                f(x+hx*pos_arr(ii),y+hy*pos_arr(jj))));
            [N(ii,jj) D(ii,jj)] = rat(W(ii,jj));
        end
    end
    K = max(D,[],'all');
    B =  N./D*K;
end

 function [dgdt] = diff_t(g0,t,dt,n)
        g(t) = g0;
        %     differentiate by t
        switch n
            case 2
                indices = -1:1;
                W = [-1/2,0,1/2];
            case 4
                indices = -2:2;
                W = [1/12, -2/3, 0, 2/3, -1/12];
            case 6
                indices = -3:3;
                W = [-1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60];
            case 8
                indices = -4:4;
                a = [1/280 -4/105, 1/5, -4/5];
                W = [a, 0, -flip(a)];
        end
    %     u = zeros(1,length(indices));
        for kk=1:length(indices)
            u(kk) = g(t+indices(kk)*dt);
        end
        dgdt = u*W'/dt;
 end
    
%   function [N,D] = getCoeffs(f_name,n)
% %             df_dx(x,y) = diff_t(f,x,hx);
%             df_dy(x,y) = diff_t(f,y,hy);
% %             df_dx2(x,y) = diff_t(df_dx,x,hx);
% %             df_dy2(x,y) = diff_t(df_dy,y,hy);
%             df_dx_dy(x,y) = diff_t(df_dy,x,hx);
%             [N,D] = rat(diff(df_dx_dy(x,y)*hx*hy, f_name));
%     end  