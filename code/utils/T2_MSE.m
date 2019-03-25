function [ T2_map, RMSE_map ] = T2_MSE( img_stack, TE_v, img_mask, lb, ub, n_p, plot_flag )

% T2 fitting using multiple spin-echoes

% Input: img_stack = row x col x TE image stack (is magnitude images)
%        TE_v = vector of echo times (time unit defines output time
%        unit)
%        img_mask = row x col mask to limit T2 fitting
%        lb = lower bounds of parameters
%        ub = upper bounds of parameters
%        plot_flag = if 1 then plot fit results
%        n_p = number of parameters in model
% Output: T2_map = row x col T1 map (time unit defined by TI_v above)
%         RMSE_map = row x col RMSE of fit

          
[n_row, n_col, n_TE] = size( img_stack );

if n_TE ~= numel(TE_v)
    error('T2_MSE: I expect the same number of time points as images')
end

lb = lb(1:n_p);
ub = ub(1:n_p);

% fit water T2
T2_map = zeros( n_row, n_col );
RMSE_map = T2_map;
options = optimset('Display','off');

for ii = 1:n_row;
    for jj = 1:n_col;
        
        if img_mask(ii,jj) == 1
            
            my_sig_v = abs( squeeze( img_stack(ii,jj,:) ) );
            
            if n_p == 3;
                
                my_fun = @(x) ( x(1) * exp(-TE_v(:)./ x(2) ) + x(3) ) - my_sig_v;
                
            elseif n_p == 2
            
                my_fun = @(x) ( x(1) * exp(-TE_v(:)./ x(2) ) ) - my_sig_v;
                
            end
            
            r = my_sig_v(1) / my_sig_v(2);
            x0_2 = ( TE_v(1) - r * TE_v(2) ) / ( 1 - r );
            if isnan(x0_2); x0_2 = 1; end
            x0 = [max(my_sig_v) x0_2 0];
            x0 = x0(1:n_p);
            [my_solution, res_norm_sq] = lsqnonlin( my_fun, x0, lb, ub, options);
            
            T2 = my_solution(2);
            T2_map(ii,jj) = T2;
            RMSE_map(ii,jj) = sqrt( res_norm_sq/n_TE );
            
            if plot_flag == 1
                figure(99);
                plot(TE_v,my_fun(my_solution) + my_sig_v);
                hold on
                plot(TE_v,my_sig_v,'x')
                hold off
                xlabel('TE')
                ylabel('signal')
                title( sprintf( 'T2 estimate at row %d col %d is %.0f',ii,jj,T2) );
                drawnow
            end
            
        end
        
    end
end



end

