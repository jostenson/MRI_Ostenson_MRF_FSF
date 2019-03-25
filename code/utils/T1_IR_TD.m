function [ T1_map, RMSE_map ] = T1_IR_TD( img_stack, TI_v, TD, img_mask, lb, ub, plot_flag )

% T1 fitting for inversion recovery with time delay

% Input: img_stack = row x col x TI image stack (is magnitude images)
%        TI_v = vector of inversion times (time unit defines output time
%        unit; minimum unit is ms)
%        TD = delay time
%        img_mask = row x col mask to limit T1 fitting
%        lb = lower bound of parameters
%        ub = upper bound of parameters
%        plot_flag = if 1 then plot fit results
% Output: T1_map = row x col T1 map (time unit defined by TI_v above)
%         RMSE_map = row x col RMSE of fit

          
[n_row, n_col, n_TI] = size( img_stack );

if n_TI ~= numel(TI_v)
    error('T1_IR: I expect the same number of time points as images')
end

% fit water T1
T1_map = zeros( n_row, n_col );
RMSE_map = T1_map;
options = optimset('Display','off');

for ii = 1:n_row;
    for jj = 1:n_col;
        
        if img_mask(ii,jj) == 1
            
            my_sig_v = abs( squeeze( img_stack(ii,jj,:) ) );
            
            my_fun = @(x) abs( x(1) * ( x(2) * ( 1 - exp( -TD/x(3) ) ) * exp( -TI_v(:)./x(3) ) + 1 - exp( -TI_v(:)./x(3) ) ) ) - my_sig_v;
            
            ds = diff(my_sig_v);
            dTI = diff(TI_v);
            T10 = abs(dTI(1))*max(my_sig_v)*2/(abs(ds(1)));
            x0 = [max(my_sig_v) 2 T10];
            [my_solution, res_norm_sq] = lsqnonlin( my_fun, x0, lb, ub, options);
            
            T1 = my_solution(3);
            T1_map(ii,jj) = T1;
            RMSE_map(ii,jj) = sqrt( res_norm_sq/n_TI );
            
            if plot_flag == 1
                figure(99);
                plot(TI_v,my_fun(my_solution) + my_sig_v);
                hold on
                plot(TI_v,my_sig_v,'x')
                hold off
                xlabel('TI')
                ylabel('signal')
                title( sprintf( 'T1 estimate at row %d col %d is %.0f',ii,jj,T1) );
                drawnow
            else
                fprintf( 'T1 estimate at row %d col %d is %.0f\n',ii,jj,T1);
            end
            
        end
        
    end
end



end

