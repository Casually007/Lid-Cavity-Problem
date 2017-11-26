% 2D lid-cavity problem

n_x = 100;
n_y = 100;
xlength = 1.0;
ylength = 1.0;
t_max = 25.1;
Re = 10000;

del_x = xlength/n_x;
del_y = ylength/n_y;

% Other variables
lambd = (1/del_x^2 + 1/del_y^2)*2;
Beta = 1.7; 
gamma = 0.9; 
tau = 0.5;
eps = 5e-6;
iter_max = 500;

% Pre allocate  u,v,p etc...
% at t = 0
u = zeros(n_x + 2, n_y + 2);
v = zeros(n_x + 2, n_y + 2);
pn = zeros(n_x + 2, n_y + 2);

F = zeros(n_x + 2, n_y + 2);
G = zeros(n_x + 2, n_y + 2);
S_ij = zeros(n_x + 2, n_y + 2);

x = zeros(n_x, n_y);
y = zeros(n_x, n_y);

% Visualization
u_avg = zeros(n_x + 1, n_y + 1);
v_avg = zeros(n_x + 1, n_y + 1);
u_mag = zeros(n_x + 1, n_y + 1);
u_avg(:, n_y + 1) = 1;
X = linspace(0, 1, n_x + 1);
Y = linspace(0, 1, n_y + 1);
fig_n = 5;
delimiterIn = ' ';
headerlinesIn = 1;
res_new= 0;
Max_div_new = 0;
count_fig = 0;
del_t = tau*(Re/lambd);
t = 0;
n = 1;

% Main loop
while (t < t_max)
   

  % North & East wall
    for i = 3:n_x
        u(i, n_y + 2) = 2 - u(i, n_y + 1);
        v(i, n_y + 1) = 0;
        u(n_x + 1, i) = 0;
        v(n_x + 2, i) = -v(n_x + 1, i);
    end  
  % South & West wall
    for i = 2:n_x+1
        u(1, i) = 0;
        v(1, i) = -v(2, i);
        v(i, 1) = 0;
        u(i, 1) = -u(i, 2);
        F(1, i) = u(1, i);
        G(i, 1) = v(i, 1);
    end
     
  % North-east corner
    
    v(n_x + 1, n_y + 1) = 0;
    u(n_x + 1, n_y + 1) = 0;
        
  % North-west corner
    
    v(2, n_y + 1) = 0;
    u(1, n_y + 1) = 0;
    u(2, n_y + 2) = 2 - u(2, n_y + 1);
      
  % South-east corner generate a valid file identifier.

    u(n_x + 1, 2) = 0;
    v(n_x + 1, 1) = 0;
    v(n_x + 2, 2) = -v(n_x + 1, 2);
     
  % Time step control  
    max_u = 0;
    max_v = 0;
    for i = 2:n_x + 1
        for j = 2:n_y + 1
            if max_u < abs(u(i,j))
                max_u = abs(u(i,j));
            end
            if max_v < abs(v(i,j))
                max_v = abs(v(i,j));
            end
        end
    end
    if (max_u >= max_v)
        max = max_u;
    else
        max = max_v;
    end
    if (max ~= 0 && del_t > tau*(del_x/max))
        del_t = tau*(del_x/max);
    end
    
  % Visualization
    if(floor(t/5) == n)
        n = n + 1;
        disp(t)
        count_fig = count_fig + 1;
    % Stream function
        q = zeros(n_x + 1, n_y + 1);
        omg = zeros(n_x - 1,n_y - 1);
        for j = 2:n_x
             for i = 2:n_y
                 omg(i-1, j-1) = (v(i + 1,j) - v(i,j))/del_x  -...
                                        (u(i,j + 1) - u(i,j))/del_y ;
             end
        end
        count_q = 0;
        while(true)
             for j = 2:n_x
                for i = 2:n_y
                    q(i, j) =  (Beta)*...
                          (omg(i - 1, j - 1)/lambd + (1/(lambd*del_x^2))...
                          *(q(i + 1, j) + q(i - 1, j)) + ...
                          (1/(lambd*del_y^2))*(q(i, j + 1) + ...
                           q(i, j - 1))) +(1 - Beta)*q(i, j);
                end
             end
             res_q = 0;
        
             for i = 2:n_x
                  for j = 2:n_y
                     temp = abs(-q(i, j) + omg(i - 1, j - 1)/lambd + ...
                         (1/(lambd*del_x^2))*(q(i + 1, j) + ...
                         q(i - 1, j)) + (1/(lambd*del_y^2))* ...
                         (q(i, j + 1) + q(i, j - 1)));
                     res_q = res_q + temp^2;
                  end
             end
            count_q = count_q + 1;
            res_q = (res_q/(n_x*n_y))^0.5;
%             disp(res_q);
            if (res_q <= eps || count_q >= iter_max)
                break;
            end
        end
        
    % Avg. Velocity        
         for j = 2:n_y
            for i = 2:n_x
                u_avg(i , j) = (u(i,j) + u(i,j + 1))/2;
                v_avg(i , j) = (v(i,j) + v(i + 1,j))/2;
            end
         end
     % Streamlines
         figure(fig_n*(count_fig - 1) + 1)
         contourf(X,Y,q',50,'k-')
         title([' Streamlines for Re = ' num2str(Re) ' at t = ' ...
             num2str(t) ' sec ']);
         plotname = sprintf('Streamlines_@_Re_%d_@_%.2f_sec.eps',Re,t);
         saveas(gca, plotname,'epsc');

   
    % Velocity field
         u_mag = sqrt(u_avg.^2 + v_avg.^2);
         figure(fig_n*(count_fig - 1) + 2)
         quiver(X,Y,(u_avg./u_mag)',(v_avg./u_mag)',.4,'k-')
         title([' Normalized velocity field for Re = ' num2str(Re) ...
             ' at t = ' num2str(t) ' sec ']);
         xlim([0 1]);
         ylim([0 1]);
         plotname = sprintf('Velocity_field_@_Re_%d_@_%.2f_sec.eps',Re,t);
         saveas(gca, plotname,'epsc');        

     % Velocity magnitude 
         figure(fig_n*(count_fig - 1) + 3)
         contour(X,Y,u_mag',20,'k-')
         title([' Velocity magnitude for Re = ' num2str(Re) ' at t = '...
             num2str(t) ' sec ']);
         plotname = sprintf('Velocity_magnitude_@_Re_%d_@_%.2f_sec.eps'...
             ,Re,t);
         saveas(gca, plotname,'epsc');           
         
    % Top wall shear 
        figure(fig_n*(count_fig - 1) + 4)
        scatter(X,(u_avg( :, n_y + 1) - u_avg( :, n_y))./(del_y), 3, 'k');
        title([' Top wall shear at ' num2str(Re) ' at t = ' num2str(t)...
            ' sec ']);
        xlabel ('x-location from left corner');
        ylabel ('Wall shear (Non-dimensional)');
        plotname = sprintf('Wall_shear_@_Re_%d_@_%.2f_sec.eps',Re,t);
        saveas(gca, plotname,'epsc');   
    
    % Midsection velocities @ x = 0.5
        figure(fig_n*(count_fig - 1) + 5)
        scatter(u_avg(floor((n_y + 1)/2), :), X, 3, 'k');
        title([' u (0.5,y) @ t = ' num2str(t) ' sec ']);
        xlabel ('u (0.5,y)');
        ylabel ('y');    
        plotname = sprintf('u(0.5,y)_@_Re_%d_@_%.2f_sec.eps',Re,t);
        saveas(gca, plotname,'epsc');      
    end
    
  
    % Main loop  
    for j = 2:n_y+1
        for i = 2:n_x+1
            
             u_e = (u(i, j) + u(i + 1, j))/2 ;
             u_w = (u(i, j) + u(i - 1, j))/2 ;
             u_n = (u(i, j) + u(i, j + 1))/2 ;
             u_s = (u(i, j) + u(i, j - 1))/2 ;
                
             v_e = (v(i, j) + v(i + 1, j))/2 ;
             v_w = (v(i, j) + v(i - 1, j))/2 ;
             v_n = (v(i, j) + v(i, j + 1))/2 ;
             v_s = (v(i, j) + v(i, j - 1))/2 ;
                

                
             CONV_u_ij = (u_e^2 - u_w^2)/del_x + gamma*...
                       (((u(i, j) - u(i + 1, j))/2)*abs(u_e) -...
                       ((u(i - 1, j) - u(i, j))/2)*abs(u_w))/del_x +...
                       (u_n*v_n - u_s*v_s)/del_y +...
                       gamma*(((u(i, j) - u(i, j + 1))/2)*abs(v_n) -...
                       ((u(i, j - 1) - u(i, j))/2)*abs(v_s))/del_y;
                            
                        
             CONV_v_ij = (v_n^2 - v_s^2)/del_y + gamma*...
                       (abs(v_n)*((v(i, j) - v(i, j + 1))/2) - ...
                       abs(v_s)*((v(i, j - 1) - v(i, j))/2))/del_y +...
                       (u_e*v_e - u_w*v_w)/del_x +...
                       gamma*(abs(u_e)*((v(i, j) - v(i + 1, j))/2) -...
                       abs(u_w)*((v(i - 1, j) - v(i, j))/2))/del_x;
                        
                
             DIFF_u_ij = (1 / Re)*( ...
                       (u(i + 1, j) - 2*u(i, j) + u(i - 1, j))/del_x^2 +...
                       (u(i, j + 1) - 2*u(i, j) + u(i, j - 1))/del_y^2);
               
             DIFF_v_ij = (1 / Re)*( ...
                       (v(i + 1, j) - 2*v(i, j) + v(i - 1, j))/del_x^2 +...
                       (v(i, j + 1) - 2*v(i, j) + v(i, j - 1))/del_y^2);
                        
             F(i, j) = u(i, j) + del_t*( DIFF_u_ij - CONV_u_ij); 
             G(i, j) = v(i, j) + del_t*( DIFF_v_ij - CONV_v_ij);
                            
             S_ij(i,j) = (1/del_t)*((F(i, j) - F(i - 1, j))/del_x + ...
                         (G(i, j) - G(i, j - 1))/del_y);

         end
     end
     count = 0;
     
   % Pressure   
     while(true)
         for j = 2:n_y+1
             for i = 2:n_x+1                
               pn(i, j) = (Beta)*...
                          (-S_ij(i, j)/lambd + (1/(lambd*del_x^2))* ...
                          (pn(i + 1, j) + pn(i - 1, j)) + ...
                          (1/(lambd*del_y^2))*(pn(i, j + 1) + ...
                           pn(i, j - 1))) +(1 - Beta)*pn(i, j);
             end
         end
        
         for i = 2:n_x+1
                pn(i, 1) = pn(i, 2);
                pn(i, n_y + 2) = pn(i, n_y + 1);
                pn(1, i) = pn(2, i);
                pn(n_x + 2, i) = pn(n_x + 1, i);
         end
        
         res_new = 0;
        
         for i = 2:n_x+1
             for j = 2:n_y+1
                 temp = abs(pn(i, j) + S_ij(i, j)/lambd - ...
                       (1/(lambd*del_x^2))*(pn(i + 1, j) + ...
                       pn(i - 1, j)) - (1/(lambd*del_y^2))* ...
                       (pn(i, j + 1) + pn(i, j - 1)));
                 res_new = res_new + temp^2;
             end
         end
         count = count + 1;
         res_new = (res_new/(n_x*n_y))^0.5;
         if (res_new <= eps || count >= iter_max)
             break;
         end
     end
     
   % Velocity update 
     for j = 2:n_y + 1
         for i = 2:n_x 
             u(i, j) = F(i, j) - (del_t/del_x)* ...
                       (pn(i + 1, j) - pn(i, j));
         end
     end
        
     for j = 2:n_y 
         for i = 2:n_x + 1
             v(i, j) = G(i, j) - (del_t/del_y)* ...
                       (pn(i , j + 1) - pn(i, j));
         end
     end
     
  % Value storage
     t = t + del_t;
end