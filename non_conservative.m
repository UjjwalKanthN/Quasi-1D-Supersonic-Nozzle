function [ sim_time_nc, th_mach_nc, th_press_nc, th_temp_nc, th_rho_nc, mach_no_nc, t_nc, rho_nc, p_nc, m_dot_nc ] = non_conservative( n, nt, x, dx, c, a, gamma, throat )

% Function for Non-Conservative Form

% Calculate Initial profiles
rho_nc = 1 - 0.3146*x; % Density
t_nc = 1 - 0.2314*x; % Temperature
v = (0.1 + 1.09*x).*t_nc.^0.5; % Velocity

% time steps
dt = min(c*dx./(t_nc.^0.5+v));

% Outer Time Loop
tic
for k = 1:nt
   
    rho_old = rho_nc;
    v_old = v;
    t_old = t_nc;
    
   % Predictor
    for j = 2:n-1
        dvdx = (v(j+1)-v(j))/dx;
        drhodx = (rho_nc(j+1) - rho_nc(j))/dx;
        dtdx = (t_nc(j+1)-t_nc(j))/dx;
        dlogadx = (log(a(j+1))-log(a(j)))/dx;
        % Continuity eq
        drhodt_p(j) = -rho_nc(j)*dvdx - rho_nc(j)*v(j)*dlogadx - v(j)*drhodx;
        % Momentum eq
        dvdt_p(j) = -v(j)*dvdx - (1/gamma)*(dtdx + (t_nc(j)/rho_nc(j)*drhodx));
        % Energy eq
        dtdt_p(j) = -v(j)*dtdx - (gamma-1)*t_nc(j)*(dvdx + v(j)*dlogadx)
        % Solution Update
        v(j) = v(j) + dvdt_p(j)*dt;
        t_nc(j) = t_nc(j) + dtdt_p(j)*dt;
        rho_nc(j) = rho_nc(j) + drhodt_p(j)*dt;
    end
    
    % Corrector
    for j = 2:n-1
        
        dvdx = (v(j)-v(j-1))/dx;
        drhodx = (rho_nc(j) - rho_nc(j-1))/dx;
        dtdx = (t_nc(j)-t_nc(j-1))/dx;
        dlogadx = (log(a(j))-log(a(j-1)))/dx;
        % Continuity eq
        drhodt_c(j) = -rho_nc(j)*dvdx - rho_nc(j)*v(j)*dlogadx - v(j)*drhodx;
        % Momentum eq
        dvdt_c(j) = -v(j)*dvdx - (1/gamma)*(dtdx + (t_nc(j)/rho_nc(j)*drhodx));
        % Energy eq
        dtdt_c(j) = -v(j)*dtdx - (gamma-1)*t_nc(j)*(dvdx + v(j)*dlogadx);
        
    end
    
    % Average Time Derivative
    dvdt = 0.5*(dvdt_p + dvdt_c);
    drhodt = 0.5*(drhodt_p + drhodt_c);
    dtdt = 0.5*(dtdt_p + dtdt_c);
    
    % Final Solution Update
    for i = 2:n-1
        v(i) = v_old(i) + dvdt(i)*dt;
        rho_nc(i) = rho_old(i) + drhodt(i)*dt;
        t_nc(i) = t_old(i) + dtdt(i)*dt;
    end
    
    % Apply boundary conditions
    % Inlet
    v(1) = 2*v(2) - v(3);
    
    % Outlet
    v(n) = 2*v(n-1) - v(n-2);
    rho_nc(n) = 2*rho_nc(n-1) - rho_nc(n-2);
    t_nc(n) = 2*t_nc(n-1) - t_nc(n-2);
    
    % Calculating Mass Flow Rate
    p_nc = rho_old.*t_old;
    mach_no_nc = v_old./(t_old).^0.5;
    m_dot_nc = rho_old.*a.*v_old;
    
    % Calculating values at throat
    th_mach_nc(k) = mach_no_nc(throat);
    th_press_nc(k) = p_nc(throat);
    th_rho_nc(k) = rho_nc(throat);
    th_temp_nc(k) = t_nc(throat);
    

% Plotting Mass Flow at different time steps
if n == 31
    
    mass_flow_nc = figure(2);
    hold on
    if k == 50
        plot(x, m_dot_nc, 'r')
    elseif k == 200
        plot(x, m_dot_nc, 'b')
    elseif k == 300
        plot(x, m_dot_nc, 'g')
    elseif k == 500
        plot(x, m_dot_nc, 'm')
        title('Non-Dimensional Mass Flow Distribution at Different Time Steps for Non-Conservative Form')
        legend('50 Deltat', '200 Deltat', '300 Deltat', '500 Deltat')
        xlabel('Non-Dimensional Length of Nozzle (x/l)')
        ylabel('Mass Flow Rate Ratio')
        grid minor
    end
end

end
sim_time_nc = toc
end