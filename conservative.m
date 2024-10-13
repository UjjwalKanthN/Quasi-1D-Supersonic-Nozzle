function [ sim_time_c, th_mach_c, th_press_c, th_temp_c, th_rho_c, mach_no_c, t_c, rho_c, p_c, m_dot_c ] = conservative( n, nt, x, dx, c, a, gamma, throat )

% Function for Conservative Form

% Calculating Initial Profiles
for i = 1:n
    if (x(i) >= 0 && x(i) <= 0.5)
        rho_c(i) = 1;
        t_c(i) = 1;
    
    elseif (x(i) >= 0.5 && x(i) <= 1.5)
        rho_c(i) = 1 - 0.366*(x(i) - 0.5);
        t_c(i) = 1 - 0.167*(x(i) - 0.5);
    
    elseif (x(i) >= 1.5 && x(i) <= 3.5)
        rho_c(i) = 0.634 - 0.3879*(x(i) - 1.5);
        t_c(i) = 0.833 - 0.3507*(x(i) - 1.5);
    end
end

v_c = 0.59./(rho_c.*a);
p_c = rho_c.*t_c;

% time steps
dt = min(c*dx./(sqrt(t_c) + v_c));

% Defining Solution Vectors
U1 = rho_c.*a;
U2 = rho_c.*a.*v_c;
U3 = rho_c.*a.*((t_c/(gamma - 1)) + (gamma/2)*(v_c.^2));

% Outer Time Loop
tic
for k = 1:nt
    
    % Saving old values
    rho_old = rho_c;
    v_old = v_c;
    t_old = t_c;
    
    U1_old = U1;
    U2_old = U2;
    U3_old = U3;
    
    % Flux vectors
    F1 = U2;
    F2 = (U2.^2./U1) + ((gamma - 1)/gamma)*(U3 - (gamma*U2.^2)./(2*U1));
    F3 = gamma*(U2.*U3./U1) - (gamma*(gamma - 1)/2)*(U2.^3./U1.^2);
    ctr = 1;
    % Predictor
    for i = 2:n-1
        
        % Source term
        J2 = (1/gamma)*rho_c(i).*t_c(i).*((a(i + 1) - a(i))/dx);
        
        dU1dt_p(i) = -(F1(i+1) - F1(i))/dx;
        dU2dt_p(i) = (-(F2(i+1) - F2(i))/dx) + J2;
        dU3dt_p(i) = -(F3(i+1) - F3(i))/dx;
        
        %Solution Update
        U1(i) = U1(i) + dU1dt_p(i)*dt;
        U2(i) = U2(i) + dU2dt_p(i)*dt;
        U3(i) = U3(i) + dU3dt_p(i)*dt;
        
        rho_c(i) = U1(i)./a(i);
        v_c(i) = U2(i)./U1(i);
        t_c(i) = (gamma - 1)*((U3(i)./U1(i)) - (gamma/2)*(v_c(i)^2));
    end
    for i = 2:n-1
        F1(i) = U2(i);
        F2(i) = (U2(i).^2./U1(i)) + ((gamma - 1)/gamma)*(U3(i) - (gamma/2)*((U2(i).^2)./U1(i)));
        F3(i) = gamma*(U2(i).*U3(i)./U1(i)) - ((gamma*(gamma - 1))/2)*(U2(i).^3./U1(i).^2);
        
        ctr = ctr + 1;
    end
    
    % Corrector
    for i = 2:n-1
        % Source term
        J2 = (1/gamma).*rho_c(i).*t_c(i).*((a(i) - a(i - 1))/dx);
        
        dU1dt_c(i) = -(F1(i) - F1(i-1))/dx;
        dU2dt_c(i) = (-(F2(i) - F2(i-1))/dx) + J2;
        dU3dt_c(i) = -(F3(i) - F3(i-1))/dx;
    end
    
    % Average Time Derivative
    dU1dt = 0.5.*(dU1dt_p + dU1dt_c);
    dU2dt = 0.5.*(dU2dt_p + dU2dt_c);
    dU3dt = 0.5.*(dU3dt_p + dU3dt_c);
    
    % Final Soution Update
    for i = 2:n-1
        U1(i) = U1_old(i) + dU1dt(i)*dt;
        U2(i) = U2_old(i) + dU2dt(i)*dt;
        U3(i) = U3_old(i) + dU3dt(i)*dt;
    end
    
    % Apply Boundary Conditions
    % Inlet
    U2(1) = 2*U2(2) - U2(3);
    
    %Outlet
    U1(n) = 2*(U1(n-1)) - U1(n-2);
    U2(n) = 2*(U2(n-1)) - U2(n-2);
    U3(n) = 2*(U3(n-1)) - U3(n-2);
    
    % Flow field variables update
    rho_c = U1./a;
    v_c = U2./U1;
    t_c = (gamma - 1)*((U3./U1) - (gamma/2)*v_c.^2);
    
    % Calculating Mass Flow Rate
    p_c = rho_old.*t_old;
    mach_no_c = v_old./(t_old).^0.5;
    m_dot_c = rho_old.*a.*v_old;
    
    % Calculating values at throat
    th_mach_c(k) = mach_no_c(throat);
    th_press_c(k) = p_c(throat);
    th_rho_c(k) = rho_c(throat);
    th_temp_c(k) = t_c(throat);


% Plotting Mass Flow at different time steps
if n == 31
    
    mass_flow_c = figure(1);
    hold on
    if k == 100
        plot(x, m_dot_c, 'r')
    elseif k == 200
        plot(x, m_dot_c, 'b')
    elseif k == 300
        plot(x, m_dot_c, 'g')
    elseif k == 550
        plot(x, m_dot_c, 'm')
        title('Non-Dimensional Mass Flow Distribution at Different Time Steps for Conservative Form')
        legend('100 Deltat', '200 Deltat', '300 Deltat', '550 Deltat')
        ylim([0.5 0.7])
        xlabel('Non-Dimensional Length of Nozzle (x/l)')
        ylabel('Mass Flow Rate Ratio')
        grid minor
    end
end

end
sim_time_c = toc
end