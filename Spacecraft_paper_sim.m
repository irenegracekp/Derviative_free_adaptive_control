global J r i omega_1 omega_2 omega_3 torque;
J = diag([114 86 87]);
r = [0,0,0,0,0,0]';
i = 1;
omega_1 = 0;
omega_2 = 0;
omega_3 = 0;
torque = 0;

tspan = [0,100];

x0 = zeros(48,1) ;
x0(1:6) = [-0.1 0.5 1 0 0 0];

[t_lc_nlp, x_lc_nlp] = ode45(@LC_NonlinearModel, tspan, x0);


t = tiledlayout(3,1); % Requires R2019b or later
nexttile
plot(t_lc_nlp, x_lc_nlp(:,1));
legend('x_{1}(t)');
grid on
nexttile
plot(t_lc_nlp, x_lc_nlp(:,2));
legend('x_{2}(t)');
ylabel('Magnitude');
grid on
nexttile
plot(t_lc_nlp, x_lc_nlp(:,3));
legend('x_{3}(t)');
xlabel('t');
grid on

t.Padding = 'compact';
t.TileSpacing = 'compact';

function xdot = LC_NonlinearModel(t,x)
    global J r i omega_1 omega_2 omega_3 torque;
    alpha = 1;
    alpham = 1;
    wn = 0.02;
    zeta = 1;
    
    if (i<300)
        r = [0,0,0,0,0,0]';
    elseif (i<600)
         r = [0.2,0.2,0.2,0,0,0]';
    else 
        r = [-0.2,-0.2,-0.2,0,0,0]';
    end
    
    xm = r; 
    
    um = [0 0 0]' ;
    Am = [zeros(3) eye(3);
          -wn^2*eye(3) -2*zeta*wn*eye(3);] ;
    Bm = [zeros(3); wn^2*eye(3);];
    Cm = [alpham*eye(3) eye(3)] ;
    
    ym = (Cm*xm)' ;
    
    sigma(1:3,1) = x(1:3);
    sigma_dot(1:3,1) = x(4:6);
             
    Skew_sigma = [0         -sigma(3)  sigma(2);...
                 sigma(3)  0          -sigma(1);...
                 -sigma(2) sigma(1)   0       ;];
             
    Skew_sigma_dot = [0         -sigma_dot(3)  sigma_dot(2);...
                 sigma_dot(3)  0          -sigma_dot(1);...
                 -sigma_dot(2) sigma_dot(1)   0       ;];
          
    T_sigma =  (1/4)*( (( 1 - sigma'*sigma)*eye(3)) + (2*Skew_sigma) + (2*(sigma*sigma')));
    
    H_sigma = inv(T_sigma')*J*inv(T_sigma) ;
    T_dot_sigma = 1/4 * ( -2*sigma'*sigma_dot*eye(3) + 2*Skew_sigma_dot + 4*sigma*sigma_dot' ) ;
    L = J*inv(T_sigma)*sigma_dot ;
    L_skew = [0 -L(3) L(2);
               L(3) 0 -L(1);
              -L(2) L(1) 0;] ;
    C_sigmasigmadot = -inv(T_sigma')*( J*inv(T_sigma)* T_dot_sigma*inv(T_sigma) + (L_skew)*inv(T_sigma) ) ;
    
    A = [ zeros(3) eye(3);
          zeros(3) -inv(H_sigma)*C_sigmasigmadot;] ;
    B = [zeros(3); inv(H_sigma)];
    C = [alpha*eye(3) eye(3)]; 
    
    y = (C*x(1:6))';
    ey = (ym - y)' ;
    
    ke = [x(13:15)';x(16:18)';x(19:21)';];
    kx = [x(22:27)';x(28:33)';x(34:39)';];
    ku = [x(40:42)';x(43:45)';x(46:48)';];
    
    ke_t = ke + ey*ey'*10^5*eye(3);
    kx_t = kx + ey*xm'*10^3*eye(6);
    ku_t = ku + ey*um'*10^3*eye(3);
    
    u = ke_t*ey + kx_t*xm + ku_t*um ;
  
    torques = T_sigma'*u;
    torque(i) = torques(1);
     omega = inv(T_sigma)*sigma_dot ; 
    omega_1(i) = omega(1);
    omega_2(i) = omega(2);
    omega_3(i) = omega(2);
    i = i+1;
    
    
    ke_dot = ey*ey'*10^5*eye(3);
    kx_dot = ey*xm'*10^3*eye(6);
    ku_dot = ey*um'*10^3*eye(3);
    
    
    xdot(1:3) = sigma_dot ;
    xdot(4:6) = ((-inv(H_sigma)*C_sigmasigmadot)*sigma_dot + inv(H_sigma)*u)' ;
    
    xdot(7:9) = x(10:12);
    xdot(10:12) =  -wn^2*eye(3)*x(7:9) + -2*zeta*wn*eye(3) * x(10:12) + wn^2*eye(3)*um ;
    
    xdot(13:15) = ke_dot(1,:); xdot(16:18) = ke_dot(2,:); xdot(19:21) = ke_dot(3,:);
    xdot(22:27) = kx_dot(1,:); xdot(28:33) = kx_dot(2,:); xdot(34:39) = kx_dot(3,:);
    xdot(40:42) = ku_dot(1,:); xdot(43:45) = ku_dot(2,:); xdot(46:48) = ku_dot(3,:);
    
    xdot = xdot';
end



