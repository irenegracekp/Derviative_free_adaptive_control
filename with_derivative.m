global J r W i omegad_1 omegad_2 omegad_3;
J = diag([114 86 87]);
r = [0,0,0,0,0,0]';
W = zeros(6,3) ;
i = 1;
omegad_1 = 0;
omegad_2 = 0;
omegad_3 = 0;

tspan = [0,1000];

x0 = zeros(66,1) ;
x0(1:6) = [-0.1 0.5 0.8 0 0 0];
x0(49:66) = ones(1,18);

[t_lc_nlp, x_lc_nlp] = ode45(@LC_NonlinearModel, tspan, x0);

t = tiledlayout(3,1); % Requires R2019b or later
nexttile
hold on
plot(t_lc_nlp, x_lc_nlp(:,1));
plot(tdf_lc_nlp, xdf_lc_nlp(:,1));
plot(R_1);
legend('Sigma_{1}(t)');
hold off
grid on
nexttile
hold on
plot(t_lc_nlp, x_lc_nlp(:,2));
plot(tdf_lc_nlp, xdf_lc_nlp(:,2));
plot(R_2);
legend('Sigma_{2}(t)');
ylabel('Magnitude');
hold off
grid on
nexttile
hold on
plot(t_lc_nlp, x_lc_nlp(:,3));
plot(tdf_lc_nlp, xdf_lc_nlp(:,3));
legend('Sigma_{3}(t)');
plot(R_3);
hold off
xlabel('t');
grid on

t.Padding = 'compact';
t.TileSpacing = 'compact';

function xdot = LC_NonlinearModel(t,x)
    global J r W i omegad_1 omegad_2 omegad_3;
    alpha = 1;
    alpham = 1;
    wn = 0.02;
    zeta = 0.7;
    
    if (t<200)
        r = [0,0,0,0,0,0]';
    elseif (t<600)
         r = [0.2,0.2,0.2,0,0,0]';
    else 
        r = [-0.2,-0.2,-0.2,0,0,0]';
    end
    
    if (t<500)
        J = diag([114 86 87]);
    else 
        J = diag([350 86 87]);
    end
    xm = r; 
    
    um = [0 0 0]' ;
    Am = [zeros(3) eye(3);
          -wn^2*eye(3) -2*zeta*wn*eye(3);] ;
    Bm = [zeros(3); wn^2*eye(3);];
    Cm = [alpham*eye(3) eye(3)] ;
    
    P = lyap(Am',Am,eye(6)) ;
    
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
    ex = (x(1:6)-xm);
    
    ke = [x(13:15)';x(16:18)';x(19:21)';];
    kx = [x(22:27)';x(28:33)';x(34:39)';];
    ku = [x(40:42)';x(43:45)';x(46:48)';];
    
    W = [x(49:54) x(55:60) x(61:66)];
    
    ke_t = ke + ey*ey'*10^5*eye(3);
    kx_t = kx + ey*xm'*10^3*eye(6);
    ku_t = ku + ey*um'*10^3*eye(3);
    
    
    %uncertainity
    
    deltax = W' * beta(sigma,sigma_dot) ;
    deltaxt = B*deltax ;
    
    phi_1 = 0.1 ;
    
    Wdash = 1.5*beta(sigma,sigma_dot)*ex'*P*B ;
       
    k1 = B\(Am - A);
    k2 = B\Bm;
    
    u = 0.55*k1*ex + k2*ey -  W'*beta(sigma,sigma_dot);
    
    omegad = inv(T_sigma)*sigma_dot ; 
    omegad_1(i) = omegad(1);
    omegad_2(i) = omegad(2);
    omegad_3(i) = omegad(2);
    i = i+1;
    
    
    ke_dot = ey*ey'*10^5*eye(3);
    kx_dot = ey*xm'*10^3*eye(6);
    ku_dot = ey*um'*10^3*eye(3);
    
    xdot(1:3) = sigma_dot + deltaxt(1:3) ;
    xdot(4:6) = ((-inv(H_sigma)*C_sigmasigmadot)*sigma_dot + inv(H_sigma)*u)' + deltaxt(4:6)';
    
    xdot(7:9) = x(10:12);
    xdot(10:12) =  -wn^2*eye(3)*x(7:9) + -2*zeta*wn*eye(3) * x(10:12) + wn^2*eye(3)*um ;
    
    xdot(13:15) = ke_dot(1,:); xdot(16:18) = ke_dot(2,:); xdot(19:21) = ke_dot(3,:);
    xdot(22:27) = kx_dot(1,:); xdot(28:33) = kx_dot(2,:); xdot(34:39) = kx_dot(3,:);
    xdot(40:42) = ku_dot(1,:); xdot(43:45) = ku_dot(2,:); xdot(46:48) = ku_dot(3,:);
    xdot(49:54) = Wdash(:,1) ;
    xdot(55:60) = Wdash(:,2) ;
    xdot(61:66) = Wdash(:,3) ;
    xdot = xdot';
   
end

function ret = beta(sigma,sigma_dot) 

    ret(1) = 1- exp(-sigma(1)) / ( 1 + exp(-sigma(1)) ) ; 
    ret(2) = 1- exp(-sigma(2)) / ( 1 + exp(-sigma(2)) ) ;
    ret(3) = 1- exp(-sigma(3)) / ( 1 + exp(-sigma(3)) ) ;
    ret(4) = 1- exp(-sigma_dot(1)) / ( 1 + exp(-sigma_dot(1)) ) ;
    ret(5) = 1- exp(-sigma_dot(2)) / ( 1 + exp(-sigma_dot(2)) ) ;
    ret(6) = 1- exp(-sigma_dot(3)) / ( 1 + exp(-sigma_dot(3)) ) ;
    ret = ret';
end
