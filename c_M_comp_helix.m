% In this program, we compute the value of c_M for the helical polygon
% case: 
% We follow the same method as inthe paper of Patxi-Vega 2017 by solving equation (60) in
%  We solve the system by using numerical method  RK4. 
% Date: November 18, 2018

 
tic 
% --------------------------------------
% PARAMETERS
% --------------------------------------
% clear ;

M = 6 ; 
L = 100 ; 
ds = 1e-3 ; 
s = 0 : ds : L ; 
% s = -s ; 
% ds = -ds ; 
N = length(s) ; 
theta = pi / 5; 
% theta = 0 ; 
b = tan(theta/2) / tan(pi/M) ;
% b = 0 ;
a = sqrt(1-b^2) ; 
rho = 2 * asin(a*sin(pi/M)) ; 
c0 = sqrt(-2*log(cos(rho/2))/pi) ; 

A1 = exp(-c0^2*pi/2) ; 
% A2 = 0.0279981163264246732081285568329364799690135848296272 ;
% A2 = 0.2484930367622096407477195713677620388301578808756232; %  M=5, theta = pi/3
A2 = 0.2861514123131838643236976725862559887127499315927088 ; % c0 = , M=6, theta=pi/5 
% A2 = 0.394884119206045929486816155514483564257765726094998190987 ; % M = 5 , theta = pi / 30 

% A2 = 0.504325759989162712692027090261232631491835472587705334513 ; % M =3, theta = pi/3

% A2 = 0.395784879281939 ; % M = 5 , b= 0 ; 
A3 = sqrt(1-exp(-c0^2*pi)-A2^2); 

% u_rot_alg = pi * c0^2 / sqrt(exp(pi*c0^2)-1) ; 


% ----------------------------------------
% Initial data for the third component X_3
% ----------------------------------------

R1 = [ 1 0 0 ; 
       0 A2/sqrt(A2^2+A3^2) A3/sqrt(A2^2+A3^2) ; 
       0 -A3/sqrt(A2^2+A3^2) A2/sqrt(A2^2+A3^2) ];

R3 = [cos(pi/M) sin(pi/M) 0 ; -sin(pi/M) cos(pi/M) 0 ; 0 0 1] ; 

R2 = [1+cos(pi/M)^2*(a*cos(theta/2)-1) cos(theta/2) * sin(pi/M-rho/2) b*cos(theta/2) ; 
        cos(theta/2)*sin(pi/M-rho/2) 1+sin(pi/M)^2*(a*cos(theta/2)-1) -b*tan(rho/2)/a ; 
        -b*cos(theta/2)     b*tan(rho/2)/a  a*cos(theta/2)  ] ; 
    
R = R2\R3*R1 ;     

%  time t = t_{1,20}
N_t = 151200*4^2 ; 
t = linspace(0,2*pi/M^2,N_t+1); 
t1 = t(1:N_t/20) ; 
X0 = [ 0 ; 0 ; 2*c0] ;                  % solution at t = 1 
X_rot0 = [-a*pi/M ; -a*(pi/M) / (tan(pi/M)) ; 0 ] + R * X0 * sqrt(t1); 

% Xrot1= (X_rot0* sqrt(t1)).'; 

Xrot = (-sqrt(Xrot1(:,1).^2+Xrot1(:,2).^2) + 1i *Xrot1(:,3)); 

return; 

X_rot1 = R * [1 ; 0 ; 0] ; % First derivative of X
X_rot2 = c0* R * [0; 1 ; 0] ; % Second derivative 


u1 = 0 ;                             % H_rot(0) 
u2 = X_rot0(3)  ;                    % H_rot'(0) 
u3 = X_rot1(3)  ;                    % H_rot''(0) 
u4 = X_rot2(3) ;                     % H_rot'''(0) 

% % -----------------------------------------
% % Initial data for the third component X_2
% % -----------------------------------------
% u1 = 0 ;                                % H_rot(0) 
% u2 = -pi /(M*tan(pi/M)) + 2 * A3 * c0 * cos(pi/M) / sqrt(A2^2+A3^2) ;    % H_rot'(0) 
% u3 = -sin(pi/M) ;                                % H_rot''(0) 
% u4 = A2 * c0 * cos(pi/M) / sqrt(A2^2+A3^2) ;       % H_rot'''(0) 
 
return ; 

u = [u1 ; u2 ; u3 ; u4] ;

u_n = u ; 
% RK4 method
u_N(:,1) = u_n ; 


for n = 0 : N-2
    
    s_n = n*ds ; 
    A = [ u_n(2) ; u_n(3) ; u_n(4) ; s_n*u_n(2)/4 - (c0^2+s_n^2/4)*u_n(3) ] ;
    uA = u_n + ds * 0.5 * A ; 
    
    s_n1 = s_n + 0.5 * ds ; 
    B = [ uA(2) ; uA(3) ; uA(4) ; s_n1*uA(2)/4 - (c0^2+s_n1^2/4)*uA(3) ] ; 
    uB = u_n + ds * 0.5 * B ; 
    
    C = [ uB(2) ; uB(3) ; uB(4) ; s_n1*uB(2)/4 - (c0^2+s_n1^2/4)*uB(3) ] ; 
    uC = u_n + ds * C ; 
    
    s_n2 = s_n + ds ; 
    D = [ uC(2) ; uC(3) ; uC(4) ; s_n2*uC(2)/4 - (c0^2+s_n2^2/4)*uC(3) ] ; 
    
    u_n = u_n + ds * (A + 2*B + 2*C + D) / 6 ; 
    u_N(:,n+2) = u_n ; 
    
end

toc 