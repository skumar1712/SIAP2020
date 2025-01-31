% The following program computes a rotation matrix which rotates vector 
% Ap = [A1 A2 A3] and Am = [A1 -A2 A3] in such a way that they match with
% the M - corner problem' initial tangent vectors 
% The algorithm uses only two matrices 

% Date: November 21, 2018 
close all; 
clear 

M = 5 ; 
theta = pi / 3; 

b = tan(theta/2) / tan(pi/M) ;

a = sqrt(1-b^2) ; 
rho = 2 * asin(a*sin(pi/M)) ; 
c0 = sqrt(-2*log(cos(rho/2))/pi) ; 

A1 = exp(-c0^2*pi/2) ; 
A2 = 0.2484930367622096407477195713677620388301578808756232; %  M=5, theta = pi/3
A3 = sqrt(1-exp(-c0^2*pi)-A2^2); 

%  ------------------------------------------
% The tangent vector for 1-corner problem 
Am =[A1 ; -A2 ; -A3] ; 
Ap =[A1 ; A2 ; A3] ; 

%  ------------------------------------------
% The tangent vector for M-corner problem 
Tm = [a * cos(2*pi/M) ; -a*sin(2*pi/M) ; b] ; 
Tp = [a ; 0;  b] ; 

ang_R1 = acos(sum(Ap.*Tp) / (norm(Ap)*norm(Tp))) ; 
u = cross(Ap,Tp) / norm(cross(Ap,Tp)) ; 

u_x = u(1) ; u_y = u(2) ; u_z = u(3) ; 

% ------------------------------------------------
% R1 rotates the vectors so that R1 * Am = Tm 
R1 = [cos(ang_R1)+u_x^2*(1-cos(ang_R1)) u_x*u_y*(1-cos(ang_R1))-u_z*sin(ang_R1) u_x*u_z*(1-cos(ang_R1))+u_y*sin(ang_R1); ...
u_y*u_x*(1-cos(ang_R1))+u_z*sin(ang_R1) cos(ang_R1) + u_y^2 * (1-cos(ang_R1)) u_y*u_z*(1-cos(ang_R1))-u_x*sin(ang_R1);...
u_z*u_x*(1-cos(ang_R1))-u_y*sin(ang_R1) u_z*u_y*(1-cos(ang_R1))+u_x*sin(ang_R1) cos(ang_R1)+u_z^2*(1-cos(ang_R1))]  ;

Am1 = R1 * Am ; 
Ap1 = R1 * Ap ; 

% ------------------------------------------------
% R2 a rotation about Tm in the plane orthogonal to Tm 
% For this we compute the angle of rotation by projection the vectors 

v = Tp ; % or v = Ap1
v_x = v(1) ; v_y = v(2) ; v_z = v(3) ; 

Am2 = cross(v,cross(Am1,v)) ; % Ap2 = Ap1 - v * sum(Ap1 .* v) ; 
Tm2 = cross(v,cross(Tm,v)) ;  % Tp2 = Tp - v * sum(Tp .* v) ; 

ang_R2 =  -acos(sum(Am2.*Tm2) / (norm(Am2)*norm(Tm2))); 

R2 = [cos(ang_R2)+v_x^2*(1-cos(ang_R2)) v_x*v_y*(1-cos(ang_R2))-v_z*sin(ang_R2) v_x*v_z*(1-cos(ang_R2))+v_y*sin(ang_R2); ...
v_y*v_x*(1-cos(ang_R2))+v_z*sin(ang_R2) cos(ang_R2) + v_y^2 * (1-cos(ang_R2)) v_y*v_z*(1-cos(ang_R2))-v_x*sin(ang_R2);...
v_z*v_x*(1-cos(ang_R2))-v_y*sin(ang_R2) v_z*v_y*(1-cos(ang_R2))+v_x*sin(ang_R2) cos(ang_R2)+v_z^2*(1-cos(ang_R2))]  ;

R = R2 * R1 ; 

RAm = R*Am ; 
RAp = R*Ap ; 





