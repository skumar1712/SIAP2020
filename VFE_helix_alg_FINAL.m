% -----------------------------------------
% In this program, we construct the solution for VFE at rational times
% ALGEBRAICALLY. In order to do so, we compute the angle theta_q 
% as a sum of the Gauss sums and theta_0/q
% The tangent vectors T(s,t) are computed from the Rotation matrix from angle rho
% and axis defined through theta_q 
% The curve X(s,t) is thus obtained from T(s,t)
% After computing X(s,t),T(s,t) we obtain the correct rotations so that the
% Algebraic solution matches the numerical solution except for the vertical
% height. 
% X(s,t) should be rotated such that the curve is orthogonal to z axis and
% the second rotation is still undecided
% Euclidean case
% Date: August 21, 2018

% The difference between the previous code is that now we consider N =
% M*q*r tangent vector values where 'r' corresponds to theta = pi/r so that
% the resutling curve X has a shift right from T. 
% -----------------------------------------

clear 

tic 

M = 6 ; 
p1 = 1 ; q1 = 8 ; 
N = 2^5 * M ; 
% SFT =
% for p1 = 0 : q1 
    
if p1 == 0 
    q = 1 ;
    p = 0 ; 
end

if ne(gcd(p1,q1),1) 
    p = p1 / gcd(p1,q1) ;
    q = q1 / gcd(p1,q1) ;
else 
    p = p1 ; 
    q = q1 ;
end

% ---------  file name 
Matfile = ['Alg_sol_VFE_M6_theta_pi5_q' num2str(q)] ; 

% --------------------------------------------------------
%     Initial values of angle rho and theta at time t = 0 
% --------------------------------------------------------   

    theta_0 = pi/5 ; 
%     if theta_0 >= 2*pi/M
%         Error('b>=1') 
%     end
    b = tan(theta_0/2) / tan(pi/M) ; 
    
%     b = 0.4 ; 
%     theta_0 = 2*atan(b*tan(pi/M)) ; 
    
    a = sqrt(1-b^2) ; 
    rho_0 = 2 * asin(a*sin(pi/M)) ; 
    
    c0 = sqrt(-2*log(cos(rho_0/2))/pi) ; 
    
% Computation of "theta_m" using Gauss sums and angle "rho" 

    for m = 0 : q-1
        if mod(q,2) == 1
            t(m+1) = -1i*log(Gauss_sum(-p,m,q) / sqrt(q)) ; 
            rho(m+1) = 2*acos(cos(rho_0/2)^(1/q)) ;
        elseif mod(q/2-m,2) == 0 
            t(m+1) = -1i*log(Gauss_sum(-p,m,q) / sqrt(2*q)) ; 
            rho(m+1) = 2*acos(cos(rho_0/2)^(2/q)) ;
        else 
            t(m+1) = 0 ;
            rho(m+1) = 0;
        end
    end

    for k = 0 : M-1
        for m = 0 : q-1
            th(q*k+m+1) = theta_0 *(k+ m/q) + t(m+1) ;
            rho_q(q*k+m+1) = rho(m+1) ; 
        end
    end
    

th = real(th) ;
th = th + p * theta_0^2 / (2*pi*q) ; 

t = real(t) ; 
j = 1 ; 


% Computation of Rotation matrix 
% % Adding the addtional phase term 
% th = th+(theta_0^2/(2*pi))*(p/q); 

for j = 1 : M*q 
    R(:,:,j) = [cos(rho_q(j))     sin(rho_q(j))*cos(th(j))     sin(rho_q(j))*sin(th(j)) ; 
        -sin(rho_q(j))*cos(th(j))    cos(rho_q(j))*cos(th(j))^2+sin(th(j))^2      (cos(rho_q(j))-1)*cos(th(j))*sin(th(j)) ; 
        -sin(rho_q(j))*sin(th(j))   (cos(rho_q(j))-1)*cos(th(j))*sin(th(j))   cos(rho_q(j))*sin(th(j))^2+cos(th(j))^2] ; 
end

% Product of M*q matrices 
% Prd_R = eye(3) ; 
% for j = 1 : M*q 
%     Prd_R = R(:,:,j) * Prd_R ; 
% end 

% Computation of Frenet frame {T,e1,e2} and tangent vectors 
Te1e2 = eye(3) ; 
T(1,:) = Te1e2(1,:) ; 

for k = 1 : M*q 
    Te1e2 = R(:,:,k) * Te1e2 ;
    T(k+1,:) = Te1e2(1,:) ; 
end 

T1 = T ; % Tangent vectors for plotting a closed curve 
T = T(1:end-1,:) ; 
% ******************************************************************
% UNLIKE WHEN b=0, FIRST TANGENT VECTOR IS AT T(0-) INSTEAD OF T(0+)
% ******************************************************************

% ----------------------------------------------
% Next, we incorporate SHIFT in our code, which we do by making r copies of
% tangent vectors 
% ----------------------------------------------
% X(1,:) = [ 0 0 0 ] ; 
% 
% if r>0 && p >0  % when b>0 , p>0 
%     for j = 1 : length(T) 
%         TT((j-1)*r + 1 : j*r ,:) = ones(r,1)*T(j,:) ; 
%     end
% 
% % Since SHIFT depends on p and q, therefore
%     
%     p = mod(p,r) ; 
%     if p == 0 , p = r ; end 
%     TT1 = TT(r-p+1:end,:) ; 
%     TT1(end+1:end+r-p,:) = TT(1:r-p,:) ; 
% % Integrate X_s = T 
%     
%     h = 2* pi /(M*q*r) ; 
% 
%     for j = 1 : length(TT1) 
%         X(j+1,:) = X(j,:) + h * TT1(j,:) ;  
%     end
%     
% else
% %     when b =0 or p = 0 ;
%     h = 2* pi /(M*q) ; 
% 
%     for j = 1 : length(T) 
%         X(j+1,:) = X(j,:) + h * T(j,:) ;  
%     end 
% % last vertex is needed to perform the first rotation (see below)! 
% end

% ---------------------------------------------
% Constructing X from the algorithm written in the draft
% It works for all kind of theta_0 or b values
% ---------------------------------------------
X(1,:) = [ 0 0 0 ] ; 

% Compute spq
spq = (2*theta_0 * p / (M*q)) ; 

rem_spq = 2*pi/(M*q) - spq ;
X(2,:) = X(1,:) + spq * T(1,:) ; 

for j = 2 : (M*q)
    X(j+1,:) = X(j,:) + (2*pi/(M*q)) * T(j,:) ;
end

j=M*q+1;
X(j+1,:) = X(j,:) + rem_spq * T(1,:) ;

if spq==0
    X = X(2:end,:);
end
% X= X(1:end-1,:) ; 
% ------------------------------------------------
% Rotation so that X satisfies X(2*pi)-X(0) = 2*pi*b 
% ------------------------------------------------

v1 = X(end,:) - X(1,:) ; 

% Computation of angle with z axis and rotation axis

ang_R1 = acos(v1(3)/norm(v1)) ; % angle
u = cross(v1,[0 0 1]) / norm(cross(v1,[0 0 1])) ; % axis
u_x = u(1) ; u_y = u(2) ; u_z = u(3) ; 

R1 = [cos(ang_R1)+u_x^2*(1-cos(ang_R1)) u_x*u_y*(1-cos(ang_R1))-u_z*sin(ang_R1) u_x*u_z*(1-cos(ang_R1))+u_y*sin(ang_R1); ...
u_y*u_x*(1-cos(ang_R1))+u_z*sin(ang_R1) cos(ang_R1) + u_y^2 * (1-cos(ang_R1)) u_y*u_z*(1-cos(ang_R1))-u_x*sin(ang_R1);...
u_z*u_x*(1-cos(ang_R1))-u_y*sin(ang_R1) u_z*u_y*(1-cos(ang_R1))+u_x*sin(ang_R1) cos(ang_R1)+u_z^2*(1-cos(ang_R1))]  ;

RX = (R1 * X.').' ; 
RT = (R1 * T.').' ; 

RT1= (R1 * T1.').' ; 

% --------------------------------------------------------------------
% Since the first two components of X(s,0) have a symmetry, so we align 
% X_alg accordingly.
% The symmetry is that X(0)-X(2*pi/M) is a multiple of (1,0,0).
% And, to incorporate rotational shift we align the solution at t=2*pi/M^2 
% also in the same way.
% --------------------------------------------------------------------
%  Second rotation 

if p==0 || (p==1 && q==1)
    
    Rx1 = [RX(p+1,1:2) 0] ;
    Rx2 = [RX(p+1+q,1:2) 0] ; 
    Rx = Rx2 - Rx1 ; 
    ang2 = acos( sum(Rx .* [1 0 0 ] / norm(Rx)) ) ; 
    v = cross(Rx,[1 0 0]) / (norm(cross(Rx,[1 0 0]))) ; 
if v(3) < 0 
    ang2 = -ang2 ; 
end 
R2 = [cos(ang2) -sin(ang2) 0; sin(ang2) cos(ang2) 0; 0 0 1 ] ; 
else
    R2 = eye(3);
end

% The third rotation corresponding to the Rotational shift 
if p == 1 && q == 1
%     ang3 = 0.204576960792854 ;  % M=5, theta = pi/3 ; 
    ang3 =  0.106569551402475 ; % b=0.4 
    R3 = [cos(ang3) -sin(ang3) 0; sin(ang3) cos(ang3) 0; 0 0 1 ] ; 
else 
    R3 = eye(3) ;
end

R2X = (R3*R2 * RX.').' ; 
R2T = (R3*R2 * RT.').' ; 

% Mean of vertices ONLY

R2X(:,1:2) = R2X(:,1:2) - mean(R2X(2:end-1,1:2));
c_M_alg = -2*log(cos(rho_0/2))/(pi*tan(pi/M)/M) ; 

R2X(:,3) = R2X(:,3) + c_M_alg*2*pi/M^2 ; 

% save(Matfile, 'T','X','M' ,'p', 'q','theta_0','rho_0','c0','RX','RT','-v7.3' ) ;

% ssqM = 2 * pi * (0:q*M) / (q * M);
% index1 = floor((q * M) * (0:N-1) / N) + 1;
% index2 = floor((q * M) * (0:N-1) / N) + 2;
% ss = 2 * pi * (0:N-1) / N;
% % ss1 = ssqM(index1);
% % ss2 = ssqM(index2);
% 
% i1 = ismember(index1,1) ;    
% n_side1 = length(i1(ne(i1,0))) ; % computing the length of first side
% % i1 = ismember(index1,1) ;    
% % n_side1 = length(i1(ne(i1,0))) ; % computing the length of last side
% n1_sft = round(2*theta_0*p*n_side1/(2*pi)) ;
% n2_sft = n_side1 - n1_sft ; 
% 
% index11 = index1(n2_sft+1:end) ;
% index22 = index11+1; 
% ss1 = ssqM(index1);
% ss2 = ssqM(index2);       
% 
% XXX = R2X(2:end,:) ;
% 
% for m = n_side1+1 : N
%     XXXX(m, :) = ((ss(m-n1_sft+1) - ss2(m)) * XXX(index1(m), :) ...
%         + (ss1(m) - ss(m-n1_sft+1)) * XXX(index2(m), :)) / (ss1(m) - ss2(m));
% %     TTTT(m, :) = TTT(index11(m), :);
% end


% for m = 1 : length(index11)
%     XXXX(m, :) = ((ss(m) - ss2(m)) * XXX(index11(m), :) ...
%         + (ss1(m) - ss(m)) * XXX(index22(m), :)) / (ss1(m) - ss2(m));
% %     TTTT(m, :) = TTT(index11(m), :);
% end

toc 

