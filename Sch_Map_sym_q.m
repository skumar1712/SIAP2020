% In this program we solve VFE and Schrodinger Map equation for a regular
% helical polygon in the Euclidean plane. We make use of the Spatial symmetry of T(s,t) and 
% reduce the compuation cost by reducing the number of elements to N/M
% Date: May 18, 2018 
% Sch_map_sym_q.m : the file deals with N mutliple of q 
% Last modified: August 30, 2018 (for adding the 'shift') 

close all;
clear  

tic
% ------------------------------------
% Parameters
% ------------------------------------
    M = 3 ; 
    n = 0 ; 
 
    tpr = 1/M^2; 
    mlt = tpr * M^2 + 0; 
    
    q1 = 2^2*3^2*5 ; % q per period  
    q = (2^2*3^2*5) * mlt ; p = 1 ;  
    q = 307 * mlt ; p = 1 ; 
    N = 2*M * q / mlt; % N should be even 
    l = 1 ; 
	N_t = 2^4*3^2*5*q ; % q = 307
%     N_t = 2^5*3*5*q * 5 ; 
    theta = pi / 2 ; 
    b =  tan(theta/2) / tan(pi/M) ;  
%     b = 0 ; 
    choice  = 1 ;           % choice = 1 : Euclidean, choice = 2 : Hyperbolic 
    L =   2*pi ; 
    s =  (0 : 1/(N) : 1-1/(N) ) * L ; s = s.' ; 
    h = s(2) - s(1) ;       % spatial step size
    
    if choice == 1 
        T_run =   mlt * 2 * pi / M^2  ; 
    else 
        T_run = 2* pi / M^2 ; 
    end

    dt = T_run / N_t ;          % temporal step size 
    
    SFT = 2 * theta / M ;       % shift at the end of one period
    SFT =  N* theta  / (M*pi) ;  % in terms of indices 
%     SFT_pq = SFT * p / q ; 
    
% ------------------------------------------------
% Stability condition   

    if dt/h^2 > 0.29
        error('Step size is too large') 
    end
    
%     Wavenumbers:
    
    k1 = 2i * pi * [0 : (N)/2-1 -(N)/2 : -1 ] / L ; k1 = k1.' ; 
    
%   wavenumbers for N/M elements 

    kk1 = 2* pi * [0 : (N)/(2*M)-1 -(N)/(2*M) : -1 ] / L ; kk1 = kk1.' ;
    k2 = M*kk1+1 ; 
    k3 = M* kk1 ; 
% ------------------------------------------------     
% Initial data 
% ------------------------------------------------
    z0 = zeros(N,1) ; 
    for j = 0 : M-1
        z0( j* N/M + (1:N/M) ) = exp(2i*pi* j / M) ; 
    end 
    j = 0 : 1 ; 
    x_0 = - (1i * pi * exp(1i * pi * (2 * j - 1 ) / M )) ./ (M * sin(pi / M) )   ;  % Initial data
     
    X00 = linspace(x_0(1), x_0(2), N/M + 1);
    Z0 = zeros(N,1);

    Z0(1:N/M) = X00(1:N/M);
    for m = 1: M-1
        Z0(m * N/M + (1:N/M)) = exp(2i * pi * m / M) * Z0(1:N/M); 
    end

    if choice == 1
        T30 = b * ones(N,1) ; T10 = sqrt(1-b^2) * real(z0) ; T20 = sqrt(1-b^2) * imag(z0) ; 
        X30 = b * s ; X10 = sqrt(1-b^2) * real(Z0); X20 = sqrt(1-b^2) * imag(Z0); 
    else
        T10 = b * ones(N,1) ; T20 = sqrt(b^2-1) * real(z0) ; T30 = sqrt(b^2-1) * imag(z0) ; 
        X10 = b * s ; X20 = sqrt(b^2-1) * real(Z0); X30 = sqrt(b^2-1) * imag(Z0); 
    end
    
% Storing the initial data
    Xfull = zeros(l*N,q+1) ; XX3full = zeros(l*N,q+1) ; 
    Tfull = zeros(l*N,q+1) ; TT3full = zeros(l*N,q+1) ; 
    
    T1 = T10; T2 = T20; T3 = T30 ; 
    X1 = X10 ; X2 = X20 ; X3 = X30 ; 
    
    T = T1 + 1i*T2 ; X = X1 + 1i*X2 ; 
    
    T = T(1:N/M) ; X = X(1:N/M) ; 
    T3 = T3(1:N/M) ; X3 = X3(1:N/M) ; 
    
    XX = [X ; exp(2i*pi/M)*X; exp(4i*pi/M)*X] ;
    XX3 = [X3 ; (2*pi*b/M)+X3 ; (4*pi*b/M)+X3] ;
    TT = [T ; exp(2i*pi/M)*T; exp(4i*pi/M)*T] ;
    TT3 = [T3 ; T3 ; T3] ;
    
    Xfull(:,1) = XX ; XX3full(:,1) = XX3 ; 
    Tfull(:,1) = TT; TT3full(:,1) = TT3 ; 
    z1(1,:) = [X1(1) X2(1) X3(1)] ;
    Xv(1,:) = [real(X(1)) imag(X(1)) X3(1)] ; 
    
% X_mean 

    X_mean = zeros(N_t,1) ; X3_mean = zeros(N_t,1) ; 
    X_mean(1) = sum(X)/(N/M) ;
    X3_mean(1) = sum(X3)/(N/M); 

   
% Rotation R 
    j = 0 : N/M-1 ;
    R = exp(-2i*pi*j/N) ; R = R.' ; 

 
% ------------------------------------------------
% Evolution of X and T
% ------------------------------------------------
n = 0 ; % the number of current period 
for r = 1 : N_t
    
    if isnan(T1(1))
        r
        r * dt
        error('Error!');
    end
    
    Tr = T .* R ; Tr_ =  fft(Tr) ; Tr_(abs(Tr_)/N < eps) = 0 ; Ts = ifft( Tr_ * 1i .* k2 ) .* conj(R) ; Tss = ifft( Tr_ .* (1i * k2 ).^2 ) .* conj(R) ; 
    T3_ =  fft(T3) ; T3_(abs(T3_) < eps) = 0 ; T3s = real(ifft( T3_ * 1i .* k3 )) ; T3ss = real(ifft( T3_ .* (1i * k3).^2 )) ; 
    AT = 1i * (-T3ss .* T  + T3 .* Tss) ; 
    AT3 = real(T) .* imag(Tss) - imag(T) .* real(Tss) ; 
    AX = 1i * (-T3s .* T  + T3 .* Ts) ; 
    AX3 = real(T) .* imag(Ts) - imag(T) .* real(Ts) ; 
    Taux = T + 0.5 * dt * AT ; 
    T3aux = T3 + 0.5 * dt * AT3 ; 
    
    Tr = Taux .* R ; Tr_ = fft(Tr) ; Tr_(abs(Tr_)/N < eps) = 0 ; Ts = ifft( Tr_ * 1i .* k2 ) .* conj(R) ; Tss = ifft( Tr_ .* (1i * k2 ).^2 ) .* conj(R) ; 
    T3_ =  fft(T3aux) ; T3_(abs(T3_) < eps) = 0 ; T3s = real(ifft( T3_ * 1i .* k3 )) ; T3ss = real(ifft( T3_ .* (1i * k3).^2 )) ; 
    BT = 1i * (-T3ss .* Taux  + T3aux .* Tss) ; 
    BT3 = real(Taux) .* imag(Tss) - imag(Taux) .* real(Tss) ; 
    BX = 1i * (-T3s .* Taux  + T3aux .* Ts) ; 
    BX3 = real(Taux) .* imag(Ts) - imag(Taux) .* real(Ts) ; 
    Taux = T + 0.5 * dt * BT ; 
    T3aux = T3 + 0.5 * dt * BT3 ; 
    
    Tr = Taux .* R ; Tr_ = fft(Tr) ; Tr_(abs(Tr_)/N < eps) = 0 ; Ts = ifft( Tr_ * 1i .* k2 ) .* conj(R) ; Tss = ifft( Tr_ .* (1i * k2 ).^2 ) .* conj(R) ; 
    T3_ = fft(T3aux) ; T3_(abs(T3_) < eps) = 0 ; T3s = real(ifft( T3_ * 1i .* k3 )) ; T3ss = real(ifft( T3_ .* (1i * k3).^2 )) ; 
    CT = 1i * (-T3ss .* Taux  + T3aux .* Tss) ; 
    CT3 = real(Taux) .* imag(Tss) - imag(Taux) .* real(Tss) ; 
    CX = 1i * (-T3s .* Taux  + T3aux .* Ts) ; 
    CX3 = real(Taux) .* imag(Ts) - imag(Taux) .* real(Ts) ; 
    Taux = T +  dt * CT ; 
    T3aux = T3 +  dt * CT3 ; 
    
    Tr = Taux .* R ; Tr_ = fft(Tr) ; Tr_(abs(Tr_)/N < eps) = 0 ; Ts = ifft( Tr_ * 1i .* k2 ) .* conj(R) ; Tss = ifft( Tr_ .* (1i * k2).^2 ) .* conj(R) ; 
    T3_ = fft(T3aux) ; T3_(abs(T3_) < eps) = 0 ; T3s = real(ifft( T3_ * 1i .* k3 )) ; T3ss = real(ifft( T3_ .* (1i * k3).^2 )) ; 
    DT = 1i * (-T3ss .* Taux  + T3aux .* Tss) ; 
    DT3 = real(Taux) .* imag(Tss) - imag(Taux) .* real(Tss) ; 
    DX = 1i * (-T3s .* Taux  + T3aux .* Ts) ; 
    DX3 = real(Taux) .* imag(Ts) - imag(Taux) .* real(Ts) ; 
    
    T = T + dt * (AT + 2 * BT + 2 * CT + DT) / 6 ; 
    T3 = T3 + dt * (AT3 + 2 * BT3 + 2 * CT3 + DT3) / 6 ; 
    X = X + dt * (AX + 2 * BX + 2 * CX + DX) / 6 ; 
    X3 = X3 + dt * (AX3 + 2 * BX3 + 2 * CX3 + DX3) / 6 ; 

    if choice == 1 
        Tnorm = sqrt(real(T).^2+imag(T).^2 + T3.^2) ; 
    else 
        Tnorm = sqrt(T1.^2-T2.^2-T3.^2) ; 
    end 
    T = T ./ Tnorm ; T3 = T3 ./ Tnorm ; 
   
    if mod(r,N_t*p/q) == 0 
        
        if mod(p,q1) == 0  % to keep the count of the period 
            n = n+1 ; 
        end
        
        Xfull(:,p+1) = XX ;  XX3full(:,p+1) = XX3 ; 
        Tfull(:,p+1) = TT ;  TT3full(:,p+1) = TT3 ; 
%         dist_pq = N / (M*q1) ; 
%         p1 = mod(p,q1) ; 
%         SFT_pq = n*SFT + SFT * p1 / q1 ;  % T_SFT is total shift before the current period 
%         SFT_pq = round(mod(SFT_pq,dist_pq) )  
%         
%         Xv(p+1,:) = [real(X(SFT_pq+1)) imag(X(SFT_pq+1)) X3(SFT_pq+1)] ; 
        p = p + 1 ; 
        
     
    end 
   
   z1(r+1,:) = [real(X(1)) imag(X(1)) X3(1)] ;
   X3_mean(r+1) = mean(X3);    
   
end

% T3full(:,:) = zeros(l*N,q+1) ; X3full = zeros(l*N,q+1); 
% TTfull = zeros(l*N,q+1) ; XXfull = zeros(l*N,q+1); 
% 
% for k = 1 : q+1
%     
%     for m = 0 : l*M-1
%         TTfull(m * N/M + (1:N/M),k) = exp(2i * pi * m / M) * Tfull(:,k); 
%         T3full(m * N/M + (1:N/M),k) =  TT3full(:,k); 
%    
%         XXfull(m * N/M + (1:N/M),k) = exp(2i * pi * m / M) * Xfull(:,k); 
%         X3full(m * N/M + (1:N/M),k) = 2*pi*b*m/M + XX3full(:,k); 
%     end
%     
% end
% 
% T1full = real(TTfull) ; T2full = imag(TTfull) ;
% X1full = real(XXfull) ; X2full = imag(XXfull) ;

% Calculation for h and c_M
% h_T = M * sum(X3) ; 
c_M = (X3_mean(end)-X3_mean(1)) / (T_run) ; 


% if choice == 1 
%     save Sch_Map_fft_sym_M3_q101_b0_3prd.mat X1full X2full X3full T1full T2full T3full X_mean X3_mean T_run z1 L M dt b c_M  
% else 
%     save Sch_Map_fft_M5_2N8_dt_hyp_eps_1.9.mat X1full X2full X3full T1full T2full T3full X1_mean X2_mean X3_mean z L  M dt b 
% end 

toc



