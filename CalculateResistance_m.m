function [Z] = CalculateResistance(p_Q,b,n,c,B_A,rho,nu,m_old,L,d)
%UNTITLED Calculates the resistance of all branches
%   Detailed explanation goes here

B_S = eye(size(B_A));

%Calculating the resistance of all �ste
R_e=zeros(n-1,1);
z_a=zeros(n-1,1);
lambda = zeros((n-1),1);

for i=1:(n-1) 
    R_e(i) = (4/(pi*d(i)*nu(i)))*(m_old(i)/rho(i)); 
    
    if R_e(i) < 2320 
        lambda(i) = 64/R_e(i);
    elseif 2320<=R_e(i) < 10^5 
        lambda(i) = 0.3164*R_e(i).^(-0.25);
    else 
        lambda(i) = 0.0032 + (0.221*R_e(i).^(-0.237));
    end
    
    z_a(i) = lambda(i) *(8*rho(i)/pi.^2)*(L(i)/d(i).^5)*(m_old(i)/rho(i)); 
end

%Sehnenwiderst�nde Z_S
Z_A=zeros(n-1); % To calculate the resistance of the Sehnen a resistance-Matrix of �ste is needed
for i=1:(n-1)
        Z_A(i,i) = z_a(i);
end

Zaehler = B_A*Z_A*m_old(1:(n-1))./rho(1:(n-1))+p_Q(1:(n-1)); %fehlen Druckquellen und gesteuerte Druckabf�lle: -B_A*P_q(1:(n-1)) +B_A*P_st_old(1:(n-1))
Nenner = B_S*m_old((n):b);

Z_s = zeros(c,1); 
for i=1:c
   Z_s(i)= Zaehler(i)/Nenner(i);  
end

%Adding the resistance of �ste and Sehnen to get a matrix with the resistance of all branches
Z_as =[z_a; Z_s];
Z=zeros(b);

for i=1:b
    Z(i,i)=Z_as(i);
end
end

