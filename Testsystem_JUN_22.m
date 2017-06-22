%% (SI-)UNITS 

% pressure p in Pa
% Temperature T in Kelvin
% length L in m
% diameter d in m
% density rho in kg/m^3
% kinematic viscosity nu in ??
% volume-vlow v in ??
% mass-flow ? in ??

%% This values depend on the system, must be changed every time

b = 6;      % branches
n = 4;      % nodes
c = 3;      % components (Maschen)

knotenindex = [-1 1 0 0; 0 -1 1 0; 0 -1 0 1 ] ; %column: knot, row: volumeflow
B_A = [-1,0,-1;0,1,-1;1,1,0]; % Äste

T1=275; % temperature at node 1
p1=1*10^5; % pressure at node 1
L = [10;10;10]; % length of tubes 
d = [2;2;2]; % diameter of tubes 


%% Initial values

p_mittel = [300000;290000;290000]; 
T_mittel = [222;222;222]; 
rho_x = CalculateDensity(p_mittel,T_mittel);
x_old = [12;1;1;133;1;1;1;1;20;2;2;2;1;1;1;1;1;1]; % Initial vector, we need to find proper starting values

%%

B = zeros(b,c);
B_S = eye(size(B_A)); % Sehnen
B = [B_A,B_S]; % Maschenmatrix
v_old = x_old((b+1):(2*b));


%%


while true

    %% Calculating the resistance
nu = CalculateViscosity( p_mittel,T_mittel,rho_x);
Z = CalculateResistance(b,n,c,B_A,rho_x,nu,v_old,L,d);

    %% Creating system of equations
Z_M = B*Z*B';
E_Z = eye(b);
E_K = eye(n-1); % size(E_K)
V_U = zeros(b);
Z_V = zeros(b);
null_one = zeros(b);
null_two = zeros(c);
null_three = zeros(c,b);
middle_A = [E_K,-B_A';null_two,Z_M];
side_A = [null_three;B];
A = [E_Z,-Z,-E_Z;null_one,middle_A,side_A;-V_U,-Z_V,E_Z];

    %% Creating Solutions-vector
p = zeros(b,1);
p(1) = 400000;
v_rhs_middle = zeros(n-1,1) ;
rhs_middle = [v_rhs_middle;B*p] ;
rhs_end = zeros(b,1);
r_k = [-p;rhs_middle;rhs_end];

    %% Calculation of the new values 
x = linsolve(A,r_k);
x_new = x; 


    %% Calculating the pressure at all nodes
kn= [x(1:(n-1))];

p_k = CalculateNodepressure( kn, p1, n, knotenindex );




    %% calculating the average pressure in all tubes:
p_mittel= zeros(n-1,1);
    for i = 1:(n-1)
    p_mittel(i)= (p_k(i)+p_k(i+1))/2; %(wenn wir Druckverluste durch Bauteile in p_Q schreiben, dann: (p_k(1)+p_k(i+1)+p_Q(1)/2 )
    end
  for i = 1:(n-1)
        for k = 1:n 
            if knotenindex((i), k) == -1
           
            p_mittel(i) = (p_k(k)+p_k(i+1))/2;  %%hier noch falsch, wenn Quelle vorhanden
            
            end
        end  
  end
  
    %% Calculating the temperature at all nodes

T_node = CalculateTemperature( T1,  n, knotenindex, p_k);

    %% Calculating the average Temperature in all tubes
T_mittel= zeros(n-1,1);
     for i = 1:(n-1)
        for k = 1:n 
            if knotenindex((i), k) == -1
           
            T_mittel(i) = (T_node(k)+T_node(i+1))/2; % wenn Bauteil vorhanden falsch; Formel nur für reine Rohre gültig
            
            end
        end  
     end

%% Calculating the average density in all tubes
    
rho_x = CalculateDensity(p_mittel,T_mittel);
    
%% Check if values are within the tolerance
if  abs(x_old - x_new) <= 0.000001 %perhaps not only compare v, -> whole solution vector -> done.
   break;
end
    x_old = x_new;
    v_old = x_new((b+1):(2*b));
end

%% Erggebnisse ausgeben:
x_new