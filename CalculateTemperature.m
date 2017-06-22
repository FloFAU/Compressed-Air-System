function [ T_node ] = CalculateTemperature( T1,  number_of_nodes, knotenindex, p_node)
% CalculateTemperature Calculates the temperature at every node (in K)
% assumption: reversibel adiabat

T_node = zeros(number_of_nodes,1);

T_node (1) = T1; %in K

cp = 1.005;  %kJ/(kg*K) Luft
cv = 0.718;  %kJ/(kg*K) Luft

n=cp/cv; %assumption: reversibel adiabat

    for i = 2:number_of_nodes
        for k = 1:number_of_nodes 
            if knotenindex((i-1), k) == -1
   
                T_node(i) = abs(T_node(k)*(p_node(i)/p_node(k))^((n-1)/n)); %Vorsicht bei der Nummerierung der Ströme und Knoten -> Fehler wenn Strom 3 aus Knoten 4 fließen würde
            
            end
        end  
    end
end

