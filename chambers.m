function Sigma = chambers(Hcal,Dcal)
    
    % This function calculated the chambers resulting from 
    % a cutting arrangement Hcal over domain Dcal.

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    % BLOCK 1: INITIALIZE VECTORS/MATRICES
    
    nc = size(Hcal,1);
    Sigma = [];

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    
    % BLOCK 2: CREATE FEASIBILITY PROBLEM STRUCT
    
    feasprob.f = [];    
    feasprob.lb = Dcal(:,1)+1e-3;
    feasprob.ub = Dcal(:,2)-1e-3;
    feasprob.solver = 'linprog';
    feasprob.options = optimoptions('linprog','display','off');
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    % BLOCK 3: CEHCK FEASIBILITY OF ALL COMBINATIONS
        
    for j=1:2.^nc
       sigma = bin2dec(dec2bin(j-1,nc)');
       feasprob.Aineq = Hcal;        
       feasprob.bineq = ones(nc,1);
       feasprob.bineq(sigma>0) = -feasprob.bineq(sigma>0);
       feasprob.Aineq(sigma>0,:) = -feasprob.Aineq(sigma>0,:);
       feasprob.bineq = feasprob.bineq-1e-3;
       [~,~,flag,~]=linprog(feasprob);
       if flag>0
           Sigma = [Sigma,sigma];
       end
    end
   
end

