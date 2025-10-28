function J = regapprox(Phi,Dcal,G)

    % This function builds a region partition, 
    % then uses fmincon as an example solver to choose all affine modes (Jcal, Kcal) 
    % to best fit nldyn on samples G, with continuity constraints.
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    % BLOCK 1: INITIALIZE / EXTRACT CUTTING HYPERPLANES
    
    options = optimoptions('fmincon',...
        'display','off','Algorithm','interior-point');
        %'EnableFeasibilityMode',true,'SubproblemAlgorithm','cg');
    d = size(G,1);

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    % BLOCK 2: ANALYZE CUTTING ARRANGEMENT
    
    Hcal = hyperplanes(Phi,Dcal);
    
    try
        Sigma = chambers(Hcal,Dcal);
        [Ccal,A] = regions(Hcal,Sigma);
        P = length(A);
        
        % SOLVE LOW-LEVEL PROBLEM

        % This is an example of solving with MATLAB fmincon
        % and inforcing the continuity constraint
        
        x0 = rand((d+1).*P,1);
        [~,J] = fmincon(@(pw) approxerr(pw,Ccal,G),x0,[],[],...
            [],[],[],[],@(pw) continuity(pw,Hcal,A),options);
    catch ME
        J = Inf;
    end
    
end