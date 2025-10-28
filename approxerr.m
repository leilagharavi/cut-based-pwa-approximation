function errn = approxerr(dv,Ccal,G)

    % This function computes the approximation error between 
    % a PWA model and a nonlinear dynamic model, measured over a set of points G.

    % dv: a vector containing model parameters
    % Ccal: calligrafic C matrix (partitioning)
    % G: matrix including the sample points (grid)
    % errn: the norm (magnitude) of the error vector

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    % BLOCK 1: INITIALIZE / EXTRACT AFFINE MODES
    
    err = zeros(size(G,2),1); % vector holding errors for each test point                  
    P = Ccal(end,end); % number of affine partitions                         
    Jcal = reshape(dv(1:end-P),P,[]); % linear coefficients of local modes
    Kcal = dv(end-P+1:end); % offset elements of local modes
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    % BLOCK 2: CALCULATE ERRORS LOOP
    
    for i=1:size(G,2)
        f = pwaapprox(Jcal,Kcal,Ccal,G(:,i));
        F = nldyn(G(:,i));
        err(i) = (f-F);
    end

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    % BLOCK 3: CALCULATE AGGREGATE ERROR

    errn = norm(err,2);
    
end