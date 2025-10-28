function f = pwaapprox(Jcal,Kcal,Ccal,x)
    
    % This function evaluates a PWA model at an input point x.

    % Jcal: calligrafic J matrix (linear coefficients of local modes)
    % Kcal: calligrafic K matrix (offset elements of local modes)
    % Ccal: calligrafic C matrix (partitioning)
    % x: input point

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    % BLOCK 1: EXTRACTING NUMBER OF REGIONS/PIECES
    
    % each region p defines a polytopic region 

    P = size(Jcal,1);
    cond = zeros(P,1);
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    % BLOCK 2: FIND THE REGION > ASSIGN CORRESPONDING MODE
    
    for p=1:P
        cond(p) = all(Ccal(Ccal(:,end)==p,1:end-1)*[x;1]<0);        
    end

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    % BLOCK 3: EVALUATE CORRESPONDING AFFINE MODE 

    f = Jcal(find(cond,1),:)*x+Kcal(find(cond,1)); 
    
end