function [Ccal,A] = regions(Hcal,Sigma)
    
    % This function turns the Sigma info into adjacency matrix A and constraint matrix Ccal.

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    % BLOCK 1: INITIALIZE VECTORS/MATRICES
    
    nc = size(Hcal,1);
    P = size(Sigma,2);
    A = zeros(P);
    Ccal = [];
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    % BLOCK 2: GENERATE A MATRIX
    
    for p=1:P
       for q=p+1:P
          delta = Sigma(:,p)-Sigma(:,q);
          if norm(delta,1)==1
              A(p,q) = find(delta);
          end
       end
    end
    A = A+A';
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    % BLOCK 3: GENERATE Ccal MATRIX
        
    nidx = zeros(P,1);
    for p=1:P
        idx = A(A(:,p)>0,p);
        nidx(p) = length(idx);
        Ccalp = [Hcal(idx,:),-ones(size(idx))];
        Ccalp(Sigma(idx,p)>0,:) = -Ccalp(Sigma(idx,p)>0,:);
        Ccalp(Sigma(idx,p)>0,end) = Ccalp(Sigma(idx,p)>0,end)-eps;
        Ccal = [Ccal;[Ccalp,p.*ones(size(idx))]];
    end
    
end

