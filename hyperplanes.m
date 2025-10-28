function Hcal = hyperplanes(Phiv,Dcal)
    
    % This function builds the hyperplane normals used to partition the input space.

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    % BLOCK 1: INITIALIZE VECTORS/MATRICES

    d = size(Dcal,1);
    nc = round(length(Phiv)./d./(d-1));
    Phi = reshape(Phiv,[d-1,d,nc]);
    rho = norm(max(abs(Dcal),[],2),2);%min(abs(Dcal),[],'all');
    Hcal = zeros(nc,d);
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    % BLOCK 2: FIND HYPERPLANES
        
    for c=1:nc    
        X = zeros(d);
        phi = squeeze(Phi(:,:,c));
        for i=1:d
            for j=1:d-1
                X(i,j) = rho.*cos(phi(j,i)).*prod(sin(phi(1:j-1,i)));            
            end
            X(i,d) = rho.*prod(sin(phi(1:d-1,i)));                        
        end
        Hcal(c,:) = (X\ones(d,1))';    
    end
   
end

