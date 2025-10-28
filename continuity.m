function [C,Ceq] = continuity(dv,Hcal,A)
    
    % This is a constraint function for continuous PWA approximation

    % dv: stacked PWA parameters
    % Hcal: calligraphic H matrix (hyperplane arrangement)
    % A: region adjacency matrix

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    % BLOCK 1: INITIALIZE / EXTRACT AFFINE MODES
    
    P    = length(A);
    Jcal = reshape(dv(1:end-P), P, []);   % P x d
    Kcal = dv(end-P+1:end);               % P x 1
    d    = size(Jcal, 2);
    
    C = [];
    Ceq = []; % we do not have equality constraints
    tol = 1e-3; % tolerance for numerical stability
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    % BLOCK 2: CHECK EDGES (COLLECTING CONSTRAINTS PER SHARED EDGE)

    for i=1:P
        for j=i+1:P
            edge = A(i,j);
            if edge > 0
               
                diff_vec = [Jcal(i,:)-Jcal(j,:), Kcal(i)-Kcal(j)];  % 1 x (d+1)
                hb = Hcal(edge, :);  % 1 x (d+1)
                nd = norm(diff_vec);
                nh = norm(hb);

                if nd == 0 || nh == 0
                    ang = 0;  % if either is zero, treat as aligned (no penalty)
                else
                    % angle via dot only, clamped for numerical safety
                    c = dot(diff_vec, hb) / (nd * nh);
                    c = max(-1, min(1, c));
                    ang = acos(c);
                end
                
                % inequality: ang - tol <= 0
                C(end+1,1) = ang - tol; 

            end
        end       
    end
end

