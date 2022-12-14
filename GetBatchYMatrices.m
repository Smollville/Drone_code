function [Fb,Gb,Qb,Rb,F,G,H] = GetBatchYMatrices(A,B,C,N,P,Q,R)
% Returns matrices for batch computation of state sequence and cost 
% functional

    F = [];
    G = [];
    H = [];
    Qb = [];
    Rb = [];
    for n = 0:N
        % state matrices
        F = [F ; A^n];
        Gi = [];
        for m = 0:N-1
            ni = n-m-1;
            Gaux = A^ni*B;
            if ni < 0
                Gaux = Gaux*0;
            end
            Gi = [Gi , Gaux];
        end
        G = [G ; Gi];
        
        H = blkdiag(H,C);

        % cost matrices
        if exist('P','var') && exist('Q','var') && exist('R','var')
            if n < N
                Qb = blkdiag(Qb,Q);
                Rb = blkdiag(Rb,R);
            else
                Qb = blkdiag(Qb,P);
            end
        end
    end

    % compute output matrices
    Fb = H*F;
    Gb = H*G;

end

