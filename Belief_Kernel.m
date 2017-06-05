function [ M ] = Belief_Kernel( TK, OB, n, N )

% Variables for Uniform Quantization of belief
d = [-0.01 1/n:1/n:1]; %define boundaries for belief quantizing


inc=100000;     % set increment for Monte Carlo integration
M = zeros(n,n,N);

for pi0 = 0:1/inc:1     % iterate over initial belief
    b0 = [1-pi0 pi0];   % put initial belief in vector form
    for u0=1:N          % iterate over previous control
        for j=1:n
            if pi0 > d(j) && pi0 <= d(j+1)
                D = j;  % define Departure bin
            end
        end
        for y1=1:2 
            
            b = (b0 * TK(:,2,u0) * OB(2,y1)) /(b0 * TK(:,:,u0) * OB(:,y1));
            
            for A=1:n
                if b > d(A) && b <= d(A+1)      %find Arrival bin
                    M(D,A,u0) = M(D,A,u0) + (b0 * TK(:,:,u0) * OB(:,y1));
                end
            end
        end
    end
end

% Take the average to normalize the Monte Carlo integration
for i=2:n
        M(i,:,:) = M(i,:,:)./(inc/n);
end
M(1,:,:) = M(1,:,:)./(inc/n+1);


end

