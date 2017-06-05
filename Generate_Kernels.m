function [ K, Q ] = Generate_Kernels( n, measure_variance, iterative_variance, percent_mass )
%%gernerates transition kernels OB and TK

x = [-1/n:-1/n:(1-n)/n];
p = 2*normcdf(x,0,iterative_variance);

Q = zeros(n,n);
K = zeros(n,n,n);
K(1,1,1) = 1;
for i=2:n
    for j=2:n
        if j==i
            K(i,j,1) = 1 - p(1);
        end
        if j<i
            K(i,j,1) = p(i-j)-p(i-j+1);
        end
    end
    K(i,1,1) = p(i-1);
end


p = normcdf(x,0,measure_variance);
Q(1,1) = 1 - sum(p)/2;
Q(n,n) = 1 - sum(p)/2;
Q(1,n) = p(n-1)/2;
Q(n,1) = Q(1,n);

for i=2:n-1
    for j=2:n-1
        if j==i
            Q(i,j) = 1 - p(1);
        end
        if j~=i
            Q(i,j) = (p(abs(i-j))-p(abs(i-j)+1))/2;
        end
        Q(1,j) = p(j-1)/2;
        Q(n,j) = p(n-j)/2;
    end
    Q(i,n) = p(n-i)/2;
    Q(i,1) = p(i-1)/2;
end

for i=1:n-1
    for k=1:n
        K(k,1,i+1) = (1-percent_mass)*K(k,1,i);
        for j=2:n-1
            K(k,j,i+1) = (1-percent_mass)*K(k,j,i);
            K(k,j,i+1) = K(k,j,i+1)+(percent_mass)*K(k,j-1,i);
        end
        K(k,n,i+1) = K(k,n,i)+((percent_mass)*K(k,n-1,i));
    end
end

end

