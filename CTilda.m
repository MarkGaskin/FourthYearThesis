%{
Partially Observed Markov Decision Process MATLAB Model

Quantization Section

Developed by: M. Galal, M. Gaskin, I. Harbell, D. Kao

This MATLAB script takes an infinite model and quantizes it
%}

clear all;
close all;
clc;



policy_indicator = 1;       % 0 indicates threshold, 1 indicates all policies



%% Variable Declaration

T = 100;                % number of time steps
E= [0 10 100];              % set Expense cost
R = [0 1 10];              %set Repair cost
C = 0;


u = zeros(1,T);     % define control random variable u
x = zeros(1,T);     % define true state random variable x
y = zeros(1,T);     % define measured state random variable y  


u_best = zeros(1,T);    %define the RVs corresponding to best_policy
x_best = zeros(1,T);     
y_best = zeros(1,T);

measure_variance = 5/4;     % error in measurement 
iterative_variance = 1/4;   % breaking on the next step chance
fix_shift = 0.9;     % fixing success
n=2;   %number of belief bins

N=2;    %number of true states 

beta = 0.99;

% Variables for Uniform Quantization of belief
% d defines the boundaries of the bin, and v assigns a value to each bin
d = [-0.01 1/n:1/n:1];  

%Non uniform quantization
%d = [0 0.7 0.75 0.8 1];             


for i = 1:n
    v(i) = (d(i+1)+d(i))/2;
end

%Generate Probabilistic Kernels
[TK,OB] = Generate_Kernels(N, measure_variance, iterative_variance, fix_shift);

%Generate all possible policies, and place in a single matrix
if policy_indicator == 1
    policy = policy_matrix(n);
end

if policy_indicator == 0
    policy = threshold_policy(n);
end


% Calculate number of policies based on policy_indicator
policy_count = n;

if policy_indicator == 1
    policy_count = 2^n;
end

CostMatrix = zeros(N,n);

P_gamma = zeros(n,n);

%BK = Belief_Kernel(TK, OB, n, N);

for i = 1:n
    for j = 1:N
        for k=1:N
            CostMatrix(j,i) =  (1-v(i))*E(k) + R(j);
        end
    end
end


% policy_count=2;


J = zeros(n, policy_count);

for k = 1:policy_count
    for j = 1:n
        C_tilda(j) = CostMatrix((policy(k,j))+1,j);
    end
    for j = 1:n
            P_gamma(j,:) = BK(j,:,policy(k,j)+1);
    end
    
    J(:,k) = inv(eye(n)-beta*P_gamma) * transpose(C_tilda);
    
    G(1,k) = mean(J(:,k));
    
    if k == 1
        Jmin = mean(J(:,k));
        best_policy = policy(k,:);
    end
    
    if mean(J(:,k)) < Jmin
        Jmin = mean(J(:,k));
        best_policy = policy(k,:);
     end
    
end

[SortedG, SortedG_index] = sort(G);



