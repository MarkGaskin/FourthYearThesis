%{
Partially Observed Markov Decision Process MATLAB Model

Quantization Section

Developed by: M. Galal, M. Gaskin, I. Harbell, D. Kao

This MATLAB script quantizes the three state model and tests policies
%}

clear all;
close all;
clc;

policy_indicator = 0;       % 0 indicates threshold, 1 indicates all policies

% Variable Declaration

T = 10000;            % number of time steps
E= [0 1 10];              % set Expense cost
R = [0 5];              %set Repair cost
C = 0;


u = zeros(1,T);     % define control random variable u
x = zeros(1,T);     % define true state random variable x
y = zeros(1,T);     % define measured state random variable y  


u_best = zeros(1,T);    %define the RVs corresponding to best_policy
x_best = zeros(1,T);     
y_best = zeros(1,T);

measure_variance = 5/4;     % error in measurement 
iterative_variance = 1/4;   % breaking on the next step chance
percent_mass = 0.9;     % fixing success
n=64;   %number of belief bins

N=3;    %number of true states 

b_count = zeros(1,n);
b = zeros(2,T);
b(2,1)= 0;      %define initial belief
b(3,1) = 1;


% d defines the boundaries of the bin, and v assigns a value to each bin
d = [0:1/n:1];  
            
for i = 1:n
    v(i) = (d(i+1)+d(i))/2;
end

%Generate Probabilistic Kernels
[TK,OB] = Generate_Kernels(N, measure_variance, iterative_variance, percent_mass);

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

Cmin=T*E(N);               % define C min as largest possible cost initially
Cost = zeros(1,policy_count);     % define Cost vector to track costs across the different policies
sumx = zeros(1,policy_count);

    
for k = 1:policy_count        % iterate through all admissible policies
    
    C=0;                % Reset cost to zero
    
    for t=2:T               % iterate through T time steps
        
        % random numbers used for the random variables' probabilities
        rand_meas = rand(1);
        rand_fail = rand(1);
        rand_break = rand(1);
        
        % take measurement of current state
        for i=1:N
            if rand_meas < sum(OB(x(t)+1,1:i)) && rand_meas >= sum(OB(x(t)+1,1:(i-1)))
                y(t) = i-1;
            end
        end
        
        
        % update belief
        for i = 1:N
            b(i,t) = (transpose(b(:,t-1)) * TK(:,i,u(t-1)+1) * OB(i,y(t)+1)) /(transpose(b(:,t-1)) * TK(:,:,u(t-1)+1) * OB(:,y(t)+1));
        end
        
        % Code for equal triangle quantization
%         if b(1,t) > 1/2
%             current_bin = 1;
%         elseif b(2,t) > 1/2
%             current_bin = 2;
%         elseif b(3,t) > 1/2
%             current_bin = 3;
%         else
%             current_bin = 4;
%         end

        %Code for threshold straight line quantization
        for i = 1:n
            if b(3,t) > ((d(i)/(1-d(i)))*b(1,t))
                current_bin = i;
            end
        end

     
    % use strategy to determine control u
    u(t) = policy(k, current_bin);
    
    % determine whether fix was successful or not
    for i=1:N
        if rand_break < sum(TK(x(t)+1,1:i,u(t)+1)) && rand_break >= sum(TK(x(t)+1,1:(i-1),u(t)+1))
            x(t+1) = i-1;
        end
    end
    
    % update cost function
    
    C = C + R(1+u(t));      % add repair cost
    
    
    C = C + E(3-x(t));      % add expense cost

end

Cost(k) = C;           % put C into Cost vector
sumx(k) = sum(x);      % track how often the system worked      

%determine best policy
if C < Cmin
        Cmin = C;
        u_best = u;
        y_best = y;
        x_best = x;
        b_best = b;
        best_policy = policy(k,:);
end
end

% Sort costs from lowest to highest along with their associated policies
[CostMatrix(1,:),CostMatrix(2,:)] = sort(Cost);
if policy_indicator == 1    
    for i=1:policy_count
        CostMatrix(2,i)=str2num(dec2bin(CostMatrix(2,i)-1,n));
    end
end  

fprintf('Best policy: %.2f %.2f %.2f %.2f \n', best_policy)
fprintf('The system was working for %d time steps out of %d at a cost of %d \n', sum(x_best),T,Cmin)
