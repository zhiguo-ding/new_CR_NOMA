clear all 

ct=100000;
M = 32;%number of users
R=1.5;
eps = 2^R-1;

%all the power levels, P_k
maxK = 1000;
Pinv = eps;
for i = 2 : maxK
    Pinv = [Pinv eps*(1+sum(Pinv(1:i-1)))];
end
for n = 1 : maxK
    eta(n) = sum(Pinv(1:n)); %for analysis
end

snrvec = [0: 5 : 20];
for k = 1: length(snrvec)
    P = 10^(snrvec(k)/10);
    sum2=0; sum3 = 0; sum4=0;
    for i = 1: ct

        h0 = complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1));
        h0 = abs(h0).^2;
        I = P*h0/eps-1;
       
        if I<Pinv(1) %P(I<=eta_1)
            K=0;
        else
            ind_K = 1; %whether to go ahead for the while loop
            temp_k = 1; % for this else, at least K=1
                while ind_K
                    temp_k = temp_k +1;% move on to the next 
                    if sum(Pinv(1:temp_k))>I %too many secondary users
                        ind_K=0;
                    end
                end
                K = temp_k-1; %back to one, since we moved on to the next
        end
        Kvectemp(i) = K;    
        
        if Kvectemp(i) == 0
            sum2 = sum2+1;
        end
    end
    Kvec(k) = mean(Kvectemp);
    %analysis
    K_est(k) = 0*(exp(-eps/P ) -exp(-(eps^2+eps)/P ) );
    for n=1: maxK-1
        K_est(k) = K_est(k) + ...
            n*(exp(-(eps+eps*eta(n))/P ) -exp(-(eps+eps*eta(n+1))/P ) );
    end
    
end
plot(snrvec,Kvec,snrvec,K_est)

% % illustration one
% P=100000000000000;
% for n =1 : maxK/10
%     KP(n) = exp(-(eps+eps*eta(n))/P ) -exp(-(eps+eps*eta(n+1))/P );
% end
% plot(KP)

% %test
% n=20;
% eta(n+1)-eta(n) - eps*(1+eps)^n
 