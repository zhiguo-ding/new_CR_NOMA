clear all 

ct=10000;
M = 8;%number of users
R=1;
eps = 2^R-1;

%all the power levels, P_k
maxK = 1000;
Pinv = eps;
for i = 2 : maxK
    Pinv = [Pinv ;eps*(1+sum(Pinv(1:i-1)))];
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
        K = min(K,M); %capped by M
        
        hm = complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1));
        hm = abs(hm).^2;
        %random schedule
        hran = hm(1:K);
        outage = (hran*P>=Pinv(1:K));
        rate_ranx(i) = sum(outage)*R;
        
        %scheduling
        hm = sort(hm,'descend');
        hsel = hm(1:K); %here hsel(1) is the largest
        hsel = flip(hsel); %here hsel(1) is the weakest
        outage = (hsel*P>=Pinv(1:K));
        rate_selx(i) = sum(outage)*R;

    end
    rate_ran(k) = mean(rate_ranx);
    rate_sel(k) = mean(rate_selx);
    
    %analysis
    sumx1=0;
    for n = 1 : M-1       
        sumx1 = sumx1 + (exp(-(eps+eps*eta(n))/P ) -exp(-(eps+eps*eta(n+1))/P ) )*R...
            *sum(exp(-Pinv(1:n)/P));
    end
    sumx2=0;
    for n = M : maxK -1  
        sumx2 = sumx2 + (exp(-(eps+eps*eta(n))/P ) -exp(-(eps+eps*eta(n+1))/P ) );
    end        
    ana_ran(k) = sumx1 + sumx2*R*sum(exp(-Pinv(1:M)/P));
    
end
plot(snrvec,rate_ran,snrvec,ana_ran,snrvec,rate_sel)