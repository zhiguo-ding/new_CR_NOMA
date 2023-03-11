clear all 

M = 8;%M+1 is the size of NOMA group
N = 8; %number of time slots in one frame
R=1;
eps = 2^R-1;
T=2; %duraing of each time slot

%all the power levels, P_k
maxK = 1001;
Pinv = eps;
for i = 2 : maxK
    Pinv = [Pinv ;eps*(1+sum(Pinv(1:i-1)))];
end
for n = 1 : maxK
    eta(n) = sum(Pinv(1:n)); %for analysis
end

P = 10; %transmit SNR, or the transmit power budget
ct = 1000;% number of super-frames

snrdb =  [0: 5 :20];
for isnr = 1: length(snrdb)
    P = 10^(snrdb(isnr)/10);
    H = complex(sqrt(0.5)*randn(M+1,(M+1)*ct),sqrt(0.5)*randn(M+1,(M+1)*ct)); %M+1 rows for user; NM*ct columns for time slots
    H = abs(H).^2; 

    last_success=0; %if ct is large, we can count from zero
    yj = 0;
    last_td =0;
    yj_td = 0;

    sum3 = 0; xxx=0; sum4=0;
    for i = 1 : (M+1)*ct
        m = mod(i,M+1);%index of the frame in one super-frame
        if m==1 %first frame

            if log2(1+P*H(1,i))>=R %success
                yj(end) = yj(end) + (i-last_success)*N*T; %increase AoI
                last_success = i;% update the location of the last success           
                yj = [yj 0]; %need to restart the counting for yj

                %for tdma
                yj_td(end) = yj_td(end) + (i-last_td)*N*T; %increase AoI
                last_td = i;% update the location of the last success           
                yj_td = [yj_td 0]; %need to restart the counting for yj
            end
        else %in mth frame, P_{M-m+2} is allocated to user 1
            if m==0
                m = M+1;
            end
            Passigned = Pinv(M-m+2);
            I = P*H(m,i)/eps-1; %H(m,i) is the mth user's channel gain

            %find how many power levels the primary user can support
            K=0;
            if I<Pinv(1) %P(I<=eta_1)
                K=0;
            else
                ind_K = 1; %whether to go ahead for the while loop
                temp_k = 1; % for this else, at least K=1
                while ind_K & (temp_k<M)
                    temp_k = temp_k +1;% move on to the next 
                    if sum(Pinv(1:temp_k))>I %too many secondary users
                        ind_K=0;
                        temp_k = temp_k -1;
                    end
                end
                K = temp_k;
                xxx = [xxx K];
            end

            if K>=M-m+2 & Passigned/H(1,i)<=P %success
                yj(end) = yj(end) + (i-last_success)*N*T; %increase AoI
                last_success = i;% update the location of the last success
                yj = [yj 0]; %need to restart the counting for yj
                if m==5
                    sum3 = sum3+1;
                end
            end    
        end
    end
    aoi_sim(isnr) = T+ mean(yj.^2)/mean(yj);
    aoi_td(isnr) = T+ mean(yj_td.^2)/mean(yj_td);

    %analysis
    for m = 2: M+1
        sum1 = 0;
        for k = M-m+2: maxK-1
            sum1 = sum1 + exp(-(eps+eps*eta(k))/P ) -exp(-(eps+eps*eta(k+1))/P );
        end    
        PEjm(m) = sum1*exp(-Pinv(M-m+2)/P);
    end
    PEjm(1) = exp(-eps/P);

    Pxj1 = prod(1-PEjm);

    %EY
    tempp = 1-PEjm;
    sum1 = 0;
    for q = 0 : 1000
        for m = 1 : M+1
            for n = 1 : M+1
                if m==M+1 & n==1
                    dfdf=0;
                end            
                sum1 = sum1 + ((M+1-m)*N*T+q*(M+1)*N*T+n*N*T)...
                    *PEjm(m)*PEjm(n)*(Pxj1)^q*...
                    prod(tempp(m+1:M+1))*prod(tempp(1:n-1));
                %for m=M+1, tempp(m+1:M+1) is empty, and prod(tempp(m+1:M+1))=1
            end
        end
    end
    sum2 = 0;
    for m = 1 : M
        for n = m+1 : M+1
            sum2 = sum2 + (n-m)*N*T*PEjm(m)*PEjm(n)*prod(tempp(m+1:n-1));
        end
    end               
    EY = sum1 + sum2;

    %EY^2
    sum1 = 0;
    for q = 0 : 1000
        for m = 1 : M+1
            for n = 1 : M+1
                sum1 = sum1 + ((M+1-m)*N*T+q*(M+1)*N*T+n*N*T)^2*PEjm(m)*PEjm(n)*(Pxj1)^q*...
                    prod(tempp(m+1:M+1))*prod(tempp(1:n-1));
            end
        end
    end
    sum2 = 0;
    for m = 1 : M
        for n = m+1 : M+1
            sum2 = sum2 + ((n-m)*N*T)^2*PEjm(m)*PEjm(n)*...
                prod(tempp(m+1:n-1));
        end
    end               
    EY2 = sum1 + sum2;

    aoi_ana(isnr) = T + EY2/EY;
    
    %TDMA
    td_ana(isnr) = T+N*T*(M+1)*(2*exp(eps/P)-1);
%     P
%     p = exp(-eps/P);
%     for n = 1: 1000
%         td_ana(isnr) = td_ana(isnr) +n*N*T*(M+1)*p*(1-p)^(n-1);
%     end
end

plot(snrdb,aoi_td,snrdb,td_ana, snrdb, aoi_sim,snrdb, aoi_ana)
%[EY/(M+1) mean(yj)]
