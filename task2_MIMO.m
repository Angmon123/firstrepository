%% task2 MIMO ZF and MMSE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

Nt=1;  %transmitteer antenna number
Nr=2;  %receiver antenna number
N=10;
L=10000;
EbN0=0:2:20; %%SNR sequence 
M=4;  %QPSK
x=randi(M-1,N*L,Nt);

%% 2*2 antenna ZF
[ber1,ber2]=ZF(Nt, Nr, N, L, EbN0, M, x)  ;
figure(1)
semilogy(EbN0,ber1,'-k*',EbN0,ber2,'-kd');  %BER range of 10^-5
legend('ZF without perfect channel knowledge','ZF with perfect channel knowledge');
xlabel('EbN0(dB)')
ylabel('BER')
title('ZF equalization 2*2 antenna');  

%% 2*5 antenna ZF
Nt=2;  %transmitteer antenna number
Nr=5;  %receiver antenna number
[ber1,ber2]=ZF(Nt, Nr, N, L, EbN0, M, x)  ;
figure(2)
semilogy(EbN0,ber1,'-k*',EbN0,ber2,'-kd');  %BER range of 10^-5
legend('ZF without perfect channel knowledge','ZF with perfect channel knowledge');
xlabel('EbN0(dB)')
ylabel('BER')
title('ZF equalization 2*5 antenna');  

%% 2*10 antenna ZF
Nt=2;  %transmitteer antenna number
Nr=10;  %receiver antenna number
[ber1,ber2]=ZF(Nt, Nr, N, L, EbN0, M, x)  ;
figure(3)
semilogy(EbN0,ber1,'-k*',EbN0,ber2,'-kd');  %BER range of 10^-5
legend('ZF without perfect channel knowledge','ZF with perfect channel knowledge');
xlabel('EbN0(dB)')
ylabel('BER')
title('ZF equalization 2*10 antenna');  

%% 2*2 antenna MMSE
Nt=2;  %transmitteer antenna number
Nr=2;  %receiver antenna number
[ber1,ber2]=MMSEhandling(Nt, Nr, N, L, EbN0, M, x)  ;
figure(4)
semilogy(EbN0,ber1,'-k*',EbN0,ber2,'-kd');  %BER range of 10^-5
legend('MMSE without perfect channel knowledge','MMSE with perfect channel knowledge');
xlabel('EbN0(dB)')
ylabel('BER')
title('MMSE 2*2 antenna');

%% 2*5 antenna MMSE
Nt=2;  %transmitteer antenna number
Nr=5;  %receiver antenna number
[ber1,ber2]=MMSEhandling(Nt, Nr, N, L, EbN0, M, x)  ;
figure(5)
semilogy(EbN0,ber1,'-k*',EbN0,ber2,'-kd');  %BER range of 10^-5
legend('MMSE without perfect channel knowledge','MMSE with perfect channel knowledge');
xlabel('EbN0(dB)')
ylabel('BER')
title('MMSE equalization 2*5 antenna');

%% 2*10 antenna MMSE
Nt=2;  %transmitteer antenna number
Nr=10;  %receiver antenna number
[ber1,ber2]=MMSEhandling(Nt, Nr, N, L, EbN0, M, x)  ;
figure(6)
semilogy(EbN0,ber1,'-k*',EbN0,ber2,'-kd');  %BER range of 10^-5
legend('MMSE without perfect channel knowledge','MMSE with perfect channel knowledge');
xlabel('EbN0(dB)')
ylabel('BER')
title('MMSE equalization 2*10 antenna');
        
%% function        
function [ber1,ber2]=ZF(Nt, Nr, N, L, EbN0, M, x)

s=pskmod(x,M,pi/4);  %% QPSK modulation; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(EbN0)  %%Loop all EbN0  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    s1=[];
    s2=[];
    s3=[];
    for k=1:L
        h=randn(Nt,Nr)+ 1i*randn(Nt,Nr); %% rayleigh fading channel;
        h=h./sqrt(2); % channel nomilization
        [q1,r1]=qr(h');
        r=r1(1:Nt,:)';
        q=q1(:,1:Nt)';  % Q,R decomposed;
        sigma1=sqrt(1/(10.^(EbN0(i)/10))); %% each antenna gaussian noise standard deviation
        n=sigma1*(randn(N,Nr) + 1i*randn(N,Nr));
        
        y=s((k-1)*N+1:k*N,:)*h*q'+n*q'; %% signal through Rayleigh channel of AWGN
        
        y1= y*inv(r);             % ZF without perfect channel knowledge
        s1=[s1;pskdemod(y1,M,pi/4)];
        
        y(:,Nt)=y(:,Nt)./(r(Nt,Nt));
        y1(:,Nt)=pskdemod(y(:,Nt),M,pi/4);
        y(:,Nt)=pskmod(y1(:,Nt),M,pi/4); %re-module
        
        y2=y;
        y3=y1;
        for m=Nt-1:-1:1 %loop all layers
            for n=m+1:Nt
                y(:,m)=y(:,m)-r(n,m).*y(:,n);  
                y2(:,m)=y2(:,m)-r(n,m).*s((k-1)*N+1:k*N,n);   % ZF with perfect channel knowledge
            end
            y1(:,m)=y2(:,m)./r(m,m);
            y2(:,m)=y2(:,m)./r(m,m);
          
            y1(:,m)=pskdemod(y(:,m),M,pi/4);
            y3(:,m)=pskdemod(y2(:,m),M,pi/4);
            y(:,m)=pskdemod(y1(:,m),M,pi/4);
            y2(:,m)=pskdemod(y3(:,m),M,pi/4);
        end
        s3=[s3;y3];
    end 
    
    [temp,ber1(i)]=biterr(x,s1,log2(M));
    [temp,ber2(i)]=biterr(x,s3,log2(M));
             
end
   
end        
              


function [ber1,ber2]=MMSEhandling(Nt, Nr, N, L, EbN0, M, x)

s=pskmod(x,M,pi/4);  %% QPSK modulation; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(EbN0)  %%Loop all EbN0  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s1=[];
    s2=[];
    for k=1:L
        h=randn(Nt,Nr)+ 1i*randn(Nt,Nr); %% rayleigh fading channel;
        h=h./sqrt(2); % channel nomilization
        sigma1=sqrt(1/(10.^(EbN0(i)/10))); %% each antenna gaussian noise standard deviation
        n=sigma1*(randn(N,Nr) + 1i*randn(N,Nr));
        w=h'*inv(h*h'+2*sigma1.^2*diag(ones(1,Nt))); %optimal solution
        
        y=s((k-1)*N+1:k*N,:)*h+n; %% signal through Rayleigh channel of AWGN
        tempy=y;
        y1=y*w;
        temp1=pskdemod(y1,M,pi/4); %demodulation；
        s1=[s1;temp1];

        temp2(:,Nt)=temp1(:,Nt);
        tempy=tempy-s((k-1)*N + 1:k*N,Nt)*h(Nt,:);
        
        h=h(1:Nt-1,:);
        for m=Nt-1:-1:1
            w=h'*inv(h*h'+2*sigma1.^2*diag(ones(1,m))); % w matrix update
            tempy2=tempy*w;
            temp2(:,m)=pskdemod(tempy2(:,m),M,pi/4);
            tempy=tempy-s((k-1)*N+1:k*N,m)*h(m,:);  %% MMSE with  with perfect channel knowledge
            h=h(1:m-1,:); 
        end
        
        s2 = [s2;temp2];
        
    end
    [~,ber1(i)]=biterr(x,s1,log2(M));
    [~,ber2(i)]=biterr(x,s2,log2(M));
end
end












        