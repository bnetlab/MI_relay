%{
THIS FUNCTION CALCULATE MUTUAL INFORMATION OF RELAY CHANNEL
NO DELAY CONSIDERED HERE
RELEASE DISTRIBUTION IS GAUSSIAN WITH var O.234
INPUT: V, sigma,d (direct channel parameters)
ASSUMTIONS: indirect channel parameters v,sigma are same as direct channel
            The length of indirect channel is 1/2 of direct channel
OUTPUT: Mutual Information (bit)
To run: F=main();
        F.MI_relay(V,sigma,d, ratio)
%}

function Fn = main
    %This function is not actual code
    % it's a wrapper : To call subfunction from outside
    Fn.MI_relay = @MI_relay;
    Fn.P_tr= @P_tr;
    Fn.P_tmr=@P_tmr;
    Fn.P_di=@P_di;
    Fn.P_id=@P_id;
    Fn.muinfo=@muinfo;
    Fn.var_mean=@var_mean;
end

function [EntC,var] =MI_relay(V,sigma,d, ratio)
    % The main function
    % Source distribution : 
    X=[0:0.1:5];
    Pa1 = normpdf(X,2.5,sqrt(0.1656));
 %   Pa1 = normpdf(X,2.5,sqrt(0.25));
    Pa=zeros(1,50);
    Pa2=[Pa Pa1] ; 
    %conditional dist calculation
    %conditional distrinution PTR
    PTR=P_tr(V,sigma,d, ratio);
    sum(PTR)
    plot(PTR, 'LineWidth',2)
    hold on
    %conditional distrinution PTMR   
    PTMR=P_tmr(V,sigma,d, ratio);
    sum(PTMR)
    plot(PTMR, 'LineWidth',2)
    %conditional distrinution PDI   
    PDI=P_di(V,sigma,d, ratio);
    sum(PDI)
    plot(PDI, 'LineWidth',2)
%     conditional distrinution PID 
    % PID is just mirror reflection of PDI about middle point 
    % because del tau is independent of first arrival time (t)
    PID=[fliplr(PDI(102:end)),PDI(101),fliplr(PDI(1:100))];
    sum(PID)
    plot(PID, 'LineWidth',2)
    %Total conditional distibution P
    P=0.25*(PTMR+PTR+PDI+PID);
    plot(P, 'LineWidth',2)
    sum(P)
    legend('PTR','PTMR','PDI','PID','P');
    % caluculate MI
    disp('MI TR :')
    muinfo(PTR,Pa2)
    disp('MI TMR')
    muinfo(PTMR,Pa2) 
    disp('MI DI')
    muinfo(PDI,Pa2)  
    disp('MI ID')
    muinfo(PID,Pa2)
    disp('MI of system :')
    EntC=muinfo(P,Pa2);
    %calculate varience
    var=var_mean(P);   
end

function out = mylog2(in)
  % calculate log2 considering 0*log2(0)=0
  out = log2(in);
  out(~in) = 0;
end

function EntC=muinfo(QYR,Pa2)
    % calculate mutual Information
    % QYR length -10sec to 10 sec 
    % Pa2 length -5 sec to 5 sec
    %calculate h(Y|X)
    EntC1=0;
    for i=1:100
        EntC1=EntC1-0.1*Pa2(i)*0.1*sum(QYR(i:100+i).*mylog2(QYR(i:100+i)));
    end    
    % calculate p(Y)
    for i=1:101
        Q12=QYR(i:i+100);
        M10(i)= 0.1* sum(fliplr(Q12).*Pa2);
    end	
	% claculate h(Y)
    EntC2=0.1*sum(-(M10(M10>0).*(mylog2(M10(M10>0)))));    
	% calculate MI
    EntC=EntC2-EntC1;
end

function QYR=P_tr(V,sigma,d, ratio)
    % THIS PROGRAMME CALCULATE CONDITIONAL DIST OF DIRECT CHANNEL
    mu= d/V;
    lambda=d^2/sigma^2;
    %conditional distrinution
    i=1;
    disp('Calculating conditional distribution PTR')
    for j=-10:0.1:10
    QY(i)=integral(@(x)pdf('InverseGaussian',x,mu,lambda).*pdf('InverseGaussian',x+j,mu,lambda),0,50);
    i=i+1;
    end
    QYR=QY;
end

function QYR=P_tmr(V,sigma,d, ratio)
    % THIS PROGRAMME CALCULATE CONDITIONAL DIST OF CASCADE CHANNEL
    d=d*ratio/2;
    mu= d/V;
    lambda=d^2/sigma^2;
    %conditional distrinution
    i=1;
    disp('Calculating conditional distribution PTMR')
    for j=-50:0.1:50
    QY(i)=integral(@(x)pdf('InverseGaussian',x,mu,lambda).*pdf('InverseGaussian',x+j,mu,lambda),0,50);
    i=i+1;
    end
    QYR=QY;    
    %looping for second channel  
    for ii=1
        X=[-500:1:500];
        ind=1;
        clear Q10
        for value = -100:1:0
            count=1;
            for i=1:1001 %%checking for sum of time deviation in both path is del t
                for j=1:1001
                    if (X(i)+X(j)==value)
                        Xa(count)=i;
                        Ya(count)=j;
                        count=count+1;
                    end
                end
            end
            Q10(ind)=sum(QYR(Xa).*QY(Ya));
            ind=ind+1;
            clear Xa
            clear Ya
        end
        Q10= 0.1*[Q10 fliplr(Q10(1:100))];
        trapz(Q10);
        QYR=Q10;
    end
end

function QYR=P_di(V,sigma,d, ratio)
    % THIS PROGRAMME CALCULATE CONDITIONAL DIST p_di (direct,indirect)
    d=d*ratio/2;
    mu= d/V;
    lambda=d^2/sigma^2;
    %QD time required to travel direct channel
    x=0:0.1:50;
    QD=pdf('InverseGaussian',x,mu,lambda);
    %QI time required to travel indirect channel
    i=1;
    for j=0:0.1:50
    QI(i)=integral(@(x)pdf('InverseGaussian',j-x,mu,lambda).*pdf('InverseGaussian',x,mu,lambda),0,50);
    i=i+1;
    end
    %conditional distrinution
    disp('Calculating conditional distribution PDI and PID')
    X=[0:1:500];
    ind=1;
    clear Q10
    for value = -100:1:100
       count=1;
       for i=1:501 %%checking for sum of time deviation in both path is del t
            for j=1:501
                if (X(i)-X(j)==value)
                    Xa(count)=i;
                    Ya(count)=j;
                    count=count+1;
                end
            end
        end
        QYR(ind)=0.1*sum(QD(Xa).*QI(Ya));
        ind=ind+1;
        clear Xa
        clear Ya
    end
    %QYR=0.1*QYR;
    %QYR=QYR(401:601);
end

% function QYR=P_id(V,sigma,d)
% 
% end

function [var,avg]=var_mean(Y)
    % calculate output variance
     X=[-10:0.1:10];
     avg=0.1*trapz(X.*Y);
     X1=(X-avg).^2;
     var=0.1*trapz(X1.*Y);
end











