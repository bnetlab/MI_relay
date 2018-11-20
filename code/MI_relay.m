%{
THIS FUNCTION CALCULATE MUTUAL INFORMATION OF RELAY CHANNEL
NO DELAY CONSIDERED HERE
RELEASE DISTRIBUTION IS GAUSSIAN WITH var O.234
INPUT: V, sigma,d (direct channel parameters)
ASSUMTIONS: indirect channel parameters v,sigma are same as direct channel
            The length of indirect channel is 1/2 of direct channel
OUTPUT: Mutual Information (bit)
To run: EntC(V,sigma,d)
%}
function EntC =MI_relay(V,sigma,d)
    % The main function
    % Source distribution : 
    X=[0:0.1:5];
    Pa1 = normpdf(X,2.5,sqrt(0.1656));
    Pa=zeros(1,50);
    Pa2=[Pa Pa1] ; 
    %conditional dist calculation
    %conditional distrinution PTR
    PTR=P_tr(V,sigma,d);
    %conditional distrinution PTMR   
    PTMR=P_tmr(V,sigma,d/2);
    %conditional distrinution PTMR   

    %conditional distrinution PTMR   
    
    %Total conditional distibution P
    P=0.5*(PTMR+PTR);
    % caluculate MI 
    EntC=muinfo(P,Pa2)
    %calculate varience
    var=var_mean(P)
    
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


function QYR=P_tmr(V,sigma,d)
    % THIS PROGRAMME CALCULATE CONDITIONAL DIST OF CASCADE CHANNEL
    mu= d/V;
    lambda=d^2/sigma^2;
    %conditional distrinution
    i=1;
    disp('Calculating conditional distribution PTMR')
    for j=-10:0.1:10
    QY(i)=integral(@(x)pdf('InverseGaussian',x,mu,lambda).*pdf('InverseGaussian',x+j,mu,lambda),0,50);
    i=i+1;
    end
    QYR=QY;    
    %looping for second channel  
    for ii=1
        X=[-100:1:100];
        ind=1;
        clear Q10
        for value = -100:1:0
            count=1;
            for i=1:201 %%checking for sum of time deviation in both path is del t
                for j=1:201
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

function QYR=P_tr(V,sigma,d)
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

function QYR=P_di(V,sigma,d)
    % THIS PROGRAMME CALCULATE CONDITIONAL DIST p_di (direct,indirect)
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

function QYR=P_id(V,sigma,d)
    % THIS PROGRAMME CALCULATE CONDITIONAL DIST p_id (indirect,direct)
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

function [var,avg]=var_mean(Y)
    % calculate output variance
     X=[-10:0.1:10];
     avg=0.1*trapz(X.*Y);
     X1=(X-avg).^2;
     var=0.1*trapz(X1.*Y);
end











