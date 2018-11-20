% mean and sd of discrete random variable
function [mean_avg,sd] = mean_sd_discrete(p,x)
mean_avg=sum(p.*x)
sd=sqrt(sum(x.^2.*p)-mean_avg.^2)
