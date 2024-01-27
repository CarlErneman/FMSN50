clc
clear
close all
load powercurve_D240
warning off


%wblrnd - Weibull random variables
%wblpdf - Weibull probability rensity function
%wblcdf - Weibull cumulative distribution function

lambda = [11.7 10.7 10.1 8.8 8.6 8.9 8.6 8.9 10 10.9 11.7 11.7];
k = [2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0];


N=10000;
step_size = 100;

%Predefined output vectors
tau1 =zeros(N/step_size,12);
lower_bound1 = zeros(N/step_size,12);
upper_bound1 = zeros(N/step_size,12);
conf_width1 =zeros(100,12);


%Crude monte-carlo
for month = 1:12

%To draw compare different sample sizes
    for samples = 100:step_size:N
        draw1 = wblrnd(lambda(month), k(month), 1, samples);
        draw_power1 = P(draw1);

        %Expected value
        tau1(samples/step_size, month) = mean(draw_power1);

        %confindence interval
        standard_dev = std(P(draw_power1));
        upper_bound1(samples/step_size, month) = tau1(samples/step_size, month)+abs(norminv(0.995))*std(P(draw1))/(sqrt(samples));
        lower_bound1(samples/step_size, month) = tau1(samples/step_size, month)-abs(norminv(0.995))*std(P(draw1))/(sqrt(samples));
        conf_width1(samples/step_size, month)=upper_bound1(samples/step_size, month)-lower_bound1(samples/step_size, month);
    end
end

month_plot = 1;

figure(1);
hold on
plot(100:step_size:N,tau1(:,month_plot),'LineWidth',2.5,'color','r')
plot(100:step_size:N,upper_bound1(:,month_plot),'--','LineWidth',1.5,'color','g')
plot(100:step_size:N,lower_bound1(:,month_plot),'--','LineWidth',1.5,'color','g')
hold off

%conf_width1 = upper_bound1(100, month_plot) - lower_bound1(100, month_plot)


%Truncated with standard monte-carlo, only when wind generates power
a=4;
b=30;
x=rand(1,N);

%Predefined output vectors
tau2 =zeros(N/step_size,12);
lower_bound2 = zeros(N/step_size,12);
upper_bound2 = zeros(N/step_size,12);
conf_width2 =zeros(100,12);

Fab = @(ab, month) wblinv(ab, lambda(month), k(month));
Inv = @(U, month) wblinv(U, lambda(month), k(month));




const1 = [5.8 6.5 6.5 6.5 6.5 6.5 6.5 6.5 6.5 6.5 6.5 6.5];
const2 = [3 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5];

