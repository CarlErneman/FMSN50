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
        upper_bound1(samples/step_size, month) = tau1(samples/step_size, month)+abs(norminv(0.995))*std(draw_power1)/(sqrt(samples));
        lower_bound1(samples/step_size, month) = tau1(samples/step_size, month)-abs(norminv(0.995))*std(draw_power1)/(sqrt(samples));
        conf_width1(samples/step_size, month)=upper_bound1(samples/step_size, month)-lower_bound1(samples/step_size, month);
    end
end

month_plot = 1;

conf_width1(100, month_plot)

figure(1);
hold on
title("Crude Monte-Carlo")
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
conf_width2 = zeros(100,12);

Fab = @(ab, month) wblcdf(ab, lambda(month), k(month));
Inv = @(U, Fa, Fb) wblinv(U*(Fb-Fa) + Fa);

for month = 1:12
%To draw compare different sample sizes
    Fa = Fab(a, month);
    Fb = Fab(b, month);

    for samples = 100:step_size:N
        %Inverse of calculated in 1b) with a scaing factor (Fb-Fa)
        draw2 = Inv(rand(1,samples), Fa, Fb)*(Fb-Fa); 


        draw_power2 = P(draw2);

        %Expected value
        tau2(samples/step_size, month) = mean(draw_power2);

        %confindence interval
        standard_dev = std(P(draw_power1));
        upper_bound2(samples/step_size, month) = tau2(samples/step_size, month)+abs(norminv(0.995))*std(draw_power2)/(sqrt(samples));
        lower_bound2(samples/step_size, month) = tau2(samples/step_size, month)-abs(norminv(0.995))*std(draw_power2)/(sqrt(samples));
        conf_width2(samples/step_size, month)=upper_bound2(samples/step_size, month)-lower_bound2(samples/step_size, month);
    end
end

conf_width2(100, month_plot)

figure(2);
hold on
title("Truncated with Crude Monte-Carlo")
plot(100:step_size:N,tau1(:,month_plot),'LineWidth',2.5,'color','r')
plot(100:step_size:N,upper_bound1(:,month_plot),'--','LineWidth',1.5,'color','g')
plot(100:step_size:N,lower_bound1(:,month_plot),'--','LineWidth',1.5,'color','g')
hold off

