function f_t = posterior_t(t,lambda, tau)
    %POSTERIOR_T Summary of this function goes here
    %   Detailed explanation goes here
    d = length(t)-1;
    t_diff = zeros(d,1);
    insid = zeros(d, 1);
    
    %Count the number of incidents during the intervall and calculate t_diff
    for i = 1:d
        t_diff(i) = t(i+1)-t(i);
        insid(i) = sum(tau >= t(i) & tau < t(i+1));
    end
    %Taken log of the expression, to prevent overflow in floting point
    %numbers. Otherwise we get NaN in acceptance
    f_t = exp(sum(insid.*log(lambda) + log(t_diff) - lambda.*t_diff));
end