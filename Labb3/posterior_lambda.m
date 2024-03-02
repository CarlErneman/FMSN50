function l_p = posterior_lambda(theta, t, tau) 
    % t_{i+1} - t_i
    t_difference = t(2:end) - t(1:end-1);
    d = length(t) - 1;
    n_i = zeros(1, d);
    %compute the number of disasters (n_i(tau)) in the given sub-interval  
    %to eventually compute lambda posterior
    for i = 1:d
        n_i(i) = sum((t(i) <= tau) & (tau < t(i+1)));
    end
    %posterior f(lambda|tau,t,theta), taking inverse of second parameter
    %due to MATLABs notations. 
    l_p = gamrnd(2 + n_i', 1./(t_difference' + theta));
end