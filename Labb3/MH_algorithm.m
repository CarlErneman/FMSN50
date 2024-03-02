function [accepted_proposals,t] = MH_algorithm(t,lambda, tau, rho)
    %MH_ALGORITHM Summary of this function goes here
    %   Detailed explanation goes here
    
    accepted_proposals = zeros(1, length(t) - 2);
    
    
    for i = 2:length(t) - 1
        t_star = 0;
        
        %Generate t_star and make sure it is in the interval [t(i), t(i+1))
        while (~(t_star > t(i-1) && t_star < t(i+1)))
            R = rho(1) * (t(i + 1) - t(i - 1));
            epsil = 2*R*rand(1,1)-R;
            t_star =  t(i) + epsil;
        end
    
        %Generate corresponding alpha function for the t_star using the
        %posterior_t function.
        alpha = min(1, posterior_t([t(1:i-1) t_star t(i+1:end)], lambda, tau)/posterior_t(t, lambda, tau));
    
        %Draw uniform dist(0,1) for accepting or regecting the proposal
        if rand(1) <= alpha
            %Accept
            t(i) = t_star;
            accepted_proposals(i-1) = accepted_proposals(i-1)+1;
        end
    end
end

