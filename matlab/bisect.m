function [ rho, lambda] = bisect( rho_old, rho_min, V, df, dg, ne)

lambda1 = 10^(-10);
lambda2 = 10^10;
n = 0.1;
rho=rho_old*0;

while (lambda2-lambda1)/(lambda1+lambda2) > 10^(-5)
    
    % calculates lambda mid
    lambda_mid = (lambda1+lambda2)/2;   
    
    % calculates rho using lambda_mid
    for e = 1:ne
        
        Be = -df(e)/(lambda_mid*dg(e));
        
        if rho_old(e) * Be^n <= rho_min
            rho(e,1) = rho_min;
            
        else
            if rho_old(e) * Be^n >= 1
                rho(e,1) = 1;
            else
                rho(e,1) = rho_old(e)*Be^n;
            end
        end
        
    end

    % calculates g constraint 
    g = rho'*dg - V;
    
    if g > 0
        lambda1 = lambda_mid;
    else
        lambda2 = lambda_mid;
    end
end

lambda = lambda_mid;

end

