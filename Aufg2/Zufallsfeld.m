function [CovMa,CovSq] = Zufallsfeld(mu_stat, sig_stat, lc, k)
    nx = length(mu_stat);
    
    R = ones(nx,nx);
    for i = 1:nx
        for j=i+1:nx
            R(i,j) = exp( -( (k(i,1)-k(j,1))^2 + (k(i,2)-k(j,2))^2 ) /(2*lc^2));
            R(j,i) = R(i,j);
        end
    end

    CovMa = R*sig_stat^2;
        
    CovSq = sqrtm(CovMa);

end