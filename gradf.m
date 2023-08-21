function df = gradf(Ke, K, u, x, b, OpKnoten)

    for i=1:length(Ke)
        %Gradienten bestimmen
        e = zeros(length(u)-20,1);
        e((2*OpKnoten)-10) = 1; %e ist für die zu optimierende Verschiebung 1 sonst 0
        lambda = K\-e;
        lambda = [zeros(10,1);lambda;zeros(10,1)];
        dKdA = Ke(:,:,i)/x(i);
        dAdr = 2*pi*sqrt(x(i)/pi);
        dKdx = dKdA*dAdr;
        u_e = [u(2*b(i,1)-1); u(2*b(i,1)); u(2*b(i,2)-1); u(2*b(i,2))];
        lambda_e = [lambda(2*b(i,1)-1); lambda(2*b(i,1)); lambda(2*b(i,2)-1); lambda(2*b(i,2))];
        dukdx(1,i) = lambda_e' * dKdx * u_e;
    end
    
    df = dukdx;
end