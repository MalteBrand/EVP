function df = gradf2(Ke, K, u, x, b, OpKnoten)
    if nargin<7
        nz = 1;
    end

    %Gradienten bestimmen
    nu = length(u);
    e = zeros(nu-20,1);
    e((2*OpKnoten)-11) = 1; %e ist für die zu optimierende Verschiebung 1 sonst 0

    lambda = K^2\-e;
    lambda = [zeros(10,1);lambda;zeros(10,1)];
    
    len = length(Ke);
    dukdx = zeros(1,len);
    for i=1:length(Ke)

        dKdA = Ke(:,:,i)/x(i);
        dKdx = dKdA;
        u_e = [u(2*b(i,1)-1); u(2*b(i,1)); u(2*b(i,2)-1); u(2*b(i,2))];
        lambda_e = [lambda(2*b(i,1)-1); lambda(2*b(i,1)); lambda(2*b(i,2)-1); lambda(2*b(i,2))];
        dukdx(1,i) = 2*lambda_e' * dKdx * u_e.^2;

    end
    
    df = dukdx;

end