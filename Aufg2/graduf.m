function dudf = graduf(u,K,OpKnoten,n,m)
    %Ableitung von u nach F
    nu = length(u);
    ek = zeros(nu-20,1); %ek = 1 an Knoten 30 in x-Richtung wo horizontale Kraft angreift
    ek(2*((floor(n/2)+1)*m)-11) = 1;
    lambda = K\-ek;
    
    ej = zeros(nu-20,1);
    ej((2*OpKnoten)-11) = 1; %ej = 1 an Knoten 26 x-Richtung wo horizontale Kraft ausgewertet wird
    dudf(1,1) = -lambda'*ej;
    
    ej = zeros(nu-20,1);
    ej((2*OpKnoten)-10) = 1; %ej = 1 an Knoten 26 y-Richtung wo horizontale Kraft ausgewertet wird
    dudf(2,1) = -lambda'*ej;
end