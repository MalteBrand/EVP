function  dfini = gradfini(k,b,E,EAs,BCs,loads,OpKnoten)   
    %Berechnung der Ableitungen über finite Differenzen
    dx = 0.0000001;
    len = length(EAs);
    for i=1:length(b)
        z = zeros(len,1);
        z(i) = E*dx;
        u_pv =  trussFEM2D.solve(k,b,EAs + z,BCs,loads); %dx muss noch in u plus dx eingebaut werden
        u_mv =  trussFEM2D.solve(k,b,EAs - z,BCs,loads);
        u_p = u_pv(2 * OpKnoten);
        u_m = u_mv(2 * OpKnoten);
        dfini(:,i) = -(u_p-u_m)/(2*dx);
        %ddfini(:,i) = (u_pv-2*u+u_mv)/dx^2;
    end

end