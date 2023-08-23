function Aufg1()

    %Initialparameter
    [n,m,rmin,rmax,dx,dy,r0,A,E,Fex,k,b,nob,EAs,rs,BCs,loads] = preprocess();
    
    %%a) Erstellen eines Referenzentwurf
    u0 = trussFEM2D.solve(k,b,EAs,BCs,loads);
    trussFEM2D.plotTruss2D(k,b,rs,2,u0,1);
    
    %%b) 
    OpKnoten = (floor(n/2)*m) + 1; %Berechnungvorschrift für Punkt P bzw. Knoten 26 
    x0 = A*ones(nob,1);
    V0 = trussFEM2D.mass(k,b,A);
    
    fhf = @(x) ziel(x, E, k, b, BCs, loads, OpKnoten);
    fhgtb = @(x) nebenBedingungen4toolbox(x, rmin, rmax, V0, k, b, E, BCs, loads);
    
    %Optimierung mit SQP
    options = optimoptions('fmincon','Algorithm','sqp','Display','iter','StepTolerance',0.000000000001);
    x_opt = fmincon(fhf,x0,[],[],[],[],[],[],fhgtb,options);
    
    A = x_opt;
    EAs = E*A;
    
    %Optimierte Verschiebung
    u = trussFEM2D.solve(k,b,EAs,BCs,loads);
    
    rs = sqrt(A/pi);
    
    %Optimierter Plot ohne Last
    trussFEM2D.plotTruss2D(k,b,rs,3);
    %Optimierter Plot mit Last
    trussFEM2D.plotTruss2D(k,b,rs,4,u,1);
    
    function [f,df,ddf] = ziel(x, E, k, b, BCs, loads, OpKnoten)
    
        EAs = E*x;
        
        %Berechnung der Verschiebung
        [u,~,Ke,K] = trussFEM2D.solve(k,b,EAs,BCs,loads); 
        u = -u; %Verschiebung kann positiv werden
        %u = abs(u);
        %u = u.^2;
    
        %Es soll die y-Verschiebung in Knoten 26 optimiert werden
        f = u(2 * OpKnoten);
        df = gradf(Ke, K, u, x, b, OpKnoten);
        ddf = zeros(174,174);
    
        %dfini = gradfini(k, b, E, EAs, BCs, loads, OpKnoten);
        % valid = dfini(1,1:10) - df(1,1:10);
    
    end
    
    function g = ungl_bed(x, rmin, rmax, k, b, E, BCs, loads)
        len = length(x);
        g(1:len)                = -(x - pi*rmin^2); %r > rmin
        g( (len+1):(2*len) )    = -(pi*rmax^2 - x); %r < rmax
        EAs = E*x;
        u = trussFEM2D.solve(k,b,EAs,BCs,loads);
        g(2*len+1) = u(2 * OpKnoten); % u < 0
    end
    function h = gl_bed(x,V0,k,b)   
        
        h = (trussFEM2D.mass(k,b,x) - V0); % V = V0 
    
    end
    
    function [nb,nbeq] = nebenBedingungen4toolbox(x, rmin, rmax, V0, k, b, E, BCs, loads)
        nb = ungl_bed(x, rmin, rmax, k, b, E, BCs, loads);
        nbeq = gl_bed(x, V0, k, b);
    end

end