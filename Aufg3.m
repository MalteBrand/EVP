function Aufg3(anzRealMC)

    [n,m,rmin,rmax,dx,dy,r0,A,E,Fex,k,b,nob,EAs,rs,BCs,loads] = preprocess();
    
    OpKnoten = (floor(n/2)*m) + 1;

    fhf = @(x) ziel(x, E, k, b, BCs, loads, OpKnoten);
    fhgtb = @(x) nebenBedingungen4toolbox(x, rmin, rmax, V0, k, b, E, BCs, loads);
    
    %Optimierung mit SQP
    options = optimoptions('fmincon','Algorithm','sqp','Display','iter','StepTolerance',0.000000000001);
    x_opt = fmincon(fhf,x0,[],[],[],[],[],[],fhgtb,options);

    function [f,df] = ziel(x, E, k, b, BCs, loads, OpKnoten,anzRealMC)
        
        [mu_u,sigma2_u,dmu_u,dsigma2_u] = MC_FEM_RDO(x, E, k, b, BCs, loads, OpKnoten,anzRealMC);

        f = mu_u + 3* sqrt(sigma2_u);
        df = dmu_u + 3* 1/(2*sqrt(sigma2_u)) * dsigma2_u;
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

    function [nb,nbeq,dnb,dnbeq] = nebenBedingungen4toolbox(x)
        nb = ungl_bed(x, rmin, rmax, k, b, E, BCs, loads);
        nbeq = gl_bed(x, V0, k, b);
    end

    function [mu_g,sigma2_g,dmu_g,dsigma2_g] = MC_FEM_RDO(x, E, k, b, BCs, loads, OpKnoten,anzRealMC)
        
        % deterministische Eingangsgrößen
        EAs = E*x;

        % stochastische Eingangsgrößen
        mu_stat = k; %Erwartungswert entspricht den unverschobenen Knotenpositionen
        sig_stat = 0.03; % 10% der kleinsten Stablänge
        lc = 1.5; %Korrelationslänge

        [CovMa,CovSq] = Zufallsfeld(mu_stat, sig_stat, lc, k);
        

        tol=1.0;  % bd rand
        % Schleife über Anzahl Realisierungen

        rng default;
        g = zeros(AnzRea,2);
        %dg = zeros(4,AnzRea);

        % Schleife über Anzahl Realisierungen
        for i=1:anzRealMC

            z = [k(1:5,:);
            randn(length(k)-10,2);
            k(51:55,:)];
            % Erzeugen von Zufallszahlen
            x_MC = CovSq * z + mu_X;

            % bd_MC = rand(4,1)*tol*2 + [b1;d1;b2;d2]-tol;  % bd rand

            % Auswerten der Zielfunktion
            [u,~,Ke,K] = trussFEM2D.solve(x_MC,b,EAs,BCs,loads);
            g1(1,i) = u(OpKnoten-1);
            g2(1,i) = u(OpKnoten);
            dg1(i,:) = gradf(Ke, K, g(i,1), x, b, OpKnoten);
            dg2(i,:) = gradf(Ke, K, g(i,2), x, b, OpKnoten,2);

        end
        g = [g1; g2];
    	dg = [dg1;dg2];
        % Mittelwert, Varianz 
        mu_g1 = mean(g1);
        mu_g2 = mean(g2);
        sigma2_g1 = var(g1);
        sigma2_g2 = var(g2);
        dmu_g1 = mean(dg1,1);
        dmu_g2 = mean(dg2,1);
        dmu_g = [dmu_g1; dmu_g2];

        dsigma2_g = 2/(AnzRea-1) * (g(:,1)*dg1 - 2*mu_g1*dmu_g2 + g2*dg2 - 2*mu_g1*dmu_g2);
    end
end