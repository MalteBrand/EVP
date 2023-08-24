function Aufg3(anzRealMC)
    anzRealMC = 100;
    [n,m,rmin,rmax,dx,dy,r0,A,E,Fex,k,b,nob,EAs,rs,BCs,loads] = preprocess();
    
    rng default;
    global z;

    z = [k(1:5,:);
    randn(length(k)-10,2);
    k(51:55,:)];

    % stochastische Eingangsgrößen
    mu_stat = k; %Erwartungswert entspricht den unverschobenen Knotenpositionen
    sig_stat = 0.03; % 10% der kleinsten Stablänge
    lc = 1.5; %Korrelationslänge

    [CovMa,CovSq] = Zufallsfeld(mu_stat, sig_stat, lc, k);

    OpKnoten = (floor(n/2)*m) + 1; %Berechnungvorschrift für Punkt P bzw. Knoten 26 
    x0 = A*ones(nob,1);
    V0 = trussFEM2D.mass(k,b,A);

    fhf = @(x) ziel(x, E, k, b, BCs, loads, OpKnoten,anzRealMC, CovSq, mu_stat);
    fhgtb = @(x) nebenBedingungen4toolbox(x, rmin, rmax, V0, k, b, E, BCs, loads, CovSq, mu_stat);
    
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

    function [f,df] = ziel(x, E, k, b, BCs, loads, OpKnoten,anzRealMC, CovSq, mu_stat)
        
        [mu_u,sigma2_u,dmu_u,dsigma2_u] = MC_FEM_RDO(x, E, k, b, BCs, loads, OpKnoten,anzRealMC, CovSq, mu_stat);
        w1 = 0.7;
        f = w1*( mu_u(1) + 3*sqrt( sigma2_u(1) ) )+ (1-w1)*( mu_u(2) + 3 * sqrt(sigma2_u(2)) ) ;

        df = w1*( dmu_u(1,:) + 3* 1/(2*sqrt( sigma2_u(1) )) * dsigma2_u(1,:) ) + (1-w1)*( dmu_u(2,:) + 3* 1/(2*sqrt( sigma2_u(2))) * dsigma2_u(2,:) );

    end

    function [g, dg] = ungl_bed(x, rmin, rmax, k, b, E, BCs, loads, CovSq, mu_stat)
        len = length(x);
        g(1:len)                = -(x - pi*rmin^2); %r > rmin
        g( (len+1):(2*len) )    = -(pi*rmax^2 - x); %r < rmax
        EAs = E*x;
        
        %[mu_u,sigma2_u,dmu_u,dsigma2_u] = MC_FEM_RDO(x, E, k, b, BCs, loads, OpKnoten,anzRealMC, CovSq, mu_stat);
        %g(2*len+1) = mu_u(2); % u < 0

        dg = zeros((2*len + 1), len);

        dg((1:len),:) = -1* eye(len);
        dg( (len+1):(2*len),: ) = eye(len);
        %dg(2*len+1,:) = dmu_u(2,:);

    end

    function h = gl_bed(x,V0,k,b)   
        len = length(x);

        h = (trussFEM2D.mass(k,b,x) - V0); % V = V0 
        
        x1=k(b(:,1),1);
        y1=k(b(:,1),2);
        x2=k(b(:,2),1);
        y2=k(b(:,2),2);

        dh = zeros(1,174);

        for i=1:len
            dh(1,i) = (sqrt((x2(i)-x1(i))^2 + (y2(i)-y1(i))^2));
        end
    
    end

    function [nb,nbeq,dnb,dnbeq] = nebenBedingungen4toolbox(x, rmin, rmax, V0, k, b, E, BCs, loads, CovSq, mu_stat)
        nb = ungl_bed(x, rmin, rmax, k, b, E, BCs, loads, CovSq, mu_stat);
        nbeq = gl_bed(x, V0, k, b);
    end

    function [mu_g,sigma2_g] = FOSM_FEM_RDO(x, E, k, b, BCs, loads, OpKnoten,anzRealMC, CovSq, mu_stat)
        
        % deterministische Eingangsgrößen
        EAs = E*x;

        %rng default;

        dmu_g = zeros(2,174);
        sigma2_g = zeros(2,1);
        dsigma2_g = zeros(2,174);

        [u,~,Ke,K] = trussFEM2D.solve(k,b,EAs,BCs,loads);
        mu_g = -u( (2*OpKnoten-1):2*OpKnoten );
        dmu_g(1,:) = gradf(Ke, K, u, x, b, OpKnoten,2);
        dmu_g(2,:) = gradf(Ke, K, u, x, b, OpKnoten);

        sigma2_g = zeros(2,1);
        dsigma2_g = zeros(2,174);
        dx = 0.00001;
        for i=1:length(k)
            %Finite Differenzen
            c = zeros(length(k),2);
            c(i,1) = dx;
            u_p =  trussFEM2D.solve((k+c),b,EAs,BCs,loads); %dx muss noch in u plus dx eingebaut werden
            u_m =  trussFEM2D.solve((k-c),b,EAs,BCs,loads);
            u_p = -u_p(2*OpKnoten -1);
            u_m = -u_m(2*OpKnoten -1);
            dgdxi = (u_p-u_m)/(2*dx);
            %Finite Differenzen
            c = zeros(length(k),2);
            c(i,2) = dx;
            u_p =  trussFEM2D.solve((k+c),b,EAs,BCs,loads); %dx muss noch in u plus dx eingebaut werden
            u_m =  trussFEM2D.solve((k-c),b,EAs,BCs,loads);
            u_p = -u_p(2*OpKnoten);
            u_m = -u_m(2*OpKnoten);
            dgdyi = (u_p-u_m)/(2*dx);
        
        
            for j=1:length(k)
                %Finite Differenzen
                c = zeros(length(k),2);
                c(j,1) = dx;
                u_p =  trussFEM2D.solve((k+c),b,EAs,BCs,loads); %dx muss noch in u plus dx eingebaut werden
                u_m =  trussFEM2D.solve((k-c),b,EAs,BCs,loads);
                u_p = -u_p(2*OpKnoten -1);
                u_m = -u_m(2*OpKnoten -1);
                dgdxj = (u_p-u_m)/(2*dx);
                %Finite Differenzen
                c = zeros(length(k),2);
                c(j,2) = dx;
                u_p =  trussFEM2D.solve((k+c),b,EAs,BCs,loads); %dx muss noch in u plus dx eingebaut werden
                u_m =  trussFEM2D.solve((k-c),b,EAs,BCs,loads);
                u_p = -u_p(2*OpKnoten);
                u_m = -u_m(2*OpKnoten);
                dgdyj = (u_p-u_m)/(2*dx);
        
                sigma2_g = sigma2_g + [dgdxi*dgdxj*CovMa(i,j); dgdyi*dgdyj*CovMa(i,j)];
            end
        end

    end

    function [mu_g,sigma2_g,dmu_g,dsigma2_g] = MC_FEM_RDO(x, E, k, b, BCs, loads, OpKnoten,anzRealMC, CovSq, mu_stat)
        
        % deterministische Eingangsgrößen
        EAs = E*x;

        %rng default;
        g1 = zeros(anzRealMC,1);
        g2 = zeros(anzRealMC,1);
        dg1 = zeros(anzRealMC,174);
        dg2 = zeros(anzRealMC,174);
        dsigma2_g = zeros(2,174);

        % Schleife über Anzahl Realisierungen
        parfor i=1:anzRealMC

            % z = [k(1:5,:);
            % randn(length(k)-10,2);
            % k(51:55,:)];

            % Erzeugen von Zufallszahlen
            x_MC = CovSq * z + mu_stat;

            % Auswerten der Zielfunktion
            [u,~,Ke,K] = trussFEM2D.solve(x_MC,b,EAs,BCs,loads);
            u = -u;
            g1(i,1) = u(OpKnoten-1);
            g2(i,1) = u(OpKnoten);
            dg1(i,:) = gradf(Ke, K, u, x, b, OpKnoten,2);
            dg2(i,:) = gradf(Ke, K, u, x, b, OpKnoten);
            
            dsigma2_g = dsigma2_g + [g1(i,1)*dg1(i,:); g2(i,1)*dg2(i,:)];
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

        mu_g = [mu_g1; mu_g2];
        dmu_g = [dmu_g1; dmu_g2];
        sigma2_g = [sigma2_g1; sigma2_g2];
        dsigma2_g = 2/(anzRealMC-1) * dsigma2_g - [2*mu_g1*dmu_g1; 2*mu_g2*dmu_g2];

    end
end