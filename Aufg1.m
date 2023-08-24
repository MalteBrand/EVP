function Aufg1()
    
    % global diff_gradabl;
    % diff_gradabl = [0,0,0];

    %Initialparameter
    [n,m,rmin,rmax,dx,dy,r0,A,E,Fex,k,b,nob,EAs,rs,BCs,loads] = preprocess();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%a) Erstellen eines Referenzentwurf
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u0 = trussFEM2D.solve(k,b,EAs,BCs,loads);
    trussFEM2D.plotTruss2D(k,b,rs,2,u0,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%b) Optimieren des Fachwerks nach Stabquerschnitt
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    % ndiff = length(diff_gradabl);
    % tdiff = 1:1:ndiff;
    % figure(5);
    % grid on
    % hold on
    % plot(tdiff, diff_gradabl(:, 1), 'g')
    % plot(tdiff, diff_gradabl(:, 2), 'b')
    % plot(tdiff, diff_gradabl(:, 3), 'r--')
    % xlabel('Anzahl der Iterationen')
    % ylabel('Wert der Ableitung')
    % legend('Differenz','Explizite Ableitung','Finite Differenzen')
    % title('Ableitung in Stab 40')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%c) Optimieren des Fachwerks nach Stabradius
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    OpKnoten = (floor(n/2)*m) + 1; %Berechnungvorschrift für Punkt P bzw. Knoten 26 
    x0 = r0*ones(nob,1);
    V0 = trussFEM2D.mass(k,b,(pi*r0.^2));
    
    fhfc = @(x) zielc(x, E, k, b, BCs, loads, OpKnoten);
    fhgtbc = @(x) nebenBedingungen4toolboxc(x, rmin, rmax, V0, k, b, E, BCs, loads);
    
    %Optimierung mit SQP
    options = optimoptions('fmincon','Algorithm','sqp','Display','iter','StepTolerance',0.000000000001);
    x_opt = fmincon(fhfc,x0,[],[],[],[],[],[],fhgtbc,options);
    
    r_opt = x_opt;
    EAs = E*(pi*r_opt.^2);
    
    %Optimierte Verschiebung
    u = trussFEM2D.solve(k,b,EAs,BCs,loads);
    
    %Optimierter Plot ohne Last
    trussFEM2D.plotTruss2D(k,b,r_opt,6);
    %Optimierter Plot mit Last
    trussFEM2D.plotTruss2D(k,b,r_opt,7,u,1);

    function [f,df,ddf] = ziel(x, E, k, b, BCs, loads, OpKnoten)
    
        EAs = E*x;
        
        %Berechnung der Verschiebung
        [u,~,Ke,K] = trussFEM2D.solve(k,b,EAs,BCs,loads); 
        u = -u; %Verschiebung kann positiv werden
    
        %Es soll die y-Verschiebung in Knoten 26 optimiert werden
        f = u(2 * OpKnoten);
        df = gradf(Ke, K, u, x, b, OpKnoten);
        ddf = zeros(174,174);
    
        % dfini = gradfini(k, b, E, EAs, BCs, loads, OpKnoten);
        % valid = abs(dfini - df);
        % [ndiff,~] = size(diff_gradabl);
        % diff_gradabl(ndiff+1,:) = [valid(1,40) , df(1,40) , dfini(1,40)];
        
    end
    
    function [g,dg] = ungl_bed(x, rmin, rmax, k, b, E, BCs, loads)
        len = length(x);
        g(1:len)                = -(x - pi*rmin^2); %r > rmin
        g( (len+1):(2*len) )    = -(pi*rmax^2 - x); %r < rmax
        EAs = E*x;
        [u,~,Ke,K] = trussFEM2D.solve(k,b,EAs,BCs,loads);
        g(2*len+1) = u(2 * OpKnoten); % u < 0

        dg((1:len),:) = -1* eye(len);
        dg( (len+1):(2*len),: ) = eye(len);
        dg(2*len+1,:) = gradf(Ke, K, u, x, b, OpKnoten);

    end
    function [h, dh] = gl_bed(x,V0,k,b)   
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
    
    function [nb,nbeq] = nebenBedingungen4toolbox(x, rmin, rmax, V0, k, b, E, BCs, loads)
        nb = ungl_bed(x, rmin, rmax, k, b, E, BCs, loads);
        nbeq = gl_bed(x, V0, k, b);
    end

    function [f,df] = zielc(x, E, k, b, BCs, loads, OpKnoten)
    
        EAs = E*(pi*x.^2);
        
        %Berechnung der Verschiebung
        [u,~,Ke,K] = trussFEM2D.solve(k,b,EAs,BCs,loads); 
        u = -u; %Verschiebung kann positiv werden

        %Es soll die y-Verschiebung in Knoten 26 optimiert werden
        f = u(2 * OpKnoten);
        df = gradfr(Ke, K, u, x, b, OpKnoten);
    
        % dfini = gradfini(k, b, E, EAs, BCs, loads, OpKnoten);
        % valid = dfini(1,1:10) - df(1,1:10);

    end
    
    function [g,dg] = ungl_bedc(x, rmin, rmax, k, b, E, BCs, loads)
        len = length(x);

        g(1:len)                = -(x - rmin); %r > rmin
        g( (len+1):(2*len) )    = -(rmax - x); %r < rmax

        EAs = E*(pi*x.^2);
        [u,~,Ke,K] = trussFEM2D.solve(k,b,EAs,BCs,loads);
        g(2*len+1) = u(2 * OpKnoten); % u < 0

        dg = zeros((2*len + 1), len);

        dg((1:len),:) = -1* eye(len);
        dg( (len+1):(2*len),: ) = eye(len);
        dg(2*len+1,:) = gradfr(Ke, K, u, x, b, OpKnoten);
    end

    function [h, dh] = gl_bedc(x,V0,k,b)   
        len = length(x);

        h = (trussFEM2D.mass(k,b,(pi*x.^2)) - V0); % V = V0
        
        x1=k(b(:,1),1);
        y1=k(b(:,1),2);
        x2=k(b(:,2),1);
        y2=k(b(:,2),2);

        dh = zeros(1,174);

        for i=1:len
            dh(1,i) = 2*pi*x(i)*(sqrt((x2(i)-x1(i))^2 + (y2(i)-y1(i))^2));
        end

    end
    
    function [nb,nbeq] = nebenBedingungen4toolboxc(x, rmin, rmax, V0, k, b, E, BCs, loads)
        nb = ungl_bedc(x, rmin, rmax, k, b, E, BCs, loads);
        nbeq = gl_bedc(x, V0, k, b);
    end

end