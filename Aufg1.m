function x_optb =  Aufg1()

    addpath(fullfile(".", "Aufg1"))

    global diff_gradabl;
    diff_gradabl = [0,0,0];

    %Initialparameter
    [n,m,rmin,rmax,dx,dy,r0,A,E,Fex,k,b,nob,EAs,rs,BCs,loads] = preprocess();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%a) Erstellen eines Referenzentwurf
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u0 = trussFEM2D.solve(k,b,EAs,BCs,loads);
    %trussFEM2D.plotTruss2D(k,b,rs,1);
    trussFEM2D.plotTruss2D(k,b,rs,1,u0,1);
    title('Referenz')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%b) Optimieren des Fachwerks nach Stabquerschnitt
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    OpKnoten = (floor(n/2)*m) + 1; %Berechnungvorschrift für Punkt P bzw. Knoten 26 
    x0 = A*ones(nob,1);
    V0 = trussFEM2D.mass(k,b,A);
    
    fhf = @(x) ziel(x, E, k, b, BCs, loads, OpKnoten);
    fhgtb = @(x) nebenBedingungen4toolbox(x, rmin, rmax, V0, k, b, E, BCs, loads);
    
    %Optimierung mit SQP
    options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'Algorithm','sqp','Display','iter','StepTolerance',0.000000000001);
    %options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'CheckGradients',true,'Algorithm','sqp','Display','iter','StepTolerance',0.000000000001);
    x_optb = fmincon(fhf,x0,[],[],[],[],[],[],fhgtb,options);
    
    A = x_optb;
    EAs = E*A;
    
    %Optimierte Verschiebung
    u1 = trussFEM2D.solve(k,b,EAs,BCs,loads);
    
    rs = sqrt(A/pi);
    
    %Optimierter Plot ohne Last
    %trussFEM2D.plotTruss2D(k,b,rs,3);
    %Optimierter Plot mit Last
    trussFEM2D.plotTruss2D(k,b,rs,2,u1,1);
    title('Optimierung nach Querschnitt')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%c) Optimieren des Fachwerks nach Stabradius
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    OpKnoten = (floor(n/2)*m) + 1; %Berechnungvorschrift für Punkt P bzw. Knoten 26 
    x0 = r0*ones(nob,1);
    V0 = trussFEM2D.mass(k,b,(pi*r0.^2));

    fhfc = @(x) zielc(x, E, k, b, BCs, loads, OpKnoten);
    fhgtbc = @(x) nebenBedingungen4toolboxc(x, rmin, rmax, V0, k, b, E, BCs, loads);

    %Optimierung mit SQP
    options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'CheckGradients',false,'Algorithm','sqp','Display','iter','StepTolerance',0.0000000000001);
    %options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'CheckGradients',true,'Algorithm','sqp','Display','iter','StepTolerance',0.000000000001);
    x_opt = fmincon(fhfc,x0,[],[],[],[],[],[],fhgtbc,options);

    r_opt = x_opt;
    EAs = E*(pi*r_opt.^2);

    %Optimierte Verschiebung
    u2 = trussFEM2D.solve(k,b,EAs,BCs,loads);

    %Optimierter Plot ohne Last
    %trussFEM2D.plotTruss2D(k,b,r_opt,5);
    %Optimierter Plot mit Last
    trussFEM2D.plotTruss2D(k,b,r_opt,3,u2,1);
    title('Optimierung nach Radius')
    
    x = 0:0.4:4;
    u0p = [u0(2),u0(12),u0(22),u0(32),u0(42),u0(52),u0(62),u0(72),u0(82),u0(92),u0(102)];
    u1p = [u1(2),u1(12),u1(22),u1(32),u1(42),u1(52),u1(62),u1(72),u1(82),u1(92),u1(102)];
    u2p = [u2(2),u2(12),u2(22),u2(32),u2(42),u2(52),u2(62),u2(72),u2(82),u2(92),u2(102)];
    figure(4);
    grid on
    hold on
    plot(x,u0p,'r')
    plot(x,u1p,'b')
    plot(x,u2p,'g')
    xlabel('X Positionen')
    ylabel('Verschiebung v')
    legend('Referenz','Optimiert Querschnitt b)','Optimiert Radius c)')
    title('Verschiebungen der untersten Knotenreihe')    

    %Plot für die Differenz
    ndiff = length(diff_gradabl);
    tdiff = 1:1:ndiff;
    figure(5);
    grid on
    hold on
    %plot(tdiff, diff_gradabl(:, 1), 'g')
    plot(tdiff, diff_gradabl(:, 2), 'b')
    plot(tdiff, diff_gradabl(:, 3), 'r--')
    xlabel('Anzahl der Iterationen')
    ylabel('Wert der Ableitung')
    legend('Differenz','Explizite Ableitung','Finite Differenzen')
    title('Ableitung in Stab 16')
    figure(6);
    grid on
    hold on
    plot(tdiff, diff_gradabl(:, 1), 'b')
    xlabel('Anzahl der Iterationen')
    ylabel('Differenz')
    %legend('Differenz','Explizite Ableitung','Finite Differenzen')
    title('Differenz zwischen Explizit/Finite Differenzen')

    function [f,df,ddf] = ziel(x, E, k, b, BCs, loads, OpKnoten)
    
        EAs = E*x;
        
        %Berechnung der Verschiebung
        [u,~,Ke,K] = trussFEM2D.solve(k,b,EAs,BCs,loads); 
        u = -u; %Verschiebung kann positiv werden
    
        %Es soll die y-Verschiebung in Knoten 26 optimiert werden
        f = u(2 * OpKnoten);
        df = gradf(Ke, K, u, x, b, OpKnoten);
        ddf = zeros(174,174);
    
        dfini = gradfini(k, b, E, EAs, BCs, loads, OpKnoten);
        valid = abs(dfini - df);
        [ndiff,~] = size(diff_gradabl);
        diff_gradabl(ndiff+1,:) = [valid(1,16) , df(1,16) , dfini(1,16)];
        
    end
    
    function [g,dg] = ungl_bed(x, rmin, rmax, k, b, E, BCs, loads)
        len = length(x);
        g((1:len))                = -(x - pi*rmin^2); %r > rmin
        g(((len+1):(2*len)))    = -(pi*rmax^2 - x); %r < rmax
        EAs = E*x;
        [u,~,Ke,K] = trussFEM2D.solve(k,b,EAs,BCs,loads);
        g((2*len+1)) = u(2 * OpKnoten); % u < 0

        dg(:,(1:len)) = -1* eye(len);
        dg(:, ((len+1):(2*len)) ) = eye(len);
        dg(:,(2*len+1)) = gradf(Ke, K, u, x, b, OpKnoten);

    end
    function [h, dh] = gl_bed(x,V0,k,b)   
        len = length(x);

        h = (trussFEM2D.mass(k,b,x) - V0); % V = V0 
        
        x1=k(b(:,1),1);
        y1=k(b(:,1),2);
        x2=k(b(:,2),1);
        y2=k(b(:,2),2);

        dh = zeros(174,1);

        for i=1:len
            dh(i,1) = (sqrt((x2(i)-x1(i))^2 + (y2(i)-y1(i))^2));
        end

    end
    
    function [nb,nbeq,dnb,dnbeq] = nebenBedingungen4toolbox(x, rmin, rmax, V0, k, b, E, BCs, loads)
        [nb,dnb] = ungl_bed(x, rmin, rmax, k, b, E, BCs, loads);
        [nbeq,dnbeq] = gl_bed(x, V0, k, b);
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

        dg = zeros(len, (2*len + 1));

        dg(:,(1:len)) = -1* eye(len);
        dg( :,(len+1):(2*len) ) = eye(len);
        dg(:,(2*len+1)) = gradfr(Ke, K, u, x, b, OpKnoten);
    end

    function [h, dh] = gl_bedc(x,V0,k,b)   
        len = length(x);

        h = (trussFEM2D.mass(k,b,(pi*x.^2)) - V0); % V = V0
        
        x1=k(b(:,1),1);
        y1=k(b(:,1),2);
        x2=k(b(:,2),1);
        y2=k(b(:,2),2);

        dh = zeros(174,1);

        for i=1:len
            dh(i,1) = 2*pi*x(i)*(sqrt((x2(i)-x1(i))^2 + (y2(i)-y1(i))^2));
        end

    end
    
    function [nb,nbeq,dnb,dnbeq] = nebenBedingungen4toolboxc(x, rmin, rmax, V0, k, b, E, BCs, loads)
        [nb,dnb] = ungl_bedc(x, rmin, rmax, k, b, E, BCs, loads);
        [nbeq,dnbeq] = gl_bedc(x, V0, k, b);
    end

end