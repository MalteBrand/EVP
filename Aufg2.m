function Aufg2(anzRealMC)
anzRealMC = 10000;

addpath(fullfile(".", "Verteilungen"))

[n,m,rmin,rmax,dx,dy,r0,A,E,Fex,k,b,nob,EAs,rs,BCs,loads] = preprocess();

OpKnoten = (floor(n/2)*m) + 1;

%% Aufgabenteil a) Monte-Carlo

loads = [(floor(n/2)+1)*m, 1, 15;
        (floor(n/2)+1)*m, 2, Fex];
u = trussFEM2D.solve(k,b,EAs,BCs,loads); %Verschiebung
u = u(OpKnoten-1:OpKnoten);
%Schleife über die Anzahl der Realisierungen
for i=1:anzRealMC
    x(i) = sig_X .* randn + mu_X; %Erzeugen einer normalverteilten Zufallszahl
    loads = [(floor(n/2)+1)*m, 1, x(i);
             (floor(n/2)+1)*m, 2, Fex]; %Zufallszahl in den Kraftvektor stecken
    u = trussFEM2D.solve(k,b,EAs,BCs,loads); %Verschiebungsw
    g(i,1) = u(OpKnoten-1);
    g(i,2) = u(OpKnoten);
    
end

%%Auswertung
mu_MC = mean(g )
var_MC = var(g )
g = sort(g );
g(round(0.01*anzRealMC))

x = g;
EX = mu_MC;
VarX = var_MC;

[a,b] = EX_VarX_to_GumbelDistribution(EX, VarX);
fh = @(x) cdf_GumbelDistribution(x,a,b);
[accept,dmax] = K_S_test_fh(x,fh)

[mue,sigma] = EX_VarX_to_LogNormalDistribution(EX, VarX);
fh = @(x) cdf_LogNormalDistribution(x,mue,sigma);
[accept,dmax] = K_S_test_fh(x,fh)

fh = @(x) cdf_NormalDistribution(x,EX, VarX);
[accept,dmax] = K_S_test_fh(x,fh)

au = min(x); bu = max(x);
fh = @(x) cdf_UniformDistribution(x,au,bu);
[accept,dmax] = K_S_test_fh(x,fh)

%% Aufgabenteil a) FOSM
loads = [(floor(n/2)+1)*m, 1, mu_X;
         (floor(n/2)+1)*m, 2, Fex];
u = trussFEM2D.solve(k,b,EAs,BCs,loads);
mu_g = u( (OpKnoten-1):OpKnoten );
loads = [(floor(n/2)+1)*m, 1, 1];
dudx = trussFEM2D.solve(k,b,EAs,BCs,loads);
sig_g = (dudx( (OpKnoten-1):OpKnoten )).^2 * cov_X;

%%Aufgabenteil b) Monte-Carlo
mu = k; %Erwartungswert entspricht den unverschobenen Knotenpositionen
sig_stat = 0.03; % 10% der kleinsten Stablänge
lc = 1.5; %Korrelationslänge

loads = [(floor(n/2)+1)*m, 2, Fex]; %Karftvektor zurücksetzen

[CovMa,CovSq] = Zufallsfeld(mu_stat, sig_stat, lc, k);

%Schleife über die Anzahl der Realisierungen
for i=1:anzRealMC
    z = [k(1:5,:);
        randn(length(k)-10,2);
        k(51:55,:)];
    x = CovSq * z + mu; %Erzeugen einer normalverteilten Zufallszahl
    loads = [(floor(n/2)+1)*m, 2, Fex]; %Zufallszahl in den Kraftvektor stecken
    u = trussFEM2D.solve(x,b,EAs,BCs,loads); %Verschiebungsw
    g(i,1) = u(OpKnoten-1);
    g(i,2) = u(OpKnoten);
    
end

%%Auswertung
mean(g )
var(g )
g = sort(g );
g(round(0.01*anzRealMC))


%%Aufgabenteil b) FOSM
u = trussFEM2D.solve(k,b,EAs,BCs,loads);
mu_g = u( (OpKnoten-1):OpKnoten );

sig_g = [0 0];
dx = 0.00001;
for i=1:length(k)
    %Finite Differenzen
    c = zeros(length(k),2);
    c(i,1) = dx;
    u_p =  trussFEM2D.solve((k+c),b,EAs,BCs,loads); %dx muss noch in u plus dx eingebaut werden
    c(i,1) = -dx;
    u_m =  trussFEM2D.solve((k+c),b,EAs,BCs,loads);
    u_p = u_p(2 * OpKnoten -1);
    u_m = u_m(2 * OpKnoten -1);
    dgdxi = (u_p-u_m)/(2*dx);
    %Finite Differenzen
    c = zeros(length(k),2);
    c(i,2) = dx;
    u_p =  trussFEM2D.solve((k+c),b,EAs,BCs,loads); %dx muss noch in u plus dx eingebaut werden
    c(i,2) = -dx;
    u_m =  trussFEM2D.solve((k+c),b,EAs,BCs,loads);
    u_p = u_p(2 * OpKnoten);
    u_m = u_m(2 * OpKnoten);
    dgdyi = (u_p-u_m)/(2*dx);


    for j=1:length(k)
        %Finite Differenzen
        c = zeros(length(k),2);
        c(j,1) = dx;
        u_p =  trussFEM2D.solve((k+c),b,EAs,BCs,loads); %dx muss noch in u plus dx eingebaut werden
        c(j,1) = -dx;
        u_m =  trussFEM2D.solve((k+c),b,EAs,BCs,loads);
        u_p = u_p(2 * OpKnoten -1);
        u_m = u_m(2 * OpKnoten -1);
        dgdxj = (u_p-u_m)/(2*dx);
        %Finite Differenzen
        c = zeros(length(k),2);
        c(j,2) = dx;
        u_p =  trussFEM2D.solve((k+c),b,EAs,BCs,loads); %dx muss noch in u plus dx eingebaut werden
        c(j,2) = -dx;
        u_m =  trussFEM2D.solve((k+c),b,EAs,BCs,loads);
        u_p = u_p(2 * OpKnoten);
        u_m = u_m(2 * OpKnoten);
        dgdyj = (u_p-u_m)/(2*dx);

        sig_g = sig_g + [dgdxi*dgdxj*CovMa(i,j) dgdyi*dgdyj*CovMa(i,j)];
    end
end

end