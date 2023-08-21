function Aufg2(anzRealMC)

%% Aufgabenteil a) Monte-Carlo
% Parameters
mu_X = 0; %Erwartungswert von Kraft F in kN
sig_X = 15; %Standardabweichung in kN
cov_X = sig_X ^ 2; % Covarianz nur Skalar, da nur eine Eingangsgröße

OpKnoten = (floor(n/2)*m) + 1; % Startwert für die Verschiebung in Punkt P bzw. Knoten 26 

n=11; m=5; %n Number nodes Horizonal %m Number of Nodes Vertical 
dx = 0.4;  dy = 0.3; %Abstand der einzelnen Knoten in x und y Richtung
rmin = 0.0001; %min Querschnittsradius
rmax = 0.01; %max Querschnittsradius
r0 = 0.005; %Startraduis
A = r0^2*pi; %Querschnittsfläche
E = 210 *10^6; %E-Modul
Fex = -100; %Äußere Kraft
% create geometry & plot initial configuration
[k,b] = trussFEM2D.truss_preProcess2D(n,m,dx,dy); %k Positionen der Knoten %b Welche Knoten miteinander Verbunden sind
nob=size(b,1);
EAs = E*A*ones(nob,1);
rs = r0*ones(nob,1);
trussFEM2D.plotTruss2D(k,b,rs,1);
% Loads and Boundary conditons
BCs = [1,1; 1,2;   m,1; m,2];
ii=1;
for i=[1:m,n*m-m+1:n*m]
    BCs(ii,:) = [i,1];
    ii=ii+1;
    BCs(ii,:) = [i,2];
    ii=ii+1;
end

%Schleife über die Anzahl der Realisierungen
for i=1:anzRealMC
    x(i) = sig_X .* randn + mu_X; %Erzeugen einer normalverteilten Zufallszahl
    loads = [(floor(n/2)+1)*m, 1, x(i), (floor(n/2)+1)*m, 2, Fex]; %Zufallszahl in den Kraftvektor stecken
    u = trussFEM2D.solve(k,b,EAs,BCs,loads); %Verschiebungsw
    g(i,1) = u(OpKnoten-1);
    g(i,2) = u(OpKnoten);
    
end

%%Auswertung
mean(g )
var(g )
g = sort(g )
g(round(0.01*anzRealMC))

%% Aufgabenteil a) FOSM
loads = [(floor(n/2)+1)*m, 1, mu_X, (floor(n/2)+1)*m, 2, Fex];
u = trussFEM2D.solve(k,b,EAs,BCs,loads);
mu_g = u( (OpKnoten-1):OpKnoten );
loads = [(floor(n/2)+1)*m, 1, 1, (floor(n/2)+1)*m, 2, 0];
dudx = trussFEM2D.solve(k,b,EAs,BCs,loads);
sig_g = dudx( (OpKnoten-1):OpKnoten ) * cov_X;

%%Aufgabenteil b)
mu_stat = 2; %??
sig_stat = 0; %??
lc = 1.5; %Korrelationslänge

x = 0:0.4:4; %Länge 4m % 0.4m Knotenabstand

nx = length(x);

mu = ones(nx,1)*mu_stat;

R = ones(nx,nx);
for i = 1:nx
    for j=i+1:nx
        R(i,j) = exp( -(x(i)-x(j))^2 /lc^2);
        R(j,i) = R(i,j);
    end
end
CovMa = R*sig_stat^2;

CovSq = sqrtm(CovMa);

figure;hold on;
for i=1:5
    yr = CovSq * randn(nx,1) + mu;
    plot(x,yr,'Marker','o')
end


end