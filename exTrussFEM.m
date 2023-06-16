%% Berechnung einer Referenzverschiebung
% Parameters
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
loads = [(floor(n/2)+1)*m, 2, Fex];
% compute deformation Ku=F and plot deformed configuration
u0 = trussFEM2D.solve(k,b,EAs,BCs,loads);
trussFEM2D.plotTruss2D(k,b,rs,2,u0);

%% Startwert für die Verschiebung in Punkt P bzw. Knoten 26 
OpKnoten = 26;
x0 = r0^2*pi*ones(nob,1);

fhf = @(x, E, k, b, BCs, loads) ziel(x, E, k, b, BCs, loads);
fhg = @(x, rmin, rmax) ungl_bed(x, rmin, rmax);
fhh = @(x, V0) gl_bed(x, V0);

%Optimierung mit SQP

function [f,df,ddf] = ziel(x, E, k, b, BCs, loads)
    
    %Parameter
    
    %Berechnung der Verschiebung in Knoten 26 für gegebene Querschnitte
    EAs = E*x;
    u = trussFEM2D.solve(k,b,EAs,BCs,loads);
    f = u(2 * OpKnoten);
    
    %
    dx = 0.1;
    u_p = trussFEM2D.solve(k,b,EAs,BCs,loads); %dx muss noch in u plus dx eingebaut werden
    u_m = trussFEM2D.solve(k,b,EAs,BCs,loads);
    u_p = u_p(2 * OpKnoten);
    u_m = u_m(2 * OpKnoten);
    
    %Berechnung der Ableitungen über finite Differenzen
    df = -(u_p-u_m)/(2*dx);
    ddf = -(u_p-2*u+u_m)/dx^2; 
end

function g = ungl_bed(x, rmin, rmax)
    n = length(x);
    g(1:n) = x - pi*rmin^2; 
    g(n+1:end) = pi*rmax^2 - x;
end
function h = gl_bed(x, V0)
    % Stablängen l muss hier noch berechnet werden
    
    h = sum(x*l) - V0; 
end 