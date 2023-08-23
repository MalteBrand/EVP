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
trussFEM2D.plotTruss2D(k,b,rs,2,u0,1);

%% Aufgabe 1 b)
OpKnoten = (floor(n/2)*m) + 1; % Startwert für die Verschiebung in Punkt P bzw. Knoten 26 
A0 = r0^2*pi*ones(nob,1); %Anfangsquerschnitte
rmin = 0.0001; %min Querschnittsradius
rmax = 0.01; %max Querschnittsradius
x0 = A0;
V0 = 0;
len_b = length(b);
l_stb = zeros(len_b,1);
for i=1:len_b
    l_stb(i) = sqrt( ( k(b(i,1),1)-k(b(i,2),1) )^2 + ( k(b(i,1),2)-k(b(i,2),2) )^2 );
    V0 = V0 + l_stb(i)*x0(i);
end

fhf = @(x) ziel(x, E, k, b, BCs, loads, OpKnoten);
fhgtb = @(x) nebenBedingungen4toolbox(x, rmin, rmax, l_stb, V0);

%Optimierung mit SQP
options = optimoptions('fmincon','Algorithm','sqp','Display','iter');
x_opt = fmincon(fhf,x0,[],[],[],[],[],[],fhgtb,options);

A = x_opt;
EAs = E*A;
u = trussFEM2D.solve(k,b,EAs,BCs,loads);

rs = sqrt(A(1:len_b)/pi);

trussFEM2D.plotTruss2D(k,b,rs,3);
trussFEM2D.plotTruss2D(k,b,rs,4,u,1);

function [f,df,ddf] = ziel(x, E, k, b, BCs, loads, OpKnoten)
    len = length(x);
    A = x;

    EAs = E*A;
    
    %Berechnung der Verschiebung
    [u,force,Ke,K] = trussFEM2D.solve(k,b,EAs,BCs,loads); %Da die Verschiebung auch positv werden kann, wird der Betrag der Verschiebung minimiert    
    u = -u;
    %u = u.^2;
    %Es soll die y-Verschiebung in Knoten 26 optimiert werden
    f = u(2 * OpKnoten);
    df = gradf(Ke, K, u, x, b, OpKnoten);
    ddf = zeros(174,174);
    
    % %Berechnung der Ableitungen über finite Differenzen
    % dx = 0.0000001;
    % for i=1:length(b)
    %     z = zeros(len,1);
    %     z(i) = E*dx;
    %     u_pv =  trussFEM2D.solve(k,b,EAs + z,BCs,loads); %dx muss noch in u plus dx eingebaut werden
    %     u_mv =  trussFEM2D.solve(k,b,EAs - z,BCs,loads);
    %     u_p = u_pv(2 * OpKnoten);
    %     u_m = u_mv(2 * OpKnoten);
    %     dfini(:,i) = -(u_p-u_m)/(2*dx);
    %     ddfini(:,i) = (u_pv-2*u+u_mv)/dx^2;
    % end
    % %Berechnung der Ableitungen über finite Differenzen
    % valid = dfini(1,1:10) - df(1,1:10);

end

function g = ungl_bed(x, rmin, rmax)
    len = length(x);
    g(1:len)                = x - pi*rmin^2; 
    g( (len+1):(2*len) )    = pi*rmax^2 - x;
    %g(2*len+1) = %%%%%Fortsetzen u<=0
end
function h = gl_bed(x, l_stb, V0)
    
    V = sum(l_stb .* x);
    
    h = V - V0; 
end

function [nb,nbeq] = nebenBedingungen4toolbox(x, rmin, rmax, l_stb, V0)
    nb = ungl_bed(x, rmin, rmax);
    nbeq = gl_bed(x, l_stb, V0);
end