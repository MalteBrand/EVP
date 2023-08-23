function [n,m,rmin,rmax,dx,dy,r0,A,E,Fex,k,b,nob,EAs,rs,BCs,loads] = preprocess()

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
    %trussFEM2D.plotTruss2D(k,b,rs,1);
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

end