function Aufg2(anzRealMC)
anzRealMC = 10000;

addpath(fullfile(".", "Verteilungen"))

[n,m,rmin,rmax,dx,dy,r0,A,E,Fex,k,b,nob,EAs,rs,BCs,loads] = preprocess();

OpKnoten = (floor(n/2)*m) + 1;

mu_X = 0; %Erwartungwert von Fx in kN
sig_X = 15; %Standardabweichung von FX in kN
cov_X = sig_X^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Aufgabenteil a) Monte-Carlo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %Schleife über die Anzahl der Realisierungen
% for i=1:anzRealMC
% 
%     x(i) = sig_X * randn + mu_X; %Erzeugen einer normalverteilten Zufallszahl
%     loads = [(floor(n/2)+1)*m, 1, x(i);
%              (floor(n/2)+1)*m, 2, Fex]; %Zufallszahl in den Kraftvektor stecken
%     u = trussFEM2D.solve(k,b,EAs,BCs,loads); %Verschiebungsw
%     g(i,1:2) = u(2*OpKnoten-1:2*OpKnoten)';
% 
% 
% end
% 
% %%Auswertung
% mu_MCa = mean(g)
% var_MCa = var(g)
% sig_MCa = std(g)
% g = sort(g);
% g(round(0.01*anzRealMC));
% 
% x = g(:,1);
% EX = mu_MCa(1,1);
% VarX = var_MCa(1,1);
% 
% figure();
% histfit(x);
% title('Histogramm der Verschiebung von Punkt P in x-Richtung (MC)');
% xlabel('Verschiebung von Punkt P in x-Richtung[m]');
% ylabel('Anzahl');
% 
% [a,c] = EX_VarX_to_GumbelDistribution(EX, VarX);
% fh = @(x) cdf_GumbelDistribution(x,a,c);
% [acceptGumbelaX,dmax1X] = K_S_test_fh(x,fh);
% 
% % [mue,sigma] = EX_VarX_to_LogNormalDistribution(EX, VarX);
% % fh = @(x) cdf_LogNormalDistribution(x,mue,sigma);
% % [acceptLogNormala,dmax2X] = K_S_test_fh(x,fh)
% 
% fh = @(x) cdf_NormalDistribution(x,EX, VarX);
% [acceptNormala,dmax3X] = K_S_test_fh(x,fh);
% 
% au = min(x); bu = max(x);
% fh = @(x) cdf_UniformDistribution(x,au,bu);
% [acceptUniforma,dmax4X] = K_S_test_fh(x,fh);
% 
% x = g(:,2);
% EX = mu_MCa(1,2);
% VarX = var_MCa(1,2);
% 
% figure();
% histfit(x);
% title('Histogramm der Verschiebung von Punkt P in y-Richtung (MC)');
% xlabel('Verschiebung von Punkt P in y-Richtung[m]');
% ylabel('Anzahl');
% 
% [a,c] = EX_VarX_to_GumbelDistribution(EX, VarX);
% fh = @(x) cdf_GumbelDistribution(x,a,c);
% [acceptGumbelaX,dmax1X] = K_S_test_fh(x,fh);
% 
% % [mue,sigma] = EX_VarX_to_LogNormalDistribution(EX, VarX);
% % fh = @(x) cdf_LogNormalDistribution(x,mue,sigma);
% % [acceptLogNormala,dmax2X] = K_S_test_fh(x,fh)
% 
% fh = @(x) cdf_NormalDistribution(x,EX, VarX);
% [acceptNormala,dmax3X] = K_S_test_fh(x,fh);
% 
% au = min(x); bu = max(x);
% fh = @(x) cdf_UniformDistribution(x,au,bu);
% [acceptUniforma,dmax4X] = K_S_test_fh(x,fh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Aufgabenteil a) FOSM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loads = [(floor(n/2)+1)*m, 1, mu_X;
         (floor(n/2)+1)*m, 2, Fex];

[u,~,~,K] = trussFEM2D.solve(k,b,EAs,BCs,loads);
dudf = graduf(u,K,OpKnoten,n,m);

%Validieren mit Finite Differenzen
dx = 0.0000000001;

loads = [(floor(n/2)+1)*m, 1, mu_X + dx;
         (floor(n/2)+1)*m, 2, Fex];
u_pv = trussFEM2D.solve(k,b,EAs,BCs,loads);
u_p = u_pv(2*OpKnoten-1:2*OpKnoten);

loads = [(floor(n/2)+1)*m, 1, mu_X - dx;
         (floor(n/2)+1)*m, 2, Fex];
u_mv = trussFEM2D.solve(k,b,EAs,BCs,loads);
u_m = u_mv(2*OpKnoten-1:2*OpKnoten);

dudfini = (u_p-u_m)/(2*dx);

valid_grad = abs(dudf - dudfini) 

mu_FOSMa = u( (2*OpKnoten-1):2*OpKnoten )
var_FOSMa = [dudf(1)^2 * cov_X; dudf(2)^2 * cov_X]
sig_FOSMa = [sqrt(var_FOSMa(1)); sqrt(var_FOSMa(2))]


x_pdf = -1*mu_MCa(1,1) : 0.001*mu_MCa(1,1) : 3*mu_MCa(1,1); % Zur Darstellung des x Achsen Abschnittes

% Monte Carlo
p_MC = pdf('Normal', x_pdf, mu_MCa(1,1), sig_MCa(1,1));
p_FOSM = pdf('Normal', x_pdf, mu_FOSMa(1,1), sig_FOSMa(1,1));
figure();
plot(x_pdf, p_MC, 'r');
plot(x_pdf, p_FOSM, 'b');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Aufgabenteil b) Monte-Carlo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu_stat = k; %Erwartungswert entspricht den unverschobenen Knotenpositionen
sig_stat = 0.03; % 10% der kleinsten Stablänge
lc = 1.5; %Korrelationslänge

loads = [(floor(n/2)+1)*m, 2, Fex]; %Karftvektor zurücksetzen

[CovMa,CovSq] = Zufallsfeld(mu_stat, sig_stat, lc, k);

%Schleife über die Anzahl der Realisierungen
for i=1:anzRealMC
    z = [k(1:5,:);
        randn(length(k)-10,2);
        k(51:55,:)];
    x = CovSq * z + mu_stat; %Erzeugen einer Gauß-normalverteilten Zufallszahl
    u = trussFEM2D.solve(x,b,EAs,BCs,loads); %Verschiebung
    g(i,1) = u(2*OpKnoten-1);
    g(i,2) = u(2*OpKnoten);
    
end

%%Auswertung
mu_MCb = mean(g )
var_MCb = var(g )
sig_MCb = std(g )

g = sort(g );
g(round(0.01*anzRealMC))

figure();
histfit(g(:,1));
figure();
histfit(g(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Aufgabenteil b) FOSM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = trussFEM2D.solve(k,b,EAs,BCs,loads);
mu_FOSMb = u( (2*OpKnoten-1):2*OpKnoten );

var_FOSMb = [0 0];
dx = 0.00001;
for i=1:length(k)
    %Finite Differenzen
    c = zeros(length(k),2);
    c(i,1) = dx;
    u_p =  trussFEM2D.solve((k+c),b,EAs,BCs,loads); %dx muss noch in u plus dx eingebaut werden
    c(i,1) = -dx;
    u_m =  trussFEM2D.solve((k+c),b,EAs,BCs,loads);
    u_p = u_p(2*OpKnoten -1);
    u_m = u_m(2*OpKnoten -1);
    dgdxi = (u_p-u_m)/(2*dx);
    %Finite Differenzen
    c = zeros(length(k),2);
    c(i,2) = dx;
    u_p =  trussFEM2D.solve((k+c),b,EAs,BCs,loads); %dx muss noch in u plus dx eingebaut werden
    c(i,2) = -dx;
    u_m =  trussFEM2D.solve((k+c),b,EAs,BCs,loads);
    u_p = u_p(2*OpKnoten);
    u_m = u_m(2*OpKnoten);
    dgdyi = (u_p-u_m)/(2*dx);


    for j=1:length(k)
        %Finite Differenzen
        c = zeros(length(k),2);
        c(j,1) = dx;
        u_p =  trussFEM2D.solve((k+c),b,EAs,BCs,loads); %dx muss noch in u plus dx eingebaut werden
        c(j,1) = -dx;
        u_m =  trussFEM2D.solve((k+c),b,EAs,BCs,loads);
        u_p = u_p(2*OpKnoten -1);
        u_m = u_m(2*OpKnoten -1);
        dgdxj = (u_p-u_m)/(2*dx);
        %Finite Differenzen
        c = zeros(length(k),2);
        c(j,2) = dx;
        u_p =  trussFEM2D.solve((k+c),b,EAs,BCs,loads); %dx muss noch in u plus dx eingebaut werden
        c(j,2) = -dx;
        u_m =  trussFEM2D.solve((k+c),b,EAs,BCs,loads);
        u_p = u_p(2*OpKnoten);
        u_m = u_m(2*OpKnoten);
        dgdyj = (u_p-u_m)/(2*dx);

        var_FOSMb = var_FOSMb + [dgdxi*dgdxj*CovMa(i,j) dgdyi*dgdyj*CovMa(i,j)];
    end
end

sig_FOSMb = [sqrt(var_FOSMb(1)); sqrt(var_FOSMb(2))]
end