function Aufg2(anzRealMC)

addpath(fullfile(".", "Aufg2\Verteilungen"))
addpath(fullfile(".", "Aufg1"))
addpath(fullfile(".", "Aufg2"))

[n,m,rmin,rmax,dx,dy,r0,A,E,Fex,k,b,nob,EAs,rs,BCs,loads] = preprocess();

OpKnoten = (floor(n/2)*m) + 1;

mu_X = 0; %Erwartungwert von Fx in kN
sig_X = 15; %Standardabweichung von FX in kN
cov_X = sig_X^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Aufgabenteil a) Monte-Carlo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Schleife über die Anzahl der Realisierungen
for i=1:anzRealMC

    x(i) = sig_X * randn + mu_X; %Erzeugen einer normalverteilten Zufallszahl
    loads = [(floor(n/2)+1)*m, 1, x(i);
             (floor(n/2)+1)*m, 2, Fex]; %Zufallszahl in den Kraftvektor stecken
    u = trussFEM2D.solve(k,b,EAs,BCs,loads); %Verschiebungsw
    g(i,1:2) = u(2*OpKnoten-1:2*OpKnoten)';


end

%%Auswertung
mu_MCa = mean(g)
var_MCa = var(g)
sig_MCa = std(g)
g = sort(g);
g(round(0.01*anzRealMC));

x = g(:,1);
EX = mu_MCa(1,1);
VarX = var_MCa(1,1);

figure(7)
histfit(x);
title('Histogramm der Verschiebung von Punkt P in X-Richtung a)');
xlabel('Verschiebung von Punkt P in X-Richtung[m]');
ylabel('Anzahl');

[a,c] = EX_VarX_to_GumbelDistribution(EX, VarX);
fh = @(x) cdf_GumbelDistribution(x,a,c);
[acceptGumbelaX,dmax1X] = K_S_test_fh(x,fh);

% [mue,sigma] = EX_VarX_to_LogNormalDistribution(EX, VarX);
% fh = @(x) cdf_LogNormalDistribution(x,mue,sigma);
% [acceptLogNormala,dmax2X] = K_S_test_fh(x,fh)

fh = @(x) cdf_NormalDistribution(x,EX, VarX);
[acceptNormalaX,dmax3X] = K_S_test_fh(x,fh);

au = min(x); bu = max(x);
fh = @(x) cdf_UniformDistribution(x,au,bu);
[acceptUniformaX,dmax4X] = K_S_test_fh(x,fh);

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

y_grada = [dudf(1,1), dudfini(1,1), valid_grad(1,1); dudf(2,1), dudfini(2,1), valid_grad(2,1)] 

fig = uifigure;
uitable(fig,'Data', y_grada, 'ColumnName', {'Explizit', 'Finite Differenzen', 'Differenz'});
title('Validierung Explizit/Finite Differenzen 2a)');

x_pdf = -4*sig_MCa(1,1) : 0.0001*sig_MCa(1,1) : 4*sig_MCa(1,1); % Zur Darstellung des x Achsen Abschnittes
%x_pdf_FOSMa = -1*mu_FOSMa(1,1) : 0.001*mu_FOSMa(1,1) : 3*mu_FOSMa(1,1);
p_MC = pdf('Normal', x_pdf, mu_MCa(1,1), sig_MCa(1,1));
p_FOSM = pdf('Normal', x_pdf, mu_FOSMa(1,1), sig_FOSMa(1,1));
figure(9)
hold on
plot(x_pdf, p_MC, 'b');
plot(x_pdf, p_FOSM, 'r--');
title('Wahrscheinlichkeitsdichtefunktionen Punkt P in X-Richtung a)');
legend('Monte-Carlo','FOSM')
xlabel('Verschiebung in X-Richtung');
ylabel('Wahrscheinlichkeit');

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
    z = [zeros(5,2);
         randn(length(k)-10,2);
         zeros(5,2)];
    %z = randn(length(k),2);
    x = CovSq * z + mu_stat; %Erzeugen einer Gauß-normalverteilten Zufallszahl
    % Keine Verschiebung der Lagerpunkte
    x = [k(1:5,:);
         x(6:50,:);
         k(51:55,:)];
    u = trussFEM2D.solve(x,b,EAs,BCs,loads); %Verschiebung
    g1(i) = u(2*OpKnoten-1);
    g2(i) = u(2*OpKnoten);
    
end

%%Auswertung
mu_g1 = mean(g1);
mu_g2 = mean(g2);
mu_MCb = [mu_g1;mu_g2]
var_g1 = var(g1);
var_g2 = var(g2);
var_MCb = [var_g1;var_g2]
sig_g1 = std(g1);
sig_g2 = std(g2);
sig_MCb = [sig_g1;sig_g2]

x = g1;
EX = mu_g1;
VarX = var_g1;

figure(10)
histfit(x);
title('Histogramm der Verschiebung von Punkt P in X-Richtung b)');
xlabel('Verschiebung von Punkt P in X-Richtung[m]');
ylabel('Anzahl');

[a,c] = EX_VarX_to_GumbelDistribution(EX, VarX);
fh = @(x) cdf_GumbelDistribution(x,a,c);
[acceptGumbelbX,dmax1X] = K_S_test_fh(x,fh);

% [mue,sigma] = EX_VarX_to_LogNormalDistribution(EX, VarX);
% fh = @(x) cdf_LogNormalDistribution(x,mue,sigma);
% [acceptLogNormala,dmax2X] = K_S_test_fh(x,fh)

fh = @(x) cdf_NormalDistribution(x,EX, VarX);
[acceptNormalbX,dmax3X] = K_S_test_fh(x,fh);

au = min(x); bu = max(x);
fh = @(x) cdf_UniformDistribution(x,au,bu);
[acceptUniformbX,dmax4X] = K_S_test_fh(x,fh);

x = g2;
EX = mu_g2;
VarX = var_g2;

figure(11)
histfit(x);
title('Histogramm der Verschiebung von Punkt P in Y-Richtung b)');
xlabel('Verschiebung von Punkt P in Y-Richtung[m]');
ylabel('Anzahl');

[a,c] = EX_VarX_to_GumbelDistribution(EX, VarX);
fh = @(x) cdf_GumbelDistribution(x,a,c);
[acceptGumbelbY,dmax1X] = K_S_test_fh(x,fh);

% [mue,sigma] = EX_VarX_to_LogNormalDistribution(EX, VarX);
% fh = @(x) cdf_LogNormalDistribution(x,mue,sigma);
% [acceptLogNormala,dmax2X] = K_S_test_fh(x,fh)

fh = @(x) cdf_NormalDistribution(x,EX, VarX);
[acceptNormalbY,dmax3X] = K_S_test_fh(x,fh);

au = min(x); bu = max(x);
fh = @(x) cdf_UniformDistribution(x,au,bu);
[acceptUniformbY,dmax4X] = K_S_test_fh(x,fh);


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
    c = [zeros(5,2);
         c(6:50,:);
         zeros(5,2)];
    u_p =  trussFEM2D.solve((k+c),b,EAs,BCs,loads); %dx muss noch in u plus dx eingebaut werden
    c(i,1) = -dx;
    c = [zeros(5,2);
         c(6:50,:);
         zeros(5,2)];
    u_m =  trussFEM2D.solve((k+c),b,EAs,BCs,loads);
    u_p = u_p(2*OpKnoten -1);
    u_m = u_m(2*OpKnoten -1);
    dgdxi = (u_p-u_m)/(2*dx);
    %Finite Differenzen
    c = zeros(length(k),2);
    c(i,2) = dx;
    c = [zeros(5,2);
         c(6:50,:);
         zeros(5,2)];
    u_p =  trussFEM2D.solve((k+c),b,EAs,BCs,loads); %dx muss noch in u plus dx eingebaut werden
    c(i,2) = -dx;
    c = [zeros(5,2);
         c(6:50,:);
         zeros(5,2)];
    u_m =  trussFEM2D.solve((k+c),b,EAs,BCs,loads);
    u_p = u_p(2*OpKnoten);
    u_m = u_m(2*OpKnoten);
    dgdyi = (u_p-u_m)/(2*dx);


    for j=1:length(k)
        %Finite Differenzen
        c = zeros(length(k),2);
        c(j,1) = dx;
        c = [zeros(5,2);
             c(6:50,:);
             zeros(5,2)];
        u_p =  trussFEM2D.solve((k+c),b,EAs,BCs,loads); %dx muss noch in u plus dx eingebaut werden
        c(j,1) = -dx;
        c = [zeros(5,2);
             c(6:50,:);
             zeros(5,2)];
        u_m =  trussFEM2D.solve((k+c),b,EAs,BCs,loads);
        u_p = u_p(2*OpKnoten -1);
        u_m = u_m(2*OpKnoten -1);
        dgdxj = (u_p-u_m)/(2*dx);
        %Finite Differenzen
        c = zeros(length(k),2);
        c(j,2) = dx;
        c = [zeros(5,2);
             c(6:50,:);
             zeros(5,2)];
        u_p =  trussFEM2D.solve((k+c),b,EAs,BCs,loads); %dx muss noch in u plus dx eingebaut werden
        c(j,2) = -dx;
        c = [zeros(5,2);
             c(6:50,:);
             zeros(5,2)];
        u_m =  trussFEM2D.solve((k+c),b,EAs,BCs,loads);
        u_p = u_p(2*OpKnoten);
        u_m = u_m(2*OpKnoten);
        dgdyj = (u_p-u_m)/(2*dx);

        var_FOSMb = var_FOSMb + [dgdxi*dgdxj*CovMa(i,j) dgdyi*dgdyj*CovMa(i,j)];
    end
end

sig_FOSMb = [sqrt(var_FOSMb(1)); sqrt(var_FOSMb(2))]

x_pdfb = (-4*sig_MCb(1,1) + mu_MCb(1,1)) : 0.0001*sig_MCb(1,1) : (4*sig_MCb(1,1) + mu_MCb(1,1)); % Zur Darstellung des x Achsen Abschnittes
%x_pdf_FOSMa = -1*mu_FOSMa(1,1) : 0.001*mu_FOSMa(1,1) : 3*mu_FOSMa(1,1);
p_MCb = pdf('Normal', x_pdfb, mu_MCb(1,1), sig_MCb(1,1));
p_FOSMb = pdf('Normal', x_pdfb, mu_FOSMb(1,1), sig_FOSMb(1,1));
figure(12)
hold on
plot(x_pdfb, p_MCb, 'b');
plot(x_pdfb, p_FOSMb, 'r--');
title('Wahrscheinlichkeitsdichtefunktionen Punkt P in X-Richtung b)');
legend('Monte-Carlo','FOSM')
xlabel('Verschiebung in X-Richtung');
ylabel('Wahrscheinlichkeit');

x_pdfb = (-4*sig_MCb(2,1) + mu_MCb(2,1)) : 0.0001*sig_MCb(2,1) : (4*sig_MCb(2,1) + mu_MCb(2,1)); % Zur Darstellung des x Achsen Abschnittes
%x_pdf_FOSMa = -1*mu_FOSMa(1,1) : 0.001*mu_FOSMa(1,1) : 3*mu_FOSMa(1,1);
p_MCb = pdf('Normal', x_pdfb, mu_MCb(2,1), sig_MCb(2,1));
p_FOSMb = pdf('Normal', x_pdfb, mu_FOSMb(2,1), sig_FOSMb(2,1));
figure(13)
hold on
plot(x_pdfb, p_MCb, 'b');
plot(x_pdfb, p_FOSMb, 'r--');
title('Wahrscheinlichkeitsdichtefunktionen Punkt P in Y-Richtung b)');
legend('Monte-Carlo','FOSM')
xlabel('Verschiebung in Y-Richtung');
ylabel('Wahrscheinlichkeit');

end