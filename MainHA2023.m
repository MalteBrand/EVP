%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Entwurfsoptimierung und probabilistische Verfahren in der Strukturmechanik
% Studenten: Malte Brand (52314),  Amelie Gildhoff (21599139)
% Abgabedatum: 27.08.2023 (SoSe 23)
% Prüfer: Prof. Kriegesmann
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc, close all;
anzRealMC = 1000;  % Anzahl Realisierungen für Monte Carlo Simulation

Aufg1()

Aufg2(anzRealMC)

%Am besten parallel-processing starten (parpool() in die Command-line)
Aufg3(anzRealMC)
