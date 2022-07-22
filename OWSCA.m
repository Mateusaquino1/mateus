% Optimal Wavelet Signal Compression - OWSCA

clear; close all; clc

global wname level 

%% Carregando sinal e definindo vari�vis e espa�o de busca

load T.txt
load SIG.txt

sig = SIG;

% Verificando o n�mero m�ximo de decomposi��es poss�vel
% lmax = wmaxlev(length(SIG),wtype);
% lmax = log2(length(SIG));         
         
% Definindo o tipo de wavelet e o n�mero de n�veis da decomposi��o
wname = 'bior1.3';  % Composi��o das caracter�sticas dos filtros da DWT
level = 12;

% Determina��o do Espa�o de busca, i.e., dos limiares inferior e superior 
% para a vari�vel Threshold

x0 =  [...
    44.321453044773079;...
    35.765000000013970;...
    26.449329150258563;...
    10.372500000026776;...
    3.901461665169336;...
    3.829999999958091;...
    2.525255092128646;...
    1.514375000144355;...
    15.607614427921362;...
    21.460000000195578;...
    36.350370883010328;...
    0.000000000000000...
    ];  

wLL = x0*0.005; % Limite inferior
wUU = x0*100; % Limite superior

%% Otimiza��o dos par�metros (limares) para compress�o

Nvars = 12; % N�mero de vari�veis -> igual a level

% ------------------------- Algoritmo gen�tico -------------------------- %

% Nger = 20; % N�mero de Gera��es
% Tolfun = 10^(-12); % Termination tolerance on fitness function value
% Tolcon = 10^(-12); % Termination tolerance on constraints
% 
% options = gaoptimset('MutationFcn',@mutationadaptfeasible,'PlotFcns',...
%           {@gaplotbestf,@gaplotbestindiv,@gaplotexpectation,@gaplotstopping},...
%           'Generations',Nger,'TolFun',Tolfun,'TolCon',Tolcon);
% 
% % Determina��o dos limiares �timos para compress�o     
% [x,f,exitflag] = ga(@(thr)objfun(SIG,T,thr,wUU,wLL),Nvars,[],[],[],[],...
%                           wLL,wUU,[],options);

% ----------------------- Enxame de Part�culas -------------------------- %                     
Nswarm = 100; % Tamanho do enxame  
options = optimoptions('particleswarm','SwarmSize',Nswarm,'Display','iter');
[x,fval,exitflag] = particleswarm(@(thr)objfun(SIG,T,thr,wUU,wLL),Nvars,wLL,wUU,options);

% -------------------- Refinamento da otimiza��o ------------------------ %
% optNLP = optimset('Algorithm','interior-point','GradObj','off','GradConstr','off',...
%                   'DerivativeCheck','off','Display','iter','TolX',10^(-10),...
%                   'TolFun',10^(-10),'TolCon',10^(-10),'MaxFunEval', 5e10);                    
%                       
[opt_thres,fval] = fmincon(@(thr)objfun(SIG,T,thr,wUU,wLL),x,[],[],[],[],wLL, wUU,[]);

figure(1)
bar(opt_thres);
xlabel('Wavelet threshold')
ylabel('Threshold value')

%% Decomposi��o do sinal pela DWT e compress�o de acordo com os limiares
% encontrados no processo de otimiza��o

% Regra de compress�o: coeficientes com valores menores que os limiares s�o
% descartados, ou seja, se d(j,k) < opt_thres(k), d(j,k) = 0

thr = opt_thres;
sorh = 'h';  % hard threshold (� um tipo de compress�o)

% Decomposi��o DWT
sigDEN = cmddenoise(SIG,wname,level,sorh,NaN,thr);

% a6 = sigDEN;
a6 = round(sigDEN);
Tr=[];
Er=[];
i=1;
k = 1;

while i<=(length(T)-1)
	j=i;
	I=[];
    bol =1;
    while bol
        if j<(length(T)-1)
        I=[I;j];
		j=j+1;
            bol = (a6(j)==a6(i));
        else
            bol=0;
        end
    
    end
    if isempty(I)
        break
    else
	I=round(mean(I));
    end
	i=j;
    if i==6000
        i=j;
    end
    
	Tr=[Tr;T(I)];
	Er=[Er;sigDEN(I)];
end

figure(2)
plot(T,SIG)
hold on
stairs(Tr,Er,'-r','LineWidth',2) % stairs(X,Y) plots the elements in Y at the locations specified by X
title('Energy vector')
xlabel('Interactions')
ylabel('Energy')
legend('Original','Compressed')

%% C�lculo de figuras de m�rito

% Interpolando sinal comprimido 
Er_int = interp1(1:numel(Er),Er,linspace(1,numel(Er),numel(SIG)))';

% Calculando coeficiente de Correla��o entre os sinais
C = abs(corr2(sig,Er_int)); % Quanto mais pr�ximo de 1, melhor a wavelet

% Root Mean Squared Error
rmse = sqrt(mean((SIG - Er_int).^2))/100; % Quanto mais pr�ximo de 0, melhor a wavelet