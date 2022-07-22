function f = objfun(SIG,T,dec,wUU,wLL)

global wname level 

thro = dec(1:12);
nP = 100;
w = 0.5;

thr = thro;
% Denoising parameters.
%----------------------
% meth = 'sqtwolog';
% scal_or_alfa = one;
sorh = 'h';    % Specified soft or hard thresholding
% thr = thro.*(wUU - wLL)  ;
% thrSettings =  [...
%     52.990582182123035 ; ...
%     48.074999999997090 ; ...
%     45.986689514451427 ; ...
%     73.385000000009313 ; ...
%     90.276322754070861 ; ...
%     137.981250000011640 ; ...
%     176.686539182090200 ; ...
%     0.000000000000000 ; ...
%     0.000000000000000 ; ...
%     0.000000000000000   ...
%     ];

% Denoise using CMDDENOISE.
%--------------------------


sigDEN = cmddenoise(SIG,wname,level,sorh,NaN,thr);


a6 = round(sigDEN);
Tr=[];
Er=[];
Eo=[];
i=1;

k = 1;
% for i=2:length(a6)
%     
%     grad(i-1) = ( a6(i) - a6(i-1) ) / a6(i) ;
%     
%     if abs(grad(i-1)) < tol;
%        Tr(k) =  T(i);
%        Er(k) = a6(i);
%        Eo(k) = SIG(i);
%        k = k+1;
%     end
%     
% end

% 
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
    Eo=[Eo;SIG(I)];
end

% length(Er)



c = length(Er) - nP;

f = w*( norm(Er - Eo) ) + (1-w)*abs(c)^2;