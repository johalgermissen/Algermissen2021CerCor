function [cope,varcope,tstat,F,Fdof] = ols(data,des,tc,Ftests)
% [COPE,VARCOPE,TSTAT,F,Fdof]=ols(DATA,DES,TC,Ftests)
% DATA IS T x V (time x # outcome variables)
% DES IS T x EV (design matrix)
% TC IS NCONTRASTS x EV  (contrast matrix): Use eye() for unspecific contrasts
% Ftests IS CELL ARRAY, NFTESTS*1; each cell contains vector that indexes contrast matrix
%
% see ols_examples.m
%
% TB 2004/ LH 2016/ JA 2019

if nargin < 3
    tc = eye(size(des,2));
end

if(size(data,1)~=size(des,1)) % different # time points
    error('OLS::DATA and DES have different number of time points');
elseif(size(des,2)~=size(tc,2)) % different # EVs
    error('OLS:: DES and TC have different number of evs')
end

pdes    = pinv(des); % invert design matrix
prevar  = diag(tc*pdes*pdes'*tc'); % contrast x design matrix, make square matrix
R       = eye(size(des,1)) - des*pdes; % error over pseudo-inverse
tR      = trace(R); % variances
pe      = pdes*data; % parameter estimates (betas)
cope    = tc*pe; % specify contrasts
if(nargout>1) % if asked for standard errors:
    res     = data-des*pe; % residual error
    sigsq   = sum(res.*res/tR); % mean squared error
    varcope = prevar*sigsq; % 
    if(nargout>2) % if asked for T-values:
        tstat = cope./sqrt(varcope); % T-values
    end
    
    if (nargout>3) % if asked for F-test:
        %see appendix of Henson and Penny, 2005
        for i = 1:size(Ftests,1)
            corth = eye(size(des,2)) - tc(Ftests{i},:)'*pinv(tc(Ftests{i},:)'); 
            % orthogonal contrast to C:
            desorth = des*corth; %design matrix of reduced model
            residR = eye(size(des,1)) - des*pinv(des); %residual forming matrix of full model
            residO = eye(size(des,1)) - desorth*pinv(desorth); %residual forming matrix of reduced model
            projectM = residO-residR; %projection matrix M due to X1
            for j = 1:size(data,2)
                F(i,j) = ((pe(:,j)'*des'*projectM*des*pe(:,j))./(data(:,j)'*residR*data(:,j)))*...
                    (size(des,1)-rank(des))/rank(projectM*des); %F-statistic
            end
            % degrees of freedom for this F-statistic:
            Fdof(i,1) = rank(projectM*des);
            Fdof(i,2) = (size(des,1)-rank(des));
        end
    end
end

end % end of function.