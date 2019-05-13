function [out,misc] = loglike_norm_psi(theta,para)
% HB likelihood with a Gaussian prior model
%  out - output
%  misc - flag for different type of error
%       = 0 (normal), = 1 (SIG non pos. def.), = 2 (lik. calc. error)

theta = theta(:);
Nd = 4;

Svec = theta(Nd+1:2*Nd).^2;
detS = prod(Svec);


if detS<=0
    out = -1e10;
    misc = 1;
    
else
    out = 0;
    for i = 1:length(para.lnEv)
        MU = para.theta{i} - repmat(theta(1:Nd)',para.Ns(i),1);
        S  = repmat(Svec',para.Ns(i),1);
        
        % log of sum of log values (avoid underflow)
        % for varying prior
 
        tmp = -0.5*Nd*log(2*pi) - 0.5*log(detS) - 0.5*sum( MU.*((MU./S)), 2) - para.pri{i};
        % for constant prior
%         tmp = -Nd/2*log(2*pi) - lndSIGd2 - sum(MU.*(SIG\MU')',2)/2 - para.pri;
        
        tmp0 = max(tmp);
        
        out = out + para.lnEv(i) - log(para.Ns(i)) + log(sum(exp(tmp - tmp0))) + tmp0;
        
    end

    if isnan(out)
        out = -1e10;
        misc = 2;

    else
        misc = 0;
    end

end

end