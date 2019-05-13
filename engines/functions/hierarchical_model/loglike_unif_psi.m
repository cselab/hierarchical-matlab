function [out,misc] = loglike_unif_psi( theta, para )
% HB likelihood with a uniform prior model

out = 0;

priorHB = @(x)  unifpdf( x(:,1), theta(1), theta(1) + theta(2) )...
             .* unifpdf( x(:,2), theta(3), theta(3) + theta(4) )...
             .* unifpdf( x(:,3), theta(5), theta(5) + theta(6) )...
             .* unifpdf( x(:,4), theta(7), theta(7) + theta(8) ) ;

         
for i = 1:length(para.lnEv)
    
    p = priorHB( para.theta{i} );

    out = out + para.lnEv(i) - log(para.Ns(i)) + log(sum( p./exp(para.pri{i}) ));

end


if( isfield(para,'CMA') && para.CMA)
    out = -out;
end

if isnan(out)
    out = -1e10;
    
    if( isfield(para,'CMA') && para.CMA)
        out = 1e10;
    end
    
    misc = 2;

else
    misc = 0;
end
    
    
