function [out,misc] = posterior_theta_unif( theta, para )
% HB likelihood with a uniform prior model


priorHB = @(x)  unifpdf( x(:,1), para.psi(:,1), para.psi(:,1) + para.psi(:,2))...
             .* unifpdf( x(:,2), para.psi(:,3), para.psi(:,3) + para.psi(:,4))...
             .* unifpdf( x(:,3), para.psi(:,5), para.psi(:,5) + para.psi(:,6))...
             .* unifpdf( x(:,4), para.psi(:,7), para.psi(:,7) + para.psi(:,8));

         

    
p = priorHB( theta );

[ ll, ~] = loglike(theta, para);

out = para.lnEv + log(para.Ns_th) - log(para.Ns_psi) + ll + log(sum( p./para.denom ));


if isnan(out)
    out = -1e10;
    misc = 2;

else
    misc = 0;
end
    
    
