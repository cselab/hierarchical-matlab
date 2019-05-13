function [out,misc] = loglike(theta, para)



data = para.data;

Np = data.Np;

std_data = theta(end);

sse = sum( ( data.y - data.modelfun(data.x,theta(1:Np)) ).^2 );

out = -0.5*data.Nd*log(2*pi) - data.Nd*log(std_data) - 0.5*sse/(std_data.^2) ;



if( isfield(para,'CMA') && para.CMA)
    out = -out;
end


misc.test = 1;

end