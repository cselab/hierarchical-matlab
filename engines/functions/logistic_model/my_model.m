function y = my_model(x,theta,FLAG)


switch(FLAG)
    case 1
        % logistic function
        f = exp(theta(3)*x);
        y = ( theta(1)*theta(2)*f )./( theta(1) + theta(2)*(f-1) );

    case 2

        % Gompertz function
        lg = log(theta(2)/theta(1));
        y = theta(1) * exp(  lg*exp(-theta(3)*x) );
end
