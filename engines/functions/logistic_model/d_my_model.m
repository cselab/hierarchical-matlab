function dy = d_my_model(x,theta,FLAG)

dy = zeros(3,1);

switch(FLAG)
    case 1
        % logistic function
        f = exp(theta(3)*x);
        dy(1) = (theta(2)^2*f*(f - 1))/(theta(1) - theta(2) + theta(2)*f)^2;
        dy(2) = (theta(1)^2*f)/(theta(1) - theta(2) + theta(2)*f)^2;
        dy(3) = (theta(1)*theta(2)*x*f*(theta(1) - theta(2)))/(theta(1) - theta(2) + theta(2)*f)^2;

    case 2

        % Gompertz function
        dy = [];
end
