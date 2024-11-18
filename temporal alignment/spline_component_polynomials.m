function b = spline_component_polynomials(xaxis,n_ctrl_pts)
%b-spline interpolation of control points is found by
%f*b where f ix n_ctrl_pts x 1 function value


    x0 = min(xaxis);
    delta = (max(xaxis)-min(xaxis))/(n_ctrl_pts - 1);
    normx = (xaxis - x0)/delta;
    b = zeros(n_ctrl_pts, length(xaxis));

    for qq = (0):(n_ctrl_pts + 1)
        i = min(max(1,qq), n_ctrl_pts);
        s = normx - qq + 1;
        b(i,:) = b(i,:) + chanCubicPolynomial(s);
    end

end


function b = chanCubicPolynomial(s)
    b = (2+s).^3/6 .*(s >= -2 & s < -1) + (2/3 - s.^2 - s.^3/2).*(s >= -1 & s < 0) + (2/3 - s.^2 + s.^3/2).*(s >= 0 & s < 1) + (2-s).^3/6 .*(s >= 1 & s < 2);
end