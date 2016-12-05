function dy = rigid(t,y)
u = sin(t);
dy = zeros(3,1);    % a column vector
dy(1) = y(2);
dy(2) = y(3);
dy(3) = (-cos(y(1))^2 * y(3) - y(2) - exp(-y(1))*u);
end