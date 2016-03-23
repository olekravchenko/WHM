function Dy = model0(x, y)

Dy = y(:);

Dy(1) = y(2);
Dy(2) = - 2 * y(1) * y(2);