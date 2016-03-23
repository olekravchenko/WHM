function Dy = model(x, y, Beta2)

Dy = y(:);

Dy(1) = y(2);
Dy(2) = y(1) * (Beta2 - 2 * y(2));