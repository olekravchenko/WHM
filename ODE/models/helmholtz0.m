function Dy = helmholtz0(x, y, par)

Dy = y;

k0 = par(1);
alpha = par(2);

Dy(1) = y(2);
Dy(2) = - (k0^2 + alpha * y(1)^2) * y(1);