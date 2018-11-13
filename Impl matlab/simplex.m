clear all

%Carga de datos del problema y definición de variables
A = load();
c = load();
b = load();

[m,n] = size(A);

%Fase I del simplex

%Definición de los parámetros para la fase I
AI = [A eye(m)];
cI = [zeros(1, n) ones(1, n)];
vnI = 1:n;
vbI = (n+1):(n+m);
xbI = b;
iout= 0;
bland = true;
niter=1;

while (iout == 0)
	fprintf("Iteración número: %4.d", niter)
	[vb2, vn2, xb2, z2, iout] = simplexP_iter(c, A, b, vb, vn, xb, z, bland);
	niter++;
end

if iout == 2
end

if iout == 3
end

%Si no nos encontramos en ninguno de esos dos casos, iout = 1 y tenemos una solución ópitma
%Fase II del simplex

if z == 0
	iout = 0;
	z = c(vb)'*xb;
	while (iout == 0)
	
	end

end




