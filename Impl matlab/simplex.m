clear all

%Carga de datos del problema y definición de variables
A = importdata('a.asc');
c = importdata('c.asc');
b = importdata('b.asc');

[m,n] = size(A);

%Fase I del simplex

%Definición de los parámetros para la fase I
AI = [A eye(m)];
cI = [zeros(1, n) ones(1, m)];
vnI = 1:n;
vbI = (n+1):(n+m);
xbI = b;
zI = cI(vbI)*b; 
iout= 0;
bland = true;
niter=1;

if bland == true
fprintf('Inicio del simplex primal con regla de Bland \n');
end

if bland == false
fprintf('Inicio del simplex primal con costos reducidos más negativos \n');
end

fprintf('Fase I \n');

while (iout == 0)
	fprintf("Iteración número: %4.d ", niter)
	[vbI, vnI, xbI, zI, iout] = simplexP_iter(cI, AI, b, vbI, vnI, xbI, zI, bland);
	niter = niter + 1;
end

if iout == 2
fprintf('Problema ilimitado');
end

if iout == 3
fprintf('La solución básica factible es degenerada');
end

%Si no nos encontramos en ninguno de esos dos casos, iout = 1 y tenemos una solución ópitma
%Fase II del simplex

if abs(zI) > 1e-5
	fprintf('Problema infactible');
end

if abs(zI) < 1e-5
	fprintf('\n Fase II');
	%Actualización de variables para comenzar la primera iteracion de la fase II
	iout = 0;
	xb = xbI;
	vb = vbI;
	vn = [];
	%Eliminación de variables artificiales de la fase I
	for i=1:n
		if vnI(i) <= n
			vn = [vn vnI(i)];
		end
	end
	z = c(vb)*xb;	
	while (iout == 0)
		fprintf("Iteración número: %4.d ", niter);
		[vb, vn, xb, z, iout] = simplexP_iter(c, A, b, vb, vn, xb, z, bland);
		niter = niter + 1;
	end
	if iout == 1
		fprintf("\n Solución óptima encontrada, iteración %4.d",niter-1);
		fprintf("z = %.4f",z );
	end
	if iout == 2
		fprintf('Problema ilimitado');
	end

	if iout == 3
		fprintf('La solución básica factible es degenerada');
	end
end




