clear all

%Carga de datos del problema y definición de variables
A = importdata('a2.asc');
c = importdata('c2.asc');
b = importdata('b2.asc');

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
bland = false;
niter=1;

if bland == true
display("Inicio del simplex primal con regla de Bland");
end

if bland == false
display("Inicio del simplex primal con costos reducidos más negativos");
end

display("Fase I");

while (iout == 0)
	fprintf("Iteración número: %4.d", niter)
	[vb, vn, xb, z, iout] = simplexP_iter(cI, AI, b, vbI, vnI, xbI, zI, bland);
	niter = niter + 1;
end

if iout == 2
display("Problema ilimitado");
end

if iout == 3
display("La solución básica factible es degenerada");
end

%Si no nos encontramos en ninguno de esos dos casos, iout = 1 y tenemos una solución ópitma
%Fase II del simplex

if z ~= 0
	display("Problema infactible");
end

if z == 0
	display("Fase II");
	%Actualización de variables para comenzar la primera iteracion de la fase II
	iout = 0;
	xb = xbI;
	vb = vbI;
	vn = []
	%Eliminación de variables artificiales de la fase I
	for i=1:n
		if vn2(i) <= n
			vn = [vn vn2(i)]
		end
	end
	z = c(vb)'*xb;	
	while (iout == 0)
		fprintf("Iteración número: %4.d", niter);
		[vb, vn, xb, z, iout] = simplexP_iter(c, A, b, vb, vn, xb, z, bland);
		niter = niter + 1;
	end
	if iout == 1
		fprintf("Solución óptima encontrada, iteración %d",niter-1);
		fprintf("z = %4.f",z );
	end
	if iout == 2
		display("Problema ilimitado");
	end

	if iout == 3
		display("La solución básica factible es degenerada");
	end
end




