%Función que hace una iteración del simplex
%Tanto p como q son las posiciones en vb y vn

function [vb2, vn2, xb2, z2, iout] = simplexP_iter(c, A, b, vb, vn, xb, z, bland)

%Informacion de tamaños
[m,n] = size(A);

%Cálculos auxiliares para costos reducidos
cn = c(vn);
cb = c(vb);
An = A(:, vn);
B = A(:, vb);
Bi = inv(B);

%Cálculo del vector de costos reducidos
r = cn - cb*Bi*An;
if (min(r) >= 0) iout = 1; return; end %SBF óptima

%Elección de variable no básica q entrante

if (bland == true)
[indice, q] = max(vn);
	for i = 1:n-m
		if(r(i) < 0)
			if vn(i) < vn(q)
			q = i;
			end
		end
	end

end

%Si se utiliza el costo reducido más negativo,
%me quedo con el elemento más pequeño del vector de costos reducidos.

if (bland == false)
[costo, q] = min(r);
end

%Cálculo de la DBF de descenso
Aq = A(:, vn(q));
db = -Bi*Aq;
if (min(db) >= 0) iout = 2; return; end %Problema ilimitado

%Cálculo del paso y selección de la variable básica de salida p
for i = 1:m
	pas = [];
	if db(i) < 0
		theta2 = -xb(i)/db(i);
		pas = [pas theta2];
		theta = min(pas);
		if theta2 <= theta
			p = i;
			end
	end
end 

if (theta == 0) iout = 3; return; end %SBF degenerada

%Actualización y cambios en las bases
aux = vb(p);
vb(p) = vn(q);
vn(q) = aux;
xb = xb + theta*db;
xb(p) = theta;	%ahora en el lugar p se encuentra el valor de x_q
z = z + theta*r(q);
iout = 0;

%printf("Iteracion 1: variable de salida %d, variable de entrada %d", a, b)  

