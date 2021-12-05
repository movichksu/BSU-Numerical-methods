clear;
clc;

#Определяем исходные данные:
n = 7; N = n**2;
eps = 1.e-5;
A = gallery('poisson', n); # матрица Пуассона
E = eye(size(A)); # Единичная матрица
f = rand(N, 1);
numb = 1000;

# Общая схема неявного итерационного процесса с переобусловливателем:
function [Err] = commonMethod(B, A, f)
    numb = 1000;
    D = inv(B) * A;
    g = inv(B) * f;
    x = zeros(size(D)(1), 1);
    r = D * x - g; 
    pogreshnost = norm(r) / norm(g);
    k = 0; #счётчик 
    Err = zeros(numb, 1); #массив из погрешностей

    while pogreshnost > eps && k < numb
        tau = (r' * r) / (( A * r)' * r);
        x = x - tau * r;
        r = D * x - g;
        pogreshnost = norm(r) / norm(g);
        k = k + 1;
        Err(k) = pogreshnost;
    endwhile 
endfunction


# Метод Якоби:
B = diag(diag(A));
err_yakobi = commonMethod(B, A, f);

# метод Зейделя
B = triu(A);
err_zeidelya = commonMethod(B, A, f);

# попеременно треугольный метод
R1 = triu(A) - diag(diag(A)) / 2;
R2 = tril(A) - diag(diag(A)) / 2;
B = (E + R1 / 2) * (E + R2 / 2);
err_triag = commonMethod(B, A, f);

x = zeros(numb ,1);
for i = 1:numb
    x(i) = i;
endfor

# Строим график 
semilogy(x, err_yakobi, x, err_zeidelya, x, err_triag);
xlim([0,80]);
legend ("Jacoby", "Zeidel", "triag");