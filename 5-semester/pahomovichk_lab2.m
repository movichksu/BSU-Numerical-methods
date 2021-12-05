clc;#очищает командное окно

v = 9;
a = randi([-9,9]);
b = a + randi([2, 9]);

list = [4, 6, 8, 10, 12]; # 
nev = zeros(size(list)(2)); #dim_sz = размер ( a , dim )
                            #Вернуть вектор-строку с размером (количеством элементов) каждого измерения для объекта a .
                            #Если задан второй аргумент dim , вернуть размер соответствующего измерения.            

#Метод квадратного корня
function [x] = SquareRoot(A, f)
    
    S = zeros(size(A));
    D = zeros(size(A));
    n = size(A)(1);
    x = zeros(n, 1);
    y = zeros(n, 1);
    
    D(1, 1) = sign(A(1, 1));
    S(1, 1) = sqrt(A(1, 1));

    
    for i = 2:n #получаем диагональные элеименты
        D(i, i) = sign(A(i, i) - sum((S(1:i - 1, i).^2)' * D(1:i - 1, 1:i - 1))); #(5)
        S(i, i) = sqrt(abs(A(i, i) - sum((S(1:i - 1, i).^2)' * D(1:i - 1, 1:i - 1)))); #(6)
        
        for j = i + 1:n #(7)
            buf_sum = 0;
            for k = 1:i - 1
                buf_sum = buf_sum + S(k, i) * D(k, k) * S(k, j);
            endfor
            S(i, j) = (A(i, j) - buf_sum) / (S(i, i) * D(i, i)); #(7)
        endfor
    endfor
     
     #В обратном ходе, решая систему S*Dy=f  с нижней треугольной матрицей S*D, 
     #найдем вектор y:
    y(1) = f(1) / (S(1, 1) * D(1, 1)); #(8)
    for i = 2:n
        buf_sum = 0;
        for j = 1:i - 1
            buf_sum = buf_sum + S(j, i) * y(j) * D(j, j);
        endfor
        y(i) = (f(i) - buf_sum) / (S(i, i) * D(i, i));
    endfor
    
    #Решая систему Sx=y с верхней треугольной матрицей S, вычислим искомое решение x
    x(n) = y(n) / S(n, n); #(9)
    for i = n - 1:-1:1
        buf_sum = 0;
        for j = i + 1:n
            buf_sum = buf_sum + S(i, j) * x(j);
        endfor    
        x(i) = (y(i) - buf_sum) / S(i, i);
    endfor
endfunction 

for z = 1:size(list)(2)
    n = list(z);
    A = zeros(n);
    f = zeros(n, 1);
    for i = 1:n
        f(i) = a + (b - a) * rand;
        for j = 1:n
            if (i == j)
                A(i, j) = 1 / (i + j + v);
            endif
        endfor
    endfor
    
    x = SquareRoot(A, f)
    A \ f
    
    nev(z) = norm(f - A * x) / norm(f);
    nev(z)
endfor

function Plot(nev, n)
    plot(n, nev, '-');
    xlabel("matrix dimension");
    ylabel("norma nevyazki");
    title("Plot norma nevyazki acc. to matrix dimension"); 
end
Plot(nev, list)