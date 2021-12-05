clear;
clc;

delt = 10**(-6); #точность

a = randi([-10,10]);
b = a + randi([-10, 10]);
 #разреженна€ матрицa ѕуассона:
n = 10;
A = zeros(n);
f = zeros(n ,1);
N = n**2;
eps = 1.e-3;
A = gallery('poisson', n);
err = zeros(1, 3);
n = N;
for i = 1:n
    f(i) = a + (b - a) * rand;
end
B = diag(diag(A));
  #1 методом скорейшего спуска
x = zeros(n, 1);
r = f - A * x; #¬ычисл€ем вектор нев€зки
for i = 1:n
    w(i,1) = r(i,1) / B(i,i); #вычисл€ем вектор поправки
endfor
t = dot(w, r) / dot(A * w, w);
x = x + t * w; #¬ычисл€ем очередное приближение
k = 1;
err(1, k) = norm(r) / norm(f);
while norm(r) / norm(f) >= delt    #¬ычислени€ продолжаем до тех пор, 
                                    #пока относительна€ норма нев€зки
                                    #не превосходит заданную точность 

    k = k + 1;
    r = f - A * x;
    for i = 1:n
        w(i,1) = r(i,1) / B(i,i);
    endfor
    t = dot(w, r) / dot(A * w, w); #формула3
    x = x + t * w; #приближение
    err(1, k) = norm(r) / norm(f);
endwhile
 
#2 методом минимальных нев€зок.
x = zeros(n, 1);

r = f - A * x;
for i = 1:n
    w(i,1) = r(i,1) / B(i,i);
endfor
t = dot(A * w, r) / dot(A * w, A * w);
x = x + t * w;
k = 1;
err(2, k) = norm(r) / norm(f);
while norm(r) / norm(f) >= delt    
    k = k + 1;
    r = f - A * x;
    for i = 1:n
        w(i,1) = r(i,1) / B(i,i);
    endfor
    t = dot(A * w, r) / dot(A * w, A * w); #формула 4
    x = x + t * w;
    err(2, k) = norm(r) / norm(f);
endwhile

#3 методом минимальных поправок
x = zeros(n, 1);
r = f - A * x;
for i = 1:n
    w(i,1) = r(i,1) / B(i,i);
endfor
t = dot(A* w, w) / dot(inv(B) * A * w, A * w); #формула 5
x = x + t * w;
k = 1;
err(3, k) = norm(r) / norm(f);
while norm(r) / norm(f) >= delt    
    k = k + 1;
    r = f - A * x;
    for i = 1:n
        w(i,1) = r(i,1) / B(i,i);
    endfor
    t = dot(A * w, w) / dot(inv(B) * A * w, A * w);
    x = x + t * w;
    err(3, k) = norm(r) / norm(f);
endwhile

x = zeros(size(err)(2) ,1);
for i = 1:size(x)
    x(i) = i;
endfor
#графики убывани€ относительной нормы нев€зки
plot(x, err(1, :));
 # отличие методов в разном подчете t