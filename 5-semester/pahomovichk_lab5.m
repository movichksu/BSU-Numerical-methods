clear;
clc;

delt = 10**(-6); #��������

a = randi([-10,10]);
b = a + randi([-10, 10]);
 #����������� ������a ��������:
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
  #1 ������� ���������� ������
x = zeros(n, 1);
r = f - A * x; #��������� ������ �������
for i = 1:n
    w(i,1) = r(i,1) / B(i,i); #��������� ������ ��������
endfor
t = dot(w, r) / dot(A * w, w);
x = x + t * w; #��������� ��������� �����������
k = 1;
err(1, k) = norm(r) / norm(f);
while norm(r) / norm(f) >= delt    #���������� ���������� �� ��� ���, 
                                    #���� ������������� ����� �������
                                    #�� ����������� �������� �������� 

    k = k + 1;
    r = f - A * x;
    for i = 1:n
        w(i,1) = r(i,1) / B(i,i);
    endfor
    t = dot(w, r) / dot(A * w, w); #�������3
    x = x + t * w; #�����������
    err(1, k) = norm(r) / norm(f);
endwhile
 
#2 ������� ����������� �������.
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
    t = dot(A * w, r) / dot(A * w, A * w); #������� 4
    x = x + t * w;
    err(2, k) = norm(r) / norm(f);
endwhile

#3 ������� ����������� ��������
x = zeros(n, 1);
r = f - A * x;
for i = 1:n
    w(i,1) = r(i,1) / B(i,i);
endfor
t = dot(A* w, w) / dot(inv(B) * A * w, A * w); #������� 5
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
#������� �������� ������������� ����� �������
plot(x, err(1, :));
 # ������� ������� � ������ ������� t