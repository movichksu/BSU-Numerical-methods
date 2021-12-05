clear;
clc;

#���������� �������� ������:
a = randi([-9,9]);
b = a + randi([2, 9]);
v = 9;
n = 10;
eps = 0.001;
E = eye(n); #��������� �������
A = zeros(n); #������� n � n
f = zeros(n ,1); #������� ������� m � n

#��������� ������� �
for i = 1:n
    f(i) = a + (b - a) * rand;
    for j = 1:n
        if i == j
            A(i, j) = 100 + v;
        else 
            A(i, j) = 1 / (i + j + v);
        endif
    endfor
endfor

#����� �������� ��������� ���
tau = [1 / (2 * norm(A)), 1 / (4 * norm(A)), 1 / (8 * norm(A))]; 

#����������� ����������
pogreshnost = zeros(length(tau)); 

#�������� ������ ������ ������� �������� ������� ������� ���
for i = 1:length(tau)
    H = E - tau(i) * A;
    phi = tau(i) * f(i);
    q = norm(H);
    xs = phi;
    xn = H * xs + phi;
    k = 0; #������� ��������(�������)
    while q / (1-q) * norm(xs - xn) >= eps
        xs = xn;
        xn = H * xs + phi;
        k = k + 1;
        pogreshnost(i,k) = q / (1 - q) * norm(xs - xn) ;
    endwhile
endfor

x = zeros(length(pogreshnost),1);
for i = 1:length(pogreshnost)
    x(i) = i;
endfor


#������ ������ �������� ����������� �� ������ �������� ��� ��������� ����������
plot(x, pogreshnost(1, :), x, pogreshnost(2, :),x , pogreshnost(3, :)); 
legend ("tau=1/(2|A|)", "tau=1/(4|A|)", "tau=1/(8|A|)");

