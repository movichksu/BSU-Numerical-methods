clear;
clc;

#����� ��������� ������
n = 10;
v = 9;
for i = 1:n
  a = v + 10 / (10 + i);
  b = v + 10 / i;
  c = v + (1 + i) / 10;
endfor

#1.	�������� ������� A
A = zeros(n);
for i = 1:n
  A(i, i) = a;
  
  if i <= n - 1
    A(i, i + 1) = b;
    A(i + 1, i) = b;
  endif
  
  if i <= n - 2
    A(i, i + 2) = c;
    A(i + 2, i) = c;
  endif
endfor

A

copy_A = A;
S = eye(n);

#2.	���������� ������� ���������a.
#�������  �������� ������� �
for i = 1:n - 1
  M_inv = eye(n);
  if A(n - i + 1, n - i) != 0
    for j = 1:n
        M_inv(n - i, j) = A(n - i + 1, j);
    endfor
  endif
  
  #����� ������� ���� ������ ������������ ������� ������� A,
  #��������� ������ ������� ����� ��� ������ ����������.
  A = M_inv * A * inv(M_inv);
  S = S * inv(M_inv);
endfor

M_inv
A

#3.	���������� ����������� �������� ������� F 
#(��� �� ����������� �������� ������� A).
eq = (-1)**n * [1 (-1 * A(1, :))];

eig_values = roots(eq);
eig_values

#4. ���������� ������������ �������
y = zeros(n, n);
eig_vectors = zeros(n, n);

for i = 1:n
  for j = 1:n
    y(j, i) = eig_values(i) ** (n - j);
  endfor
  eig_vectors(:, i) = S * y(:, i);
endfor

y


eig_vectors
#5.	�������� ����������� ��������. 
display("Check eigenvalues: ");

[vectors, D] = eigs(copy_A, n);
values = diag(D);
values
vectors

#6.	��������, ��� ��� ����������� ������� ������� A,
# ���������� � ���� ������������, ������������� �������� �� ������������
# ���������: �.�. ��������� ����� �������� 
for i = 1:n
  display(norm(copy_A * eig_vectors(:, i) - eye(n) * eig_values(i) * eig_vectors(:, i)));
endfor 

