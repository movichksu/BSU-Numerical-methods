clear;
clc;

# ���������� ����������� ������
n = 10;
v = 9;
eps = 0.5 * 10**(-3);
for i = 1:n
  a = v + 10 / (10 + i);
  b = v + 10 / i;
  c = v + (1 + i) / 10;
endfor

# ����� ������� A
z = ones(n, 1);
y = ones(n, 1);
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
  z(i) = rand();
  y(i) = rand();
endfor
display(A);

# ��������� �����
y = A * z;
lambda = y ./ z;
max_lambda = max(lambda);
min_lambda = min(lambda);
while abs(max_lambda - min_lambda) > eps  # �������� ����������� ����� �������� ������ ��������
  z = y / norm(y);                        # 1) n-������ ������������� ������
  y = A * z;                              # 2)
  lambda = y ./ z;                        #��������� ��������� ��������� ��������
  max_lambda = max(lambda);               #������� ������������ � ����������� ��������� ������
  min_lambda = min(lambda);
endwhile
x1 = z;


# ����� ��������� ������������
y = A * z;
lambda_n = dot(y, z) / dot(z, z);
lambda_s = 0;
while abs(lambda_s - lambda_n) > eps      # �������� ����������� ����� �������� ������ ��������
  z = y / norm(y);                        # 1) n-������ ������������� ������ 
  y = A * z;                              # 2) ���������� ��������� ������ ����������� � max �� ������ ������������ �������
  lambda_s = lambda_n;
  lambda_n = dot(y, z) / dot(z, z);
endwhile
x2 = z;

# ������� ����������
max(eigs(A, n)) 
[v,D] = eig(A)
v

display(max_lambda);
display(x1);
display(lambda_n);
display(x2);

