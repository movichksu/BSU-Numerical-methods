clear;
clc;

#���������� �������� ������:
n = 7; N = n**2;
eps = 1.e-5;
A = gallery('poisson', n); # ������� ��������
E = eye(size(A)); # ��������� �������
f = rand(N, 1);
numb = 1000;

# ����� ����� �������� ������������� �������� � ��������������������:
function [Err] = commonMethod(B, A, f)
    numb = 1000;
    D = inv(B) * A;
    g = inv(B) * f;
    x = zeros(size(D)(1), 1);
    r = D * x - g; 
    pogreshnost = norm(r) / norm(g);
    k = 0; #������� 
    Err = zeros(numb, 1); #������ �� ������������

    while pogreshnost > eps && k < numb
        tau = (r' * r) / (( A * r)' * r);
        x = x - tau * r;
        r = D * x - g;
        pogreshnost = norm(r) / norm(g);
        k = k + 1;
        Err(k) = pogreshnost;
    endwhile 
endfunction


# ����� �����:
B = diag(diag(A));
err_yakobi = commonMethod(B, A, f);

# ����� �������
B = triu(A);
err_zeidelya = commonMethod(B, A, f);

# ����������� ����������� �����
R1 = triu(A) - diag(diag(A)) / 2;
R2 = tril(A) - diag(diag(A)) / 2;
B = (E + R1 / 2) * (E + R2 / 2);
err_triag = commonMethod(B, A, f);

x = zeros(numb ,1);
for i = 1:numb
    x(i) = i;
endfor

# ������ ������ 
semilogy(x, err_yakobi, x, err_zeidelya, x, err_triag);
xlim([0,80]);
legend ("Jacoby", "Zeidel", "triag");