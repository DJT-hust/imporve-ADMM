clear
close all
clc

tic

%% 初始化
data

for t = 2:T + 1 %为了时间统一角标

    F_n = (E_n(:, :, t - 1) - E_n(:, :, t - 1)') / 2;
    % F_n_avg=sum(F_n)
    E_n_k = zeros(N, Nir + 1);
    E_nm_k = zeros(N, N, Nir + 1);

    %% 原问题

    for n = 1:N

        Deltamu_t = zeros(N, Nir);
        Deltadelta_t = zeros(N, Nir);
        M_set = find(Connection_mat(:, n));
        N_omega_n = size(M_set, 2); %n的连接数

        %平均值
        F_n_avg = sum(F_n(n, :)) / N_omega_n;
        price_avg = sum(price_tensor(n, :, t - 1)) / N_omega_n;
        E_n_avg = sum(E_n(n, :, t - 1)) / N_omega_n;
        V_n = (price_avg + rou * F_n_avg + eta * E_n_avg - b(n, t)) / (a(n, t) + (rou + eta) / N_omega_n);

        for k = 1:Nir
            W = zeros(N_omega_n, 1);
            E_n_k(n, k + 1) = clip(V_n + Deltadelta_t(n, k) / (a(n, t) + (rou + eta) / N_omega_n), En_underline(n, t), En_overline(n, t)); %(4.19)

            for index = 1:N_omega_n
                W(index) = (price(n, M(index), t - 1) + rou * F_n(n, M_set(index)) + eta * E_n(n, M_set(index), t - 1) - b(n, t)) / (rou + eta);

                if ismember(n, Consumer)
                    E_nm_k(n, M_set(index), k + 1) = clip(W(index) + (Deltamu_t(n, k) - a(n, t) * E_n_k(n, k + 1)) / (rou + eta), 0, E_n_k(n, k + 1)); %(4.21)
                else
                    E_nm_k(n, M_set(index), k + 1) = clip(W(index) + (Deltamu_t(n, k) - a(n, t) * E_n_k(n, k + 1)) / (rou + eta), E_n_k(n, k + 1), 0); %(4.22)
                end

            end

        end

        Deltamu_t(n, k + 1) = (a(n, t) + (rou + eta) / N_omega_n) * (E_n_k(n, k + 1) - V_n) - Deltadelta_t(n, k);
        Deltadelta_t(n, k + 1) = (rou + eta) * (sum(E_nm_k(n, :, k + 1)) / N_omega_n - mean(W)) - Deltamu_t(n, k + 1) + a(n, t) * E_n_k(n, k + 1);

    end

    Deltamu(:, t - 1) = Deltamu_t(:, Nir + 1);
    Deltadelta(:, t - 1) = Deltadelta_t(:, Nir + 1);

    %% 价格更新
    price_tensor(:, :, t) = price_tensor(:, :, t - 1) - rou / 2 * (E_nm_k(:, :, Nir + 1) + E_nm_k(:, :, Nir + 1)');

end

toc
