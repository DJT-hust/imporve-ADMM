clear
close all
clc

tic

%% 初始化
data

for t = 2:T+1 %为了时间统一角标
    t
    F_nm = (E_nm(:, :, t - 1) - E_nm(:, :, t - 1)')/ 2;
    % F_n_avg=sum(F_n)
    E_n_k = 1e-5*ones(N, 2);
    E_nm_k = zeros(N, N, 2);

    %% 原问题

    for n = 1:N

        Deltamu_t = zeros(N, 1);
        Deltadelta_t = zeros(N, 1);
        M_set = find(Connection_mat(:, n));
        N_omega_n = size(M_set, 1); %n的连接数

        %平均值
        F_n_avg = sum(F_nm(n, :)) / N_omega_n;
        price_avg = sum(price_tensor(n, :, t - 1)) / N_omega_n;
        E_n_avg = sum(E_nm(n,:, t - 1)) / N_omega_n;
        V_n = (price_avg + rou * F_n_avg + eta * E_n_avg - b(n, t)) / (a(n, t) + (rou + eta) / N_omega_n);
        k=1;
        while abs(E_n_k(n,k)-sum(E_nm_k(n,:,k)))>1e-7
            W = zeros(N_omega_n, 1);
%             E_n_k(n, k + 1) = clip(V_n + Deltadelta_t(n, k) / (a(n, t) + (rou + eta) / N_omega_n), En_underline(n, t), En_overline(n, t)); %(4.19)
            if n==3
                    E_n_k(n, k + 1) = clip(V_n + Deltadelta_t(n, k) / (a(n, t) + (rou + eta) / N_omega_n), max(En_underline(n, t),E_n(n,t)-0.5),min( En_overline(n, t),E_n(n,t)+0.5)); %(4.19)
            elseif n==1
                    E_n_k(n, k + 1) = clip(V_n + Deltadelta_t(n, k) / (a(n, t) + (rou + eta) / N_omega_n), max(En_underline(n, t),-30-sum(E_n(n,2:t-1))), En_overline(n, t)); %(4.19)
            else
                    E_n_k(n, k + 1) = clip(V_n + Deltadelta_t(n, k) / (a(n, t) + (rou + eta) / N_omega_n), En_underline(n, t), En_overline(n, t)); %(4.19)
            end
            for index = 1:N_omega_n
                W(index) = (price_tensor(n, M_set(index), t - 1) + rou * F_nm(n, M_set(index)) + eta * E_nm(n, M_set(index), t - 1) - b(n, t)) / (rou + eta);

                if ismember(n, Producer)
                    E_nm_k(n, M_set(index), k + 1) = clip(W(index) + (Deltamu_t(n, k) - a(n, t) * E_n_k(n, k + 1)) / (rou + eta), 0, E_n_k(n, k + 1)); %(4.21)
                else
                    E_nm_k(n, M_set(index), k + 1) = clip(W(index) + (Deltamu_t(n, k) - a(n, t) * E_n_k(n, k + 1)) / (rou + eta), E_n_k(n, k + 1), 0); %(4.22)
                end

            end
            Deltamu_t(n, k + 1) = (a(n, t) + (rou + eta) / N_omega_n) * (E_n_k(n, k + 1) - V_n) - Deltadelta_t(n, k);
            Deltadelta_t(n, k + 1) = (rou + eta) * (sum(E_nm_k(n, :, k + 1)) / N_omega_n - mean(W)) - Deltamu_t(n, k + 1) + a(n, t) * E_n_k(n, k + 1);
            k=k+1;
        end
        % for k = 1:Nir
        %     W = zeros(N_omega_n, 1);
        %     E_n_k(n, k + 1) = clip(V_n + Deltadelta_t(n, k) / (a(n, t) + (rou + eta) / N_omega_n), En_underline(n, t), En_overline(n, t)); %(4.19)

        %     for index = 1:N_omega_n
        %         W(index) = (price_tensor(n, M_set(index), t - 1) + rou * F_nm(n, M_set(index)) + eta * E_nm(n, M_set(index), t - 1) - b(n, t)) / (rou + eta);

        %         if ismember(n, Consumer)
        %             E_nm_k(n, M_set(index), k + 1) = clip(W(index) + (Deltamu_t(n, k) - a(n, t) * E_n_k(n, k + 1)) / (rou + eta), 0, E_n_k(n, k + 1)); %(4.21)
        %         else
        %             E_nm_k(n, M_set(index), k + 1) = clip(W(index) + (Deltamu_t(n, k) - a(n, t) * E_n_k(n, k + 1)) / (rou + eta), E_n_k(n, k + 1), 0); %(4.22)
        %         end

        %     end

        % end
        

        % Deltamu_t(n, k + 1) = (a(n, t) + (rou + eta) / N_omega_n) * (E_n_k(n, k + 1) - V_n) - Deltadelta_t(n, k);
        % Deltadelta_t(n, k + 1) = (rou + eta) * (sum(E_nm_k(n, :, k + 1)) / N_omega_n - mean(W)) - Deltamu_t(n, k + 1) + a(n, t) * E_n_k(n, k + 1);

    end

    Deltamu(:, t - 1) = Deltamu_t(:, k);
    Deltadelta(:, t - 1) = Deltadelta_t(:, k);

    %% 价格更新
    price_tensor(:, :, t) = price_tensor(:, :, t - 1) - rou / 2 * (E_nm_k(:, :, k) + E_nm_k(:, :, k)');

    %% 映射更新交易量
    E_nm_l = zeros(N, N, 2);
    E_nm_l(:, :, 1) = E_nm_k(:, :, k);
    E_n_l = ones(N, 2);
    F_nm_l = zeros(N, N, 1);
    
    OD=[];
    for n=1:N
        OD=[OD ; find_connect(n,Connection_mat)];
    end
    l=1;
    N_OD=size(OD,1);
    while Jud_Enm(OD,E_nm_l(:,:,l),1e-7) && l<10000
        for index=1:N_OD
            n=OD(index,1);
            m=OD(index,2);
            F_nm_l(n, m, l) = (E_nm_l(n, m, l) - E_nm_l(m, n, l)) / 2;
            E_nm_l(n, m, l + 1) = F_nm_l(n, m, l);
        end
        for index=1:N_OD
            n=OD(index,1);
            m=OD(index,2);
            E_n_l(n,l+1)=sum(E_nm_l(n,:,l+1));
%             E_proj = max(min(E_n_l(n, l + 1), En_overline(n, t)), En_underline(n, t));
            if n==1
                E_proj = max(min(E_n_l(n, l + 1), En_overline(n, t)), max(En_underline(n, t),-30-sum(E_n(n,2:t-1))));
            elseif n==3
                E_proj = max(min(E_n_l(n, l + 1), min(En_overline(n, t),E_n(n,t-1)+0.5)), max(En_underline(n, t),E_n(n,t-1)-0.5));
            else
                E_proj = max(min(E_n_l(n, l + 1), En_overline(n, t)), En_underline(n, t));
            end
            E_nm_l(n, m, l + 1) = E_nm_l(n, m, l + 1) .* E_proj / E_n_l(n, l + 1);
        end
        l=l+1;
    end
    % for index=1:N_OD
    %     n=OD(index,1);
    %     m=OD(index,2);
    %     l=1;
    %     while abs(E_nm_l(n,m,l)+E_nm_l(m,n,l))>1e-3
    %         F_nm_l(n, m, l) = (E_nm_l(n, m, l) - E_nm_l(m, n, l)) / 2;
    %         E_nm_l(n, m, l + 1) = F_nm_l(n, m, l);
    %         E_n_l(n,l+1)=sum(E_nm_l(n,:,l+1));
    %         E_proj = max(min(E_n_l(n, l + 1), En_underline(n, t)), En_overline(n, t));
    %         E_nm_l(n, :, l + 1) = E_nm_l(n, :, l + 1) .* E_proj(n) / E_n_l(n, l + 1);
    %         l=l+1;
    %     end
    % end

    E_n(:,t)=E_n_l(:,l);
    E_nm(:,:,t)=E_nm_l(:,:,l);

end

toc
