T = 10;
N = 3;
rou = sqrt(T);
eta = sqrt(T);
Nir = 10000; %迭代次数
Deltamu = zeros(N, T);
Deltadelta = zeros(N, T);
E_n = zeros(N, T+1);
E_nm=zeros(N,N,T+1);
%1为用户，2为CG，3为RE；2，3分别与1独立相连
Connection_mat = [0 1 1; 1 0 0; 1 0 0];
Consumer = [1];
Producer = [2 3];
% a=[0.0144 0.0210 0.01];
% b=[6.4149 15.0413 5];

% En_underline=[-5.3252 0 1.3021];
% En_overline=[-2.4533 4.9014 1.3021];
a = [0.0144	0.021	0.01;
    0.0144	0.021	0.01;
    0.0177	0.021	0.01;
    0.0119	0.021	0.01;
    0.0145	0.021	0.01;
    0.0171	0.021	0.01;
    0.0128	0.021	0.01;
    0.0166	0.021	0.01;
    0.0112	0.021	0.01;
    0.0196	0.021	0.01;
    0.0159	0.021	0.01;
    ];
a = a';
b = [6.4149	15.0413	5;
    6.9078	15.0413	5;
    8.976	15.0413	5;
    7.4488	15.0413	5;
    8.2316	15.0413	5;
    8.7734	15.0413	5;
    8.3985	15.0413	5;
    5.8131	15.0413	5;
    7.4918	15.0413	5;
    6.7019	15.0413	5;
    6.1191	15.0413	5;
    ];
b = b';
% En_underline = [-5.3252	0	1.3021;
%             - 5.3252	0	1.3021;
%             - 5.3252	0	1.3067;
%             - 5.3252	0	1.3067;
%             - 5.3252	0	1.2739;
%             - 5.3252	0	1.2985;
%             - 5.3252	0	1.2748;
%             - 5.3252	0	1.3012;
%             - 5.3252	0	1.2758;
%             - 5.3252	0	1.3018;
%             - 5.3252	0	1.3018;
%             ];
En_underline = [-5.3252	1e-4	1.3021;
            - 5.3252	1e-4	1.3021;
            - 5.3252	1e-4	1.3067;
            - 5.3252	1e-4	1.3067;
            - 5.3252	1e-4	1.2739;
            - 5.3252	1e-4	1.2985;
            - 5.3252	1e-4	1.2748;
            - 5.3252	1e-4	1.3012;
            - 5.3252	1e-4	1.2758;
            - 5.3252	1e-4	1.3018;
            - 5.3252	1e-4	1.3018;
            ];
En_underline = En_underline';
En_overline = [-2.4533	4.9014	1.3021;
            - 2.4533	4.9014	1.3021;
            - 2.4533	4.9014	1.3067;
            - 2.4533	4.9014	1.3067;
            - 2.4533	4.9014	1.2739;
            - 2.4533	4.9014	1.2985;
            - 2.4533	4.9014	1.2748;
            - 2.4533	4.9014	1.3012;
            - 2.4533	4.9014	1.2758;
            - 2.4533	4.9014	1.3018;
            - 2.4533	4.9014	1.3018;
            ];
En_overline = En_overline';
price = 16*Connection_mat;
% price=[0 15 14.7; 15 0 0 ; 14.7 0 0];
price_tensor = zeros(N, N, 1);
price_tensor(:, :, 1) = price;
