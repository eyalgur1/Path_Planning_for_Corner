format compact

t = -1:0.005:1;  % time
Alpha = logspace(0,1.5,4);  % x-axis acceleration bound
Beta = Alpha;  % y-axis acceleration bound

%%
xbar = @(t,eps)-1i*sqrt(2)*eps*ellipticE(1i*asinh(t/eps),0.5);  % the function x_bar(t)
Cbar = @(eps)1-sqrt(1+eps^2);
Dbar = @(eps)1 - xbar(1,eps);
eta_x = @(t,eps)xbar(t,eps) + Dbar(eps);
eta_y = @(t,eps)sqrt(t.^2 + eps^2) + Cbar(eps);
C_deviate = @(eps) 1 - 2*xbar(1,eps);
x_acc = @(t,eps)(t*eps^2)./(sqrt(t.^2+2*eps^2).*(t.^2+eps^2).^(1.5));
y_acc = @(t,eps)(eps^2)./((t.^2+eps^2).^(1.5));

%%
legend_titles = {'$\gamma\left(t\right)$'};
plot(t, abs(t), 'black', 'LineWidth', 3); hold on
for alpha = Alpha
    opt_eps = epsilon_fun(alpha, alpha);
    C_const = Cbar(opt_eps);
    D_const = Dbar(opt_eps);
    plot(eta_x(t,opt_eps), eta_y(t,opt_eps), 'LineWidth', 2)
    legend_titles{end+1} = ['$\alpha=\beta=$ ',num2str(alpha)];
end
xlim([-1 1])
legend(legend_titles,'Location','southeast','FontSize',18,'interpreter','latex')
hold off

%%
C_deviate_val = [];
for eps_val=logspace(-3,2,200)
    C_deviate_val(end+1) = C_deviate(eps_val);
end
plot(abs(-1-C_deviate_val), 'LineWidth', 2)
xticks(0:40:200)
xticklabels({'0','0.0095','0.0966','0.9771','9.8850','100'})
hold on

B_deviate_val = [];
for eps_val=logspace(-3,2,200)
    B_deviate_val(end+1) = min(t.^2+eta_y(t,eps_val).^2);
end
plot(sqrt(B_deviate_val), 'LineWidth', 2)
xticks(0:40:200)
xticklabels({'0','0.0095','0.0966','0.9771','9.8850','100'})

int_deviate_val = [];
for eps_val=logspace(-3,2,200)
    int_deviate_val(end+1) = integral(@(t)eta_y(t,eps_val)-abs(t),-1,1);
end
plot(abs(int_deviate_val), 'LineWidth', 2)
xticks(0:40:200)
xticklabels({'0','0.0095','0.0966','0.9771','9.8850','100'})
legend('$\textrm{Deviation from }\mathcal{C}$','$\textrm{Deviation from }\mathcal{B}$','$\textrm{Difference betwee curves}$','Location','southeast','FontSize',18,'interpreter','latex')
%\int_{-1}^1\left(\gamma\left(t\right)-\eta\left(t\right)\right)dt

%%
epsilon_fun(1, 1)
alp = 2;
bet = 8;
x_acc_val = x_acc(t,epsilon_fun(alp, bet));
y_acc_val = y_acc(t,epsilon_fun(alp, bet));
plot(t, x_acc_val, 'LineWidth', 2)
hold on
plot(t, y_acc_val, 'LineWidth', 2)
legend('$x$-axis','$y$-axis','Location','southeast','FontSize',18,'interpreter','latex')
title('$\textrm{Acceleration for }\alpha=2\textrm{ and }\beta=8$','FontSize',18,'interpreter','latex')
grid on

function opt_epsilon = epsilon_fun(alpha, beta)
const1 = (3*sqrt(3*sqrt(10)-9))/((1+sqrt(10))^(1.5));
const2 = sqrt(3/(sqrt(10)-2));
const3 = 4/sqrt(189+33*sqrt(33));
phi = @(delta)(delta^2)/(sqrt(1+2*delta^2)*(1+delta^2)^(1.5)) - alpha;
if max(1/beta, const1/alpha)<=const2
    opt_epsilon = max(1/beta, const1/alpha);
elseif alpha>=const3
    opt_epsilon = 1/beta;
else
    delta_star = fzero(phi,1);
    opt_epsilon = max(1/beta, abs(delta_star));
end
end


