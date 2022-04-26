%% Look for stability and phase portrait
% The code below follows the usual steps for linear stability

%% Author
% Daniele Proverbio, 20/08/2019
% daniele.proverbio@outlook.com
% University of Luxembourg

%% Prepare env
clear; close all; clc;

%% Perform calculation

syms x            % Working with symbolic manipulation (more precise)
K = 0.1;          % Basal expression (constant, within accepted range)
c = 1:0.002:3.5;   % Max Production (control parameter)

solutions=[];   % Vector of equilibria solutions
c_vector=[];    % Vector of accepted control parameters
cmap = [];      % Vector of colors

enne = 10;

for m = 1:length(c)     % Check all c
    
    f = K + (c(m)*(x^enne))/(1+(x^enne)) - x;     % My equation (vector field)
    soly = vpasolve(f == 0, x);         % Look for roots
    f_prime = diff(f);                  % Estimate derivatives
    
    for n = 1:length(soly)
        if isreal(soly(n))              % Only real roots allowed, obviously
            solutions = [solutions,soly(n)];
            c_vector = [c_vector,c(m)];
            if(vpa(subs(f_prime,x,soly(n))) < 0 )   % Check for linear stability -> eigenvalue <> 0
                cmap = [cmap; [0,0.6,1]];             % If Stable point, color it in blue
            else
                cmap = [cmap; [1,0.6,0]];             % If unstable, color it in red
                                                    % Color shades are
                                                    % changed according to
                                                    % n ([. , 0 -> 0.6 , .])
            end
        end
    end
end

% plot
figure;
scatter(c_vector,solutions,10,cmap,'filled');
ylim([0 3])
xlim([1, 3])
ax = gca;
ax.FontSize = 18; 
set(gca,'XMinorTick','on','YMinorTick','on')
%title('Equilibria of the system','FontSize',15);
xlabel('c','FontSize',30);
ylabel('$\mathbf{\tilde{x}}$','FontSize',30,'interpreter','latex');
%saveas(gcf,'equilibria.eps');


%% Estimate Focal Width (Delta)
% NB: from now on, I call "p" what is "c" in the code, to be consistent
% with analytical notation.

% I know it's a local fold -> x_s ~ Sqrt(p)
% Delta is given by the multiplier of Sqrt(p)
% Hence I can fit x_s = a*Sqrt(p) + b next to the fold and get a = Delta
% More here https://www.quora.com/How-do-I-determine-the-vertex-focus-and-directrix-of-a-parabola

% So: 1) Identify p_c
%     2) Extrapolate points (p,x_s)=(p_c + something; x_s_c + something)
%     3) Fit p = a + b*x_s + c*x_s^2 (it's better in order to account for
%     translations etc.)
%     4) Get Delta
%     5) Repeat for every n so to have Delta = Delta(n)
%     6) Use theoretical results to link n and the behavior of EWS

% ------
% 1) Identify p_c
% the first time c is repeated in c_vector, that is the first time I have
% bistability, hence the fold.

store_index = 0;
for ii = 3 : length(c_vector)
    if mod(enne,2) == 0
        if(c_vector(ii)==c_vector(ii-1)  && store_index==0)
        store_index=ii+1;
        end
    else
        if(c_vector(ii)==c_vector(ii-1)  && store_index==0) && c_vector(ii-1)==c_vector(ii-2)
        store_index=ii+1;
        end
    end
end


% ------
% 2) Estrapolate points (p_points,x_points) in an interval of 30 points
% after p_c

if mod(enne,2)==0
    n_points = 30;
    %p_points = (c_vector(store_index:3:store_index+n_points));                  % p_points = (c_vector(store_index+4:4:store_index+120)); for odd n (and +120 instead of 90)
    p_points = [(c_vector(store_index-1:3:store_index+n_points)) (c_vector(store_index:3:store_index+n_points))]; 
    %x_points = double(solutions(store_index:3:store_index+n_points));
    x_points = [double(solutions(store_index-1:3:store_index+n_points)) double(solutions(store_index:3:store_index+n_points))];

else 
    n_points = 40;
    %p_points = (c_vector(store_index:4:store_index+n_points));
    p_points = [(c_vector(store_index-1:4:store_index+n_points)) (c_vector(store_index:4:store_index+n_points))]; 
    %x_points = double(solutions(store_index:4:store_index+n_points));
    x_points = [double(solutions(store_index-1:4:store_index+n_points)) double(solutions(store_index:4:store_index+n_points))];

end
% figure()
% scatter(x_points,p_points)


% ------
% 3) Fit

fit = polyfit(x_points,p_points,2);

% Compare data with fit
p_fit = fit(1).*x_points.*x_points + fit(2).*x_points + fit(3);

% figure()
% hold on
% plot(x_points,p_fit);
% plot(x_points,p_points);
% hold off
% 
% % Inspect the fitted parabula to some other values of x
% x_test = -5:0.01:5;
% p_test = fit(1).*x_test.*x_test + fit(2).*x_test + fit(3);
% plot(x_test,p_test);


% ------
% 4) Get the focus and the focal distance
x_focus = -fit(2)/(2*fit(1));
y_focus = (1-(fit(2)*fit(2)-4*fit(1)*fit(3)))/(4*fit(1));

syms z
x_comparison = max(double(solve(fit(1)*z^2 + fit(2)*z + fit(3) - y_focus == 0, z)));

focal_width = abs(x_comparison - x_focus);


% -------
% 5) estimate uncertainties

fit2 = polyfitn(x_points,p_points,2);
std_x_focus = sqrt((0.5*fit2.ParameterStd(1)/fit2.Coefficients(1)).^2+(fit2.ParameterStd(2) * 0.5*fit2.Coefficients(2)/(fit2.Coefficients(1)).^2).^2);
std_y_focus = sqrt((fit2.ParameterStd(2)*fit2.Coefficients(2)/(2*fit2.Coefficients(1)))^2 + (fit2.ParameterStd(1)*(1-fit2.Coefficients(2)^2)/(4*fit2.Coefficients(1)^2))^2 + fit2.ParameterStd(3)^2);



%% Focal_width vs n
% After having computed the focal width for different n, I check now their
% relationship
% Data on Quad 2.1, 12/11/19
% Updated on Quad 2.4, 1/4/22

n = [2,3,4,5,6,7,8];
FW = [0.5208,0.4270,0.3338,0.2686,0.2230,0.1885,0.1640]; % collection of focal_width
errors = FW*20/100; % after error propagation on z and y_focus and x_focus, this is a good estimate of the relative error (cf Quad 2.3, 14/12/21)

coeffs = coeffvalues(fit_fw);
R_squared = 0.9989;% from Fitting toolbox

figure(Position=[1,1,450,400])
errorbar(n,FW,errors,'o', LineWidth=1.5,color='black');
xlabel('n');
ylabel('Focal Width');
ax = gca;
ax.FontSize = 20;
ylim([0,0.8])
xlim([1,9])
%axis([1 9 0.7 inf])
txt = ['Best fit:'];
txt1 = ['FW = 2.7 / (-0.30 n -1.7 )^2']; % the values come from the fitting (Toolbox) on n_to10 and FW_to10 below
txt2 = ['Adj. R-squared = 0.9945'];
text(4.5,0.75,txt,fontsize=15)
text(4.5,0.7,txt1,fontsize=15)
text(4.5,0.64,txt2,fontsize=15)


%% NB: c'e' parecchia dipendenza dei numeri esatti su quanto lontano dalla biforcazione vado (sia usando c_vector piu' coarse, sia prendendo piu' datapoint per il fitting)
%% In ogni caso, la relazione tra i vari punti e' mantenuta (cio' che cambia e' solo unaa traslazione di tutti i punti verso il basso.
%% Magari possiamo discutere questo aspetto nel Sup Mat

n_to10 = [2,3,4,5,6,7,8,9,10];
FW_to10 = [0.5208,0.4270,0.3338,0.2686,0.2230,0.1885,0.1640,0.1446,0.1290];

% check trends of FW vs n
% refer to Quad 2.4, 5/4/2022
x = [0,0.01:30];
cofs_A = [2.723,-0.2944,-1.685];
cofs_B = [0.7818,-0.2061];
cofs_C = [0.7803,-0.2809,0.08087];
figure()
hold on
plot(x,cofs_A(1)./(cofs_A(2)*x + cofs_A(3)).^2)
plot(x,cofs_B(1)*exp(cofs_B(2)*x))
plot(x,cofs_C(1)*exp(cofs_C(2)*x)+cofs_C(3))
legend({'A',"B",'C'})



