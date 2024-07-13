function[T, Y] = calrates()
t = [0:7200];
p0=101325; %1 atm
v0=0.001; % initial volume in m^3
y1=0.45*p0*v0/(8.314*973); %mol P1V0/RT
y6=0.1*p0*v0/(8.314*973); %mol P6V0/RT
yar=y1; % CH4 and Ar have same pressure initially
n0=y1+y6+yar; % total moles initially
na=6.023*10^17; % Avo. const

s=[y1,0,0,0,0,y6,0,0,0,0,0,0,0,0,0,0];
[T, Y] = ode15s(@rates,t,s);
nt = sum(Y, 2) + yar; %total mole after initial
v = v0 * nt / n0; %volume other than v0

kf=[1.1*10^-8, 2.2*10^9, 1.5*10^11, (na)*2.2*10^-21, (na)*3.2*10^-22, (na)*1*10^-23, (na)*9.4*10^-27, (na)*2.3*10^-23, (na)*3.1*10^-21, (na)*6.9*10^-13, (na)*9.9*10^-24, 9.5*10^9, 8.7*10^10, 6.8*10^3, 8.4*10^5, 7.1*10^3, 5*10^9, 5.8*10^8, 4.4*10^-1, 1.4*10^8, 2.4*10^1, 4.7*10^6, (na)*1.9*10^-10, (na)*2.8*10^-10];
kb=[(na)*4.1*10^-10, (na)*1.5*10^-9, (na)*4.4*10^-9, (na)*1.9*10^-10, 4.1*10^9, (na)*1.9*10^-12, (na)*6.6*10^-11, 3.9*10^2, 4.7*10^4, (na)*5.7*10^-15, (na)*3.8*10^-19, (na)*3.3*10^-9, (na)*7.3*10^-9, (na)*3*10^-9, (na)*2.3*10^-10, (na)*1.4*10^-9, (na)*1.5*10^-9, (na)*3.1*10^-9, (na)*3.7*10^-12, (na)*5.1*10^-10, (na)*5*10^-10, (na)*3.9*10^-11, 1.8*10^-10, 2.7*10^-4];

kfn = 1.01 * kf; % Changed (increased) by 1%
keq = kf ./ kb;
kbn = kfn ./ keq;

plot(T,Y(:,1),'-')
hold on
title('solution using ode15s');
xlabel('time_t');
ylabel('conc_y');
x1 = v .* (-kfn(1) * (Y(:, 1) ./ v) + kbn(1) * (Y(:, 2) ./ v) .* (Y(:, 3) ./ v));
x2=0;
x3=0;
x4 = v .* (-kfn(4) * (Y(:, 1) ./ v) .* (Y(:, 6) ./ v) + kbn(4) * (Y(:, 8) ./ v) .* (Y(:, 2) ./ v));
x5 = v .* (-kfn(5) * (Y(:, 6) ./ v) .* (Y(:, 1) ./ v) + kbn(5) * (Y(:, 11) ./ v));
x6 = v .* (-kfn(6) * (Y(:, 1) ./ v) .* (Y(:, 5) ./ v) + kbn(6) * (Y(:, 10) ./ v) .* (Y(:, 2) ./ v));
x7 = v .* (-kfn(7) * (Y(:, 1) ./ v) .* (Y(:, 5) ./ v) + kbn(7) * (Y(:, 16) ./ v) .* (Y(:, 3) ./ v));
x8 = v .* (-kfn(8) * (Y(:, 1) ./ v) .* (Y(:, 5) ./ v) + kbn(8) * (Y(:, 12) ./ v));
x9 = v .* (-kfn(9) * (Y(:, 1) ./ v) .* (Y(:, 7) ./ v) + kbn(9) * (Y(:, 13) ./ v));
x10 = v .* (-kfn(10) * (Y(:, 1) ./ v) .* (Y(:, 3) ./ v) + kbn(10) * (Y(:, 2) ./ v) .* (Y(:, 4) ./ v));
x11 = v .* (-kfn(11) * (Y(:, 1) ./ v) .* (Y(:, 2) ./ v) + kbn(11) * (Y(:, 14) ./ v) .* (Y(:, 3) ./ v));

figure
plot(T(1:10), x1(1:10), 'DisplayName', 'x1', 'Color', 'black');
hold on;
plot(T(1:10), x4(1:10), 'DisplayName', 'x4', 'Color', 'magenta');
plot(T(1:10), x5(1:10), 'DisplayName', 'x5', 'Color', 'green');
plot(T(1:10), x6(1:10), 'DisplayName', 'x6', 'Color', 'yellow');
plot(T(1:10), x7(1:10), 'DisplayName', 'x7', 'Color', 'red');
plot(T(1:10), x8(1:10), 'DisplayName', 'x8', 'Color', '#FFA500');
plot(T(1:10), x9(1:10), 'DisplayName', 'x9', 'Color', 'blue');
plot(T(1:10), x10(1:10), 'DisplayName', 'x10', 'Color', "#EDB120");
plot(T(1:10), x11(1:10), 'DisplayName', 'x11', 'Color', 'cyan');
legend;
end