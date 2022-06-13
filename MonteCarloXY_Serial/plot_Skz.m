L = 6;
kz = 0:2*(pi)/L:pi*(2*L-1)/L;

for j1 = 0:1:3
    figure()
    hold on
    title(['$J_{z1} = ' num2str(j1) '$, $J_{z2} = ' num2str(j1) '$'], 'Interpreter', 'latex', 'FontSize', 18)
    set(gca, 'ColorOrder', flipud(parula(2)));
    legstr = {};
    for b = 0:1:1
        try
            corr = csvread(['corr_j' num2str(j1) '_beta' num2str(b) '.csv']);
            corr = reshape(corr, L, L, L);
            sk = real(fft(corr(:,1,1))); %real(fft(sum(sum(corr, 3), 2)));
            plot(kz/pi, log(sk), '*-')
        end
    end
    box on
    xlabel('$k_z/\pi$', 'Interpreter', 'latex', 'FontSize', 18)
    ylabel('$log(S(k_z))$', 'Interpreter', 'latex', 'FontSize', 18)
end