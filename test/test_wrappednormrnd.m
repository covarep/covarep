
% Test WND random number generator


if 1
    figure
    sigmas = exp(-7:0.01:log(6));
    estsigmas = zeros(size(sigmas));
    tic;
    for si=1:length(sigmas)
        wgr = wrappednormrnd(0, sigmas(si), 10000, 1);
        [dummy es] = wrappednormestim(wgr);
        estsigmas(si) = es;
    end
    toc;

    hold off;
    plot(sigmas, sigmas, 'k');
    hold on;
    plot(sigmas, estsigmas, 'b');
    grid on;
end

if 0
    figure
    m = pi/3;
    s = 4;
    ems = [];
    ess = [];
    while length(ess)<1001
        % Generate
        wgr = wrappednormrnd(m, s, (length(ess)+1), 1);

        % Estimate
        [em es] = wrappednormestim(wgr);

        ems = [ems, log10(abs(angle(exp(1i*m)*exp(-1i*em)))/pi)];
        ess = [ess, log10(abs((s-es)/s))];

        subplot(211);
        plot(ems, 'k');
%          ylim([-6 0]);
        title(['Log10 norm abs error of mean. m=' num2str(m) ' em=' num2str(em)]);
        subplot(212);
        plot(ess, 'k');
%          ylim([-1 1]);
        title(['Log10 relative-error of standard-deviation. s=' num2str(s) ' es=' num2str(es)]);
        pause(0.01);
    end
end


% Test accuracy of the distribution computation
if 0
    xs = -pi:0.001:pi;

    N = 10;
    s = 8;

    p = wrappednormpdf(xs, 0, s, 10000);

    plot(xs, log(p), 'k');
    hold on;
    grid on;

    pdef = wrappednormpdf(xs, 0, s, N);
    plot(xs, log(pdef), 'b');

    disp(['mean log error=' num2str(sqrt(mean(log(p)-log(pdef)).^2))]);
    disp(['uniformity with s=' num2str(s) ', (the smaller the more uniform)=' num2str(std(log(p)))]);

    ylim([-1.9 log(max(p)+0.01)]);
end
