function [] = plotstuff( out )


%plot(vt(4, :)')
subplot(4, 1, 1)
%plot(debug(:, 1:out.group_sizes(1)));
plot(squeeze(out.vart(1, :, :))');
title('Neuron spike arrival time');
%legend({'N1', 'N2', 'N3'});

subplot(4, 1, 2)
plot(squeeze(out.vart(1, :, :))');
%set(gca, 'Color', 'k')
%legend({'N1', 'N2', 'N3'});
title('variance');

subplot(2, 2, 3);
%plot(out.vt( net.group_sizes(1) + 1:end, :)')
plot(out.vt(11:12, :)');
title('output response');

subplot(2, 2, 4);
plot(out.spike_time_trace(out.spike_time_trace(:, 2) == 11, 1), mod(out.spike_time_trace(out.spike_time_trace(:, 2) == 11, 1), 500))
hold on
plot(out.spike_time_trace(out.spike_time_trace(:, 2) == 12, 1), mod(out.spike_time_trace(out.spike_time_trace(:, 2) == 12, 1), 500))
hold off
axis([-Inf Inf 5 14]);
title('output spike time');


end