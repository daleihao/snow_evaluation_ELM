M = shaperead('../tl_2017_us_state/tl_2017_us_state');
figure;
hold on
for k=1:length(M)
    plot(M(k).X,M(k).Y,'color','k','linewidth',1);
end