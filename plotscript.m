%% common things to plot for individual parameters

    xm = prctile(data,50,3);
    xl = prctile(data,100*(1-ci)./2,3);
    xr = prctile(data,100-100*(1-ci)./2,3);
    [xs,I] = sort(xm);
    %for ii = 1:numel(xm)
        l1 = plot(xm(I),(1:numel(xm))'./numel(xm),'.k');
    for ii = 1:numel(xm)
        l2 = plot([xl(I(ii)),xr(I(ii))],repmat(ii./numel(xm),[1 2]),'-k');
    end
    ksgrid = linspace(min(xm)-0.1*range(xm),max(xm)+0.1*range(xm),101);
    ks = ksdensity(xm,ksgrid);
    l4 = plot(ksgrid,ks./max(ks),'-k','LineWidth',2);