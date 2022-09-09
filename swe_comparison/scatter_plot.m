
i = 2;
a = ELM_swes(:,:,12*9+i);
b = UA_swes(:,:,12*9+i);
c = SNODAS_swes(:,:,12*9+i);

filters = a>0 & b>0 & c>0;


f = scatter(a(filters),b(filters),10,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);

corrcoef(a(filters),b(filters))
p = polyfit(a(filters),b(filters),1);
x1 = [0 1000];
y1 = polyval(p,x1);
hold on
plot(x1,y1)
hold off


%% spatial r2
Rs = nan(144, 168);
for i = 1:144
    for j = 1:168
        a = squeeze(ELM_swes(i,j,:));
        b = squeeze(UA_swes(i,j,:));
        
        filters = a>0 & b>0;
        if(sum(filters)>12)
            tmp = corrcoef(a(filters),b(filters));
            Rs(i,j) = tmp(1,2);
        end
    end
end
figure
imagesc(Rs,[-1 1])


Rs = nan(144, 168);
for i = 1:144
    for j = 1:168
        a = squeeze(ELM_swes(i,j,:));
        b = squeeze(SNODAS_swes(i,j,:));
        
        filters = a>0 & b>0;
        if(sum(filters)>12)
            tmp = corrcoef(a(filters),b(filters));
            Rs(i,j) = tmp(1,2);
        end
    end
end
figure
imagesc(Rs,[-1 1])

%% snow cover
%% spatial r2
Rs = nan(144, 168);
for i = 1:144
    for j = 1:168
        a = squeeze(ELM_fsnos(i,j,1:(19*12)));
        b = squeeze(MODSCAG_fsnos(i,j,:));
        
        filters = a>0 & b>0;
        if(sum(filters)>12)
            tmp = corrcoef(a(filters),b(filters));
            Rs(i,j) = tmp(1,2);
        end
    end
end
figure
imagesc(Rs,[-1 1])

%% spatial r2
Rs = nan(144, 168);
for i = 1:144
    for j = 1:168
        a = squeeze(MODSCAG_fsnos(i,j,1:(19*12)));
        b = squeeze(SPIRES_fsnos(i,j,:));
        
        filters = a>0 & b>0;
        if(sum(filters)>12)
            tmp = corrcoef(a(filters),b(filters));
            Rs(i,j) = tmp(1,2);
        end
    end
end
figure
imagesc(Rs,[-1 1])


%% spatial r2
Rs = nan(144, 168);
for i = 1:144
    for j = 1:168
        a = squeeze(ELM_BSAs(i,j,1:(19*12)));
        b = squeeze(MCD43_BSAs(i,j,1:(19*12)));
        snows = squeeze(ELM_fsnos(i,j,1:(19*12)));
        filters = a>0 & b>0 &  snows>0;
        if(sum(filters)>12)
            tmp = corrcoef(a(filters),b(filters));
            Rs(i,j) = tmp(1,2);
        end
    end
end
figure
imagesc(Rs,[-1 1])

%% albedo
Rs = nan(144, 168);
for i = 1:144
    for j = 1:168
        a = squeeze(ELM_snowmelts(i,j,1:(19*12)));
        b = squeeze(SNODAS_snowmelts(i,j,1:(19*12)));
        filters = a>0 & b>0;
        if(sum(filters)>12)
            tmp = corrcoef(a(filters),b(filters));
            Rs(i,j) = tmp(1,2);
        end
    end
end
figure
imagesc(Rs,[-1 1])

%% albedo
Rs = nan(144, 168);
for i = 1:144
    for j = 1:168
        a = squeeze(ELM_grainsizes(i,j,1:(19*12)));
        b = squeeze(MODSCAG_grainsizes(i,j,1:(19*12)));
        filters = a>0 & b>0;
        if(sum(filters)>12)
            tmp = corrcoef(a(filters),b(filters));
            Rs(i,j) = tmp(1,2);
        end
    end
end
figure
imagesc(Rs,[-1 1])

%% albedo
Rs = nan(144, 168);
for i = 1:144
    for j = 1:168
        a = squeeze(ELM_deltaAlbedos(i,j,1:(19*12)));
        b = squeeze(MODSCAG_deltaAlbedos(i,j,1:(19*12)));
        filters = a>0 & b>0;
        if(sum(filters)>12)
            tmp = corrcoef(a(filters),b(filters));
            Rs(i,j) = tmp(1,2);
        end
    end
end
figure
imagesc(Rs,[-1 1])