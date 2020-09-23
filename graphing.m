clear all
%load data
eValues = readmatrix('eigenValues.csv');
eVectors = readmatrix('eigenVectors.csv');

%remove the indexing from the columns 
eValues = eValues(:,2);
eVectors = eVectors(:,2:287);

%---------explained varience ratio--------%
explained = (cumsum(eValues) ./ sum(eValues)).';

x = 1:56;

plot(x,explained,'b-o')
title('Explained varience ratio')
xlabel('n')
ylabel('cumsun(eV)/sum(eV)')
xlim([0 56]) 
ylim([0.92 1.05])
saveas(gcf,'Explained varience ratio.png')


%---------Projecting spectra on the PCs--------%
%we pick the first 3 P.Cs, if we use the method recommended 
%to us in the lectures

PCs = eVectors(1:3,:);
%read data
data = readmatrix('DS19hH2_dk0_FTIR_Spectra_instant_coffee.csv');
%discard the first row
wavelengths = data(1,:);
data = data(2:end,:);

%each row represents a spectra of the bean [1-29 are for Arabica beans] 
%and [30-57 are for Robusta beans]

for i = 1:56
    for j = 1:3
       projections(i,:,j) = dot(data(i,:), PCs(j,:));
    end
end

arabica = projections(1:29,:,:);
robusta = projections(30:end,:,:);

%plot projection

scatter3(arabica(:,:,1),arabica(:,:,2),arabica(:,:,3),'b','filled')
hold on
scatter3(robusta(:,:,1),robusta(:,:,2),robusta(:,:,3),'r','filled')
title('dim reduction [3-D]')
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
view([-30 10])
legend('Arabica bean', 'Robusta bean')
saveas(gcf,'Dim reducation 3-D PCs [1,2,3].png')
hold off


%-------Finding the two PCs that best seperate the two coffee beans------%

%find 56 choose 2
comb = nchoosek(1:56,2);

%since we have 1540 combinations,  manually checking the scatter plot of
%every is not possible, to solve this I will impliment a soft metahuristic
%to remove most combination, only leaving a small number that I can 
%later check manually.

%Metahuristic: compute the projection data for that combination, then find the 
%sum euclidean distance between each arabic point and all robust points,
%then check if that distance is larger than some arbitrary value. then I 
%plot the total number of 'true' values, this will revel the top 20-ish
%combinations that can be checked manually.

n = length(comb);
dists = zeros(29,27,n);
for i = 1:n
    PC1 = eVectors(comb(i,1),:);
    PC2 = eVectors(comb(i,2),:);
    PCs = [PC1;PC2];
    
    for k = 1:56
        for j = 1:2
            pj(k,:,j) = dot(data(k,:), PCs(j,:));
        end
    end
    %discard the last two arabica points, so the number of points match
    arabica = pj(1:29,:,:);
    arabica = squeeze(arabica);
    robusta = pj(30:end,:,:);
    robusta = squeeze(robusta);
    
    %lets calculate the distance from each point in the robusta set to each
    %point in the arabic set, if it exceeds some small distance, we set the
    %value to true, we then sum the total number of distances set as true.
    %ploting a histogram of the results should show us suitable
    %combinations.
    
    for k = 1:length(arabica)
       for j = 1:length(robusta)
           dists(k,j,i) = sqrt(((arabica(k,1)-robusta(j,1))^2+...
               ((arabica(k,2)-robusta(j,2))^2)));
       end
    end
    
    
end

%lets transform dists array into a logical array, based on if the distances
%exceed some arbitrary value.
limit = 180.0;
logical_dists = dists > limit;
hist_data = squeeze(sum(sum(logical_dists)));
histogram(hist_data);
title('Histogram indicating fittness of each PC combination')
xlabel('Combination n.o')
ylabel('incidents that points were further then distance set')
saveas(gcf,'hist_goodnessOfFittness.png')

%as expected the first 20-ish combinations have the highest difference,
%this is due to the fact that the power method finds the dominate
%eigenpairs first. 

%to confirm the results lets plot combination 832 i.e. 18th and 51st
%PCs. this has very low difference in x and y directions, so we should see
%very small differences between the two types of beans.
alldata = {eVectors, data, comb(832,:),1};
plot2dscatter(alldata)

%now lets look at the combinations that socred above 0, in the test above
logical_hist = hist_data > 0;
good_comb = comb(logical_hist,:);
for i = 1:length(good_comb)
    %change the 0 to a 1, if you want the images to be saved
    %this was implimented for later use in the report.
    alldata = {eVectors, data, good_comb(i,:),0};
    plot2dscatter(alldata)
end

%from the graphs we can clearly see that PCs 1 and 3 offer the best
%seperation of the two types of beans. combination number 2.

%save plot.
 alldata = {eVectors, data, good_comb(2,:),1};
 plot2dscatter(alldata)

%-------Spectra reconstruction------%

%lets reconstruct the first arabica and Robusta spectra

arabica_real = data(1,:);
robusta_real = data(30,:);

%calculate column means

means = mean(data);

%calculate projections

PCs = eVectors(1:6,:);

for i = 1:56
    for j = 1:6
       projections(i,:,j) = dot(data(i,:), PCs(j,:));
    end
end

projections = squeeze(projections);

%get projections for each spectra 
arabica = projections(1,:);
robusta = projections(30,:);

arabica_constrct = zeros(1,286);

for i = 1:286
    for j = 1:6
        %multiply projection by PCs.
        arabica_constrct(i) = arabica_constrct(i) + arabica(j) * PCs(j,i);
    end
end

plot(wavelengths, arabica_real,'r','DisplayName','True spectra')
hold on
plot(wavelengths,arabica_constrct,'b-','DisplayName','Reconstructed spectra')
title('Arabica coffee first spectra reconstruction using 6 PCs')
xlabel('Wavelength /nm')
ylabel('Reflectance')
legend
%legend('True spectra', 'Reconstructed spectra')
saveas(gcf,'Arabica reconstruction.png')
hold off



robusta_constrct = zeros(1,286);

for i = 1:286
    for j = 1:6
        %multiply projection by PCs.
        robusta_constrct(i) = robusta_constrct(i) + robusta(j) * PCs(j,i);
    end
end

plot(wavelengths, robusta_real,'r','DisplayName','True spectra')
hold on
plot(wavelengths,robusta_constrct,'b-','DisplayName','Reconstructed spectra')
title('Robusta coffee first spectra reconstruction using 6 PCs')
xlabel('Wavelength /nm')
ylabel('Reflectance')
legend
%legend('True spectra', 'Reconstructed spectra')
hold off
saveas(gcf,'Robusta reconstruction.png')
stop = 1;

function [] = plot2dscatter(varargin)
    %extract data
    eVectors = varargin{1}{1};
    data = varargin{1}{2};
    comb = varargin{1}{3};
    save = varargin{1}{4};
    %compute the projections
    PC1 = eVectors(comb(1,1),:);
    PC2 = eVectors(comb(1,2),:);
    PCs = [PC1;PC2];
    
    for k = 1:56
        for j = 1:2
            pj(k,j) = dot(data(k,:), PCs(j,:));
        end
    end
    
    arabica = pj(1:29,:);
    robusta = pj(30:end,:);
    
    %plot data
    
    scatter(arabica(:,1),arabica(:,2),'b','filled')
    hold on
    scatter(robusta(:,1),robusta(:,2),'r','filled')
    title_text = "dim reduction [2-D] PCs [" + num2str(comb(1,1)) + ',' + ...
        num2str(comb(1,2)) + ']';
    title(title_text(1))
    xlabel('PC1')
    ylabel('PC2')
    legend('Arabica bean', 'Robusta bean')
    
    
    if save == 1
        name = "2-D_Scatter_" + num2str(comb(1,1)) + '_' + num2str(comb(1,2))...
            + ".png";
            
       saveas(gcf,name) 
    end
    
    hold off
    
    
end



