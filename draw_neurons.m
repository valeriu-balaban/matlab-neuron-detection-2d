function draw_neurons(I, Neurons)

figure;
imshow(I,[]), hold on;
axis on, xlabel x, ylabel y;
plot(Neurons(:,1),Neurons(:,2),'.','color','red', 'MarkerSize', 5);