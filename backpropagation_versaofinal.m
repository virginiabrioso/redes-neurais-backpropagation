close all, clear all, clc
TrainData = [5.1 3.5 1.4 0.2   % Iris-setosa     
               4.9 3.0 1.4 0.2   % Iris-setosa  
               4.7 3.2 1.3 0.2   % Iris-setosa   
               4.6 3.1 1.5 0.2   % Iris-setosa   
               5.0 3.6 1.4 0.2   % Iris-setosa   
               5.4 3.9 1.7 0.4   % Iris-setosa   
               4.6 3.4 1.4 0.3   % Iris-setosa   
               5.0 3.4 1.5 0.2   % Iris-setosa   
               4.4 2.9 1.4 0.2   % Iris-setosa   
               4.9 3.1 1.5 0.1   % Iris-setosa   
               5.4 3.7 1.5 0.2   % Iris-setosa   
               4.8 3.4 1.6 0.2   % Iris-setosa   
               4.8 3.0 1.4 0.1   % Iris-setosa   
               4.3 3.0 1.1 0.1   % Iris-setosa   
               5.8 4.0 1.2 0.2   % Iris-setosa   
               5.7 4.4 1.5 0.4   % Iris-setosa   
               5.4 3.9 1.3 0.4   % Iris-setosa   
               5.1 3.5 1.4 0.3   % Iris-setosa 
               5.7 3.8 1.7 0.3   % Iris-setosa 
               5.1 3.8 1.5 0.3   % Iris-setosa 
               5.4 3.4 1.7 0.2   % Iris-setosa 
               5.1 3.7 1.5 0.4   % Iris-setosa 
               4.6 3.6 1.0 0.2   % Iris-setosa 
               5.1 3.3 1.7 0.5   % Iris-setosa 
               4.8 3.4 1.9 0.2   % Iris-setosa 
               5.0 3.0 1.6 0.2   % Iris-setosa 
               5.0 3.4 1.6 0.4   % Iris-setosa 
               5.2 3.5 1.5 0.2   % Iris-setosa 
               5.2 3.4 1.4 0.2   % Iris-setosa 
               4.7 3.2 1.6 0.2   % Iris-setosa 
               5.0 3.5 1.3 0.3   % Iris-setosa 
               4.5 2.3 1.3 0.3   % Iris-setosa 
               4.4 3.2 1.3 0.2   % Iris-setosa 
               5.0 3.5 1.6 0.6   % Iris-setosa 
               5.1 3.8 1.9 0.4   % Iris-setosa 
               4.8 3.0 1.4 0.3   % Iris-setosa 
               5.1 3.8 1.6 0.2   % Iris-setosa 
               4.6 3.2 1.4 0.2   % Iris-setosa 
               5.3 3.7 1.5 0.2   % Iris-setosa 
               5.0 3.3 1.4 0.2   % Iris-setosa 
               7.0 3.2 4.7 1.4   % Iris-versicolor
               6.4 3.2 4.5 1.5   % Iris-versicolor
               6.9 3.1 4.9 1.5   % Iris-versicolor
               5.5 2.3 4.0 1.3   % Iris-versicolor
               6.5 2.8 4.6 1.5   % Iris-versicolor
               5.7 2.8 4.5 1.3   % Iris-versicolor
               6.3 3.3 4.7 1.6   % Iris-versicolor
               4.9 2.4 3.3 1.0   % Iris-versicolor
               6.6 2.9 4.6 1.3   % Iris-versicolor
               5.2 2.7 3.9 1.4   % Iris-versicolor
               5.0 2.0 3.5 1.0   % Iris-versicolor
               5.9 3.0 4.2 1.5   % Iris-versicolor
               6.0 2.2 4.0 1.0   % Iris-versicolor
               6.1 2.9 4.7 1.4   % Iris-versicolor
               5.6 2.9 3.6 1.3   % Iris-versicolor
               6.7 3.1 4.4 1.4   % Iris-versicolor
               5.6 3.0 4.5 1.5   % Iris-versicolor
               5.8 2.7 4.1 1.0   % Iris-versicolor
               6.2 2.2 4.5 1.5   % Iris-versicolor
               5.6 2.5 3.9 1.1   % Iris-versicolor
               5.9 3.2 4.8 1.8   % Iris-versicolor
               6.1 2.8 4.0 1.3   % Iris-versicolor
               6.3 2.5 4.9 1.5   % Iris-versicolor
               6.1 2.8 4.7 1.2   % Iris-versicolor
               6.4 2.9 4.3 1.3   % Iris-versicolor
               6.6 3.0 4.4 1.4   % Iris-versicolor
               6.8 2.8 4.8 1.4   % Iris-versicolor
               6.7 3.0 5.0 1.7   % Iris-versicolor 
               6.0 2.9 4.5 1.5   % Iris-versicolor 
               5.7 2.6 3.5 1.0   % Iris-versicolor
               5.5 2.4 3.8 1.1   % Iris-versicolor
               5.5 2.4 3.7 1.0   % Iris-versicolor
               5.8 2.7 3.9 1.2   % Iris-versicolor
               6.0 2.7 5.1 1.6   % Iris-versicolor
               5.4 3.0 4.5 1.5   % Iris-versicolor
               6.0 3.4 4.5 1.6   % Iris-versicolor
               6.7 3.1 4.7 1.5   % Iris-versicolor
               6.3 2.3 4.4 1.3   % Iris-versicolor
               5.6 3.0 4.1 1.3   % Iris-versicolor
               5.5 2.5 4.0 1.3   % Iris-versicolor
               6.3 3.3 6.0 2.5   % Iris-verginica
               5.8 2.7 5.1 1.9   % Iris-verginica
               7.1 3.0 5.9 2.1   % Iris-verginica
               6.3 2.9 5.6 1.8   % Iris-verginica
               6.5 3.0 5.8 2.2   % Iris-verginica
               7.6 3.0 6.6 2.1   % Iris-verginica
               4.9 2.5 4.5 1.7   % Iris-verginica
               7.3 2.9 6.3 1.8   % Iris-verginica
               6.7 2.5 5.8 1.8   % Iris-verginica
               7.2 3.6 6.1 2.5   % Iris-verginica
               6.5 3.2 5.1 2.0   % Iris-verginica
               6.4 2.7 5.3 1.9   % Iris-verginica
               6.8 3.0 5.5 2.1   % Iris-verginica
               5.7 2.5 5.0 2.0   % Iris-verginica
               5.8 2.8 5.1 2.4   % Iris-verginica
               6.4 3.2 5.3 2.3   % Iris-verginica
               6.5 3.0 5.5 1.8   % Iris-verginica
               7.7 3.8 6.7 2.2   % Iris-verginica
               7.7 2.6 6.9 2.3   % Iris-verginica
               6.0 2.2 5.0 1.5   % Iris-verginica
               6.9 3.2 5.7 2.3   % Iris-verginica
               5.6 2.8 4.9 2.0   % Iris-verginica
               7.7 2.8 6.7 2.0   % Iris-verginica
               6.3 2.7 4.9 1.8   % Iris-verginica
               6.7 3.3 5.7 2.1   % Iris-verginica
               7.2 3.2 6.0 1.8   % Iris-verginica
               6.2 2.8 4.8 1.8   % Iris-verginica
               6.1 3.0 4.9 1.8   % Iris-verginica
               6.4 2.8 5.6 2.1   % Iris-verginica
               7.2 3.0 5.8 1.6   % Iris-verginica
               7.4 2.8 6.1 1.9   % Iris-verginica
               7.9 3.8 6.4 2.0   % Iris-verginica
               6.4 2.8 5.6 2.2   % Iris-verginica
               6.3 2.8 5.1 1.5   % Iris-verginica
               6.1 2.6 5.6 1.4   % Iris-verginica
               7.7 3.0 6.1 2.3   % Iris-verginica
               6.3 3.4 5.6 2.4   % Iris-verginica
               6.4 3.1 5.5 1.8   % Iris-verginica
               6.2 3.4 5.4 2.3   % Iris-verginica
               5.9 3.0 5.1 1.8]; % Iris-verginica
  
 TestData = [  4.8 3.1 1.6 0.2   % Iris-setosa 
               5.4 3.4 1.5 0.4   % Iris-setosa 
               5.2 4.1 1.5 0.1   % Iris-setosa 
               5.5 4.2 1.4 0.2   % Iris-setosa 
               4.9 3.1 1.5 0.1   % Iris-setosa 
               5.0 3.2 1.2 0.2   % Iris-setosa 
               5.5 3.5 1.3 0.2   % Iris-setosa 
               4.9 3.1 1.5 0.1   % Iris-setosa 
               4.4 3.0 1.3 0.2   % Iris-setosa 
               5.1 3.4 1.5 0.2   % Iris-setosa 
               5.5 2.6 4.4 1.2   % Iris-versicolor
               6.1 3.0 4.6 1.4   % Iris-versicolor
               5.8 2.6 4.0 1.2   % Iris-versicolor
               5.0 2.3 3.3 1.0   % Iris-versicolor
               5.6 2.7 4.2 1.3   % Iris-versicolor
               5.7 3.0 4.2 1.2   % Iris-versicolor
               5.7 2.9 4.2 1.3   % Iris-versicolor
               6.2 2.9 4.3 1.3   % Iris-versicolor
               5.1 2.5 3.0 1.1   % Iris-versicolor
               5.7 2.8 4.1 1.3   % Iris-versicolor 
               6.0 3.0 4.8 1.8   % Iris-verginica
               6.9 3.1 5.4 2.1   % Iris-verginica
               6.7 3.1 5.6 2.4   % Iris-verginica
               6.9 3.1 5.1 2.3   % Iris-verginica
               5.8 2.7 5.1 1.9   % Iris-verginica
               6.8 3.2 5.9 2.3   % Iris-verginica
               6.7 3.3 5.7 2.5   % Iris-verginica
               6.7 3.0 5.2 2.3   % Iris-verginica
               6.3 2.5 5.0 1.9   % Iris-verginica
               6.5 3.0 5.2 2.0]; % Iris-verginica
%codificação das classes
setosa = [0 0 1]';
versicolor = [0 1 0]';
verginica = [1 0 0]';

% define targets
targets = [repmat(setosa,1,40) repmat(versicolor,1,40) repmat(verginica,1,40)];
targets = targets';

%parametros
alpha = 0.07; % taxa de aprendizagem 
epochs = 1000; %epocas
I = 4; %entradas
H = 15; %neuronios ocultos
O = 3; %saídas

%matrizes aleatorias pro peso
wH = randn(I,H);
wH2 = randn(H,O);

%normalizadas
wH = wH./sqrt(I);
wH2 = wH2./sqrt(H);

error = zeros(epochs,3);

tic;
disp('Training the Neural Network...');

for t = 1:epochs   
    X=['Epoch#',num2str(t)];
    disp(X);
    for k = 1:length(TrainData)
            % FORWARD PASS
            layer1 = logsig(TrainData(k,:)*wH); 
            output = logsig(layer1*wH2);
            
            % ERROR CALCULATION
            erro = targets(k,:) - output;
            
            % BACKWARD PASS
            dwH2=layer1'*(erro.*(output.*(1-output))); 
            dwH=((TrainData(k,:)'*(erro.*(output.*(1-output))))*wH2').*(layer1.*(1-layer1));
            
            % weights updation
            wH2=wH2+alpha.*dwH2;
            wH=wH+alpha.*dwH;  
            
    end
    
    error(t,:)= erro;
end

error = error';
sse=sum((error(:,1:epochs).^2),1);

plot(sse);
title('error square plot for training');
xlabel('no of iterations');
ylabel('error.^2');
toc;

%%%%%%%testing%%%%%%%%%

out = [];

 for i=1:size(TestData,1)
     % FORWARD PASS
            layer1 = logsig(TestData(i,:)*wH); %hidden layer
            output = logsig(layer1*wH2);
            
            out(i,:) = output;
            class(1,i) = find(out(i,:)>0.5);
                % Predicting the values for each TestData input            
 end
 
%targets2 = [repmat(3,1,10) repmat(2,1,10) repmat(1,1,10)];
%C = confusionmat(targets2,class);
%cm = confusionchart(C);
%cm.Title = 'Classificação';