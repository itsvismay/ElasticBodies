import numpy as np
import csv
#-----------------------------------------------
def sphere(x):
    return np.sum(x**2)
#-----------------------------------------------
def fitQuadSurface(X, Y):
    n = X.shape[1]
    minNumSamples = (n+1)*(n+2)/2 # minimum number of samples for quadratic fit
    Xnew = np.zeros((m,minNumSamples))
    Xnew[:,0] = np.ones(m)    
    print Xnew.shape
    # Fill in Xnew matrix:
    for i in range(n):
        Xnew[:,i+1] = X[:,i]
        Xnew[:,i+6] = X[:,i]**2
    
    count = 11    
    for i in range(n-1):
        for k in range(i+1,n):
            Xnew[:,count] = X[:,i] * X[:,k]
            count += 1
        
    coef = np.linalg.pinv(np.dot(Xnew.T, Xnew))
    coef = np.dot(coef,Xnew.T)
    coef = np.dot(coef,Y)
    return coef
#-----------------------------------------------
def predict(x, coef):
    x_new = np.zeros(coef.shape)
    x_new[0] = 1.0
    n = len(x)
    
    for i in range(n):
        x_new[i+1] = x[i]
        x_new[i+6] = x[i]**2
    count = 1+2*n  
    for i in range(n-1):
        for k in range(i+1,n):
            x_new[count] = x[i] * x[k]
            count += 1
            
    pred = np.sum(coef * x_new)
    return pred
#-----------------------------------------------   

fileID = open('OptimalLHS_05D_seed01.txt','rb')
reader = csv.reader(fileID, delimiter=',', quoting=csv.QUOTE_NONE)
data = list(reader)
X = np.array(data).astype('float')
m, n = X.shape
print X.shape
Y = map(sphere, X)

coef = fitQuadSurface(X, Y)
x = np.zeros(n) # unknown location
y_hat = predict(x, coef) # Predict the function value at an unknown location



 
    
    
    
    
    