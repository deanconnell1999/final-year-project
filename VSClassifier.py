'''
Author: Dean Connell
Description: Converts data in folders into a dataframe
and tests various classifiers to determine which is most
accurate in classifying the data for VS task.
'''

import numpy as np
import matplotlib.pyplot as plot
from scipy.io import loadmat
from scipy.io.matlab.miobase import MatReadWarning
from sklearn.cross_decomposition import PLSRegression
from statistics import mean
from sklearn.linear_model import LogisticRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.cluster import KMeans
from sklearn.model_selection import train_test_split, LeaveOneOut
from sklearn.svm import SVC
from sklearn.metrics import confusion_matrix, accuracy_score
from sklearn.neighbors import (NeighborhoodComponentsAnalysis, KNeighborsClassifier)
from sklearn.pipeline import Pipeline
from sklearn.neural_network import MLPClassifier
import warnings

warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
warnings.filterwarnings("ignore", category=MatReadWarning)


# splitting the frequency array into n segments
def segmentFrequencyBand(freqarr, n):
    average = len(freqarr) / float(n)
    output = []
    final = 0.0
    while final < len(freqarr):
        output.append([int(final), int(final + average)])
        final = final + average
    return output


# need to define our frequency bands and extract alpha band.
freqs = [0]
i = 1
while i < 280:
    freqs.append(freqs[-1]+0.25)
    i = i+1

freqbands = segmentFrequencyBand(freqs, 4)
alphaBand = [freqbands[0][0], freqbands[0][1]]

# defining powermap - names of struct elements
powermap = ["flcmap", "lcmap", "rcmap", "frcmap", "acmap"]
respnames = ["FLC", "LC", "RC", "FRC", "AC"]


# store periodogram data
def storeMetricData(fileArr, path):
    mapArr = []
    temp = np.zeros(5)
    resp = 0
    for i in range(len(fileArr)):
        filein = loadmat(path + str(fileArr[i]) + ".mat").get('patient_erp')
        while resp < 5:
            response = filein[powermap[resp]]
            response = response[0][0][alphaBand[0]:alphaBand[1]]
            temp[resp] = sum(response)/len(response)
            resp = resp+1
        mapArr.append(temp)
        temp = np.zeros(5)
        resp = 0
    return mapArr


# numbering files so we can search for them in folders specified
controlFiles = [15, 20, 22]
strokeFiles = [53, 54, 55, 56, 58, 59, 61]

# importing the files necessary for analysis
P7ConT1 = storeMetricData(controlFiles, "ERP Files\\VS\\Controls\\T1\\P7\\")
P8ConT1 = storeMetricData(controlFiles, "ERP Files\\VS\\Controls\\T1\\P8\\")
P7StrT1 = storeMetricData(strokeFiles, "ERP Files\\VS\\Patients\\T1\\P7\\")
P8StrT1 = storeMetricData(strokeFiles, "ERP Files\\VS\\Patients\\T1\\P8\\")

P7ConT2 = storeMetricData(controlFiles, "ERP Files\\VS\\Controls\\T2\\P7\\")
P8ConT2 = storeMetricData(controlFiles, "ERP Files\\VS\\Controls\\T2\\P8\\")
P7StrT2 = storeMetricData(strokeFiles, "ERP Files\\VS\\Patients\\T2\\P7\\")
P8StrT2 = storeMetricData(strokeFiles, "ERP Files\\VS\\Patients\\T2\\P8\\")


# assignment to values 1 (stroke) and 0 (control):
def createDataframe(elec1, elec2, id):
    if id == 0:
        dataframe = np.array([[elec1], [elec2], [np.zeros(len(elec1))]])
    else:
        dataframe = np.array([[elec1], [elec2], [np.ones(len(elec1))]])
    return dataframe


# creating stroke and control dataframes
controlT1 = createDataframe(P7ConT1, P8ConT1, 0)
strokeT1 = createDataframe(P7StrT1, P8StrT1, 1)
controlT2 = createDataframe(P7ConT2, P8ConT2, 0)
strokeT2 = createDataframe(P7StrT2, P8StrT2, 1)

# defining training data in useful format
x1 = np.array([np.concatenate((controlT1[0, 0], strokeT1[0, 0]), axis=0),
               np.concatenate((controlT1[1, 0], strokeT1[1, 0]), axis=0)])
y1 = np.concatenate((controlT1[2, 0], strokeT1[2, 0]), axis=0)
x2 = np.array([np.concatenate((controlT2[0, 0], strokeT2[0, 0]), axis=0),
               np.concatenate((controlT2[1, 0], strokeT2[1, 0]), axis=0)])
y2 = np.concatenate((controlT2[2, 0], strokeT2[2, 0]), axis=0)
y1 = y1.astype(int)
y2 = y2.astype(int)

# print(x1) # P7, P8 data
# print(y1) # 0 or 1 stroke or control


# function that takes training data as array and extracts response data for response at index n
def extractResponseData(trainingData, n):
    respData = np.zeros((len(trainingData[0]), 2))
    rowcount = 0
    colcount = 0
    for row in trainingData:
        for element in row:
            respData[rowcount][colcount] = element[n]
            rowcount = rowcount+1
        rowcount = 0
        colcount = colcount+1
    return respData


# response types are turned into training sets
FLCT1 = extractResponseData(x1, 0)
LCT1 = extractResponseData(x1, 1)
RCT1 = extractResponseData(x1, 2)
FRCT1 = extractResponseData(x1, 3)
ACT1 = extractResponseData(x1, 4)

FLCT2 = extractResponseData(x2, 0)
LCT2 = extractResponseData(x2, 1)
RCT2 = extractResponseData(x2, 2)
FRCT2 = extractResponseData(x2, 3)
ACT2 = extractResponseData(x2, 4)


# plots the Power Spectral Densities of one eletrode against another for both groups.
def PSDPlots(StrP7, StrP8, ConP7, ConP8, StrFiles, ConFiles):
    n = 0
    while n < 5:
        plot.subplot(2, 3, n+1)
        i = 0
        while i < len(StrP8):
            plot.scatter(StrP7[i][n], StrP8[i][n], color='red')
            plot.annotate(str(StrFiles[i]), (StrP7[i][n], StrP8[i][n]))
            i = i+1
        i = 0
        while i < len(ConP8):
            plot.scatter(ConP7[i][n], ConP8[i][n], color='blue')
            plot.annotate(str(ConFiles[i]), (ConP7[i][n], ConP8[i][n]))
            i = i+1
        plot.title("P7 vs. P8 Power (" + respnames[n] + ")")
        n = n+1
    plot.xlabel("P7 Power (dB/Hz)")
    plot.ylabel("P8 Power (dB/Hz)")
    plot.show()

# Generating scatterplots for T1 and T2
# PSDPlots(P7StrT1, P8StrT1, P7ConT1, P8ConT1, strokeFiles, controlFiles)
# PSDPlots(P7StrT2, P8StrT2, P7ConT2, P8ConT2, strokeFiles, controlFiles)


# Method 1: PLS Discriminant Analysis
def PLSDisAnalysis(xTrain, yTrain, xTest, yTest):
    print("Partial Least Squares Discriminant Analysis:")
    pls = PLSRegression(n_components=2)
    pls.fit(xTrain, yTrain)
    plspred = (pls.predict(xTest)[:, 0] > 0.5).astype('uint8')
    confusion = confusion_matrix(yTest, plspred)
    accuracy = accuracy_score(yTest, plspred)
    print("Predictions: ", plspred)
    print("Actual:      ", yTest)
    # print("Confusion Matrix:\n", confusion)
    print("Accuracy:    ", accuracy)
    print()
    return accuracy


# Method 2: K-Means Clustering
def kMeansClassification(xTrain, yTrain, xTest, yTest):
    print("k-means Classifier:")
    kmeans = KMeans(n_clusters=2, random_state=0)
    kmeans.fit(xTrain, yTrain)
    kmpred = kmeans.predict(xTest)
    confusion = confusion_matrix(yTest, kmpred)
    accuracy = accuracy_score(yTest, kmpred)
    print("Predictions: ", kmpred)
    print("Actual:      ", yTest)
    # print("Confusion Matrix:\n", confusion)
    print("Accuracy:    ", accuracy)
    print()
    return accuracy


# Method 3: SVC
def SupportVectorClassifier(xTrain, yTrain, xTest, yTest):
    print("Support Vector Classification:")
    svc = SVC(kernel='poly', gamma='auto')
    svc.fit(xTrain, yTrain)
    svcpred = svc.predict(xTest)
    confusion = confusion_matrix(yTest, svcpred)
    accuracy = accuracy_score(yTest, svcpred)
    print("Predictions: ", svcpred)
    print("Actual:      ", yTest)
    # print("Confusion Matrix:\n", confusion)
    print("Accuracy:    ", accuracy)
    print()
    return accuracy


# Method 4: Linear Discriminant Analysis
def LinearDisAnalysis(xTrain, yTrain, xTest, yTest):
    print("Linear Discriminant Analysis:")
    linreg = LinearDiscriminantAnalysis()
    linreg.fit(xTrain, yTrain)
    linpred = linreg.predict(xTest)
    confusion = confusion_matrix(yTest, linpred)
    accuracy = accuracy_score(yTest, linpred)
    print("Predictions: ", linpred)
    print("Actual:      ", yTest)
    # print("Confusion Matrix:\n", confusion)
    print("Accuracy:    ", accuracy)
    print()
    return accuracy


# Method 5: Logistic Regression
def LogisticClassifier(xTrain, yTrain, xTest, yTest):
    print("Logistic Regression:")
    logreg = LogisticRegression(random_state=0)
    logreg.fit(xTrain, yTrain)
    logpred = logreg.predict(xTest)
    confusion = confusion_matrix(yTest, logpred)
    accuracy = accuracy_score(yTest, logpred)
    print("Predictions: ", logpred)
    print("Actual:      ", yTest)
    # print("Confusion Matrix:\n", confusion)
    print("Accuracy:    ", accuracy)
    print()
    return accuracy


# Method 6: K-Nearest Neighbours
def KNN(xTrain, yTrain, xTest, yTest):
    print("k Nearest Neighbours Classification:")
    knn = KNeighborsClassifier(n_neighbors=2)
    knn.fit(xTrain, yTrain)
    knnpred = knn.predict(xTest)
    confusion = confusion_matrix(yTest, knnpred)
    accuracy = accuracy_score(yTest, knnpred)
    print("Predictions: ", knnpred)
    print("Actual:      ", yTest)
    # print("Confusion Matrix:\n", confusion)
    print("Accuracy:    ", accuracy)
    print()
    return accuracy


# Method 7: NeighbourhoodComponentsAnalysis, KNN Classifier
def NCAKNN(xTrain, yTrain, xTest, yTest):
    print("Neighbourhood Component Analysis with kNN Classification")
    nca = NeighborhoodComponentsAnalysis(random_state=42)
    knn = KNeighborsClassifier(n_neighbors=2)
    ncapipe = Pipeline([('nca', nca), ('knn', knn)])
    ncapipe.fit(xTrain, yTrain)
    ncapipepred = ncapipe.predict(xTest)
    confusion = confusion_matrix(yTest, ncapipepred)
    accuracy = accuracy_score(yTest, ncapipepred)
    print("Predictions: ", ncapipepred)
    print("Actual:      ", yTest)
    # print("Confusion Matrix:\n", confusion)
    print("Accuracy:    ", accuracy)
    print()
    return accuracy


# Method 8: Neural Networks
def NeuralNetwork(xTrain, yTrain, xTest, yTest):
    print("Neural Network (MLP):")
    nn = MLPClassifier(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=(5, 2), random_state=1)
    nn.fit(xTrain, yTrain)
    nnpred = nn.predict(xTest)
    confusion = confusion_matrix(yTest, nnpred)
    accuracy = accuracy_score(yTest, nnpred)
    print("Predictions: ", nnpred)
    print("Actual:      ", yTest)
    # print("Confusion Matrix:\n", confusion)
    print("Accuracy:    ", accuracy)
    print()
    return accuracy


def LeaveOneOutT1(tasks):
    PLS = []
    SVC = []
    KMC = []
    LiC = []
    LoC = []
    Knn = []
    NCA = []
    NN = []
    for i in range(len(tasks)):
        task = tasks[i]
        loo = LeaveOneOut()
        loo.get_n_splits(task, y1)
        for train_index, test_index in loo.split(task, y1):
            x1Train, x1Test = task[train_index], task[test_index]
            y1Train, y1Test = y1[train_index], y1[test_index]
            # calling classification functions
            PLS.append(PLSDisAnalysis(x1Train, y1Train, x1Test, y1Test))
            SVC.append(SupportVectorClassifier(x1Train, y1Train, x1Test, y1Test))
            KMC.append(kMeansClassification(x1Train, y1Train, x1Test, y1Test))
            LiC.append(LinearDisAnalysis(x1Train, y1Train, x1Test, y1Test))
            LoC.append(LogisticClassifier(x1Train, y1Train, x1Test, y1Test))
            Knn.append(KNN(x1Train, y1Train, x1Test, y1Test))
            NCA.append(NCAKNN(x1Train, y1Train, x1Test, y1Test))
            NN.append(NeuralNetwork(x1Train, y1Train, x1Test, y1Test))
        print("Leave-One-Out Accuracies:")
        print("PLS Discriminant Analysis:    ", mean(PLS))
        print("Support Vector Classifier:    ", mean(SVC))
        print("k-means Classifier:           ", mean(KMC))
        print("Linear Discriminant Analysis: ", mean(LiC))
        print("Logistic Classifier:          ", mean(LoC))
        print("k Nearest Neighbours:         ", mean(Knn))
        print("NCA with kNN:                 ", mean(NCA))
        print("Neural Network (MLP):         ", mean(NN))


# note: if calling this function, it may be beneficial
# to print confusion matrices in each classification function.
def TrainTestSplitT2(tasks, responses):
    PLS = []
    SVC = []
    KMC = []
    LiC = []
    LoC = []
    Knn = []
    NCA = []
    NN = []
    # making test and training data
    for iterations in range(30):
        # 30 iterations for accuracies to converge
        for i in range(len(tasks)):
            x1Train, x1Test, y1Train, y1Test = train_test_split(tasks[i], y1, test_size=0.1)
            y1Train = y1Train.astype(int)
            y1Test = y1Test.astype(int)
            # Calling classification functions:
            PLS.append(PLSDisAnalysis(x1Train, y1Train, responses[i], y2))
            SVC.append(SupportVectorClassifier(x1Train, y1Train, responses[i], y2))
            KMC.append(kMeansClassification(x1Train, y1Train, responses[i], y2))
            LiC.append(LinearDisAnalysis(x1Train, y1Train, responses[i], y2))
            LoC.append(LogisticClassifier(x1Train, y1Train, responses[i], y2))
            Knn.append(KNN(x1Train, y1Train, responses[i], y2))
            NCA.append(NCAKNN(x1Train, y1Train, responses[i], y2))
            NN.append(NeuralNetwork(x1Train, y1Train, responses[i], y2))
    print("Train Test Split Accuracies:")
    print("PLS Discriminant Analysis:    ", mean(PLS))
    print("Support Vector Classifier:    ", mean(SVC))
    print("k-means:                      ", mean(KMC))
    print("Linear Discriminant Analysis: ", mean(LiC))
    print("Logistic Classifier:          ", mean(LoC))
    print("k Nearest Neighbours:         ", mean(Knn))
    print("NCA and kNN:                  ", mean(NCA))
    print("Neural Network (MLP):         ", mean(NN))


tasks = [FLCT1, LCT1, RCT1, FRCT1, ACT1]
responses = [FLCT2, LCT2, RCT2, FRCT2, ACT2]

# TrainTestSplitT2(tasks, responses)
LeaveOneOutT1(tasks)
