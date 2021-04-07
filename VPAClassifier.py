'''
Author: Dean Connell
Description: Converts data in folders into a dataframe
and tests various classifiers to determine which is most
accurate in classifying the data for VPA task.
'''

import numpy as np
import matplotlib.pyplot as plot
from scipy.io import loadmat
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
    freqs.append(freqs[-1] + 0.25)
    i = i + 1

freqbands = segmentFrequencyBand(freqs, 4)
alphaBand = [freqbands[0][0], freqbands[0][1]]

# defining powermap - names of struct elements
powermap = ["fccpow", "ficpow", "tccpow", "ticpow"]
respnames = ["FCC", "FIC", "TCC", "TIC"]


# store periodogram data
def storeMetricData(fileArr, path):
    mapArr = []
    temp = np.zeros(4)
    resp = 0
    for i in range(len(fileArr)):
        filein = loadmat(path + str(fileArr[i]) + ".mat").get('patient_erp')
        while resp < 4:
            response = filein[powermap[resp]]
            response = response[0][0][alphaBand[0]:alphaBand[1]]
            temp[resp] = sum(response) / len(response)
            resp = resp + 1
        mapArr.append(temp)
        temp = np.zeros(4)
        resp = 0
    return mapArr


# numbering files so we can search for them in folders specified
controlFiles = [13, 14, 15, 16, 17, 21]
strokeFiles = [53, 54, 55, 56, 57, 58, 59, 61]
controlFilesT2 = [13, 14, 15, 21]

# importing the files necessary for analysis
PO3ConT1 = storeMetricData(controlFiles, "ERP Files\\VPA\\Controls\\T1\\PO3\\")
PO4ConT1 = storeMetricData(controlFiles, "ERP Files\\VPA\\Controls\\T1\\PO4\\")
PO3StrT1 = storeMetricData(strokeFiles, "ERP Files\\VPA\\Patients\\T1\\PO3\\")
PO4StrT1 = storeMetricData(strokeFiles, "ERP Files\\VPA\\Patients\\T1\\PO4\\")

PO3ConT2 = storeMetricData(controlFilesT2, "ERP Files\\VPA\\Controls\\T2\\PO3\\")
PO4ConT2 = storeMetricData(controlFilesT2, "ERP Files\\VPA\\Controls\\T2\\PO4\\")
PO3StrT2 = storeMetricData(strokeFiles, "ERP Files\\VPA\\Patients\\T2\\PO3\\")
PO4StrT2 = storeMetricData(strokeFiles, "ERP Files\\VPA\\Patients\\T2\\PO4\\")


# assignment to values 1 (stroke) and 0 (control):
def createDataframe(elec1, elec2, id):
    if id == 0:
        dataframe = np.array([[elec1], [elec2], [np.zeros(len(elec1))]])
    else:
        dataframe = np.array([[elec1], [elec2], [np.ones(len(elec1))]])
    return dataframe


# creating stroke and control dataframes
controlT1 = createDataframe(PO3ConT1, PO4ConT1, 0)
strokeT1 = createDataframe(PO3StrT1, PO4StrT1, 1)
controlT2 = createDataframe(PO3ConT2, PO4ConT2, 0)
strokeT2 = createDataframe(PO3StrT2, PO4StrT2, 1)

# defining training data in useful format
x1 = np.array([np.concatenate((controlT1[0, 0], strokeT1[0, 0]), axis=0),
               np.concatenate((controlT1[1, 0], strokeT1[1, 0]), axis=0)])
y1 = np.concatenate((controlT1[2, 0], strokeT1[2, 0]), axis=0)
x2 = np.array([np.concatenate((controlT2[0, 0], strokeT2[0, 0]), axis=0),
               np.concatenate((controlT2[1, 0], strokeT2[1, 0]), axis=0)])
y2 = np.concatenate((controlT2[2, 0], strokeT2[2, 0]), axis=0)
y1 = y1.astype(int)
y2 = y2.astype(int)

# print(x1) # PO3, PO4 data
# print(y1) # 0 or 1 stroke or control


# function that takes training data as array and extracts response data for response at index n
def extractResponseData(trainingData, n):
    respData = np.zeros((len(trainingData[0]), 2))
    rowcount = 0
    colcount = 0
    for row in trainingData:
        for element in row:
            respData[rowcount][colcount] = element[n]
            rowcount = rowcount + 1
        rowcount = 0
        colcount = colcount + 1
    return respData


# response types are turned into training sets
FCCT1 = extractResponseData(x1, 0)
FICT1 = extractResponseData(x1, 1)
TCCT1 = extractResponseData(x1, 2)
TICT1 = extractResponseData(x1, 3)

FCCT2 = extractResponseData(x2, 0)
FICT2 = extractResponseData(x2, 1)
TCCT2 = extractResponseData(x2, 2)
TICT2 = extractResponseData(x2, 3)


# plots the Power Spectral Densities of one eletrode against another for both groups.
def PSDPlots(StrPO3, StrPO4, ConPO3, ConPO4, StrFiles, ConFiles):
    n = 0
    while n < 4:
        plot.subplot(2, 2, n+1)
        i = 0
        while i < len(StrPO4):
            plot.scatter(StrPO3[i][n], StrPO4[i][n], color='red')
            plot.annotate(str(StrFiles[i]), (StrPO3[i][n], StrPO4[i][n]))
            i = i+1
        i = 0
        while i < len(ConPO4):
            plot.scatter(ConPO3[i][n], ConPO4[i][n], color='blue')
            plot.annotate(str(ConFiles[i]), (ConPO3[i][n], ConPO4[i][n]))
            i = i+1
        plot.title("PO3 vs. PO4 Power (" + respnames[n] + ")")
        n = n+1
    plot.xlabel("PO3 Power (dB/Hz)")
    plot.ylabel("PO4 Power (dB/Hz)")
    plot.show()

# Generating scatterplots for T1 and T2
# PSDPlots(PO3StrT1, PO4StrT1, PO3ConT1, PO4ConT1, strokeFiles, controlFiles)
# PSDPlots(PO3StrT2, PO4StrT2, PO3ConT2, PO4ConT2, strokeFiles, controlFiles)


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
    print("Neighbourhood Component Analysis withs kNN Classification")
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
        # print("k-means Classifier:           ", mean(KMC))
        print("Linear Discriminant Analysis: ", mean(LiC))
        print("Logistic Classifier:          ", mean(LoC))
        # print("k Nearest Neighbours:         ", mean(Knn))
        # print("NCA with kNN:                 ", mean(NCA))
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


tasks = [FCCT1, FICT1, TCCT1, TICT1]
responses = [FCCT2, FICT2, TCCT2, TICT2]

# TrainTestSplitT2(tasks, responses)
LeaveOneOutT1(tasks)
