import csv
import numpy as np
from matplotlib import pyplot as plt

def readhistorical():
    '''读取历史数据函数'''
    with open(r'/home/xitefei/group work orriginal code/anasazi-model-group4/data/target_data.csv') as f1:   #读取文件，''中是我电脑的数据路径，你们验证时要改成自己的
        reader = csv.reader(f1)                                        #创建reader
        historicaldata = []                                            #创建空列表，用于存放历史数据
        for row in reader:
            data1 = float(row[1])                                      #遍历csv文件中的每一行，将其中字符串数据转换为浮点型
            historicaldata.append(data1)                               #将数据存入列表
    return historicaldata

def readsimulate():
    '''读取模拟数据函数'''
    with open(r'/home/xitefei/group work orriginal code/anasazi-model-group4/NumberOfHousehold.csv') as f2:
        reader = csv.reader(f2)
        header_row=next(reader)
        simulatedata = []
        for row in reader:
            data2 = float(row[1])
            simulatedata.append(data2)
    return simulatedata

def readmaxCapacity():
    '''读取maxCapacity'''
    with open(r'/home/xitefei/group work orriginal code/anasazi-model-group4/NumberOfHousehold.csv') as f3:
        reader = csv.reader(f3)
        header_row=next(reader)
        maxCapacity = []
        for row in reader:
            data3 = float(row[2])
            maxCapacity.append(data3)
    return maxCapacity

def readyear():
    '''读取年份'''
    with open(r'/home/xitefei/group work orriginal code/anasazi-model-group4/data/target_data.csv') as f4:
        reader = csv.reader(f4)
        year = []
        for row in reader:
            data4 = float(row[0])
            year.append(data4)
    return year


def squaredifference(a,b):
    '''计算残差平方和'''
    v1 = np.array(a)                                                    #转换为numpy数组形式
    v2 = np.array(b)
    sum = ((v1-v2)**2).sum()
    return sum

h = readhistorical()
s = readsimulate()
y = readyear()
m = readmaxCapacity()
sum = squaredifference(h,s)
print(h)
print(s)
print(m)
print(sum)

plt.figure(figsize=(8,4))
plt.title('Anasazi Model')
plt.xlabel('year')
plt.ylabel('population')
plt.plot(y,h,c='red')
plt.plot(y,s,c='blue')
plt.plot(y,m,c='green')
plt.show()
