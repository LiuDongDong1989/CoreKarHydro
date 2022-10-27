# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 19:56:51 2021
@author: liudongdong
一、帕累托最优参数集
程序作用：实现multi-objective complex evolution (MOCOM-UA) global optimization method
输出结果：帕累托最优参数集、参数敏感性分析
模块1：参数空间初始化
模块2：计算目标函数
模块3：根据目标函数值，进行帕累托排序
模块4：基于帕累托分类的参数选择
模块5：多目标单纯形生成，多目标下山法搜索
二、参数敏感性分析
程序作用：实现RBD-FAST - Random Balance Designs Fourier Amplitude Sensitivity Test
输出结果：全局敏感性
三、编程日志
2021.1.8：完成框架搭建,完成拉丁超立方抽样代码
2021.1.9：完成模块2和模块3
2021.1.20:
    1）完成模块4和模块5, 准备debug和精简程序
    2）修改pre_error.py
    3）try模块
2021.1.24:修复了若干BUG,排序代码疑似出错
2021.1.30:
    1）修复is_dominate的逻辑错误
    2）修改cal_objectives_single代码
2021.2.20:
    1）改进后处理程序
    2）删除旧的LH方法，采用SAlib包的拉丁超立方采样(弃用，改回原来的方法)
    3）改进参数输入，用字典来表示
2021.2.22：
    1）更新参数敏感性分析代码
    2）修改cal_objectives函数
    3）统一代码中的变量名称
2021.2.25：
    1）改为三目标函数
2022.4.3:
    1)改回两目标函数
    2）增加了通用性
"""

import numpy as np
import pre_Error

class ParetoSet(object): 
        
    def __init__(self, problem): 
    #:param num_vars:参数个数
    #:param bounds:参数对应范围（list）
    #:param objNum:目标函数个数
    #:param population:种群数量
    #:param func:正演程序
        self.problem=problem
        self.num_vars=problem['num_vars']
        self.bounds=problem['bounds']
        self.objNum=problem['objNum']
        self.population=problem['population']
        self.func=problem['func']
        self.outputs=problem['outputs']
        self.treatment=problem['treatment']
                
    #Space=ParameterSpace+ObjectiveSpace+RankingSpace+SelectionProbability
        self.ParameterSpace=self.LHSample()
        self.ObjectiveSpace=np.zeros([self.population, self.objNum])
        self.RankingSpace=np.zeros([self.population, 1])
        self.SelectionProbability=np.zeros([self.population, 1])
        self.Space=np.column_stack((self.ParameterSpace, self.ObjectiveSpace, self.RankingSpace, self.SelectionProbability)) 
        self.MaxRank=0  
        
    def LHSample(self):
    #拉丁超立方抽样          
        result = np.empty([self.population, self.num_vars])
        temp = np.empty([self.population])
        d = 1.0 / self.population
        for i in range(self.num_vars):   
            for j in range(self.population):
                temp[j] = np.random.uniform(low=j * d, high=(j + 1) * d, size = 1)[0]
            np.random.shuffle(temp)            
            for j in range(self.population):
                result[j, i] = temp[j]
        
        #对样本数据进行拉伸
        b = np.array(self.bounds)
        lower_bounds = b[:, 0]
        upper_bounds = b[:, 1]
        if np.any(lower_bounds > upper_bounds):
            print ("范围出错")
            return None
        
        #sample * (upper_bound - lower_bound) + lower_bound
        np.add(np.multiply(result, (upper_bounds - lower_bounds), out=result), lower_bounds, out=result)
        
        return result             
        
    def value_valid(self, a):
    #判断取值时否在范围内  
        for i in range(self.num_vars):
            if (a[i]<self.bounds[i][0]) or (a[i]>self.bounds[i][1]): 
                flag=False
                print('反射点取值不在范围内')  
                break
            else:
                flag=True  
        return True if flag==True else False
                  
    def cal_objectives_single(self, parameters):  
    #计算单个参数空间的目标函数
        try:
            data=self.func(parameters, treatment=self.treatment)             
            print('求解成功，参数为{}'.format(parameters) )
            #主要采用Hydro_Err方法计算，sklearn_Err方法辅助
            flag=self.outputs
            if flag=='Test_func':              
                Error=pre_Error.Hydro_Err(np.array([3.01,5.01,7.04,9.09]), pred=data)
                return [Error['rmse'], Error['nse']], parameters
            elif flag=='VadoseZone':
                Error=pre_Error.Hydro_Err(true=data[0], pred=data[1]) 
                return [Error['rmse'], Error['nse']], parameters
        except RuntimeError:
            print('求解出错，参数为{}'.format(parameters) )
            return [float('inf'), -float('inf')], parameters
                                                    
    def cal_objectives(self):  
    #计算全部参数空间的目标函数
        for i, element in enumerate(self.ParameterSpace):
            self.ObjectiveSpace[i, :], self.ParameterSpace[i] =self.cal_objectives_single(element)                 
                          
    def is_dominate(self, a, b):
    #3个目标函数的情况
        '''
        判定a是否被b强支配（强支配区别于弱支配）
        '''
        flag= a[0]>b[0] and a[1]<b[1]
        return True if flag else False

    def Mul_is_dominate(self, a, b): 
    #多个目标函数的判断（更新中）
        '''
        判定a是否被b强支配
        a和b分别是字典
        '''
        flag=True #True:a被b支配
        for key in range(list(a)):
            if a[key] > b[key]: #b在该项上优于a， >表示离最优值的距离较远
                flag=True
            else:
                flag=False  
                break
        return True if flag else False
        
    def ParetoRanking(self):
    #帕累托排序
        print('执行ParetoRanking')
        self.MaxRank, index, NewRange=0, [], list(range(self.population))
        while len(NewRange)>0:
            self.MaxRank+=1
            for i in NewRange:
                for j in NewRange :
                    if self.is_dominate(self.ObjectiveSpace[i,:], self.ObjectiveSpace[j,:]):
                        flag=False #不能排在前面
                        break
                    else: 
                        flag=True #能排序在前面
                if flag==True:    #能排序在前面
                    self.RankingSpace[i]=self.MaxRank   
                    index.append(i)            
            [NewRange.remove(k) for k in index]
            index=[]           
    
    def cal_SelectionProbability(self):
    #基于帕累托分类结果的参数选择概率
        print('执行cal_SelectionProbability')
        for i in range(self.population):
            self.SelectionProbability[i,0] =(self.MaxRank - self.RankingSpace[i,0] + 1 )/( (self.MaxRank + 1)*self.population - sum( self.RankingSpace[:,0] ) )
        self.Space=np.column_stack( (self.ParameterSpace, self.ObjectiveSpace, self.RankingSpace, self.SelectionProbability) ) 
        self.Space=self.Space[self.Space[:,-2].argsort()]#按Ranking列整体排序
        
    def ComplexEvolution(self):   
    #单纯形个数=最差的点的个数
    #每个单纯形包含（self.num_vars+1)个点，self.num_vars是参数个数，self.num_vars个点取于“较好的点”，1个点取于“最差的点 ”
        print('执行ComplexEvolution')
        self.ComplexNumber=np.sum(self.Space[:,-2] == self.MaxRank )
        self.BetterSpace=self.Space[: -self.ComplexNumber, :]
        self.WorstSpace=self.Space[-self.ComplexNumber : , :]  
     #多目标单纯形生成、下山法搜索      
        Probability_normal=self.BetterSpace[:,-1]/self.BetterSpace[:,-1].sum()#概率归一化操作
        for i in range(self.ComplexNumber):
            #基于概率抽样
            index_better=np.random.choice(self.BetterSpace.shape[0], size=self.num_vars, p=Probability_normal)
            S_g=np.mean( self.BetterSpace[index_better, :], axis=0 )
            S_w=self.WorstSpace[i, :]
            
            #产生反射点和约束点参数空间
            self.ReflectionSpace=(lambda γ : γ*S_g+(1-γ)*S_w ) (2) 
            self.ContrationSpace=(lambda γ : γ*S_g+(1-γ)*S_w ) (0.5)
            #print('第{}个单纯形的反射点是{}|约束点是{}'.format(i, self.ReflectionSpace, self.ContrationSpace))
            
            #判定反射点是否被拒绝（接受：仅当在区间范围内且是不支配解的时候）   
            if self.value_valid(self.ReflectionSpace)==True:
                self.ReflectionSpace[self.num_vars : self.num_vars+self.objNum ], self.ReflectionSpace[0: self.num_vars] = self.cal_objectives_single( self.ReflectionSpace[0: self.num_vars] ) 
                for j in index_better:
                    if self.is_dominate(self.ReflectionSpace[ self.num_vars : self.num_vars+self.objNum], self.BetterSpace [ j, self.num_vars : self.num_vars+self.objNum ] ):
                        flag=False
                        break
                    else:
                        flag=True 
            else: flag=False
                        
            #真：用反射点代替“最差的点” |假：用约束点代替“最差的点” 
            if flag==True:
                #print('真：用反射点代替“最差的点”')
                self.WorstSpace[i,:]=self.ReflectionSpace 
                
            else:
                #print('假：用约束点代替“最差的点”')
                self.ContrationSpace[self.num_vars : self.num_vars+self.objNum], self.ContrationSpace[0: self.num_vars] = self.cal_objectives_single( self.ContrationSpace[0: self.num_vars] ) 
                self.WorstSpace[i,:]=self.ContrationSpace   
            
        #更新下一代种群
        self.Space=np.row_stack( (self.BetterSpace, self.WorstSpace) )
        self.ObjectiveSpace= self.Space[ :, self.num_vars : self.num_vars+self.objNum ]
                            
    def MOCOM_UA(self):
    #反演程序
        flag=0
        self.cal_objectives()     
        self.ParetoRanking()
        while self.MaxRank!=1 and flag<100000:  
            flag+=1
            self.cal_SelectionProbability()
            self.ComplexEvolution() 
            self.ParetoRanking()
            self.cal_SelectionProbability()
            template='第{}代种群的目标函数:{},最大排序{}'.format(flag, self.Space[:,-(2+self.objNum):-1], self.MaxRank)
            print(template) 
            if (flag >= 2000) and (flag % 2000==0):
                template_space='第{}代种群的数据'.format(flag)
                np.savetxt(template_space, self.Space)
        np.savetxt("ParetoSet.txt", self.Space)
        
        
        
        

