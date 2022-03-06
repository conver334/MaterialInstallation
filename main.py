""" Function:
    The software is developed for the designing of the installation of photoelectrocatalytic building materials and flexible solar cells. Photoelectrocatalytic building materials are assembled on the outer surface of the building. Photoelectrocatalytic reaction occurs with solar light to convert carbon dioxide into methanol, releasing oxygen as well. While methanol is an economical chemical product. The work efficiency of the building materials is affected by solar illumination. Flexible solar cells are also assembled on the outer surface of the building, which can convert solar radiation into electric energy, whose efficiency is also affected by solar illumination.
    However, the exterior area of a building is always limited. Meanwhile, affected by architectural modeling, solar radiation intensity of different positions on building surface is usually different, which would affect the work efficiency of photoelectrocatalytic building materials and flexible solar cells. In order to plan the installation area of the two materials reasonably and obtain the best ecological and economic benefits, we developed this software.
        
-
What you need to know:
    In the final results, areas where flexible solar cells are installed would be marked with grid lines;
    Building surface you selected should be one or a set of surfaces in type of Mesh and entered by Mesh component in grasshopper;
    The default parameters of the component are the parameters suitable for a certain building in Guangzhou area. If you want to change the region or building, refer to the annotation of each parameter to enter corresponding parameters.
    The feasible latitude is between [-66.5,66.5]. This method is not applicable at high latitudes above the polar circle.
-
Provided by Simiao Zhang
            
        
    Inputs:
        _geometry:The surface of a building where flexible solar cells and arrays of functional units are placed. Geometry for which radiation analysis will be conducted. Geometry must be either a Mesh or a list of Meshes.
        H: total horizontal radiation. The default H is set to 11.375
        hd: Horizontal scattering radiation. The default hd is set to 7.08
        rou:Ground reflectivity
        isc:Solar constant. unit MJ/m^2/day. The default hd is set to 118.11
        latitude:The latitude of the weather file location, between [-66.5,66.5], This method is not applicable at high latitudes above the polar circle. unit degree. The default hd is set to 23.16
        lowlim: Building opening time,between[0,24].lowlim<uplim. default 7
        uplim: Building closing time,between[0,24]. default 18
        E:Annual electricity consumption per capita in the community. unit kwh .default 1400
        C:Per capita daily carbon emissions. unit carbon emission/ person .default 6.2
        T:T hours one day. default 24
        yita: Solar electrode efficiency [0,1]. default 0.08
        r0: Environmental influence factor of external power supply. unit 1/kwh.default 1
        N: Days in a year [365,366]. unit day. default 365
        p:Set electricity price per kilowatt hour. unit CNY/kwh. default 0.5
        q:Annual output value per square meter of phototropic photoanode material. unit CNY.default 95.87
        a:Environmental effect weight. 1>=a>=0, a+b=1. default 0.9
        b:Economic benefit weight. 1>=b>=0, a+b=1. default 0.1
        rp:Environmental impact factor of per capita carbon footprint. default 1.13
        I:Annual average solar radiation per square meter on the light plane. default 5700   
        k:The sunlight energy received per mole of carbon dioxide consumed. unit e/td. default 36.55,
        ratio: the area utilization rate of the photocatalytic material. unit percentage. default 0.1963 
        Gaussian_mu:In a day, the flow of people is normal distribution, the expectation of Gaussian function
        Gaussian_sigma:In a day, the flow of people is normal distribution, Variance of Gaussian function
        Gaussian_peak:In a day, the flow of people is normal distribution, The peak value of Gaussian function
    
    Output:
        point1: First cornor of grids. Connect this four points output to a "Mesh" grasshopper component to preview this output seperately from the others of this component. 
        point2: Second cornor of grids
        point3: Third cornor of grids
        point4: Forth cornor of grids
        meshNum: Number of faces selected for solar installation
        mesharea: Solar panel installation area.
                
"""

from ghpythonlib.componentbase import executingcomponent as component
import Grasshopper, GhPython
import System
import Rhino
import json
import csv
import rhinoscriptsyntax as rs
import Grasshopper.Kernel as gh
import Rhino as rc
import scriptcontext as sc
import System.Guid
import math
from math import pi
#import sys
#import json
#import os
#import System.Threading.Tasks as tasks
#import System
#import time
# from Grasshopper import DataTree
# from Grasshopper.Kernel.Data import GH_Path
# import rhino_tools
inf = 2147483647
__author__ = "ZHANG Simiao"
__version__ = "2021.11.10"


cfg = {
    """
    The default setting is the parameter of Guangdong region
    Sunlight line and equatorial crossing angle(rad):(ps:must less than 45)
    """
    # H: total horizontal radiation
    'H': 11.375,
    # hd: Horizontal scattering radiation
    'hd': 7.08, 
    # Ground reflectivity
    'rou' : 0.2,
    # Solar constant    
    #1 kwh = 3.6*10^6 J
    #1 MJ = 1*10^6 J
    #1 kwh = 3.6 MJ
    #1 kw/m^2 =  86.4 MJ/m^2/day
    #isc = 1.367 kW/m^2 = 118.11 MJ/m^2  
    'isc' : 118.11, # MJ/m^2/day
    #latitude  [-66.5,66.5]
    'latitude' : 23.16,
    # Building opening time
    'lowlim': 7, 
    # Building closing time, lowlim<uplim  
    'uplim' : 18,
    # Annual electricity consumption per capita in the community
    'E': 1400,
    # Per capita daily carbon emissions
    'C':6.2,
    # t hours one day
    'T_const': 24,
    # Solar electrode efficiency  [0,1]
    'yita':0.08,
    # Environmental influence factor of external power supply
    'r0': 1,
    # Days in a year [365,366]
    'N': 365,
    # Set electricity price per kilowatt hour
    'p' : 0.5,
    # Annual output value per square meter of phototropic photoanode material
    'q': 95.87,
    #Environmental effect weight. a+b=1
    'a': 0.9,
    # Economic benefit weight
    'b':0.1,
    # Environmental impact factor of per capita carbon footprint
    'rp': 1.13,
    # Annual average solar radiation per square meter on the light plane
    'I':5700,
    #unit e/td. The sunlight energy received per mole of carbon dioxide consumed
    'k': 36.55,
    #ratio is the area utilization rate of the photocatalytic material, the unit is percentage
    'ratio': 0.1963,
    # The people flow function is approximated as a Gaussian function
    'Gaussian_mu':12,
    'Gaussian_sigma':1,
    'Gaussian_peak':150,
    'Ei':[],
    'Si':[]
} 
    
def pdf2(x,mu,sigma,c):
    return c * math.exp(-(x-mu) ** 2 / (2* sigma**2)) 
        
def sum_fun_xk(xk, func,mu,sigma,c):
    return sum([func(each,mu,sigma,c) for each in xk])
        
def integral(a, b, n, func,mu,sigma,c):
    h = (b - a)/float(n)
    xk = [a + i*h for i in range(1, n)]
    return h/2 * (func(a,mu,sigma,c) + 2 * sum_fun_xk(xk, func,mu,sigma,c) + func(b,mu,sigma,c))
def cdfd(lowlim,uperlime,c,mu,sigma):
    return integral(lowlim,uperlime,10000,pdf2,mu,sigma,c)

def savefile_csv(file_path, *args):
    f = open(file_path,'w')
    writer = csv.writer(f)
    for i in args:
        writer.writerow(i)
    f.close()

def savefile_json(file_path, **kwargs):
    filejson = json.dumps(kwargs)
    f = open(file_path, 'w')
    f.write(filejson)
    f.close()

def getei_si():
    #analysis 1: With the change of a, b, the ratio of materials and flexible solar cells to the total surface of the optimal solution
    globals().update(cfg)         
    result = cdfd(lowlim,uplim,Gaussian_peak,Gaussian_mu,Gaussian_sigma)  
    Sa = k * C * N * (uplim - lowlim) / I / ratio #m^2/person     

    n=len(Si)
    Si_sum = [0] * n             
    Si_sum[0] = Si[0]
    for i in range(1, n):
        Si_sum[i] = Si_sum[i - 1] + Si[i]
            
    Ei_sum = [0]*n      #mj
    Ei_sum[0] = Ei[0]
    for i in range(1, n):
        Ei_sum[i] = Ei_sum[i - 1]  + Ei[i]
            
    ans_T = -inf
    ans_j = 0

    red_list=[]
    white_list=[]
    anssheet={}  
    a_list=[]
    b_list=[]
    
    onesheet={}
    onesheet['k']=k
    onesheet['q']=q
    
    Tecology_list=[0]*n
    Teconomy_list=[0]*n
    temin = inf
    temax = -temin
    tecmin = inf
    tecmax = -tecmin
    for j in range(0, n):
        Tenv = r0 * (E / T_const * result - yita * Ei_sum[j]) + (rp * N * (uplim - lowlim)) / T_const * (result / (uplim - lowlim) - (Si_sum[n-1] - Si_sum[j]) / Sa)  #kw.h
        Teconomy = (-1) * p * (E / T_const * result - yita * Ei_sum[j]) + q * (Si_sum[n-1] - Si_sum[j]) #yuan
        Tecology_list[j]=Tenv
        Teconomy_list[j]=Teconomy
        temin = min(temin,Tenv)
        temax = max(temax,Tenv)
        tecmin = min(tecmin,Teconomy)
        tecmax = max(tecmax,Teconomy)
    telen = temax - temin
    teclen = tecmax - tecmin
    
    for i in range(0, n):
        Tecology_list[i] = (Tecology_list[i]-temin)/telen
        Teconomy_list[i] = (Teconomy_list[i]-tecmin)/teclen

    a=0.00
    ei_si=[]
    ei_pro=[]
    while(a<=1.001):
        b=1-a
        ans_T = -1000000000
        ans_j = 0
        
        for j in range(0, n):
            T = (-1) * a * Tecology_list[j] + b * Teconomy_list[j]
            if T > ans_T:
                ans_T = T
                ans_j = j
        a_list.append(a)
        b_list.append(b)
        ei_si.append(Ei_sum[ans_j]/Si_sum[ans_j])#平均太阳辐射量
        ei_pro.append(Si_sum[ans_j]/Si_sum[n-1])
        a+=0.01
    file_path = '/Users/zsm/Documents/project/writebook/code/ei_si.csv'
    savefile_csv(file_path, a_list, b_list, ei_si, ei_pro)
    
    return [],[],[],[]    

def gete_area_ratio():    
    #analysis 2: Environmental and ecological impact factors of unit materials and solar panels as a and b change
    globals().update(cfg)         
    result = cdfd(lowlim,uplim,Gaussian_peak,Gaussian_mu,Gaussian_sigma)  
    Sa = k * C * N * (uplim - lowlim) / I / ratio #m^2/person     
    
    n=len(Si)
    Si_sum = [0] * n             
    Si_sum[0]=Si[0]
    for i in range(1, n):
        Si_sum[i] = Si_sum[i - 1] + Si[i]
            
    Ei_sum = [0]*n      #mj
    Ei_sum[0] = Ei[0]
    for i in range(1, n):
        Ei_sum[i] = Ei_sum[i - 1]  + Ei[i]
            
    ans_T = -inf
    ans_j = 0
    

    red_list=[]
    white_list=[]
    anssheet={}  
    a_list=[]
    b_list=[]
    
    onesheet={}
    onesheet['k']=k
    onesheet['q']=q
    
    Tecology_list=[0]*n
    Teconomy_list=[0]*n
    temin = 10000000000
    temax = -temin
    tecmin = 1000000000
    tecmax = -tecmin
    for j in range(0, n):
        Tenv = r0 * (E / T_const * result - yita * Ei_sum[j]) + (rp * N * (uplim - lowlim)) / T_const * (result / (uplim - lowlim) - (Si_sum[n-1] - Si_sum[j]) / Sa)  #kw.h
        Teconomy = (-1) * p * (E / T_const * result - yita * Ei_sum[j]) + q * (Si_sum[n-1] - Si_sum[j]) #yuan
        Tecology_list[j]=Tenv
        Teconomy_list[j]=Teconomy
        temin = min(temin,Tenv)
        temax = max(temax,Tenv)
        tecmin = min(tecmin,Teconomy)
        tecmax = max(tecmax,Teconomy)
    telen = temax - temin
    teclen = tecmax - tecmin
    
    
    for i in range(0, n):
        Tecology_list[i] = (Tecology_list[i]-temin)/telen
        Teconomy_list[i] = (Teconomy_list[i]-tecmin)/teclen

    a_list=[]
    red_list=[]
    white_list=[]
    anssheet={} 
    a_list=[]
    b_list=[]
    red_list=[]
    white_list=[]
    onesheet={}
    onesheet['k']=k
    onesheet['q']=q
        
    a=0.00
    ecology_red=[]
    ecology_white=[]
    economy_red=[]
    economy_white=[]

    while(a<=1.001):
        b=1-a
        ans_T = -1000000000
        ans_j = 0
        now_ecology_red=0
        now_ecology_white=0
        now_economy_red=0
        now_economy_white=0
        for j in range(0, n):
            T = (-1) * a * Tecology_list[j] + b * Teconomy_list[j]
            if T > ans_T:
                ans_T = T
                ans_j = j
                now_ecology_red=-1*Tecology_list[j]/Ei_sum[j]
                now_ecology_white=-1*Tecology_list[j]/ (Si_sum[n-1] - Si_sum[j])
                now_economy_red=Teconomy_list[j]/Ei_sum[j]
                now_economy_white=Teconomy_list[j]/ (Si_sum[n-1] - Si_sum[j])
        a_list.append(a)
        b_list.append(b)
        ecology_red.append(now_ecology_red)
        ecology_white.append(now_ecology_white)
        economy_red.append(now_economy_red)
        economy_white.append(now_economy_white)
        a+=0.01

    onesheet['a_list']=a_list
    onesheet['b_list']=b_list
    onesheet['ecology_red']=ecology_red
    onesheet['ecology_white']=ecology_white
    onesheet['economy_red']=economy_red
    onesheet['economy_white']=economy_white
        
    file_path = '/Users/zsm/Documents/project/writebook/code/four_json.json'
    savefile_json(file_path, onesheet)  

def get_kchange():
    #analysis 3: The effect of k changes on the results
    globals().update(cfg)         
    result = cdfd(lowlim,uplim,Gaussian_peak,Gaussian_mu,Gaussian_sigma)  
    
    cha = 1e-4

    n=len(Si)
    Si_sum = [0] * n
    Si_sum[0]=Si[0]
    for i in range(1, n):
        Si_sum[i] = Si_sum[i - 1] + Si[i]
            
    Ei_sum = [0] * n
    Ei_sum[0] = Ei[0]
    for i in range(1, n):
        Ei_sum[i] = Ei_sum[i - 1]  + Ei[i]
            
    ans_T = -1000000000
    ans_j = 0
    kreg=36.55
    qreg=95.87
    chengji = kreg * qreg
    # bili=[1/1000, 1.0/5.0 , 1.0/4.0 , 1.0/3.0 , 1.0/2.0 , 1 , 2 , 3 , 4 , 5]
    bili =  [0.02, 0.04, 0.06 , 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2.0, 4.0, 6.0, 8, 10, 12, 18, 20, 40, 60 ]
    a_list=[]
    red_list=[]
    white_list=[]
    anssheet={}
    ratio = 0.1963
    nowk = 1
    heng = []
    zong = []
    path = '/Users/zsm/Documents/project/writebook/code/notregular/'
    for i in range(0,100):
        # k=kreg*bili[i]
        # q=qreg/bili[i]
        k = nowk
        heng.append(k)
        nowk = nowk + 1
        q = chengji / k
        a_list=[]
        b_list=[]
        red_list=[]
        white_list=[]
        # filename="k="+str(k)
        # worksheet = workbook.add_sheet(filename)
        Sa = k * C * N * (uplim - lowlim) / I / ratio
        onesheet={}
        onesheet['k']=k
        onesheet['q']=q
        qwq=-1
        # worksheet.write(0, 0, k) 
        # worksheet.write(1,0,'a')
        # worksheet.write(2,0,'b')
        # worksheet.write(3,0,'s_red')
        # worksheet.write(4,0,'s_white')
        # qwq=-1
        a=0
        # temin = 10000000000
        # temax = -temin
        # tecmin = 1000000000
        # tecmax = -tecmin
        Tecology_list=[0]*n
        Teconomy_list=[0]*n
        for j in range(0, n):
            Tenv = r0 * (E / T_const * result - yita * Ei_sum[j]) + (rp * N * (uplim - lowlim)) / T_const * (result / (uplim - lowlim) - (Si_sum[n-1] - Si_sum[j]) / Sa)  #kw.h
            Teconomy = (-1) * p * (E / T_const * result - yita * Ei_sum[j]) + q * (Si_sum[n-1] - Si_sum[j]) #yuan
            Tecology_list[j]=Tenv
            Teconomy_list[j]=Teconomy
            # temin = min(temin,Tenv)
            # temax = max(temax,Tenv)
            # tecmin = min(tecmin,Teconomy)
            # tecmax = max(tecmax,Teconomy)
        # telen = temax - temin
        # teclen = tecmax - tecmin
        # Normalization affects monotonicity
        # for j in range(0, n):
        #     Tecology_list[j] = (Tecology_list[j]-temin)/telen
        #     Teconomy_list[j] = (Teconomy_list[j]-tecmin)/teclen

        happen = 0
    
        while(a<=1.001):
            b=1-a
            ans_T = -1000000000
            ans_j = 0  
            for j in range(0, n):
                T = (-1) * a * Tecology_list[j] + b * Teconomy_list[j]
                if T > ans_T:
                    ans_T = T
                    ans_j = j
            # if(qwq!=ans_j):
            a_list.append(a)
            b_list.append(b)
            red_list.append(Si_sum[ans_j]/Si_sum[n-1])
            if Si_sum[ans_j]/Si_sum[n-1]>cha and happen==0:
                zong.append(a)
                happen = 1
            white_list.append(1-Si_sum[ans_j]/Si_sum[n-1])
            # qwq=ans_j
            a+=0.01
        # count=1
        # for numa in a_list:
        # onesheet['a_list']=a_list
        # onesheet['b_list']=b_list
        # onesheet['red_list']=red_list
        # onesheet['white_list']=white_list
        # count+=1
            # worksheet.write(1,count,numa)
            # worksheet.write(2,count,1-numa)
            # worksheet.write(3,count,red_list[count-1])
            # worksheet.write(4,count,1-red_list[count-1])
        # anssheet[i] = onesheet
        
        # writer.writerow(a_list)
        # writer.writerow(b_list)
        # writer.writerow(red_list)
        # writer.writerow(white_list)
        i+=1
    # savefile_json('/Users/zsm/Documents/project/writebook/code/kchange_json.json',anssheet)  
    savefile_csv(path+'k='+str(k)+'.csv','w',heng,zong)     
    return 0,0,0,0

def get_ratiochange():
    #analysis 4: The effect of ratio changes on the results
    globals().update(cfg)         
    result = cdfd(lowlim,uplim,Gaussian_peak,Gaussian_mu,Gaussian_sigma)  
    Sa = k * C * N * (uplim - lowlim) / I / ratio #m^2/person     

    n=len(Si)
    Si_sum = [0] * n
    Si_sum[0]=Si[0]
    for i in range(1, n):
        Si_sum[i] = Si_sum[i - 1] + Si[i]
            
    Ei_sum = [0] * n
    Ei_sum[0] = Ei[0]
    for i in range(1, n):
        Ei_sum[i] = Ei_sum[i - 1]  + Ei[i]
            
    ans_T = -1000000000
    ans_j = 0
    k=36.55
    q=95.87
    ratio = 0.1
    # bili=[1/93.45, 1.0/5.0 , 1.0/4.0 , 1.0/3.0 , 1.0/2.0 , 1 , 2 , 3 , 4 , 5]
    a_list=[]
    red_list=[]
    white_list=[]
    anssheet={}
    ratio = 0.1
    i = 1
    path = '/Users/zsm/Documents/project/writebook/code/'
    while ratio<1.01:
        a_list=[]
        b_list=[]
        red_list=[]
        white_list=[]
        # filename="k="+str(k)
        # worksheet = workbook.add_sheet(filename)
        Sa = k * C * N * (uplim - lowlim) / I / ratio
        onesheet={}
        onesheet['ratio']=ratio
        qwq=-1
        # worksheet.write(0, 0, k) 
        # worksheet.write(1,0,'a')
        # worksheet.write(2,0,'b')
        # worksheet.write(3,0,'s_red')
        # worksheet.write(4,0,'s_white')
        # qwq=-1
        a=0
        temin = 10000000000
        temax = -temin
        tecmin = 1000000000
        tecmax = -tecmin
        Tecology_list=[0]*n
        Teconomy_list=[0]*n
        for j in range(0, n):
            Tenv = r0 * (E / T_const * result - yita * Ei_sum[j]) + (rp * N * (uplim - lowlim)) / T_const * (result / (uplim - lowlim) - (Si_sum[n-1] - Si_sum[j]) / Sa)  #kw.h
            Teconomy = (-1) * p * (E / T_const * result - yita * Ei_sum[j]) + q * (Si_sum[n-1] - Si_sum[j]) #yuan
            Tecology_list[j]=Tenv
            Teconomy_list[j]=Teconomy
            temin = min(temin,Tenv)
            temax = max(temax,Tenv)
            tecmin = min(tecmin,Teconomy)
            tecmax = max(tecmax,Teconomy)
        telen = temax - temin
        teclen = tecmax - tecmin
        
        
        for j in range(0, n):
            Tecology_list[j] = (Tecology_list[j]-temin)/telen
            Teconomy_list[j] = (Teconomy_list[j]-tecmin)/teclen

        while(a<=1.001):
            b=1-a
            ans_T = -1000000000
            ans_j = 0  
            for j in range(0, n):
                T = (-1) * a * Tecology_list[j] + b * Teconomy_list[j]
                if T > ans_T:
                    ans_T = T
                    ans_j = j
        
            a_list.append(a)
            b_list.append(b)
            red_list.append(Si_sum[ans_j]/Si_sum[n-1])
            white_list.append(1-Si_sum[ans_j]/Si_sum[n-1])
            
            a+=0.01
        # onesheet['a_list']=a_list
        # onesheet['b_list']=b_list
        # onesheet['red_list']=red_list
        # onesheet['white_list']=white_list
        # anssheet[i] =onesheet
        savefile_csv(path+'ratio='+str(i)+'.csv',a_list,b_list,red_list,white_list)
        i+=1
        ratio += 0.1 
    # savefile_json('/Users/zsm/Documents/project/writebook/code/ratiochange_json.json',anssheet)       
    return 0,0,[],[]

def check(k,q,a,b,Si_sum,Ei_sum,n,ratio): 
    uplim=cfg['uplim']
    lowlim=cfg['lowlim']  
    E=cfg['E']            
    C=cfg['C']
    I=cfg['I']
    yita=cfg['yita']
    r0=cfg['r0']
    N=cfg['N']  
    p=cfg['p']
    T_const=24
    rp=cfg['rp']
    result =cdfd(lowlim,uplim,cfg['Gaussian_peak'],cfg['Gaussian_mu'],cfg['Gaussian_sigma'])   
    ans_T = -1000000000
    ans_j = 0  
    Sa = k * C * N * (uplim - lowlim) / I / ratio
    for j in range(0, n):
        Tecology = r0 * (E / T_const * result - yita * Ei_sum[j]) + (rp * N * (uplim - lowlim)) / T_const * (result / (uplim - lowlim) - (Si_sum[n-1] - Si_sum[j]) / Sa)
        Teconomy = (-1) * p * (E / T_const * result - yita * Ei_sum[j]) + q * (Si_sum[n-1] - Si_sum[j])
        T = (-1) * a * Tecology + b * Teconomy
        if T > ans_T:
            ans_T = T
            ans_j = j
    return ans_j    
def getZero():
    #analysis 5: Find the critical k value
    globals().update(cfg)         
    Si_sum = [0] * len(Si)
            
    Si_sum[0]=Si[0]
    n=len(Si)
    for i in range(1, n):
        Si_sum[i] = Si_sum[i - 1] + Si[i]
            
    Ei_sum = [0]*n
            
    Ei_sum[0] = Ei[0]
    for i in range(1, n):
        Ei_sum[i] = Ei_sum[i - 1]  + Ei[i]
            
    ans_T = -1000000000
    ans_j = 0
    a_list=[]
    red_list=[]
    white_list=[]
    anssheet={}
    l=0
    r=1.0/5
    cha=0.00001
    ans_mid=0
    
    while(abs(r-l)>cha):
        mid=(l+r)/2
        if(check(kreg*mid,qreg/mid,1,0,Si_sum,Ei_sum,n,ratio)>0.5):  
            r=mid
            ans_mid=mid
        else:
            l=mid
    
    return 0,0,0,l*kreg

                
    
# def myacos(xx):
#     if xx>=-1 and xx<=1 :
#          return math.acos(xx)
#     elif xx<-1: 
#         return math.pi
#     return 0
                        
# def myasin(xx):
#     if xx>=-1 and xx<=1 : return math.asin(xx)
#     elif xx<-1: return -math.pi/2
#     return math.pi/2
    
def c_tuple(td_in):
    return (float(td_in[0]), float(td_in[1]), float(td_in[2]))
def c_list(td_in):
    return [float(td_in[0]), float(td_in[1]), float(td_in[2])]

def getAzmuth(ver):
    #print('normal vector :',ver)
    k=(0,0,1)
    val1=abs(ver[2]/abs(math.sqrt(ver[0]**2+ver[1]**2+ver[2]**2)))  #0～
    s=math.acos(val1)
    if ver[0]==0 and ver[1]==0:val2=0
    else :
        val2=-ver[1]/abs(math.sqrt(ver[0]**2+ver[1]**2))
    azmuth=math.acos(val2)
    if ver[0]>0 :azmuth*=-1
    return s,azmuth

def myasin(xx):
    if xx>=-1 and xx<=1 :
        return math.asin(xx)
    # elif xx<-1: 
        # print("debug sin out of range ")
        # return -math.pi/2
        # return inf
    else : 
        # print("debug sin out of range ")
        # return math.pi/2
        return inf
def myacos(xx,name):
    if xx>=-1 and xx<=1 :
        return math.acos(xx)
    # elif xx<-1: 
        # print(name,"debug cos out of range ",xx)
        # return math.pi
        # return inf
    else:
        # print(name,"debug cos out of range ",xx) 
        # return 0
        return inf

def getRadiant(s,azmuth, latitude):
    DayOfMonth = [17,47,75,105,135,162,198,228,258,288,318,344]
    days=[31,28,31,30, 31, 30, 31, 31, 30, 31, 30, 31]
    latitude=latitude/180*pi
    sume=0
    for i in range(12):
        hudu = 360*(284+DayOfMonth[i])/365 * pi/180
        solarray=23.45*math.sin(hudu) #woc it doesn't change unit
        solarray=solarray/180*pi
        tmp2=-math.tan(latitude)*math.tan(solarray)
        hsunsetangle=myacos(tmp2,"hsunsetangle")
        a=math.sin(solarray)*(math.sin(latitude)*math.cos(s)-math.cos(latitude)*math.sin(s)*math.cos(azmuth))
        b=math.cos(solarray)*(math.cos(latitude)*math.cos(s)+math.sin(latitude)*math.sin(s)*math.cos(azmuth))
        c=math.cos(solarray)*math.sin(s)*math.sin(azmuth)
        D=math.sqrt(b**2+c**2)
        ssunsetangle=min(hsunsetangle,myacos(-a/D,"ssunsetangle")+myasin(c/D))
        wr=-min(hsunsetangle,abs(-myacos(-a/D,"wr")+myasin(c/D)))
        print("wr",wr)
        tmp1=pi/180*(ssunsetangle-wr)*math.sin(solarray)*(math.sin(latitude)*math.cos(s)-math.cos(latitude)*math.sin(s)*math.cos(azmuth))
        tmp2=math.cos(solarray)*(math.sin(ssunsetangle)-math.sin(wr))*(math.cos(latitude)*math.cos(s)+math.sin(latitude)*math.sin(s)*math.cos(azmuth))
        tmp3=(math.cos(ssunsetangle)-math.cos(wr))*math.cos(solarray)*math.sin(s)*math.sin(azmuth)
        print("lati",latitude,solarray,hsunsetangle)
        rb = (tmp1+tmp2+tmp3)/2/(math.cos(latitude)*math.cos(solarray)*math.sin(hsunsetangle)+(math.pi/180)*hsunsetangle*math.sin(latitude)*math.sin(solarray))
        h0=24/pi*cfg['isc']*(1+0.033*math.cos(360*DayOfMonth[i]/365))*(math.cos(latitude)*math.cos(solarray)*math.sin(hsunsetangle)+(2*pi*hsunsetangle/360)*math.sin(latitude)*math.sin(solarray))
        # MJ/m^2/day 
        hb=cfg['H']-cfg['hd']
        # MJ/m^2/day   1kwh = 3.6 MJ 
        hbt=hb*rb
        # MJ/m^2/day
        hdt=cfg['hd']*((hb/h0)*rb+0.5*(1-hb/h0)*(1+math.cos(s)))
        # MJ/m^2/day   
        hrt=0.5*cfg['rou']*cfg['H']*(1-math.cos(s))
        # MJ/m^2/day
        ht=hbt+hdt+hrt
        sume+=ht*days[i]
        
    return sume   #MJ/m^2 /year

def getlen(ver):
    return math.sqrt(ver[0]**2+ver[1]**2+ver[2]**2)

def hailun(aa,bb,cc):
    p=(aa+bb+cc)/2
    if aa+bb<=cc or bb+cc<=aa or aa+cc<=bb: return 0
    s=math.sqrt(p*(p-aa)*(p-bb)*(p-cc))
    return s

def getArea(x):    
    v12=[x[1][0]-x[0][0],x[1][1]-x[0][1],x[1][2]-x[0][2]]
    v13=[x[2][0]-x[0][0],x[2][1]-x[0][1],x[2][2]-x[0][2]]
    v14=[x[3][0]-x[0][0],x[3][1]-x[0][1],x[3][2]-x[0][2]]
#    v3=np.cross(v1,v2)
    v23=[x[1][0]-x[2][0],x[1][1]-x[2][1],x[1][2]-x[2][2]]
    v34=[x[2][0]-x[3][0],x[2][1]-x[3][1],x[2][2]-x[3][2]]
#    v4=np.cross(v1,v2)
    s1=hailun(getlen(v12),getlen(v13),getlen(v23))
    s2=hailun(getlen(v14),getlen(v13),getlen(v34))
#    print(s1+s2)
    return s1+s2

class Mian(object):
    def __init__( self, poin,nor,num,faceid):
        self.poin = poin
        self.nor = nor
        self.num = num
        self.faceid = faceid
        self.s,self.azmuth = getAzmuth(nor) 
        #print(self.azmuth)
        self.si=getArea(poin)
        # self.ei=getRadiant(self.s,self.azmuth, cfg['latitude'])*self.si/3.6
        self.eper=getRadiant(self.s,self.azmuth, cfg['latitude'])/3.6 #3.6 is change MJ to kwh
        self.ei=self.eper*self.si #3.6 is change MJ to kwh
        
    def __lt__(self,other):
        return self.eper > other.eper

def getAns():
    #main method
    globals().update(cfg)
    result = cdfd(lowlim,uplim,Gaussian_peak,Gaussian_mu,Gaussian_sigma)  
    
    result = 10
    Sa = k * C * N * (uplim - lowlim) / I / ratio #m^2/person     
    
    n=len(Si)
    Si_sum = [0] * n             
    Si_sum[0]= Si[0]
    
    
    for i in range(1, n):
        Si_sum[i] = Si_sum[i - 1] + Si[i]
            
    Ei_sum = [0]*n      #mj
    Ei_sum[0] = Ei[0]
    for i in range(1, n):
        Ei_sum[i] = Ei_sum[i - 1] + Ei[i]

    red_list=[]
    white_list=[]
    anssheet={}  
    a_list=[]
    b_list=[]

    Tecology_list=[0]*n
    Teconomy_list=[0]*n
    temin = inf
    temax = -inf
    tecmin = inf
    tecmax = -inf
    for j in range(0, n):
        Tenv = r0 * (E / T_const * result - yita * Ei_sum[j]) + (rp * N * (uplim - lowlim)) / T_const * (result / (uplim - lowlim) - (Si_sum[n-1] - Si_sum[j]) / Sa)  #kw.h
        Teconomy = (-1) * p * (E / T_const * result - yita * Ei_sum[j]) + q * (Si_sum[n-1] - Si_sum[j]) #yuan
        Tecology_list[j]=Tenv
        Teconomy_list[j]=Teconomy
        temin = min(temin,Tenv)
        temax = max(temax,Tenv)
        tecmin = min(tecmin,Teconomy)
        tecmax = max(tecmax,Teconomy)
    telen = temax - temin
    teclen = tecmax - tecmin
    
    for i in range(0, n):
        Tecology_list[i] = (Tecology_list[i]-temin)/telen
        Teconomy_list[i] = (Teconomy_list[i]-tecmin)/teclen
    
    ans_T = -inf
    ans_j = 0
    for j in range(0, n):
        T = (-1) * a * Tecology_list[j] + b * Teconomy_list[j]           
        if T > ans_T:
            ans_T = T
            ans_j = j

    return ans_j,Si_sum[ans_j],0,0

class MyComponent(component):
    
    def RunScript(self, _geometry, H, hd, rou, isc, latitude, uplim, lowlim, E, C, T, yita, r0, N, p, q, a, b, rp, I, k, ratio, Gaussian_mu, Gaussian_sigma, Gaussian_peak):
        analysisMesh=None
        radiationResult=None
        self.Name = "A software for the designing of the installation of photoelectrocatalytic building materials and flexible solar cells."
        self.NickName = 'MaterialInstallationSoftware'
        self.Message = 'VER 1.1.1\nNov_10_2021'
        self.Category = "GreenChemistry"
        self.SubCategory = "1 | MaterialInstallation"
        
        def RandomColor():
            red = random.randint(0,255)
            green = random.randint(0,255)
            blue = random.randint(0,255)
            return rs.coercecolor((red,green,blue))

        def init():
        #   global cfg
            if H == None : cfg['H']=11.375
            else: cfg['H']=H
                
            if hd == None: cfg['hd']=7.08
            else: cfg['hd']=hd
                
            if rou == None : cfg['rou']=0.2
            else: cfg['rou']=rou
                
            if isc == None: cfg['isc']=118.11
            else: cfg['isc']=isc
                
            if latitude == None : cfg['latitude']=23.16
            else: cfg['latitude']=latitude
                
            if uplim == None: cfg['uplim']=18
            else: cfg['uplim']=uplim
                
            if lowlim == None : cfg['lowlim']=7
            else: cfg['lowlim']=lowlim
                
            if E == None: cfg['E']=1400
            else: cfg['E']=E
                
            if C == None : cfg['C']=6.2
            else: cfg['C']=C
                
            if T == None: cfg['T']=24
            else: cfg['T']=T
                
            if yita == None : cfg['yita']=0.08
            else: cfg['yita']=yita
                
            if r0 == None: cfg['r0']=1
            else: cfg['r0']=r0
                
            if N == None : cfg['N']=365
            else: cfg['N']=N
                
            if p == None: cfg['p']=0.5
            else: cfg['p']=p
                
            if q == None: cfg['q']=95.87
            else: cfg['q']=q
                    
            if a == None: cfg['a']=0.1
            else: cfg['a']=a
                
            if b == None: cfg['b']=0.9
            else: cfg['b']=b
                
            if rp == None: cfg['rp']=1.13
            else: cfg['rp']=rp
                
            if I == None: cfg['I']=5700
            else: cfg['I']=I   
                
            if k == None: cfg['k']=36.55
            else: cfg['k']=k   

            if ratio == None: cfg['ratio']=0.1963
            else: cfg['ratio']=ratio   
            
            if Gaussian_mu == None: cfg['Gaussian_mu']=12
            else: cfg['Gaussian_mu']=Gaussian_mu
                
            if Gaussian_sigma == None: cfg['Gaussian_sigma']=1
            else: cfg['Gaussian_sigma']=Gaussian_sigma
                
            if Gaussian_peak == None: cfg['Gaussian_peak']=150
            else: cfg['Gaussian_peak']=Gaussian_peak   
            
        def main():
            init()
            if _geometry==None :num=0
            else : num = len(_geometry)
            if num != 0 and _geometry[0] != None :

                face_id = 0
                totalmesh=[]
                faces_all=[]
                point_list = []
                for i in range(num):
                    normals = rs.MeshFaceNormals(_geometry[i])#face unit normal for each face of a mesh object
                    faces_ver = rs.MeshFaces(_geometry[i],True)#return mesh point tra
                            
                    cou = 0
                    faces_all.append(faces_ver)
                    
                    j = 0  #in this piece of large surface, its triangular surface number
                    rs.EnableRedraw(False)
                            
                    while( j<len(faces_ver) ):
                        piece_face = (c_tuple(faces_ver[j]), c_tuple(faces_ver[j+1]), c_tuple(faces_ver[j+2]),c_tuple(faces_ver[j+3]))
                        point_list.append(c_list(faces_ver[j]))
                        point_list.append(c_list(faces_ver[j+1]))
                        point_list.append(c_list(faces_ver[j+2]))
                        point_list.append(c_list(faces_ver[j+3]))
                        x= Mian(piece_face,c_tuple(normals[cou]),j,face_id )
                        j += 4
                        totalmesh.append(x)
                        cou+=1
                    
                    face_id+=1   #face_id is a big surface
                totalmesh.sort()
                Ei=[]
                Si=[]

                for x in totalmesh:
                    Ei.append(x.ei)
                    Si.append(x.si)
                cfg['Ei'] = Ei
                cfg['Si'] = Si
                 
                #global analysisMesh,radiationResult,list1,list2
                #analysisMesh=faces_all
                #radiationResult=Ei

                ansmesh,ansarea,print1,print2 = getAns()
    
                minis = inf
                maxis = -inf
                sizes = len(Si)
                for x in range (sizes):
                    minis = min(minis,Si[x])
                    maxis = max(maxis,Si[x])

                print(ansmesh)
                polygon_list = []
                rs.EnableRedraw(True)
                list1=[]
                list2=[]
                list3=[]
                list4=[]
                """
                # rs.MaterialColor(material_index, color=(255,0,0))
                for i in range(0,ansmesh+1):
                    id = totalmesh[i].num
                    fd =  totalmesh[i].faceid 
                    # print_face = faces_all[fd][id], faces_all[fd][id+1], faces_all[fd][id+2], faces_all[fd][id+3]
                    list1.append( faces_all[fd][id])
                    list2.append( faces_all[fd][id+1])
                    list3.append( faces_all[fd][id+2])
                    list4.append( faces_all[fd][id+3])
                    """
                return list1,list2,list3,list4,ansmesh,ansarea
               
        (point1,point2,point3,point4, meshNum,mesharea) = main()
       
        return (point1,point2,point3,point4, meshNum, mesharea)
