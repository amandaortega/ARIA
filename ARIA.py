#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 20:22:58 2017

@author: eulerr
"""
import random
from tqdm import trange


import math
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
import scipy.io
from scipy.spatial import distance
import copy


matD31 = scipy.io.loadmat('d31.mat')

dataD31 = matD31['dados']
labelD31 = matD31['rotulo'] 



matShan = scipy.io.loadmat('shan_compound.mat')

dataShan = matShan['dados']
labelShan = matShan['rotulo'] 




def plotstate():
    global Ag
    global Ab
    fig, ax = plt.subplots()
    ax.cla()
    ax.scatter(Ag[:,0],Ag[:,1],c = 'blue', edgecolors='black',s = 30,alpha = 0.9,linewidths=1, marker = "+", )
    
    for idx,item in enumerate(Ab):
        ci = plt.Circle((item[0], item[1]), item[2], color='r',linewidth=1.5, fill=False)
        ax.scatter(item[0],item[1], c ='r', edgecolors='black',s = 30,alpha = 0.9,linewidths=1, marker = "8")
        ax.add_artist( ci )
        
       
    plt.title('Espaço de busca')
    plt.xlabel('x')
    plt.ylabel('y')
    #plt.yticks(y_major_ticks)
    plt.grid(which= 'both')
    plt.show()




# Não precisa do raio para calcular o fitness

 

def maturacao(N_MAT,MUTB):
    global Ag
    global Ab
  
    dist = 0
    a = 0.0
    b = 0.0
    fitness = np.zeros((len(Ab),1))
    
      
    for j in range(N_MAT):
        # Escolhe randomicamente um Ag e apresenta ele para todos Ab
        index = np.random.choice(Ag.shape[0],1, replace=False)[0]        
        a = Ag[index,0]
        b = Ag[index,1]
        
        
        for idx, item in enumerate(Ab):
        
            cx = item[0]
            cy = item[1] 
            
            dist = distance.euclidean([a,b],[cx,cy])
            
            
            fitness[idx,0] = 1/(math.pow(dist,2) + 0.001 )
        # Calcula o fitness de cada anticorpo para o determinado Ag
        
        # Seleciona maior fitness 
        Lambda = np.argmax(fitness)
        
        # Muta indivíduo de maior fitness na direção do Ag que ele representa       
        mut(MUTB,Lambda,index)


 
       
def mut(MUTB,Lambda,i):
    global Ag
    global Ab
    
    Ab[Lambda][0] = Ab[Lambda][0] + MUTB*random.uniform(0,1)*(Ag[i,0] - Ab[Lambda][0])
    Ab[Lambda][1] = Ab[Lambda][1] + MUTB*random.uniform(0,1)*(Ag[i,1] - Ab[Lambda][1])


  

  
def supressao():
    global Ag
    global Ab

              
    a = 0.0
    b = 0.0
    l = np.zeros((len(Ab),1))
#    fitness = np.zeros((len(Ab),1))
    supr = []

    for i in range(Ag.shape[0]):
        a = Ag[i,0]
        b = Ag[i,1]
        
      
        
        for idx, item in enumerate(Ab):
        
            cx = item[0]
            cy = item[1] 
            r =  item[3]
            dist = distance.euclidean([a,b],[cx,cy])
             
            
            if dist < r :
                l[idx] +=1 

    supr = [ idx for idx,item in enumerate(l) if item == 0 ]   # Pega Ab que não reconheceram nehum Ag
            
    if len(supr) ==len(Ab): # Evita supressão de todos os Ab

        supr.remove(random.choice(supr))

    Ab =  [v for i, v in enumerate(Ab) if i not in supr] # Suprimi Ab perdedores





# O anticorpo que vence é aquele que tem menor distancia de um AG fora do raio.

  
def clonagem(MUTB,E):
    global Ag
    global Ab
    global AbGamma
    
    
    l = [[] for i in range(len(AbGamma))]
    fitness = np.zeros((len(AbGamma),1))
    clon = []
    Ab1= []
    Ags= []

    for i in range(Ag.shape[0]): # Apresenta um Ag para todos Ab
        a = Ag[i,0]
        b = Ag[i,1]
        
        
        for idx, item in enumerate(AbGamma):
        
            cx = item[0]
            cy = item[1] 


            dist = distance.euclidean([a,b],[cx,cy])  
 
            fitness[idx,0] = 1/(math.pow(dist,2) + 0.001 )
        # Seleciona menor distância (  Ab que vence para Ag )
        Lambda = np.argmax(fitness)
        
        aux = math.sqrt((1/fitness[Lambda,0]) - 0.001 ) # Recupera a distância 
        
        if aux > AbGamma[Lambda][2] : # se a distância entre o Ab e Ag é maior que o R de Ab
            l[Lambda].append(np.array([a,b]))  # Salva Ag fora do raio de Ab             

        
    clon = [ idx for idx,i in enumerate(l) if i ]   # Pega o indice dos Ab que reconheceram  Ag > r
    
    Ags =  [ i for idx,i in enumerate(l) if i ]   # Pega a lista de Ag reconhecido fora do raio de Ab 

    Ab1 = [v for i, v in enumerate(AbGamma) if i in clon ]
 
    
    for idx,item in enumerate(Ab1):
#        aux = Ags[idx][0]
        aux = random.choice(Ags[idx])  # Escolhe aleatóriamente Ag reconhecido fora do raio de Ab 
        Ab1[idx][0] = item[0] + MUTB*random.uniform(0,1)*(aux[0] - item[0])
        Ab1[idx][1] = item[1] + MUTB*random.uniform(0,1)*(aux[1] - item[1])

    Ab = [*Ab,*Ab1]  # Concatena os Ab's clonados
    


 
  
    



def densidade(E):
    global Ag
    global Ab

    a = 0.0
    b = 0.0
    
    for idx,item in enumerate(Ab): # Zera densidade local para atualização
        Ab[idx][3] = 0 
       
    for i in range(Ag.shape[0]):
        a = Ag[i,0]
        b = Ag[i,1]
        
        
        for idx, item in enumerate(Ab):
        
            cx = item[0]
            cy = item[1] 

                        
            dist = distance.euclidean([a,b],[cx,cy])
            
            if dist < E:
               Ab[idx][3]+=1   # Atualiza valor de densidade com o # de AG
  




def atualR(r):
    global Ag
    global Ab
    global dim
    
    den = [item[3] for item in Ab]
    for idx, item in enumerate(Ab):

        if  Ab[idx][3] == 0: 
            Ab[idx][2] = r
            
        else:    
            Ab[idx][2] = r * math.pow((max(den)/Ab[idx][3]), 1/dim)
    





def suprAuto():
    supr = []
    global Ab    
  
    for idx, item in enumerate(Ab):
        a = item[0]
        b = item[1]
        r1 = item[2]
 
        for j, i in enumerate(Ab):

            if j == idx:
                pass
            
            else:
                cx = i[0]
                cy = i[1] 
                r2 = i[2]
                
                dist = distance.euclidean([a,b],[cx,cy])
                
                if dist < r1 or dist < r2:                   
                    if r1 > r2: 
                        supr.append(idx)
                    else:    
                        supr.append(j)    

    if len(supr) == len(Ab): # Evita supressão de todos os Ab
        supr.remove(random.choice(supr))

    Ab =  [v for i, v in enumerate(Ab) if i not in supr] # Suprimi Ab que se reconhecem
    





def main(pop_size,MAX_IT,N_MAT,r):
    global Ag
    global Ab
    global AbGamma
    
    random.seed(64)
    
    
    pop_size = pop_size # Tamanho população inicial
    
    MAX_IT = MAX_IT # Número máximo de interações
    N_MAT = N_MAT # Número de gerações
    
    Ab = []
    AbGamma= [] # Anticorpos pais
    
    MUTB = 1     # Taxa de mutação
    decay = 0.95 # Taxa de decaimento
    r = r     # Menor raio possível
    E = 5     # Raio que define a vizinhança para o cálculo da densidade local
    
    
    
    h={ 'supr': [],    # histórico da evolução
        'clon':  [], 
        'suprsi': [], 
        'E':     [], 
            }


    for i in range(pop_size) :
        Ab.append([random.uniform(0, 30),random.uniform(0, 30),random.uniform(r, 2*r),random.uniform(0, 10)])
    # cx,cy,r,densidade local
    
    #plotstate()
    
    for i in trange(MAX_IT):
        
        AbGamma = copy.deepcopy(Ab) # Faz cópia profunda para usar no processo de expansão clonal
        
        maturacao(N_MAT,MUTB)     
        #     plotstate()
        
        # Supressão
        aux= len(Ab)
        supressao()  #Ab's que não reconhecem nenhum Ag
        #print("Foram suprimidos {} anticorpos na {} iteração" .format((aux-len(Ab)),MAX_IT),end='\r')
        
        sup= aux-len(Ab)
           # plotstate()
        
        
        #Clonagem
        aux= len(Ab)
        clonagem(MUTB,E)   
        clon = len(Ab)-aux
        
           # plotstate()
        
        # Calcula densidade
        densidade(E)
        
        # Atualiza Raio de cada Ab
        atualR(r)
        
        
        #plotstate()
        
        # Supressão dos mais próximos entre si
         
        aux = len(Ab)
        suprAuto()
        #print("\nForam elimnados {} anticorpos que se auto-reconheceram" .format(aux-len(Ab))) 
        #plotstate()
        suprsi= aux-len(Ab)
        
        
        # Atualiza valor de E
        aux = [item[2] for item in Ab]
        E = sum(aux)/float(len(aux))   # Média dos Raios de Ab
        
        
        # Atualiza taxa de mutação
        MUTB = MUTB * decay


        h['supr'].append(sup)
        h['clon'].append(clon)
        h['suprsi'].append(suprsi)
        h['E'].append(E)


    return h
 
    

  
Ab = []
AbGamma= []
dim= dataShan.ndim
Ag= dataShan



h = main(pop_size=1,MAX_IT=15,N_MAT=240,r=0.92)
plotstate()



fig = plt.figure(figsize=(15,7))

plt.plot([None] + h['supr'], label="Ab's suprimidos")
plt.plot([None] + h['clon'], label="Ab's clonados")
plt.plot([None] + h['suprsi'], label="Ab's suprimidos (Auto Reconhecimento)")
plt.grid()
plt.xlabel("Iterações")
plt.ylabel("Número de Ag")
plt.title("Evolução Ab's")
plt.legend()

fig = plt.figure(figsize=(15,7))

plt.plot([None] + h['E'], label="Vizinhança E")
plt.grid()
plt.xlabel("Iterações")
plt.ylabel("")
plt.title('E')
plt.legend()



