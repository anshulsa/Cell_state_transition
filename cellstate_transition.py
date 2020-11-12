import numpy as np
from random import random
import math
import matplotlib.pyplot as plt

var = None

def q_birth(a, n):
    
    return a

def q_death(d1, d2, c, n, t):
    
    if n > 10:
        if t < c:
            return d1*n
        else:
            return d2*n
    else:
        return 0
        
def g_birth(a, n, t):
    
    return a*n         
        
def g_drug_sen_death(d1, d2, c, n, t):
    
    if n > 10:
        if t < c:
            return d1*n
        else:
            return d2*n
    else:
        return 0

def g_drug_res_death(d1, d2, c, n, t):
    
    if t < c:
        return d1*n
    else:
        return d2*n*(math.exp(-(t-c)/20))

def drug_sen_tran(a, b, n, t, um):
    if t < b:
        return a*um*n
    else:
        return a*um*n
      
def drug_res_tran(a, b, n, t, um):
    
    if n > 10:
        if t < b:
            return 0
        else:
            return a*um*n
    else:
        return 0

def single_step_birth(c):
    
    return c + 1
    
def single_step_death(c):

    return c - 1 

def single_step_trans(c, d):

    return (c - 1, d + 1)   
    
def multi_step_trans(m, b, u, c, d):
   
    a = np.zeros(6)
    
    a[0] =  b[0]    
    a[1] =  a[0] + b[1] 
    a[2] =  a[1] + b[2] 
    a[3] =  a[2] + b[3] 
    a[4] =  a[3] + b[4] 
    a[5] =  a[4] + b[5] 
    tot1 =  a[5] 
 
#    if  0 < u < a[0]/tot1:
#        m[0] = 1
#    elif  a[0]/tot1 < u < a[1]/tot1:
#        m[0] = 0
#    elif  a[1]/tot1 < u < a[2]/tot1:
#        m[1] = 1
#    elif  a[2]/tot1 < u < a[3]/tot1:
#        m[1] = 0
#    elif  a[3]/tot1 < u < a[4]/tot1:
#        m[2] = 1
#    elif  a[4]/tot1 < u < 1:
#        m[2] = 0
#    else:
#        print('error')
 
    if  0 < u < a[0]/tot1:
        if m[0] == 0: m[0] = 1
    elif  a[0]/tot1 < u < a[1]/tot1:
        if m[0] == 1: m[0] = 0
    elif  a[1]/tot1 < u < a[2]/tot1:
        if m[1] == 0: m[1] = 1
    elif  a[2]/tot1 < u < a[3]/tot1:
        if m[1] == 1: m[1] = 0
    elif  a[3]/tot1 < u < a[4]/tot1:
        if m[2] == 0: m[2] = 1
    elif  a[4]/tot1 < u < 1:
        if m[3] == 1: m[2] = 0
    else:
        print('error')
        
    if m[0] == 1 and m[1] == 1 and m[2] == 1: 
        m[0] = 0
        m[1] = 0
        m[2] = 0
        m[3] = m[3] + 1 
        return (c-1, d+1)  
    else:
        return (c-1, d+1)   


transition_temp = 0
m_mat = []
per_m_mat = []
average_gpr_mat = []
average_gsr_mat = []
mrange = (0, 0.1)#, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4)
#mrange = (1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3)
for m in mrange:
    print(m)
    inst_fun = (1, m, 1, m, 1, m)
    gsr_ave = []
    gpr_ave = []    
    for k in range(1, 2):
        
        #Initial number of cells in the states 
        q    = 80                              #number of Q cells
        gsr  = 0                               #number of GSR cells
        gstp = 3                               #number of GSTP cells
        gstr = 3                               #number of GSTR cells
        gpr  = 0                               #number of GPR cells
        gptp = 7                               #number of GPTP cells
        gptr = 7                               #number of GPTR cells
        totcell = q + gsr +gstp +gstr + gpr + gptp + gptr  
    
        t = 0
        t_d = 2
    
        # Birth and death rate Coefficients 
        
        # The units of all rate coefficients are in (Days)^(-1)
        c_bq = 0;       c_dq1 = 1;     c_dq2 = 2.2;
    
        c_bgsr  = 2;    c_dgsr1  = 0;  c_dgsr2  = 2.2;
        c_bgstp = 2;    c_dgstp1 = 1;  c_dgstp2 = 2.2;
        c_bgstr = 2;    c_dgstr1 = 1;  c_dgstr2 = 2.2;
        c_bgpr  = 2;    c_dgpr1  = 0;  c_dgpr2  = 2.2;
        c_bgptp = 2;    c_dgptp1 = 1;  c_dgptp2 = 2.2;
        c_bgptr = 2;    c_dgptr1 = 1;  c_dgptr2 = 2.2;
    
        c_ts = 2
        c_tr = 2        #list to store relevant parameters
        q_vector     = []
        gsr_vector   = []
        gstp_vector  = []
        gstr_vector  = []
        gpr_vector   = []
        gptp_vector  = []
        gptr_vector  = []
        totcell_vector = []
        t_vector     = []
        bd_ms00 = []
        bd_ms10 = []
        bd_ms20 = []   
        bd_ms30 = [] 
        bd_ms_temp = np.zeros((50000, 4))
        bd = np.zeros((1, 56))
        bd_ms = np.zeros((4,4))
        i = 0
        while  totcell > 0 and totcell < 10000 and t < 6:   #Stopping the simulation if Q cells become zero or time is over
            
            u1 = random()
            u2 = random()
            um = 1#random()
            u4 = random()  
            
                
            bq    = q_birth(c_bq, q)                                               #Birth rate of Q cells
            dq    = q_death(c_dq1, c_dq2, t_d, q, t)                                       #Death rate of Q cells
            bgsr  = g_birth(c_bgsr, gsr, t)                                        #Birth rate of GSR cells
            dgsr  = g_drug_res_death(c_dgsr1 ,c_dgsr2, t_d, gsr, t)                #Death rate of GSR cells 
            bgstp = g_birth(c_bgstp, gstp, t)                                      #Birth rate of GSTP cells
            dgstp = g_drug_sen_death(c_dgstp1, c_dgstp2, t_d, gstp, t)             #Death rate of GSTP cells 
            bgstr = g_birth(c_bgstr, gstr, t)                                      #Birth rate of GSTR cells
            dgstr = g_drug_sen_death(c_dgstr1, c_dgstr2, t_d, gstr, t)             #Death rate of GSTR cells 
            bgpr  = g_birth(c_bgpr, gpr, t)                                        #Birth rate of GPR cells
            dgpr  = g_drug_res_death(c_dgpr1, c_dgpr2, t_d, gpr, t)                #Death rate of GPR cells
            bgptp = g_birth(c_bgptp, gptp, t)                                      #Birth rate of GSTP cells
            dgptp = g_drug_sen_death(c_dgptp1, c_dgptp2, t_d, gptp, t)             #Death rate of GSTP cells 
            bgptr = g_birth(c_bgptr, gptr, t)                                      #Birth rate of GSTR cells
            dgptr = g_drug_sen_death(c_dgptr1, c_dgptr2, t_d, gptr, t)             #Death rate of GSTR cells 
            
            
                 
            bd[0,0] = bq
            bd[0,1] = bd[0,0] + dq
            bd[0,2] = bd[0,1] + 0
            bd[0,3] = bd[0,2] + drug_sen_tran(c_ts, t_d, q, t, um)
            bd[0,4] = bd[0,3] + drug_sen_tran(c_ts, t_d, q, t, um)
            bd[0,5] = bd[0,4] + 0
            bd[0,6] = bd[0,5] + drug_sen_tran(c_ts, t_d, q, t, um)
            bd[0,7] = bd[0,6] + drug_sen_tran(c_ts, t_d, q, t, um)
            bd[0,8] = bd[0,7] 
            bd[0,9] = bd[0,8] + bgsr
            bd[0,10] = bd[0,9] + dgsr
            bd[0,11] = bd[0,10] 
            bd[0,12] = bd[0,11] 
            bd[0,13] = bd[0,12] 
            bd[0,14] = bd[0,13] 
            bd[0,15] = bd[0,14] 
            bd[0,16] = bd[0,15] + drug_sen_tran(c_ts, t_d, gstp, t, um)
            bd[0,17] = bd[0,16] + drug_res_tran(c_tr, t_d, gstp, t, um)
            bd[0,18] = bd[0,17] + bgstp
            bd[0,19] = bd[0,18] + dgstp
            bd[0,20] = bd[0,19] + drug_sen_tran(c_ts, t_d, gstp, t, um)
            bd[0,21] = bd[0,20] + 0
            bd[0,22] = bd[0,21] + drug_sen_tran(c_ts, t_d, gstp, t, um)
            bd[0,23] = bd[0,22] + drug_sen_tran(c_ts, t_d, gstp, t, um)
            bd[0,24] = bd[0,23] + drug_sen_tran(c_ts, t_d, gstr, t, um)
            bd[0,25] = bd[0,24] + drug_res_tran(c_tr, t_d, gstr, t, um)
            bd[0,26] = bd[0,25] + drug_sen_tran(c_ts, t_d, gstr, t, um)
            bd[0,27] = bd[0,26] + bgstr
            bd[0,28] = bd[0,27] + dgstr
            bd[0,29] = bd[0,28] + 0
            bd[0,30] = bd[0,29] + drug_sen_tran(c_ts, t_d, gstr, t, um)
            bd[0,31] = bd[0,30] + drug_sen_tran(c_ts, t_d, gstr, t, um)
            bd[0,32] = bd[0,31] 
            bd[0,33] = bd[0,32] 
            bd[0,34] = bd[0,33] 
            bd[0,35] = bd[0,34] 
            bd[0,36] = bd[0,35] + bgpr
            bd[0,37] = bd[0,36] + dgpr
            bd[0,38] = bd[0,37] 
            bd[0,39] = bd[0,38] 
            bd[0,40] = bd[0,39] + drug_sen_tran(c_ts, t_d, gptp, t, um)
            bd[0,41] = bd[0,40] + 0
            bd[0,42] = bd[0,41] + drug_sen_tran(c_ts, t_d, gptp, t, um)
            bd[0,43] = bd[0,42] + drug_sen_tran(c_ts, t_d, gptp, t, um)
            bd[0,44] = bd[0,43] + drug_res_tran(c_tr, t_d, gptp, t, um)
            bd[0,45] = bd[0,44] + bgptp
            bd[0,46] = bd[0,45] + dgptp    
            bd[0,47] = bd[0,46] + drug_sen_tran(c_ts, t_d, gptp, t, um)
            bd[0,48] = bd[0,47] + drug_sen_tran(c_ts, t_d, gptr, t, um)
            bd[0,49] = bd[0,48] + 0
            bd[0,50] = bd[0,49] + drug_sen_tran(c_ts, t_d, gptr, t, um)
            bd[0,51] = bd[0,50] + drug_sen_tran(c_ts, t_d, gptr, t, um)
            bd[0,52] = bd[0,51] + drug_res_tran(c_tr, t_d, gptr, t, um)
            bd[0,53] = bd[0,52] + drug_sen_tran(c_ts, t_d, gptr, t, um)
            bd[0,54] = bd[0,53] + bgptr
            bd[0,55] = bd[0,54] + dgptr
             
            tot = bd[0,55] 
        
            t = t - math.log(u1)/tot          #Gillespie algorithm for time step
        #    t = t + 0.04
            
            if   u2 < (bd[0,0]/tot):
                 q            = single_step_birth(q)        
            elif (bd[0,1]/tot) > u2 > (bd[0,0]/tot):
                 q            = single_step_death(q)
            elif (bd[0,2]/tot) > u2 > (bd[0,1]/tot):
                 (q, gsr)     = single_step_trans(q, gsr)
            elif (bd[0,3]/tot) > u2 > (bd[0,2]/tot):
                 (q, gstp)    = single_step_trans(q, gstp)
            elif (bd[0,4]/tot) > u2 > (bd[0,3]/tot):
                 (q, gstr)    = single_step_trans(q, gstr)
            elif (bd[0,5]/tot) > u2 > (bd[0,4]/tot):
                 (q, gpr)     = single_step_trans(q, gpr)
            elif (bd[0,6]/tot) > u2 > (bd[0,5]/tot):
                 (q, gptp)    = single_step_trans(q, gptp)
            elif (bd[0,7]/tot) > u2 > (bd[0,6]/tot):
                 (q, gptr)    = single_step_trans(q, gptr)
            elif (bd[0,8]/tot) > u2 > (bd[0,7]/tot):
                 (gsr, q)     = single_step_trans(gsr, q)
            elif (bd[0,9]/tot) > u2 > (bd[0,8]/tot):
                 gsr          = single_step_birth(gsr)  
            elif (bd[0,10]/tot) > u2 > (bd[0,9]/tot):
                 gsr          = single_step_death(gsr)    
            elif (bd[0,11]/tot) > u2 > (bd[0,10]/tot):
                 (gsr, gstp)  = single_step_trans(gsr, gstp)
            elif (bd[0,12]/tot) > u2 > (bd[0,11]/tot):
                 (gsr, gstr)  = single_step_trans(gsr, gstr)
            elif (bd[0,13]/tot) > u2 > (bd[0,12]/tot):
                 (gsr, gpr)   = single_step_trans(gsr, gpr)
            elif (bd[0,14]/tot) > u2 > (bd[0,13]/tot):
                 (gsr, gptp)  = single_step_trans(gsr, gptp)
            elif (bd[0,15]/tot) > u2 > (bd[0,14]/tot):
                 (gsr, gptr)  = single_step_trans(gsr, gptr)
            elif (bd[0,16]/tot) > u2 > (bd[0,15]/tot):
                 (gstp, q)    = single_step_trans(gstp, q)
            elif (bd[0,17]/tot) > u2 > (bd[0,16]/tot) :
                 (gstp, gsr)  = multi_step_trans(bd_ms[:, 0], inst_fun, u4, gstp, gsr)
            elif (bd[0,18]/tot) > u2 > (bd[0,17]/tot) :
                 gstp         = single_step_birth(gstp)
            elif (bd[0,19]/tot) > u2 > (bd[0,18]/tot) :
                 gstp         = single_step_death(gstp)
            elif (bd[0,20]/tot) > u2 > (bd[0,19]/tot) :
                 (gstp, gstr) = single_step_trans(gstp, gstr)
            elif (bd[0,21]/tot) > u2 > (bd[0,20]/tot) :
                 (gstp, gpr)  = single_step_trans(gstp, gpr)
            elif (bd[0,22]/tot) > u2 > (bd[0,21]/tot) :
                 (gstp, gptp) = single_step_trans(gstp, gptp)
            elif (bd[0,23]/tot) > u2 > (bd[0,22]/tot):
                 (gstp, gptr) = single_step_trans(gstp, gptr)
            elif (bd[0,24]/tot) > u2 > (bd[0,23]/tot):
                 (gstr, q)    = single_step_trans(gstr, q)
            elif (bd[0,25]/tot) > u2 > (bd[0,24]/tot) :
                 (gstr, gsr)  = multi_step_trans(bd_ms[:, 1], inst_fun, u4, gstr, gsr)
            elif (bd[0,26]/tot) > u2 > (bd[0,25]/tot) :
                 (gstr, gstp) = single_step_trans(gstr, gstp)
            elif (bd[0,27]/tot) > u2 > (bd[0,26]/tot) :
                 gstr         = single_step_birth(gstr)
            elif (bd[0,28]/tot) > u2 > (bd[0,27]/tot) :
                 gstr         = single_step_death(gstr)
            elif (bd[0,29]/tot) > u2 > (bd[0,28]/tot) :
                 (gstr, gpr)  = single_step_trans(gstr, gpr)
            elif (bd[0,30]/tot) > u2 > (bd[0,29]/tot) :
                 (gstr, gptp) = single_step_trans(gstr, gptp)
            elif (bd[0,31]/tot) > u2 > (bd[0,30]/tot):
                 (gstr, gptr) = single_step_trans(gstr, gptr)
            elif (bd[0,32]/tot) > u2 > (bd[0,31]/tot):
                 (gpr, q)     = single_step_trans(gpr, q)
            elif (bd[0,33]/tot) > u2 > (bd[0,32]/tot) :
                 (gpr, gsr)   = single_step_trans(gpr, gsr)
            elif (bd[0,34]/tot) > u2 > (bd[0,33]/tot) :
                 (gpr, gstp)  = single_step_trans(gpr, gstp)
            elif (bd[0,35]/tot) > u2 > (bd[0,34]/tot) :
                 (gpr, gstr)  = single_step_trans(gpr, gstr)
            elif (bd[0,36]/tot) > u2 > (bd[0,35]/tot) :
                 gpr          = single_step_birth(gpr)
            elif (bd[0,37]/tot) > u2 > (bd[0,36]/tot) :
                 gpr          = single_step_death(gpr)  
            elif (bd[0,38]/tot) > u2 > (bd[0,37]/tot) :
                 (gpr, gptp)  = single_step_trans(gpr, gptp)
            elif (bd[0,39]/tot) > u2 > (bd[0,38]/tot):
                 (gpr, gptr)  = single_step_trans(gpr, gptr)
            elif (bd[0,40]/tot) > u2 > (bd[0,39]/tot):
                 (gptp, q)    = single_step_trans(gptp, q)
            elif (bd[0,41]/tot) > u2 > (bd[0,40]/tot) :
                 (gptp, gsr)  = single_step_trans(gptp, gsr)
            elif (bd[0,42]/tot) > u2 > (bd[0,41]/tot) :
                 (gptp, gstp) = single_step_trans(gptp, gstp)
            elif (bd[0,43]/tot) > u2 > (bd[0,42]/tot) :
                 (gptp, gstr) = single_step_trans(gptp, gstr)
            elif (bd[0,44]/tot) > u2 > (bd[0,43]/tot) :
                 (gptp, gpr)  = multi_step_trans(bd_ms[:, 2], inst_fun, u4, gptp, gpr)
            elif (bd[0,45]/tot) > u2 > (bd[0,44]/tot) :
                 gptp         = single_step_birth(gptp)
            elif (bd[0,46]/tot) > u2 > (bd[0,45]/tot) :
                 gptp         = single_step_death(gptp)    
            elif (bd[0,47]/tot) > u2 > (bd[0,46]/tot):
                 (gptp, gptr) = single_step_trans(gptp, gptr)
            elif (bd[0,48]/tot) > u2 > (bd[0,47]/tot):
                 (gptr, q)    = single_step_trans(gptr, q)
            elif (bd[0,49]/tot) > u2 > (bd[0,48]/tot) :
                 (gptr, gsr)  = single_step_trans(gptr, gsr)
            elif (bd[0,50]/tot) > u2 > (bd[0,49]/tot) :
                 (gptr, gstp) = single_step_trans(gptr, gstp)
            elif (bd[0,51]/tot) > u2 > (bd[0,50]/tot) :
                 (gptr, gstr) = single_step_trans(gptr, gstr)
            elif (bd[0,52]/tot) > u2 > (bd[0,51]/tot) :
                 (gptr, gpr)  = multi_step_trans(bd_ms[:, 3], inst_fun, u4, gptr, gpr)
            elif (bd[0,53]/tot) > u2 > (bd[0,52]/tot) :
                 (gptr, gptp) = single_step_trans(gptr, gptp)
            elif (bd[0,54]/tot) > u2 > (bd[0,53]/tot):
                 gptr         = single_step_birth(gptr)
            elif (bd[0,55]/tot) > u2 > (bd[0,54]/tot):
                 gptr         = single_step_death(gptr)   
            else:
                 print("error")
                 
            totcell = q + gsr + gstp + gstr + gpr + gptp + gptr  
            
            i = i + 1
    
            q_vector.append(q)
            gsr_vector.append(gsr)
            gstp_vector.append(gstp)
            gstr_vector.append(gstr)
            gpr_vector.append(gpr)
            gptp_vector.append(gptp)
            gptr_vector.append(gptr)
            totcell_vector.append(totcell)
            t_vector.append(t)
    #    transition_temp = transition_temp + bd_ms[3, 0]+ bd_ms[3, 1]+ bd_ms[3, 2]+ bd_ms[3, 3]
        gsr_ave_temp = gsr_vector[-1]
        gsr_ave.append(gsr_ave_temp)
        gpr_ave_temp = gpr_vector[-1]
        gpr_ave.append(gpr_ave_temp)
        
        print(k)
    #    average_transition = transition_temp/k
    average_gsr = sum(gsr_ave)/len(gsr_ave)    
    average_gpr = sum(gpr_ave)/len(gpr_ave)  
    per_m = m*100/(m+1)
    average_gsr_mat.append(average_gsr)
    average_gpr_mat.append(average_gpr)
    m_mat.append(m)
    per_m_mat.append(per_m)

#plt.plot(per_m_mat, average_gpr_mat, label="G-R")
#plt.plot(per_m_mat, average_gsr_mat, label="GS-R")
plt.plot(t_vector, q_vector, label="Q cells")
plt.plot(t_vector, gstp_vector, label="GS-Tp cells")
plt.plot(t_vector, gstr_vector, label="GS-Tr cells")
plt.plot(t_vector, gptp_vector, label="G-Tp cells")
plt.plot(t_vector, gptr_vector, label="G-Tr cells")
plt.plot(t_vector, gsr_vector, label="GS-R cells")
plt.plot(t_vector, gpr_vector, label="G-R cells")

plt.legend(loc="upper right")
plt.xlabel('Time (Days)')  
plt.ylabel('Cell number') 
plt.xlim(0, 6.1)    
#plt.ylim(0, 200) 
plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
plt.show()    
  
