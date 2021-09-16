import sys
import math
import numpy as np
from scipy.sparse import random
from scipy import stats

#######################
## CONSTRAINT VALUES ##
#######################

def generate_cc_data(n,m,omega,density,beta, seed, positive=False):
    # random instances for chance constraints
    np.random.seed(seed)

    ## Density - which varaibles are active or not
    x = random(m,n, density = density, data_rvs=np.ones)
    x = x.toarray()
    x = x.astype('int32') 
    #x = np.random.randint(0, 2, size=(m,n) )

    if positive:
        ## Linear values - vectors a
        a = np.random.randint(1, 100, size=(m,n) )
        a = np.multiply(a,x)

        ## Covariance term 
        d = np.random.randint(1, 40, size=(m,n,n) )
        x_aux = np.empty( (n,n) )

    else:
        ## Linear values - vectors a
        a = np.random.randint(-50, 50, size=(m,n) )
        a = np.multiply(a,x)

        ## Covariance term 
        d = np.random.randint(-20, 20, size=(m,n,n) )
        x_aux = np.empty( (n,n) )



    for i in range(m):
        x_aux[:] = x[i]
        x_aux = np.multiply( x_aux, x_aux.transpose() )
        d[i] = np.multiply(d[i],x_aux)

    ## LHS constant - b vector
    b = np.zeros(m)
    aux = 0

    for i in range(m):
        aux = 0
        # sum quadratic terms
        for j in range(n):
            aux += np.power( max(d[i][j][d[i][j]> 0].sum(), -d[i][j][d[i][j]< 0].sum()), 2)
        
        bmax = a[i][a[i] > 0].sum() + omega*np.sqrt(aux) 
        bmin = a[i][a[i] < 0].sum() + 0 

        b[i] = np.floor(beta*( bmax + bmin ))

    b = b.astype('int32') 

    ## Profit - c 
    c = np.random.randint(0, 100, size=n )

    return a, d, b, c


def generate_ks_data(n,m,omega,density,beta, seed):
    # Following the description in Atamturk et. al 2009

    np.random.seed(seed)

    ## Density - which varaibles are active or not
    x = random(m,n, density = density, data_rvs=np.ones)
    x = x.toarray()
    x = x.astype('int32') 
    #x = np.random.randint(0, 2, size=(m,n) )

    ## Linear values - vectors a
    a = np.random.randint(0, 100, size=(m,n) )
    a = np.multiply(a,x)

    ## Diagonal values - vector d
    d = np.zeros([m,n], dtype='int32')

    for i in range(m):
        for j in range(n):
            if a[i,j] >0:
                d[i,j] = np.random.randint(0,a[i][j])

    ## LHS constant - b vector
    b = np.zeros(m)
    for i in range(m):
        b[i] = beta*np.floor( np.sum(a[i,:]) + omega*np.linalg.norm(d[i,:]))
    b = b.astype('int32') 

    ## Profit - c 
    c = np.random.randint(0, 100, size=n )

    return a, d, b, c


################
## WRITE FILE ##
################
def write_problem(n, m, a,b, c, d, omega, beta, c_type, seed):
    folder = ""
    if c_type== "knapsack":
        folder = "random-diag-knapsack/large/"
    elif c_type== "chance":
        folder = "random-chance-constraint/medium/"

    if n >= 100:
        filename = folder + "%s_n%dm%do%db%d_%d.txt"%(c_type, n, m, omega, beta*10, seed)
    else:
        filename = folder + "%s_n0%dm%do%db%d_%d.txt"%(c_type, n, m, omega, beta*10, seed)

    f = open(filename, "w+")

    print(filename)

    # Write n, m and omega
    f.write("%d %d\n"%(n,m))
    f.write("%d\n"%omega)

    # Write profit vector c
    for i in range(n):
        if i < n-1: f.write("%d "%c[i])
        else:       f.write("%d\n"%c[i])
    
    # Write vector b
    for i in range(m):
        if i < m-1: f.write("%d "%b[i])
        else:       f.write("%d\n"%b[i])
    
    # Write matrix a
    for i in range(m):
        for j in range(n):
            if j < n-1: f.write("%d "%a[i,j])
            else:       f.write("%d\n"%a[i,j])

    # Write matrix d for knapsack
    if c_type== "knapsack" or c_type== "gub-knapsack":
        for i in range(m):
            for j in range(n):
                if j < n-1: f.write("%d "%d[i,j])
                else:       f.write("%d\n"%d[i,j])

    # Write matrix d for chance constraints    
    elif c_type== "chance":
        for i in range(m):
            for j in range(n):
                for k in range(n):
                    if k < n-1: f.write("%d "%d[i,j,k])
                    else:       f.write("%d\n"%d[i,j,k])

    f.close()

####
# Write all knapsack constraints

def write_all_knapsack():
    ##Tenting large instances
    # c_type = "knapsack"
    # beta = 0.5
    # range_n = [175, 200, 225, 250]

    # ##Large data set
    # for n in range_n:
    #     density = 2/math.sqrt(n)
    #     for m in [20]:
    #         for omega in [3]:
    #             for seed in range(5):
    #                 a, d, b, c = generate_ks_data(n, m, omega, density, beta, seed)
    #                 write_problem(n, m, a, b, c, d, omega, beta, c_type, seed)

    ##Standar data set
    c_type = "knapsack"
    beta = 0.5
    range_n = [100, 125, 150]

    for n in range_n:
        density = 2/math.sqrt(n)
        for m in [10, 20]:
            for omega in [1,3,5]:
                for seed in range(5):
                    a, d, b, c = generate_ks_data(n, m, omega, density, beta, seed)
                    write_problem(n, m, a, b, c, d, omega, beta, c_type, seed)



####
# Write all chance constraints

def write_all_chance_constraints():
    c_type = "chance"
    beta = 0.1
    #range_n = [25, 50, 75]
    range_n = [50, 75, 100, 125]

    for n in range_n:
        density = 2/math.sqrt(n)
        for m in [10, 20]:
            for omega in [1,3,5]:
                for seed in range(5):
                    a, d, b, c = generate_cc_data(n,m,omega,density,beta, seed)
                    write_problem(n, m, a, b, c, d, omega, beta, c_type, seed)



####
# Main

#Parameters
#c_type = "knapsack"    # diagonal knapsack
# c_type="chance"          # chance constraints
# n = 5
# m = 1
# density = 1
# seed = 1
# omega = 1
# beta = 0.3 #chance constraints need a beta smaller 

# if c_type == "knapsack" :
#     a, d, b, c = generate_ks_data(n,m,omega,density,beta, seed)
#     write_problem(n, m, a, b, c, d, omega, beta, c_type, seed)
# elif c_type== "chance":
#     a, d, b, c = generate_cc_data(n,m,omega,density,beta, seed, positive=False)
#     write_problem(n, m, a, b, c, d, omega, beta, c_type, seed)


##All Atamturk instances
write_all_knapsack()
write_all_chance_constraints()
