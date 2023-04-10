import numpy as np
import scipy.stats
import matplotlib.pyplot as plt


def set_gen(a1, a2, a3):
    coeff_a = []
    for i in range(len(a1)):

        for j in range(len(a2)):

            for k in range(len(a3)):
                n = [a1[i], a2[j], a3[k]]
                coeff_a.append(n)
    return(coeff_a)


def coord(X4_apparent, A, X):
    #ai = fractional contribution of corresponding intermediate specie on the total spectra (ideally 0.25, 0.5, 0.75)
    #xi = fraction of corresponding intermediate specie present in the total spectra
    a1 = A[0]
    a2 = A[1]
    a3 = A[2]
    x1 = X[0]
    x2 = X[1]
    x3 = X[2]
    
    X0_apparent = 1 - X4_apparent
    
    x4 = X4_apparent - (a1*x1 + a2*x2 + a3*x3)
    
    x0 = X0_apparent - ((1-a1)*x1 + (1-a2)*x2 + (1-a3)*x3)
    
    N_Cl = []
    if 0 <= x4 <= 1:
        if 0 <= x0 <= 1:
            N = x1 + 2*x2 + 3*x3 + 4*x4
            N_Cl.append(N)
            return(N_Cl)



def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return h
        
def final(X4_apparent, coeff_a, coeff_x):
    lst = []
    for i in range(len(coeff_a)):
        for j in range(len(coeff_x)):
            lst.append(coord(X4_apparent, coeff_a[i], coeff_x[j]))

    lst2 = [x for x in lst if x]
#print(len(lst))
#print(len(lst2))
#print(lst2)

    lst3 = [x for l in lst2 for x in l]


#print(lst3)

    mx = max(lst3)
    mn = min(lst3)

    rng = mx -mn
    
    conf_int = mean_confidence_interval(lst3, confidence=0.95)
    
    med = np.median(lst3)
    p15 = np.percentile(lst3, 15)
    p5 = np.percentile(lst3, 5)
    p95 = np.percentile(lst3, 95)
    p85 = np.percentile(lst3, 85)
    
    #plot histogram of N_Cl
    n_cl = np.array(lst3)
    plt.hist(n_cl, bins = 100)
    plt.show()
            

    #print(f'Total number of satisfying N_Cl values found = {len(lst3)}')
    #print(f'Max = {mx}', f'Min = {mn}', f'Range = {rng}')
    #print(f'Mean = {np.mean(lst3)}', f'SD = {np.std(lst3)}')
    
    lst4 = [X4_apparent, len(lst3), np.mean(lst3), np.std(lst3), mn, mx, rng, conf_int, med, p5, p15, p85, p95]
    
    return(lst4)
    #return(X4_apparent, len(lst3), np.mean(lst3), np.std(lst3), mn, mx, rng)


#writing fitted data to file
def write_to_file(data, name, delimiter='\t'): 
    # for csv file, change delimiter to ','
    if delimiter is ',':
        ext = '.csv'
    else:
        ext = '.txt'
    
    # Write the file
    with open(f'{name}{ext}', 'w') as file:
        # Headers first
        file.write(f'X4_apparent{delimiter}N_data{delimiter}Mean{delimiter}S.D.{delimiter}Min{delimiter}Max{delimiter}Range{delimiter}conf_int{delimiter}median{delimiter}5percentile{delimiter}15percentile{delimiter}85percentile{delimiter}95percentile\n')
        # Then the data
        for i in range(len(data)):
            file.write(f'{data[i][0]}{delimiter}{data[i][1]}{delimiter}{data[i][2]}{delimiter}{data[i][3]}{delimiter}{data[i][4]}{delimiter}{data[i][5]}{delimiter}{data[i][6]}{delimiter}{data[i][7]}{delimiter}{data[i][8]}{delimiter}{data[i][9]}{delimiter}{data[i][10]}{delimiter}{data[i][11]}{delimiter}{data[i][12]}\n')