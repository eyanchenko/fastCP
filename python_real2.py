import cpnet
import numpy as np
import networkx as nx
import rdata
import time

def cp_be(A):
    
    algorithm = cpnet.BE()
    G = nx.from_numpy_array(A)
    algorithm.detect(G)
    x = algorithm.get_coreness()
    
    return list(x.values())
    
def BE(A, C):
    n = len(C)
    k = sum(C)
    phat = sum(sum(A) / (n*(n-1)))
    dhat = (0.5*k*(k-1) + k*(n-k)) / (0.5*n*(n-1))

    obj = 0
    for i in range(n):
        for j in range(n):
            obj += A[i,j] * (C[i] + C[j] - C[i]*C[j])/2

    obj  = obj -  0.5*n*(n-1)*phat*dhat
    obj  = obj / (0.5*n*(n-1)*(phat*(1-phat)*dhat*(1-dhat))**(0.5))

    return obj
    
   

file_path = "/Users/ericyanchenko/Documents/Research/Srijan/CP_clean/Data/email4_adj.RData"
parsed = rdata.parser.parse_file(file_path)
converted = rdata.conversion.convert(parsed)
A = list(converted.values())[0]

start = time.time()
C = cp_be(A)
end = time.time()

print("Network: Email 4")
print("Obj fun:", round(BE(A,C),2))
print("Runtime:", round(end - start,3))
print("Coresize:", sum(C))

file_path = "/Users/ericyanchenko/Documents/Research/Srijan/CP_clean/Data/congress_adj.RData"
parsed = rdata.parser.parse_file(file_path)
converted = rdata.conversion.convert(parsed)
A = list(converted.values())[0]

start = time.time()
C = cp_be(A)
end = time.time()

print("Network: Twitter Congress")
print("Obj fun:", round(BE(A,C),2))
print("Runtime:", round(end - start,3))
print("Coresize:", sum(C))

file_path = "/Users/ericyanchenko/Documents/Research/Srijan/CP_clean/Data/hosp_adj.RData"
parsed = rdata.parser.parse_file(file_path)
converted = rdata.conversion.convert(parsed)
A = list(converted.values())[0]

start = time.time()
C = cp_be(A)
end = time.time()

print("Network: Hospital")
print("Obj fun:", round(BE(A,C),2))
print("Runtime:", round(end - start,3))
print("Coresize:", sum(C))

file_path = "/Users/ericyanchenko/Documents/Research/Srijan/CP_clean/Data/copenBT_adj.RData"
parsed = rdata.parser.parse_file(file_path)
converted = rdata.conversion.convert(parsed)
A = list(converted.values())[0]

start = time.time()
C = cp_be(A)
end = time.time()

print("Network: Copenhagen BT")
print("Obj fun:", round(BE(A,C),2))
print("Runtime:", round(end - start,3))
print("Coresize:", sum(C))

file_path = "/Users/ericyanchenko/Documents/Research/Srijan/CP_clean/Data/bio2_adj.RData"
parsed = rdata.parser.parse_file(file_path)
converted = rdata.conversion.convert(parsed)
A = list(converted.values())[0]

start = time.time()
C = cp_be(A)
end = time.time()

print("Network: Biological 2")
print("Obj fun:", round(BE(A,C),2))
print("Runtime:", round(end - start,3))
print("Coresize:", sum(C))

file_path = "/Users/ericyanchenko/Documents/Research/Srijan/CP_clean/Data/bio3_adj.RData"
parsed = rdata.parser.parse_file(file_path)
converted = rdata.conversion.convert(parsed)
A = list(converted.values())[0]

start = time.time()
C = cp_be(A)
end = time.time()

print("Network: Biological 3")
print("Obj fun:", round(BE(A,C),2))
print("Runtime:", round(end - start,3))
print("Coresize:", sum(C))
