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
    
  
file_path = "/Users/ericyanchenko/Documents/Research/Srijan/CP_clean/Data/ukfaculty_adj.RData"
parsed = rdata.parser.parse_file(file_path)
converted = rdata.conversion.convert(parsed)
A = list(converted.values())[0]

start = time.time()
C = cp_be(A)
end = time.time()

print("Network: UK Faculty")
print("Obj fun:", round(BE(A,C),2))
print("Runtime:", round(end - start,3))
print("Coresize:", sum(C))

file_path = "/Users/ericyanchenko/Documents/Research/Srijan/CP_clean/Data/BritishMP_adj.RData"
parsed = rdata.parser.parse_file(file_path)
converted = rdata.conversion.convert(parsed)
A = list(converted.values())[0]

start = time.time()
C = cp_be(A)
end = time.time()

print("Network: British MP")
print("Obj fun:", round(BE(A,C),2))
print("Runtime:", round(end - start,3))
print("Coresize:", sum(C))

file_path = "/Users/ericyanchenko/Documents/Research/Srijan/CP_clean/Data/polblogs_adj.RData"
parsed = rdata.parser.parse_file(file_path)
converted = rdata.conversion.convert(parsed)
A = list(converted.values())[0]

start = time.time()
C = cp_be(A)
end = time.time()

print("Network: Political Blogs")
print("Obj fun:", round(BE(A,C),2))
print("Runtime:", round(end - start,3))
print("Coresize:", sum(C))

file_path = "/Users/ericyanchenko/Documents/Research/Srijan/CP_clean/Data/school_adj.RData"
parsed = rdata.parser.parse_file(file_path)
converted = rdata.conversion.convert(parsed)
A = list(converted.values())[0]

start = time.time()
C = cp_be(A)
end = time.time()

print("Network: School")
print("Obj fun:", round(BE(A,C),2))
print("Runtime:", round(end - start,3))
print("Coresize:", sum(C))


file_path = "/Users/ericyanchenko/Documents/Research/Srijan/CP_clean/Data/dblp_adj.RData"
parsed = rdata.parser.parse_file(file_path)
converted = rdata.conversion.convert(parsed)
A = list(converted.values())[0]

start = time.time()
C = cp_be(A)
end = time.time()

print("Network: DBLP")
print("Obj fun:", round(BE(A,C),2))
print("Runtime:", round(end - start,3))
print("Coresize:", sum(C))

file_path = "/Users/ericyanchenko/Documents/Research/Srijan/CP_clean/Data/bio1_adj.RData"
parsed = rdata.parser.parse_file(file_path)
converted = rdata.conversion.convert(parsed)
A = list(converted.values())[0]

start = time.time()
C = cp_be(A)
end = time.time()

print("Network: Biological 1")
print("Obj fun:", round(BE(A,C),2))
print("Runtime:", round(end - start,3))
print("Coresize:", sum(C))

file_path = "/Users/ericyanchenko/Documents/Research/Srijan/CP_clean/Data/FB_adj.RData"
parsed = rdata.parser.parse_file(file_path)
converted = rdata.conversion.convert(parsed)
A = list(converted.values())[0]

start = time.time()
C = cp_be(A)
end = time.time()

print("Network: Facebook")
print("Obj fun:", round(BE(A,C),2))
print("Runtime:", round(end - start,3))
print("Coresize:", sum(C)) 
    
