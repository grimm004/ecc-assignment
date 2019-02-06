#function hammingGeneratorMatrix
#input: a number r
#output: G, the generator matrix of the (2^r-1,2^r-r-1) Hamming code
def hammingGeneratorMatrix(r):
    n = 2**r-1
    
    #construct permutation pi
    pi = []
    for i in range(r):
        pi.append(2**(r-i-1))
    for j in range(1,r):
        for k in range(2**j+1,2**(j+1)):
            pi.append(k)

    #construct rho = pi^(-1)
    rho = []
    for i in range(n):
        rho.append(pi.index(i+1))

    #construct H'
    H = []
    for i in range(r,n):
        H.append(decimalToVector(pi[i],r))

    #construct G'
    GG = [list(i) for i in zip(*H)]
    for i in range(n-r):
        GG.append(decimalToVector(2**(n-r-i-1),n-r))

    #apply rho to get Gtranpose
    G = []
    for i in range(n):
        G.append(GG[rho[i]])

    #transpose    
    G = [list(i) for i in zip(*G)]

    return G


#function decimalToVector
#input: numbers n and r (0 <= n<2**r)
#output: a string v of r bits representing n
def decimalToVector(n,r): 
    v = []
    for s in range(r):
        v.insert(0,n%2)
        n //= 2
    return v


def vectorToDecimal(v):
    return sum(v[i] * (2 ** (len(v) - i - 1)) for i in range(len(v)))


def vectorToMatrix(v):
    return [v]


def matrixToVector(m):
    return m[0]


def getH(n):
    return [[(j >> (n - i - 1)) & 1 for j in range(1, 2 ** n)] for i in range(n)]


def getN(r):
    return (2 ** r) - 1


def getK(r):
    return getN(r) - r


def validateCodeword(l):
    r = 2
    k = 1
    while True:
        k = getK(r)
        if l == k:
            break
        elif l < k:
            return (False, 0, 0)
        r += 1
    return (True, r, k)


def validateHammingcode(l):
    r = 2
    n = 3
    while True:
        n = getN(r)
        if l == n:
            break
        elif l < n:
            return (False, 0, 0)
        r += 1
    return (True, r, n)


def message(a):
    l = len(a)
    r = 2
    k = 1
    while True:
        k = getK(r)
        if not k - r < l: break
        r += 1
    l_bin = decimalToVector(l, r)
    return [(l_bin[i] if i < r else (a[i - r] if i < r + l else 0)) for i in range(k)]


def multMat(m0, m1):
    if len(m0) == 0 or len(m1) == 0 or (len(m0[0]) != len(m1)): return []
    return [[sum([(m0[i][k] * m1[k][j]) for k in range(len(m1))]) % 2 for j in range(len(m1[0]))] for i in range(len(m0))]


def zeroMat(m, n):
    return [[0 for _ in range(n)] for _ in range(m)]


def identMat(n):
    return [[1 if i == j else 0 for i in range(n)] for j in range(n)]


def transposeMat(m):
    return [[m[j][i] for j in range(len(m))] for i in range(len(m[0]))]


def hammingEncoder(m):
    valid, r, _ = validateCodeword(len(m))
    return matrixToVector(multMat(vectorToMatrix(m), hammingGeneratorMatrix(r))) if valid else []


def hammingDecoder(v):
    valid, r, _ = validateHammingcode(len(v))
    if not valid: return []
    i = vectorToDecimal(matrixToVector(multMat(vectorToMatrix(v), transposeMat(getH(r))))) - 1
    v[i] = (v[i] + 1) % 2
    return v


def messageFromCodeword(c):
    valid, r, _ = validateHammingcode(len(c))
    if not valid: return []
    return [c[(2 ** i) - 1] for i in range(r)]


def dataFromMessage(m):
    valid, r, k = validateCodeword(len(m))
    if not valid: return []

    l = vectorToDecimal(m[0:r])
    if l > k - r or 1 in [m[i] for i in range(r + l, len(m))]:
        return []

    return [m[i] for i in range(r, r + l)]


def repetitionEncoder(m, n):
    if len(m) != 1: return []
    return [m[0] for _ in range(n)]


def repetitionDecoder(v):
    if len(v) == 0: return []
    zero_count = sum(1 for b in v if b == 0)
    one_count = len(v) - zero_count
    if zero_count == one_count: return []
    return [0 if zero_count > one_count else 1]


def printMat(m):
    if len(m) == 0:
        print("NaM")
        return
    mt = "["
    for i in range(len(m)):
        mt += "["
        for j in range(len(m[0])):
            mt += "%s, " % m[i][j]
        mt = mt[:-2]
        mt += "],\n "
    mt = mt[:-3]
    mt += "]"
    print(mt)


if __name__ == '__main__':
    print(messageFromCodeword(hammingDecoder(hammingEncoder(message([0, 1, 1])))))
    """
    print(dataFromMessage(message([1])) == [1])
    print(dataFromMessage(message([0, 0, 1])) == [0, 0, 1])
    print(dataFromMessage(message([0, 1, 1, 0])) == [0, 1, 1, 0])
    print(dataFromMessage(message([1, 1, 1, 1, 0, 1])) == [1, 1, 1, 1, 0, 1])
    print(dataFromMessage(message([0, 1, 1, 0, 1])) == [0, 1, 1, 0, 1])
    for i in range(16):
        a = (i >> 3) & 1
        b = (i >> 2) & 1
        c = (i >> 1) & 1
        d = (i >> 0) & 1
        v = [a, b, c, d]
        encoded = hammingEncoder(v)
        hammingDecoder(encoded)
        #print(v, encoded)
    """
