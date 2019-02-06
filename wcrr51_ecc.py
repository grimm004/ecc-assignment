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


def printMatrix(m):
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


vectorToDecimal = lambda v: sum(v[i] * (2 ** (len(v) - i - 1)) for i in range(len(v)))
vectorToMatrix = lambda v: [v]
matrixToVector = lambda m: m[0] if len(m) > 0 else []

getH = lambda n: [[(j >> (n - i - 1)) & 1 for j in range(1, 2 ** n)] for i in range(n)]
getN = lambda r: (2 ** r) - 1
getK = lambda r: getN(r) - r

multiplyMatrix = lambda m0, m1: [] if len(m0) == 0 or len(m1) == 0 or (len(m0[0]) != len(m1)) else [[sum([(m0[i][k] * m1[k][j]) for k in range(len(m1))]) % 2 for j in range(len(m1[0]))] for i in range(len(m0))]
zeroMatrix = lambda m, n: [[0 for _ in range(n)] for _ in range(m)]
identityMatrix = lambda n: [[1 if i == j else 0 for i in range(n)] for j in range(n)]
transposeMat = lambda m: [[m[j][i] for j in range(len(m))] for i in range(len(m[0]))]


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


def hammingEncoder(m):
    valid, r, _ = validateCodeword(len(m))
    return matrixToVector(multiplyMatrix(vectorToMatrix(m), hammingGeneratorMatrix(r))) if valid else []


def hammingDecoder(v):
    v = v.copy()
    valid, r, _ = validateHammingcode(len(v))
    if not valid: return []
    i = vectorToDecimal(matrixToVector(multiplyMatrix(vectorToMatrix(v), transposeMat(getH(r))))) - 1
    if i > -1: v[i] = (v[i] + 1) % 2
    return v


def messageFromCodeword(c):
    valid, _, n = validateHammingcode(len(c))
    if not valid: return []
    m = [c[i] for i in range(n) if not (((i + 1) & i) == 0)]
    return m


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
    one_count = sum(v)
    zero_count = len(v) - one_count
    if zero_count == one_count: return []
    return [0 if zero_count > one_count else 1]


if __name__ == '__main__':
    from random import randint
    allPass = True
    for i in range(100000):
        numberToSend = randint(0, (2**16) - 1)
        nBits = 16
        originalData = decimalToVector(numberToSend, nBits)
        encodedMessage = message(originalData)
        hammingEncodedMessage = hammingEncoder(encodedMessage)

        receivedHammingMessage = hammingEncodedMessage.copy()
        receivedHammingMessage[randint(0, len(receivedHammingMessage) - 1)] ^= 1

        hammingDecodedMessage = hammingDecoder(receivedHammingMessage)
        decodedMessage = messageFromCodeword(hammingDecodedMessage)
        receivedData = dataFromMessage(decodedMessage)
        receivedNumber = vectorToDecimal(receivedData)

        allPass &= numberToSend == receivedNumber

    print("all passed" if allPass else "all failed")








































    """
    correct = True
    correct &= dataFromMessage(message([1])) == [1]
    correct &= dataFromMessage(message([0, 0, 1])) == [0, 0, 1]
    correct &= dataFromMessage(message([0, 1, 1, 0])) == [0, 1, 1, 0]
    correct &= dataFromMessage(message([1, 1, 1, 1, 0, 1])) == [1, 1, 1, 1, 0, 1]
    correct &= dataFromMessage(message([0, 1, 1, 0, 1])) == [0, 1, 1, 0, 1]
    print("message() and dataFromMessage():", "pass" if correct else "fail")
    
    correct = True
    correct &= hammingEncoder([1, 1, 1]) == []
    correct &= hammingEncoder([1, 0, 0, 0]) == [1, 1, 1, 0, 0, 0, 0]
    correct &= hammingDecoder([1, 0, 1, 1]) == []
    correct &= hammingDecoder([0, 1, 1, 0, 0, 0, 0]) == [1, 1, 1, 0, 0, 0, 0]
    print("hammingEncoder() and hammingDecoder():", "pass" if correct else "fail")
    
    correct = True
    correct &= messageFromCodeword([1, 1, 0, 1]) == []
    correct &= messageFromCodeword([1, 1, 1, 0, 0, 0, 0]) == [1, 0, 0, 0]
    print("messageFromCodeword():", "pass" if correct else "fail")
    
    correct = True
    correct &= repetitionEncoder([0], 4) == [0, 0, 0, 0]
    correct &= repetitionDecoder([1, 0, 0, 0]) == [0]
    print("repetitionEncoder() and repetitionDecoder():", "pass" if correct else "fail")
    
    for i in range(16):
        a = (i >> 0) & 1
        b = (i >> 1) & 1
        c = (i >> 2) & 1
        d = (i >> 3) & 1
        v = [a, b, c, d]
        encoded = hammingEncoder(v)
        hammingDecoder(encoded)
        #print(v, encoded)
    """
