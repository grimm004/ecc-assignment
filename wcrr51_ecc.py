#function HammingG
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


def getK(r):
    return (2 ** r) - r - 1


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


def hammingEndocder(m):
    return []


def hammingDecoder(v):
    return []


def messageFromCodeword(c):
    return []


def dataFromMessage(m):
    r = 2
    k = 1
    while True:
        k = getK(r)
        if len(m) == k:
            break
        elif len(m) < k:
            return []
        r += 1

    l = sum(m[i] * (2 ** (r - i - 1)) for i in range(r))
    if l > k - r:
        return []

    return [m[i] for i in range(r, r + l)]


def repetitionEncoder(m, n):
    return []


def repetitionDecode(v):
    return []


if __name__ == '__main__':
    print(dataFromMessage(message([1])) == [1])
    print(dataFromMessage(message([0, 0, 1])) == [0, 0, 1])
    print(dataFromMessage(message([0, 1, 1, 0])) == [0, 1, 1, 0])
    print(dataFromMessage(message([1, 1, 1, 1, 0, 1])) == [1, 1, 1, 1, 0, 1])
    print(dataFromMessage(message([0, 1, 1, 0, 1])) == [0, 1, 1, 0, 1])
    
    print(dataFromMessage([1, 0, 0, 1, 1, 0, 1, 0]))
    print(dataFromMessage([1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0]))
    print(dataFromMessage([0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0]))
    #print(dataFromMessage(message([0, 1, 1, 0, 1])))

