from collections import namedtuple
from utils.matrix import IdentityMatrix, MM, MM_GROUP, MMMatrix


class PS():

    def gen(self, pars):
        raise NotImplementedError

    def prove(self, pk, x, w):
        raise NotImplementedError

    def verify(self, pk, sk, x, proof):
        raise NotImplementedError


BPS_PK = namedtuple('BPS_PK', ['kDI', 'kID', 'kDD'])
BPS_SK = namedtuple('BPS_SK', ['k'])


class BenignPS(PS):

    def __init__(self, mscheme: MM):
        super().__init__()
        self.mscheme = mscheme

    def gen(self, D: MMMatrix):
        n, d = D.dimension
        k = self.mscheme.sample(n ** 2)
        sk = BPS_SK(k=k)

        I = self.mscheme.sample_from(IdentityMatrix(n, MM_GROUP))
        pk = BPS_PK(kDI=k.T() * (D @ I),
                    kID=k.T() * (I @ D),
                    kDD=k.T() * (D @ D)
                    )
        return pk, sk

    def prove(self, pk: BPS_PK, u, r):
        return pk.kDD * (r @ r)

    def verify(self, sk: BPS_SK, u, proof):
        return proof == sk.k.T() * (u @ u)

    def peval(self, pk: BPS_PK, u, proof, r):
        return proof + (pk.kID * (u @ r)) + (pk.kDI * (r @ u)) + (pk.kDD * (r @ r))
