from charm.toolbox.PKEnc import PKEnc
from charm.toolbox.pairinggroup import ZR, G1, G2, GT, pair

from collections import namedtuple
from functools import reduce

from utils.matrix import MM, MatrixDistribution

FFHR_SK = namedtuple('FFHR_SK', ['a', 'f', 'g', 'F', 'G'])
FFHR_PK = namedtuple(
    'FFHR_PK', ['D', 'E', 'aTD', 'fTD', 'FTD', 'gTE', 'GTE', 'GDx', 'FE'])
FFHR_CIPHER = namedtuple('FFHR_CIPHER', ['u', 'p', 'v', 'pi'])


class FFHR19(PKEnc):

    def __init__(self, mscheme: MM, k: int, dist: MatrixDistribution):
        super().__init__()
        self.mscheme = mscheme
        self.k = k
        self.dist = dist

    def keygen(self):
        a, f, g = (self.mscheme.sample(self.k+1) for i in range(3))
        F = self.mscheme.sample(self.k+1, self.k+1)
        G = self.mscheme.sample(self.k+1, self.k+2)
        D, E = (self.mscheme.sample_from(self.dist) for i in range(2))
        sk = FFHR_SK(a=a, f=f, g=g, F=F, G=G)

        aTD = a.T() * D
        Dx = D | aTD

        pk = FFHR_PK(D=D >> G1,
                     E=E >> G2,
                     aTD=aTD >> G1,
                     fTD=(f.T() * D) >> GT,
                     FTD=(F.T() * D) >> G1,
                     gTE=(g.T() * E) >> GT,
                     GDx=(G * Dx) >> G1,
                     GTE=(G.T() * E) >> G2,
                     FE=(F * E) >> G2)
        return pk, sk

    def encrypt(self, pk: FFHR_PK, msg):
        r = self.mscheme.sample(self.k)
        s = self.mscheme.sample(self.k)
        u = pk.D * r
        #msg = self.mscheme.sample(1, gtype=G1)
        p = pk.aTD * r + msg
        x = u | p
        v = pk.E * s
        pi1 = pk.fTD * r + (pk.FTD * r).pair_with(v)
        pi2 = pk.gTE * s + x.pair_with(pk.GTE * s)
        pi = pi1 + pi2
        return FFHR_CIPHER(u=u, p=p, v=v, pi=pi)

    def decrypt(self, pk, sk, c):
        x = c.u | c.p
        msg = c.p - sk.a.T() * c.u
        pi1 = (sk.F * (c.v)).T() * c.u + sk.f.T() * (c.u >> GT)
        pi2 = (sk.G * x).T() * c.v + sk.g.T() * (c.v >> GT)

        if(pi1 + pi2 == c.pi):
            return msg
        return None

    def rand(self, pk, c):
        r, s = (self.mscheme.sample(self.k) for i in range(2))
        u = c.u + pk.D * r
        p = c.p + pk.aTD * r
        x = u | p
        v = c.v + pk.E * s
        pi1 = pk.fTD * r + (pk.FTD * r).pair_with(v) + c.u.pair_with(pk.FE * s)
        pi2 = pk.gTE * s + (pk.GDx * r).pair_with(c.v) + \
            x.pair_with(pk.GTE * s)

        return FFHR_CIPHER(u=u, p=p, v=v, pi=c.pi+pi1+pi2)
