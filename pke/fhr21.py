from charm.toolbox.PKEnc import PKEnc
from charm.toolbox.pairinggroup import G1, GT

from collections import namedtuple

from utils.matrix import MM, MatrixDistribution
from ps.benign import BenignPS

FFHR21_SK = namedtuple('FFHR21_SK', ['a', 'f', 'F', 'psk'])
FFHR21_PK = namedtuple(
    'FFHR21_PK', ['D', 'aTD', 'fTD', 'FTD', 'FDx', 'ppk'])
FFHR21_CIPHER = namedtuple('FFHR21_CIPHER', ['u', 'p', 'y', 'pi'])


class FHR21(PKEnc):

    def __init__(self, mscheme: MM, n: int, d: int, dist: MatrixDistribution):
        super().__init__()
        self.mscheme = mscheme
        self.n = n
        self.d = d
        self.dist = dist

    def keygen(self):
        n, d = self.n, self.d
        a, f = (self.mscheme.sample(n) for i in range(2))
        F = self.mscheme.sample(n, n+1)
        D = self.mscheme.sample_from(self.dist)
        D1 = D >> G1
        ppk, psk = BenignPS(self.mscheme).gen(D1)
        sk = FFHR21_SK(a=a, f=f, F=F, psk=psk)

        aTD = a.T() * D
        Dx = D | aTD

        pk = FFHR21_PK(D=D1,
                       aTD=aTD >> G1,
                       fTD=(f.T() * D) >> GT,
                       FTD=(F.T() * D) >> G1,
                       FDx=(F * Dx) >> G1,
                       ppk=ppk
                       )
        return pk, sk

    def encrypt(self, pk: FFHR21_PK, msg):
        n, d = self.n, self.d
        r = self.mscheme.sample(d)
        u = pk.D * r

        pi = BenignPS(self.mscheme).prove(pk.ppk, u, r)
        p = pk.aTD * r + msg
        x = u | p
        y = pk.fTD * r + (x.pair_with(pk.FTD) * r)

        return FFHR21_CIPHER(u=u, p=p, y=y, pi=pi)

    def decrypt(self, pk, sk, c):
        x = c.u | c.p
        msg = c.p - sk.a.T() * c.u
        y = (sk.f.T() * c.u >> GT) + (sk.F * x).pair_with(c.u)
        b1 = y == c.y
        b2 = BenignPS(self.mscheme).verify(sk.psk, c.u, c.pi)
        return msg if b1 and b2 else None

    def rand(self, pk, c):
        r = self.mscheme.sample(self.d)
        x = c.u | c.p
        u_hat = c.u + pk.D * r
        p_hat = c.p + pk.aTD * r
        y_hat = (pk.fTD * r) + x.pair_with(pk.FTD * r) + \
            (pk.FDx * r).pair_with(u_hat)
        pi = BenignPS(self.mscheme).peval(pk.ppk, c.u, c.pi, r)

        return FFHR21_CIPHER(u=u_hat, p=p_hat, y=c.y + y_hat, pi=pi)
