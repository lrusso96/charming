from charm.toolbox.PKEnc import PKEnc
from charm.toolbox.ecgroup import ECGroup, ZR, G

from collections import namedtuple


CS98_PK = namedtuple('CS98_PK', ['g1', 'g2', 'c', 'd', 'h', 'H'])
CS98_SK = namedtuple('CS98_SK', ['x1', 'x2', 'y1', 'y2', 'z'])
CS98_CIPHER = namedtuple('CS98_CIPHER', ['u1', 'u2', 'e', 'v'])


class InputError(Exception):
    pass


class CS98(PKEnc):

    def __init__(self, curve):
        PKEnc.__init__(self)
        self.group = ECGroup(curve)

    def keygen(self, secparam):
        g1, g2 = (self.group.random(G) for i in range(2))
        x1, x2, y1, y2, z = (self.group.random() for i in range(5))
        c = (g1 ** x1) * (g2 ** x2)
        d = (g1 ** y1) * (g2 ** y2)
        h = g1 ** z
        pk = CS98_PK(g1=g1, g2=g2, c=c, d=d, h=h, H=self.group.hash)
        sk = CS98_SK(x1=x1, x2=x2, y1=y1, y2=y2, z=z)
        return (pk, sk)

    def encrypt(self, pk, m):
        if len(m) != 20:
            error_msg = f"Curve prime192v1 supports message of length 20 bytes only, but this message has length {len(m)}"
            raise InputError(error_msg)
        r = self.group.random()
        u1 = pk.g1 ** r
        u2 = pk.g2 ** r
        e = self.group.encode(m) * (pk.h ** r)
        alpha = self.group.hash((u1, u2, e))
        v = (pk.c ** r) * (pk.d ** (r * alpha))
        return CS98_CIPHER(u1=u1, u2=u2, e=e, v=v)

    def decrypt(self, pk, sk, c):
        alpha = pk.H((c.u1, c.u2, c.e))
        v_prime = (c.u1 ** (sk.x1 + (sk.y1 * alpha))) * \
            (c.u2 ** (sk.x2 + (sk.y2 * alpha)))
        if (c.v == v_prime):
            return self.group.decode(c.e / (c.u1 ** sk.z))
