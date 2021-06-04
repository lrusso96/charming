from utils.matrix import MM_GROUP, MM
from ps.benign import BenignPS


import unittest


class TestBenignPS(unittest.TestCase):

    def test_basic(self):

        curve = 'SS512'
        group = PairingGroup(curve)
        g1 = group.random(G1)
        g2 = g1
        gT = pair(g1, g2)
        MM_SS512 = MM_GROUP(G=group, g1=g1, g2=g2, gT=gT)

        # Define security parameters
        n = 6
        d = 2

        mscheme = MM(MM_SS512)

        ps = BenignPS(mscheme)

        D = mscheme.sample(n, d) >> G1

        pk, sk = ps.gen(D)
        r = mscheme.sample(d)
        u = D * r
        pi = ps.prove(pk, u, r)

        self.assertEqual(True, ps.verify(sk, u, pi), "Should be equal")

        r = mscheme.sample(d)
        pi2 = ps.peval(pk, u, pi, r)
        u = u + D*r
        self.assertEqual(True, ps.verify(sk, u, pi2), "Should be equal")


if __name__ == '__main__':
    unittest.main()
