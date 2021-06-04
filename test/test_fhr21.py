from charm.toolbox.pairinggroup import PairingGroup, G1, pair
from utils.matrix import D_MDDH, MM_GROUP, MM
from pke.fhr21 import FHR21


import unittest


class TestFHR21(unittest.TestCase):

    def test_basic(self):
        # create a MM_GROUP on Curve SS512
        curve = 'SS512'
        group = PairingGroup(curve)
        g1 = group.random(G1)
        g2 = g1
        gT = pair(g1, g2)
        MM_SS512 = MM_GROUP(G=group, g1=g1, g2=g2, gT=gT)

        # Define security parameters
        n = 6
        d = 2
        D6_2_MDDH = D_MDDH(n, d, MM_SS512)

        mscheme = MM(MM_SS512)
        scheme = FHR21(mscheme, n, d, D6_2_MDDH)

        pk, sk = scheme.keygen()
        msg = mscheme.sample(1, gtype=G1)
        c = scheme.encrypt(pk, msg)
        m = scheme.decrypt(pk, sk, c)
        self.assertEqual(m, msg, "Should be equal")

        c2 = scheme.rand(pk, c)
        m = scheme.decrypt(pk, sk, c2)
        self.assertEqual(m, msg, "Should be after rand too")


if __name__ == '__main__':
    unittest.main()
