from charm.toolbox.pairinggroup import PairingGroup, G1, G2, pair
from utils.matrix import DK_MDDH, MM_GROUP, MM
from pke.ffhr19 import FFHR19


import unittest


class TestFFHR19(unittest.TestCase):

    def test_basic(self):
        # create a MM_GROUP on Curve MNT159
        curve = 'MNT159'
        group = PairingGroup(curve)
        g1 = group.random(G1)
        g2 = group.random(G2)
        gT = pair(g1, g2)
        MM_SS512 = MM_GROUP(G=group, g1=g1, g2=g2, gT=gT)

        # Define security parameters
        k = 2
        D3_2_MDDH = DK_MDDH(k, MM_SS512)

        mscheme = MM(MM_SS512)
        scheme = FFHR19(mscheme, k, D3_2_MDDH)

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
